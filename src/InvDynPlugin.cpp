/*
 * Copyright (C) 2017 IIT-ADVR
 * Author:
 * email:
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#include <InvDynPlugin/InvDynPlugin.h>
#include <boost/make_shared.hpp>
#include <OpenSoT/SubTask.h>
#include <OpenSoT/constraints/GenericConstraint.h>
#include <OpenSoT/tasks/MinimizeVariable.h>
#include <OpenSoT/tasks/acceleration/Contact.h>
#include <malloc_finder/malloc_finder.h>

/* Specify that the class XBotPlugin::ForceAccExample is a XBot RT plugin with name "ForceAccExample" */
REGISTER_XBOT_PLUGIN_(XBotPlugin::InvDynPlugin)


const std::string floating_base_name = "base_link";

bool XBotPlugin::InvDynPlugin::init_control_plugin(XBot::Handle::Ptr handle)
{

    _xbot_handle = handle;

    _robot = handle->getRobotInterface();
    _logger = XBot::MatLogger::getLogger("/tmp/opensot_force_acc_example");

    _robot->getStiffness(_k0);
    _robot->getDamping(_d0);

    _imu = _robot->getImu().begin()->second;

    _model = XBot::ModelInterface::getModel(_xbot_handle->getPathToConfigFile());

    Eigen::VectorXd qhome;
    _model->getRobotState("home", qhome);
    _model->setJointPosition(qhome);
    _model->update();

    _contact_links = {"l_ankle", "r_ankle"};

    _invdyn = boost::make_shared<OpenSoT::utils::InverseDynamics>(_contact_links, *_model);

    using TypeConstraint = OpenSoT::constraints::GenericConstraint::Type;

    for(int i = 0; i < _contact_links.size(); i++){

        auto cl = _contact_links[i];
        

        _feet_cartesian.push_back(
            boost::make_shared<OpenSoT::tasks::acceleration::Cartesian>(cl + "_cartesian",
                                                                        *_model,
                                                                        cl,
                                                                        "world",
                                                                        _invdyn->getJointsAccelerationAffine())
                                 );



    }

    Eigen::VectorXd qddot_max;
    qddot_max.setConstant(_invdyn->getJointsAccelerationAffine().getOutputSize(), 50);

    auto qddot_lims = boost::make_shared<OpenSoT::constraints::GenericConstraint>("qddot_limit", 
                                                                                  _invdyn->getJointsAccelerationAffine(), 
                                                                                  qddot_max, -qddot_max,
                                                                                  TypeConstraint::CONSTRAINT);



    _postural_task = boost::make_shared<OpenSoT::tasks::acceleration::Postural>(*_model,
                                                                                _invdyn->getJointsAccelerationAffine());




    _dyn_feas = boost::make_shared<OpenSoT::constraints::acceleration::DynamicFeasibility>(
                                                            "DYN_FEAS",
                                                            *_model,
                                                            _invdyn->getJointsAccelerationAffine(),
                                                            _invdyn->getContactsWrenchAffine(),
                                                            _contact_links
                                                            );

    auto feet_cart_aggr = _feet_cartesian[0] + _feet_cartesian[1];

    _waist_task = boost::make_shared<OpenSoT::tasks::acceleration::Cartesian>("waist_task",
                                                                        *_model,
                                                                        floating_base_name,
                                                                        "world",
                                                                        _invdyn->getJointsAccelerationAffine());


    std::list<uint> pos_idx = {0,1,2};
    std::list<uint> or_idx = {3,4,5};

    _autostack = (( feet_cart_aggr + _waist_task  ) / ( _postural_task ));

    _autostack << _dyn_feas << qddot_lims;

    _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(),
                                                         _autostack->getBounds(),
                                                         1e4);



    _autostack->update(_x);

    if(!_solver->solve(_x))
    {
        Logger::error("Unable to solve!!! \n");
        return false;
    }
    
    
    /* Memory allocation */
    
    _model->getJointPosition(_q);
    _model->getJointVelocity(_qdot);
    _qddot.setZero(_model->getJointNum());
    _tau = _q;
    _x.setZero(_invdyn->getSerializer()->getSize());
    
    
    _sh_fb_pos = handle->getSharedMemory()->getSharedObject<Eigen::Vector3d>("/gazebo/floating_base_position");
    _sh_fb_vel = handle->getSharedMemory()->getSharedObject<Eigen::Vector6d>("/gazebo/floating_base_velocity");
    _sh_fb_pos.set(Eigen::Vector3d::Zero());
    _sh_fb_vel.set(Eigen::Vector6d::Zero());
    _fbest = boost::make_shared<OpenSoT::floating_base_estimation::qp_estimation>(_model, _imu, _contact_links);

    return true;
    
    XBot::Utils::MallocFinder::SetOnMalloc( XBot::Utils::MallocFinder::PrintBacktrace );
    XBot::Utils::MallocFinder::SetThrowOnMalloc(true);
    XBot::Utils::MallocFinder::SetThrowOnFree(true);
}

void XBotPlugin::InvDynPlugin::on_start(double time)
{
    _start_time = time;

    sync_model(0.0);

    for(auto ct : _feet_cartesian)
    {
        ct->resetReference();
    }

    _waist_task->resetReference();
}

void XBotPlugin::InvDynPlugin::control_loop(double time, double period)
{
    
    const bool ENABLE_FEEDBACK = true;
    
    sync_model(period);
    
    set_gains();

    XBot::Utils::MallocFinder::Enable();
    solve();
    XBot::Utils::MallocFinder::Disable();

    if(!ENABLE_FEEDBACK){
        
        integrate(period);
        _robot->setReferenceFrom(*_model, XBot::Sync::Position);
    }
    else
    {

        double k_factor = _dynreconfig._impedance_gain.load();
        double d_factor = std::sqrt(k_factor);

        _k = _k0 * k_factor;
        _d = _d0 * d_factor;

        _robot->setStiffness(_k);
        _robot->setDamping(_d);
        _robot->setReferenceFrom(*_model, XBot::Sync::Effort);
    }

    
    _logger->add("x", _x);
    
    _robot->move();
    
    

}


void XBotPlugin::InvDynPlugin::sync_model(double period)
{
        _model->syncFrom(*_robot, XBot::Sync::MotorSide, XBot::Sync::Position, XBot::Sync::Velocity);
        
        _fbest->update(period);

//         Eigen::Affine3d w_T_fb;
//         Eigen::Matrix3d w_R_fb;
//         Eigen::Vector6d fb_twist;
//         Eigen::Vector3d fb_pos;
// 
//         _sh_fb_pos.get(fb_pos);
//         _sh_fb_vel.get(fb_twist);
// 
//         _imu->getOrientation(w_R_fb);
// 
//         w_T_fb.linear() = w_R_fb;
//         w_T_fb.translation() = fb_pos;
// 
//         _model->setFloatingBaseState(w_T_fb, fb_twist);
//         _model->update();
}

void XBotPlugin::InvDynPlugin::set_gains()
{
    _postural_task->setLambda(_dynreconfig._joints_lambda.load());
    _waist_task->setLambda(_dynreconfig._waist_lambda.load());
}


void XBotPlugin::InvDynPlugin::solve()
{
    /* Update stack */
    _autostack->update(Eigen::VectorXd::Zero(0));

    /* Solve QP */
    _x.setZero(_x.size());
    if(!_solver->solve(_x))
    {
        Logger::error("Unable to solve!!! \n");
        return;
    }

    /* Compute torque */
    _invdyn->computedTorque(_x, _tau, _qddot);
    _model->setJointEffort(_tau);
}


void XBotPlugin::InvDynPlugin::integrate(double period)
{
    /* Update model */
    _model->getJointPosition(_q);
    _model->getJointVelocity(_qdot);

    _q.noalias() += 0.5*period*period*_qddot + period*_qdot;
    _qdot.noalias() += period*_qddot;

    _model->setJointPosition(_q);
    _model->setJointVelocity(_qdot);
    _model->update();
}


XBotPlugin::DynReconfigure::DynReconfigure()
{
    dynamic_reconfigure_advr::Server<QPPVM_RT_plugin::QppvmConfig>::CallbackType f;
    f = boost::bind(&DynReconfigure::cfg_callback, this, _1, _2);
    _server.setCallback(f);
    
    _impedance_gain.store(1.0);
    _joints_lambda.store(10.0);
    _waist_lambda.store(10.0);
    _stiffness_Feet_gain.store(0.0);
    _damping_Feet_gain.store(0.0);
}

void XBotPlugin::DynReconfigure::cfg_callback(QPPVM_RT_plugin::QppvmConfig& config, uint32_t level)
{
    _impedance_gain.store(config.impedance_gain);
    _stiffness_Waist_gain.store(config.stiffness_waist);
    _damping_Waist_gain.store(config.damping_waist);
    _stiffness_Feet_gain.store(config.stiffness_feet);
    _damping_Feet_gain.store(config.damping_feet);
    _joints_lambda.store(config.joints_lambda);
    _waist_lambda.store(config.waist_lambda);

    Logger::info(Logger::Severity::HIGH, "\nSetting impedance gain to %f \n", config.impedance_gain);
    Logger::info(Logger::Severity::HIGH, "Setting joints gain to %f \n", config.joints_lambda);
    Logger::info(Logger::Severity::HIGH, "Setting waist gain to %f \n", config.waist_lambda);
    Logger::info(Logger::Severity::HIGH, "Setting waist stiffness gain to %f and damping gain to %f\n", config.stiffness_waist, config.damping_waist);
    Logger::info(Logger::Severity::HIGH, "Setting feet stiffness gain to %f and damping gain to %f\n", config.stiffness_feet, config.damping_feet);
}