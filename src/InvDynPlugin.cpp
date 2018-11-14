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

    
    _model = XBot::ModelInterface::getModel(_xbot_handle->getPathToConfigFile());
    _model->initLog(_logger, 1e5);

    Eigen::VectorXd qhome;
    _model->getRobotState("home", qhome);
    _model->setJointPosition(qhome);
    _model->update();
    
    
    
    if(!_model->getRobotState("torque_offset", _tau_offset))
        _tau_offset.setZero(_model->getJointNum());
    

    //_contact_links = {"l_wrist", "r_wrist"};
    _contact_links = {"l_sole", "r_sole"};
    
    std::string anchor_link = _contact_links[1];
    std::string fb_link; _model->getFloatingBaseLink(fb_link);
    
    Eigen::Affine3d p; _model->getPose(anchor_link,fb_link,p);
    _model->setFloatingBasePose(p);
    _model->update();
    
    _invdyn = boost::make_shared<OpenSoT::utils::InverseDynamics>(_contact_links, *_model);

    using TypeConstraint = OpenSoT::constraints::GenericConstraint::Type;

    for(int i = 0; i < _contact_links.size(); i++){

        auto cl = _contact_links[i];
        

        _feet_cartesian.push_back(
            boost::make_shared<OpenSoT::tasks::acceleration::Cartesian>(cl + "_cartesian",
                                                                        *_model,
                                                                        cl,
                                                                        "world",//"DWYTorso",
                                                                        _invdyn->getJointsAccelerationAffine())
                                 );



    }
    

    Eigen::VectorXd qddot_max;
    qddot_max.setConstant(_invdyn->getJointsAccelerationAffine().getOutputSize(), 300);

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
                                                                        "Waist",
                                                                        "world",
                                                                        _invdyn->getJointsAccelerationAffine());
    
    _com_task = boost::make_shared<OpenSoT::tasks::force::CoM>(_invdyn->getContactsWrenchAffine(),
                                                               _contact_links, *_model);
    g.setZero(); g[2] = -9.81;
    Eigen::Vector3d com_ref = _com_task->getLinearReference();
    _com_task->setLinearReference(g, Eigen::Vector3d::Zero(), com_ref);
    

    std::list<uint> pos_idx = {0,1,2};
    std::list<uint> or_idx = {3,4,5};

    _autostack = ( feet_cart_aggr/  (_waist_task)  /  _postural_task );
    _autostack << _dyn_feas << qddot_lims;
//     _autostack = boost::make_shared<OpenSoT::AutoStack>(_postural_task);
//     _autostack = ((feet_cart_aggr) / ( _postural_task ));
    _autostack <<  qddot_lims;

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
    _tau = _q; _tau.setZero(_tau.size());
    _x.setZero(_invdyn->getSerializer()->getSize());
    
//     _qdot_filtered = _qdot;
//     _alpha_filter = 1.;
    
    _sh_fb_pos = handle->getSharedMemory()->getSharedObject<Eigen::Vector3d>("/gazebo/floating_base_position");
    _sh_fb_vel = handle->getSharedMemory()->getSharedObject<Eigen::Vector6d>("/gazebo/floating_base_velocity");
    _sh_fb_pos.set(Eigen::Vector3d::Zero());
    _sh_fb_vel.set(Eigen::Vector6d::Zero());
    _imu = _robot->getImu().begin()->second;
    
    //_fbest = boost::make_shared<OpenSoT::floating_base_estimation::qp_estimation>(_model, _imu, _contact_links);
    //_fbest->update(0.0);
    
    Eigen::Affine3d fb_T_world, fb_T_anchor;
    _model->getFloatingBasePose(fb_T_world);
    _model->getPose(anchor_link, fb_link, fb_T_anchor);
    _fbest_kinematics = boost::make_shared<OpenSoT::floating_base_estimation::kinematic_estimation>(
        _model, anchor_link, fb_T_world.inverse()*fb_T_anchor);

    /* Init cartesian ifc */
    YAML::Node yaml_file = YAML::LoadFile(handle->getPathToConfigFile());
    XBot::Cartesian::ProblemDescription ik_problem(yaml_file["CartesianInterface"]["problem_description"], _model);
    _ci = std::make_shared<XBot::Cartesian::CartesianInterfaceImpl>(_model, ik_problem);
    _ci->enableOtg(0.002);
    _sync_from_nrt = std::make_shared<XBot::Cartesian::Utils::SyncFromIO>("/xbotcore/cartesian_interface", handle->getSharedMemory());
    
    return true;
    
    
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
    
    Eigen::Vector3d com_ref; _model->getCOM(com_ref);
    _com_task->setLinearReference(g, Eigen::Vector3d::Zero(), com_ref);

    
    _first_sync_done = false;
    
}

void XBotPlugin::InvDynPlugin::control_loop(double time, double period)
{
    
    _model->log(_logger,time);
    
    const bool ENABLE_FEEDBACK = true;
    
    sync_model(period);
    
    set_gains();
    
    sync_cartesian_ifc(time, period);
    
    Eigen::Vector3d com_ref = _com_task->getLinearReference();
    _com_task->setLinearReference(g, Eigen::Vector3d::Zero(), com_ref);
    
    solve(period);

    if(!ENABLE_FEEDBACK){
        
        integrate(period);
        _robot->setReferenceFrom(*_model, XBot::Sync::Position);
    }
    else
    {

        double k_factor = _dynreconfig._impedance_gain.load();
        double d_factor = std::sqrt(k_factor);

        _k = _k0 * k_factor;
        _d = _d0;// * d_factor;

        _robot->setStiffness(_k);
        _robot->setDamping(_d);
        _robot->setReferenceFrom(*_model, XBot::Sync::Effort);
    }

    
    
    _robot->move();
    
    _logger->add("x", _x);
    _autostack->log(_logger);
    

}


void XBotPlugin::InvDynPlugin::sync_model(double period)
{
        _model->syncFrom(*_robot, /*XBot::Sync::MotorSide,*/ XBot::Sync::Position, XBot::Sync::Velocity);
        
        //HERE WE FILTER and we apply the filtered velocities to the model
//     setFilter(period, 10.);
//     _model->getJointVelocity(_qdot);
//     _qdot_filtered.noalias() = _alpha_filter*_qdot + (1.-_alpha_filter)*_qdot_filtered;
//     _model->setJointVelocity(_qdot_filtered);
//     _model->update();
//     _logger->add("qdot", _qdot);
//     _logger->add("qdot_filtered", _qdot_filtered);
        
       // _fbest->update(period);
       // _fbest->log(_logger);
        
        _fbest_kinematics->update(true);
//          Eigen::Affine3d fb_T_world_corrected;
//          _model->getFloatingBasePose(fb_T_world_corrected);
//          fb_T_world_corrected.translation() = _fbest_kinematics->getFloatingBasePose().translation();
//          _model->setFloatingBasePose(fb_T_world_corrected);
//          _model->update();
//         
        

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
    _postural_task->setLambda(_dynreconfig._joints_lambda.load(), _dynreconfig._joints_lambda2.load());
    _waist_task->setLambda(_dynreconfig._waist_lambda.load(), _dynreconfig._waist_lambda2.load());
    _feet_cartesian[0]->setLambda(_dynreconfig._feet_lambda.load(), _dynreconfig._feet_lambda2.load());
    _feet_cartesian[1]->setLambda(_dynreconfig._feet_lambda.load(), _dynreconfig._feet_lambda2.load());
}

void XBotPlugin::InvDynPlugin::setFilter(const double period, const double cut_off_freq)
{
    _period = period;
    _cut_off_freq = cut_off_freq;
    _alpha_filter = 1-std::exp(-2.*M_PI*cut_off_freq*period);
}


void XBotPlugin::InvDynPlugin::solve(double period)
{   
    /* Update stack */
    _autostack->update(_x);

    /* Solve QP */
    _x.setZero(_x.size());
    if(!_solver->solve(_x))
    {
        Logger::error("Unable to solve!!! \n");
        return;
    }

    /* Compute torque */
    _invdyn->computedTorque(_x, _tau, _qddot);
    _tau += _tau_offset;

    /* Ramp the torque reference continuously */
    static int iter = 0;
    const int ITER_MAX = 1000;
    if(iter < ITER_MAX)
    {
        _tau *= iter/double(ITER_MAX);
        iter++;
    }
    
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
    dynamic_reconfigure_advr::Server<InvDyn::InvDynConfig>::CallbackType f;
    f = boost::bind(&DynReconfigure::cfg_callback, this, _1, _2);
    _server.setCallback(f);
    
    _impedance_gain.store(1.0);
    _joints_lambda.store(0.0);
    _waist_lambda.store(0.);
    _feet_lambda.store(0.);
    _joints_lambda2.store(.0);
    _waist_lambda2.store(.0);
    _feet_lambda2.store(.0);
}

void XBotPlugin::DynReconfigure::cfg_callback(InvDyn::InvDynConfig& config, uint32_t level)
{
    _impedance_gain.store(config.impedance_gain);
    _joints_lambda.store(config.joints_lambda_pos);
    _waist_lambda.store(config.waist_lambda_pos);
    _feet_lambda.store(config.feet_lambda_pos);
    _joints_lambda2.store(config.joints_lambda_vel);
    _waist_lambda2.store(config.waist_lambda_vel);
    _feet_lambda2.store(config.feet_lambda_vel);

    Logger::info(Logger::Severity::HIGH, "\nSetting impedance gain to %f \n", config.impedance_gain);
    
    Logger::info(Logger::Severity::HIGH, "Setting joints gain to %f, %f \n", 
                 config.joints_lambda_pos, config.joints_lambda_vel);
    
    Logger::info(Logger::Severity::HIGH, "Setting waist gain to %f, %f \n", config.waist_lambda_pos, config.waist_lambda_vel);
    
    Logger::info(Logger::Severity::HIGH, "Setting feet gain to %f, %f \n", config.feet_lambda_pos, config.feet_lambda_vel);
}

void XBotPlugin::InvDynPlugin::sync_cartesian_ifc(double time, double period)
{
    /* Sync cartesian references from ROS */
    
    if(!_first_sync_done)
    {
        if(_sync_from_nrt->try_reset(_model, time))
        {
            _first_sync_done = true;
            XBot::Logger::info(Logger::Severity::HIGH, "Resetting NRT CI \n");
        }
    }
    
    if(_first_sync_done)
    {
        _sync_from_nrt->try_sync(time, period, _ci, _model);
    }
    

    if(!_ci->update(time, period))
    {
        XBot::Logger::error("CartesianInterface: unable to solve \n");
        return;
    }
    
    
    /* Update qppvm references */
    Eigen::Affine3d T_ref;
    if(_ci->getPoseReference(_waist_task->getDistalLink(), T_ref) )
    {
        if(_ci->getBaseLink(_waist_task->getDistalLink()) == _waist_task->getBaseLink() )
        {
            _waist_task->setReference(T_ref);
        }
    }
    
    if(_ci->getPoseReference(_feet_cartesian[0]->getDistalLink(), T_ref) )
    {
        if(_ci->getBaseLink(_feet_cartesian[0]->getDistalLink()) == _feet_cartesian[0]->getBaseLink() )
        {
            _feet_cartesian[0]->setReference(T_ref);
        }
    }
    
    if(_ci->getPoseReference(_feet_cartesian[1]->getDistalLink(), T_ref) )
    {
        if(_ci->getBaseLink(_feet_cartesian[1]->getDistalLink()) == _feet_cartesian[1]->getBaseLink() )
        {
            _feet_cartesian[1]->setReference(T_ref);
        }
    }

}




