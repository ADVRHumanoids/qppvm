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

#include <ForceAccPlugin/ForceAcc.h>
#include <boost/make_shared.hpp>
#include <OpenSoT/SubTask.h>
#include <OpenSoT/constraints/GenericConstraint.h>
#include <OpenSoT/tasks/acceleration/Contact.h>
#include <OpenSoT/constraints/force/FrictionCone.h>

/* Specify that the class XBotPlugin::ForceAccExample is a XBot RT plugin with name "ForceAccExample" */
REGISTER_XBOT_PLUGIN_(XBotPlugin::ForceAccExample)

std::string floating_base_name;

bool XBotPlugin::ForceAccExample::init_control_plugin(XBot::Handle::Ptr handle)
{
    
    _robot = handle->getRobotInterface();
    
    /* Robot-specific params (TBD param server?) */
    
    _contact_links = {"foot_fl", "foot_fr", "foot_hr", "foot_hl"};
    bool optimize_contact_torque = false;
    Eigen::MatrixXd contact_matrix(5,6);
    contact_matrix << 1, 0, 0, 0, 0, 0,
                      0, 1, 0, 0, 0, 0,
                      0, 0, 1, 0, 0, 0,
                      0, 0, 0, 1, 0, 0,
                      0, 0, 0, 0, 0, 1;
    
    if(_robot->getUrdf().name_ == "cogimon")
    {
        _contact_links = {"l_sole", "r_sole"};
        optimize_contact_torque = true;
        contact_matrix = Eigen::MatrixXd::Identity(6,6);
    }

    OpenSoT::constraints::force::FrictionCone::friction_cones friction_cones;
    
    for(auto cl : _contact_links)
    {
        friction_cones.emplace_back(cl, 0.3);
    }
    
    
    
    _logger = XBot::MatLogger::getLogger("/tmp/opensot_force_acc_example");
    
    _robot->getStiffness(_k);
    _robot->getDamping(_d);
    _k *= 0;
    _d *= 0;
    
    _imu = _robot->getImu().begin()->second;
    
    _model = XBot::ModelInterface::getModel(handle->getPathToConfigFile());
    _model_fbest = XBot::ModelInterface::getModel(handle->getPathToConfigFile());
    _model->getFloatingBaseLink(floating_base_name);
    
    _fb_estimator = std::make_shared<estimation::FloatingBaseEstimator>(_model, _imu, _contact_links, contact_matrix);
    
    Eigen::VectorXd qhome;
    _model->getRobotState("home", qhome);
    _model->setJointPosition(qhome);
    _model->update();
    
    _model->initLog(_logger, 10000);
    
    _sh_fb_pos = handle->getSharedMemory()->getSharedObject<Eigen::Vector3d>("/gazebo/floating_base_position");
    _sh_fb_vel = handle->getSharedMemory()->getSharedObject<Eigen::Vector3d>("/gazebo/floating_base_velocity");
    _sh_fb_pos.set(Eigen::Vector3d::Zero());
    _sh_fb_vel.set(Eigen::Vector3d::Zero());
    
    
    
    
    
    _wrench_value.assign(_contact_links.size(), Eigen::VectorXd::Zero(6));
    
    OpenSoT::OptvarHelper::VariableVector vars;
    vars.emplace_back("qddot", _model->getJointNum());
    
    for(auto cl : _contact_links){
        vars.emplace_back(cl, optimize_contact_torque ? 6 : 3); // put 6 for full wrench
    }
    
    OpenSoT::OptvarHelper opt(vars);
    
    _qddot = opt.getVariable("qddot");
    
    Eigen::VectorXd wrench_ub(6), wrench_lb(6);
    wrench_ub << 1000, 1000, 1000, 50, 50, 50;
    wrench_lb << -1000, -1000, 0, -50, -50, -50;
    
    std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_bounds;
    std::vector<OpenSoT::tasks::acceleration::Contact::Ptr> contact_tasks;
    OpenSoT::constraints::force::FrictionCone::Ptr friction_contraint;
    
    for(auto cl : _contact_links){
        
        _wrenches.emplace_back(opt.getVariable(cl) / 
                               OpenSoT::AffineHelper::Zero(opt.getSize(), optimize_contact_torque ? 0 : 3)
                              );
        
        _feet_cartesian.push_back(
            boost::make_shared<OpenSoT::tasks::acceleration::Cartesian>(cl + "_cartesian", 
                                                                        *_model, 
                                                                        cl, 
                                                                        "world", 
                                                                        _qddot)
                                 );
        
        contact_tasks.push_back(
             boost::make_shared<OpenSoT::tasks::acceleration::Contact>(cl + "_contact", 
                                                                        *_model, 
                                                                        cl, 
                                                                        _qddot, 
                                                                        contact_matrix
                                                                      )
            
        );
        
        _feet_cartesian.back()->setLambda(1);
        
        
        wrench_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"_bound", 
                                                                                             _wrenches.back(), 
                                                                                             wrench_ub, 
                                                                                             wrench_lb) 
                               );
        
        
        
    }
    
    friction_contraint = boost::make_shared<OpenSoT::constraints::force::FrictionCone>(_wrenches, *_model, friction_cones);
    
    
    
    
    
    _com_task = boost::make_shared<OpenSoT::tasks::force::CoM>(_wrenches, _contact_links, *_model);
    
    _postural_task = boost::make_shared<OpenSoT::tasks::acceleration::Postural>("POSTURAL", 
                                                                                *_model, 
                                                                                _qddot);
    
    _postural_task->setWeight(0.001 * Eigen::MatrixXd::Identity(_model->getActuatedJointNum(), 
                                                                _model->getActuatedJointNum()));
    
    _dyn_feas = boost::make_shared<OpenSoT::constraints::acceleration::DynamicFeasibility>("DYN_FEAS", 
                                                                                           *_model, 
                                                                                           _qddot, 
                                                                                           _wrenches,
                                                                                            _contact_links
                                                                                          );
    
    auto feet_cart_aggr = _feet_cartesian[0] + _feet_cartesian[1];
    for(int i = 2; i < _contact_links.size(); i++){
        feet_cart_aggr = feet_cart_aggr + _feet_cartesian.at(i);
    }
    
    auto feet_contact_aggr = contact_tasks[0] + contact_tasks[1];
    for(int i = 2; i < _contact_links.size(); i++){
        feet_contact_aggr = feet_contact_aggr + contact_tasks.at(i);
    }
    

    _waist_task = boost::make_shared<OpenSoT::tasks::acceleration::Cartesian>("waist_task", 
                                                                        *_model, 
                                                                        floating_base_name, 
                                                                        "world", 
                                                                        _qddot);
    
    _waist_task->setLambda(25);
                                 
    
    std::list<uint> pos_idx = {0,1,2};
    std::list<uint> or_idx = {3,4,5};
    auto waist_or = boost::make_shared<OpenSoT::SubTask>(_waist_task, or_idx);
    
    _autostack = ( ( _waist_task  ) / (  feet_cart_aggr + _postural_task  ) ) << 
        _dyn_feas;
    
        
    for(int i = 0; i < _contact_links.size(); i++){
        _autostack << wrench_bounds[i];
    }
    
    _autostack << friction_contraint;
    
    _solver = boost::make_shared<OpenSoT::solvers::QPOases_sot>(_autostack->getStack(), 
                                                                _autostack->getBounds(), 
                                                                1e4);
    
    
    return true;
}

void XBotPlugin::ForceAccExample::on_start(double time)
{

    _start_time = time;
    _qdot.setZero(_q.size());
    
    sync_model();
    
    for(auto ct : _feet_cartesian)
    {
        ct->resetReference();
    }
    
    _waist_task->resetReference();
    
    _model->getPointPosition(floating_base_name, Eigen::Vector3d::Zero(),_initial_com);
}

void XBotPlugin::ForceAccExample::control_loop(double time, double period)
{
    period = 0.001;
    
    const bool enable_torque_ctrl = true;
    const bool enable_feedback = true;
    
    if(enable_feedback){
        
        sync_model();
        
    }
    
    /* Set reference*/
    if( (time - _start_time) <= 2.0 ){
        _waist_task->setPositionReference(_initial_com - 0.1*Eigen::Vector3d::UnitZ()*(time - _start_time));
    }
    
    /* Update stack */
    _autostack->update(Eigen::VectorXd::Zero(0));
    _autostack->log(_logger);
    
    /* Solve QP */
    _x.setZero(_x.size());
    if(!_solver->solve(_x))
    {
        Logger::error("Unable to solve!!!");
        return;
    }
    
    /* Retrieve values from QP solution */
    _qddot.getValue(_x, _qddot_value);
    
    for(int i = 0; i < _wrenches.size(); i++){
        _wrenches[i].getValue(_x, _wrench_value[i]);
        _logger->add(_contact_links[i] + "_wrench", _wrench_value[i]);
    }
    
    /* Compute torques due to contacts */
    _tau_c.setZero(_model->getJointNum());
    for(int i = 0; i < _contact_links.size(); i++){
        _model->getJacobian(_contact_links[i], _Jtmp);
        _tau_c.noalias() += _Jtmp.transpose()*_wrench_value[i];
    }
    
    /* Set solution inside model */
    _model->setJointAcceleration(_qddot_value);
    _model->update();
    
    /* Compute full ID */
    _model->computeInverseDynamics(_tau);
    _tau -= _tau_c;
    _model->setJointEffort(_tau);
    
    if( !enable_feedback ){
        /* Update model */
        _model->getJointPosition(_q);
        _model->getJointVelocity(_qdot);
        
        _q.noalias() += 0.5*period*period*_qddot_value + period*_qdot;
        _qdot.noalias() += period*_qddot_value;
        
        _model->setJointPosition(_q);
        _model->setJointVelocity(_qdot);
        _model->update();
    }
    
    /* Log */
    _logger->add("tau", _tau);
    _logger->add("tau_c", _tau_c);
    _logger->add("qddot_value", _qddot_value);
    _logger->add("x", _x);
    
    /* Send commands to robot */
    if(enable_torque_ctrl){
        _robot->setStiffness(_k);
        _robot->setDamping(_d);
        _robot->setReferenceFrom(*_model, XBot::Sync::Position, XBot::Sync::Effort);
    }
    else{
        _robot->setReferenceFrom(*_model, XBot::Sync::Position, XBot::Sync::Effort);
    }
    
    _robot->move();
    _model->log(_logger, time);
    
    
    
}


void XBotPlugin::ForceAccExample::sync_model()
{
        _model->syncFrom(*_robot);

        _fb_estimator->update(0.001);
        _fb_estimator->log(_logger);
        
}

