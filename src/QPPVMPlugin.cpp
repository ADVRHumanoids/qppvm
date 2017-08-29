/*
 * Copyright (C) 2016 IIT-ADVR
 * Author: Arturo Laurenzi
 * email:  arturo.laurenzi@iit.it
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

#include <QPPVM_RT_plugin/QPPVMPlugin.h>
#include <qpOASES.hpp>

// #include <rtdk.h>
// #define DPRINTF rt_printf

REGISTER_XBOT_PLUGIN(QPPVMPlugin, demo::QPPVMPlugin)

using namespace demo;
using namespace OpenSoT::constraints::torque;
using namespace OpenSoT::tasks::torque;
using namespace OpenSoT::solvers;

QPPVMPlugin::QPPVMPlugin()
{

}

bool QPPVMPlugin::init_control_plugin(  std::string path_to_config_file,
                                        XBot::SharedMemory::Ptr shared_memory,
                                        XBot::RobotInterface::Ptr robot)
{
    _matlogger = XBot::MatLogger::getLogger("/tmp/qppvm_log");

    _robot = robot;
    _model = XBot::ModelInterface::getModel(path_to_config_file);
    _model->syncFrom(*_robot);
    
    _model->initLog(_matlogger, 30000);

    _model->getEffortLimits(_tau_max_const);
    _tau_min_const = -_tau_max_const;
    std::cout<<"tau_max_const: "<<_tau_max_const<<std::endl;

    _tau_d.resize(_model->getJointNum());
    _tau_d.setZero(_tau_d.size());

    sense();
    _model->computeNonlinearTerm(_h);
    _tau_max = _tau_max_const-_h;
    _tau_min = _tau_min_const-_h;

    _model->getRobotState("home", _q_home);
    _q = _q_home;
    _q_ref = _q;

    _k.setZero(_robot->getJointNum());
    _d.setZero(_robot->getJointNum());

    Eigen::VectorXd k0, d0;
    _model->getStiffness(k0);
    _model->getDamping(d0);


//     _k[_model->getDofIndex("j_arm2_5")] = k0[_model->getDofIndex("j_arm2_5")];
//     _k[_model->getDofIndex("j_arm2_6")] = k0[_model->getDofIndex("j_arm2_6")];
//     _k[_model->getDofIndex("j_arm2_7")] = k0[_model->getDofIndex("j_arm2_7")];
//     _k[_model->getDofIndex("j_arm1_5")] = k0[_model->getDofIndex("j_arm1_5")];
//     _k[_model->getDofIndex("j_arm1_6")] = k0[_model->getDofIndex("j_arm1_6")];
//     _k[_model->getDofIndex("j_arm1_7")] = k0[_model->getDofIndex("j_arm1_7")];
// 
//     _d[_model->getDofIndex("j_arm2_5")] = d0[_model->getDofIndex("j_arm2_5")];
//     _d[_model->getDofIndex("j_arm2_6")] = d0[_model->getDofIndex("j_arm2_6")];
//     _d[_model->getDofIndex("j_arm2_7")] = d0[_model->getDofIndex("j_arm2_7")];
//     _d[_model->getDofIndex("j_arm1_5")] = d0[_model->getDofIndex("j_arm1_5")];
//     _d[_model->getDofIndex("j_arm1_6")] = d0[_model->getDofIndex("j_arm1_6")];
//     _d[_model->getDofIndex("j_arm1_7")] = d0[_model->getDofIndex("j_arm1_7")];


    Eigen::MatrixXd _k_matrix = k0.asDiagonal();
    Eigen::MatrixXd _d_matrix = d0.asDiagonal();

    _k_matrix *= Eigen::VectorXd(_k_matrix.size()).setConstant(10.0).asDiagonal();
    _d_matrix *= Eigen::VectorXd(_d_matrix.size()).setConstant(100.0).asDiagonal();
    
    std::cout<<"_d_matrix: "<<_d_matrix<<std::endl;
    std::cout<<"_k_matrix: "<<_k_matrix<<std::endl;

    _torque_limits.reset(new TorqueLimits(_tau_max, _tau_min));
    _joint_task.reset(new JointImpedanceCtrl(_q_home, *_model));
    _joint_task->setStiffness(_k_matrix);
    _joint_task->setDamping(_d_matrix);
    std::vector<bool> active_joint_mask = _joint_task->getActiveJointsMask();
    for(unsigned int i = 0; i < active_joint_mask.size(); ++i)
        active_joint_mask[i] = false;
    active_joint_mask[_model->getDofIndex("torso_yaw")] = true;
    active_joint_mask[_model->getDofIndex("j_arm1_1")] = true;
    active_joint_mask[_model->getDofIndex("j_arm1_2")] = true;
    active_joint_mask[_model->getDofIndex("j_arm1_3")] = true;
    active_joint_mask[_model->getDofIndex("j_arm1_4")] = true;
    active_joint_mask[_model->getDofIndex("j_arm2_1")] = true;
    active_joint_mask[_model->getDofIndex("j_arm2_2")] = true;
    active_joint_mask[_model->getDofIndex("j_arm2_3")] = true;
    active_joint_mask[_model->getDofIndex("j_arm2_4")] = true;
   // _joint_task->setActiveJointsMask(active_joint_mask);
    _joint_task->useInertiaMatrix(true);
    _joint_task->update(_q);

//    _model->getJointLimits(_q_min, _q_max);
//    Eigen::VectorXd q_range = _q_max - _q_min;
//    _q_max -= q_range * 0.1;
//    _q_min += q_range * 0.1;
//    std::cout << "QMIN: " << _q_min.transpose() << std::endl;
//    std::cout << "QHOME: " << _q_home.transpose() << std::endl;
//    std::cout << "QMAX: " << _q_max.transpose() << std::endl;


    _ee_task_left.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("LEFT_ARM", 
                                                                       _q, 
                                                                       *_model, 
                                                                       _robot->chain("left_arm").getTipLinkName(),
                                                                       "world"
                                                                       ) );
    _ee_task_left_pos.reset(new OpenSoT::SubTask(_ee_task_left, OpenSoT::Indices::range(0,2)));
    
    Eigen::MatrixXd Kc(6,6); Kc.setIdentity(6,6); Kc = 100.*Kc;
    Eigen::MatrixXd Dc(6,6); Dc.setIdentity(6,6); Dc = 10.0*Dc; 
    _ee_task_left->setStiffnessDamping(Kc, Dc);
    for(unsigned int i = 0; i < active_joint_mask.size(); ++i)
        active_joint_mask[i] = false;
    active_joint_mask[_model->getDofIndex("torso_yaw")] = true;
    active_joint_mask[_model->getDofIndex("j_arm1_1")] = true;
    active_joint_mask[_model->getDofIndex("j_arm1_2")] = true;
    active_joint_mask[_model->getDofIndex("j_arm1_3")] = true;
    active_joint_mask[_model->getDofIndex("j_arm1_4")] = true;
    //_ee_task_left->setActiveJointsMask(active_joint_mask);
    _ee_task_left->useInertiaMatrix(true);
    _ee_task_left_pos->update(_q);
    
    _ee_task_right.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("RIGHT_ARM", 
                                                                       _q, 
                                                                       *_model, 
                                                                       _robot->chain("right_arm").getTipLinkName(),
                                                                       "world"
                                                                       ) );
    _ee_task_right_pos.reset(new OpenSoT::SubTask(_ee_task_right, OpenSoT::Indices::range(0,2)));
    
    _ee_task_right->setStiffnessDamping(Kc, Dc);
    for(unsigned int i = 0; i < active_joint_mask.size(); ++i)
        active_joint_mask[i] = false;
    active_joint_mask[_model->getDofIndex("torso_yaw")] = true;
    active_joint_mask[_model->getDofIndex("j_arm2_1")] = true;
    active_joint_mask[_model->getDofIndex("j_arm2_2")] = true;
    active_joint_mask[_model->getDofIndex("j_arm2_3")] = true;
    active_joint_mask[_model->getDofIndex("j_arm2_4")] = true;
    //_ee_task_right->setActiveJointsMask(active_joint_mask);
    _ee_task_right->useInertiaMatrix(true);
    _ee_task_right_pos->update(_q);
    
//    _elbow_task_left.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("LEFT_ELBOW",
//                                                                          _q_home,
//                                                                           *_model,
//                                                                           _robot->chain("left_arm").getUrdfLinks()[4]->name,
//                                                                          "world"
//                                                                         ) );
    
//    _elbow_task_right.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("RIGHT_ELBOW",
//                                                                          _q_home,
//                                                                           *_model,
//                                                                           _robot->chain("right_arm").getUrdfLinks()[4]->name,
//                                                                          "world"
//                                                                         ) );
    

//    _joint_limits.reset( new OpenSoT::constraints::torque::JointLimits(_q_home, _q_max, _q_min, *_model) );
//    _joint_limits->setGains(k0*10, d0*20);
//    _joint_limits->update(_q_home);

//     _joint_task->getConstraints().push_back(_torque_limits);
//     _joint_task->useInertiaMatrix(false);
     

    _autostack = ( (_ee_task_left_pos + _ee_task_right_pos) 
                    /*/ (_elbow_task_left + _elbow_task_right)*/ 
                    / _joint_task ) << _torque_limits;
    
    
//      QPOases_sot::Stack stack_of_tasks;
//      stack_of_tasks.push_back(_ee_task_left);
//      stack_of_tasks.push_back(_joint_task);
     

    _solver.reset(new QPOases_sot(_autostack->getStack(), _autostack->getBounds(), 1.0 ) );
    
//     qpOASES::Options opt;
//     opt.setToReliable();
//     opt.printLevel = qpOASES::PL_LOW;
//     
//     _solver->setOptions(0, opt);
//     _solver->setOptions(1, opt);

    return true;
}

void QPPVMPlugin::QPPVMControl()
{
    _tau_max = _tau_max_const-_h;
    _tau_min = _tau_min_const-_h;
    _torque_limits->setTorqueLimits(_tau_max, _tau_min);


//    _q_ref = _q_home;
//    _joint_task->setReference(_q_ref);


     //_torque_limits->update(_q);
     //_joint_task->update(_q);
     //_ee_task_left->update(_q);
     _autostack->update(_q);
    
    
     
//      std::cout << "Left task error: \n" << _ee_task_left->getSpringForce() + _ee_task_left->getDamperForce() << std::endl;
     _matlogger->add("left_error", _ee_task_left->getSpringForce() + _ee_task_left->getDamperForce());
     

    if(!_solver->solve(_tau_d)){
        _tau_d.setZero(_tau_d.size());
        std::cout<<"SOLVER ERROR!"<<std::endl;
    }
    
    _matlogger->add("left_task_error", _ee_task_left->getA()*_tau_d - _ee_task_left->getb());

    
    
    _tau_d = _tau_d + _h;
}

void QPPVMPlugin::on_start(double time)
{
    sense();
    
    _start_time = time;
    _robot->setStiffness(_k);
    _robot->setDamping(_d);
    _robot->move();

    Eigen::Affine3d left_ee_pose;
    _model->getPose(_ee_task_left->getDistalLink(), left_ee_pose);
    
    Eigen::Affine3d right_ee_pose;
    _model->getPose(_ee_task_right->getDistalLink(), right_ee_pose);
    
    _ee_task_left->setReference(left_ee_pose.matrix());
    _ee_task_right->setReference(right_ee_pose.matrix());
    _joint_task->setReference(_q);
    
    //_torque_limits->update(_q);
    //_joint_task->update(_q);
    //_ee_task_left->update(_q);
    _autostack->update(_q);
   
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Right task error: \n" << _ee_task_right->getSpringForce() + _ee_task_right->getDamperForce() << std::endl;
    std::cout << "Left task error: \n" << _ee_task_left->getSpringForce() + _ee_task_left->getDamperForce() << std::endl;
    std::cout << "Joint task error: \n" << _joint_task->getSpringForce() + _joint_task->getDampingForce() << std::endl;
    
    
}       


void QPPVMPlugin::control_loop(double time, double period)
{
    
    sense();
    _model->computeNonlinearTerm(_h);


    QPPVMControl();
    
    // set the joint effort on the model and then synchronize the effort on the robot
    _model->setJointEffort(_tau_d);

    _robot->setReferenceFrom(*_model, XBot::Sync::Effort);

    _matlogger->add("tau_desired", _tau_d);
    _matlogger->add("time_matlogger", time);

    
    _model->log(_matlogger, time);
    

    _robot->move();
}

void QPPVMPlugin::sense()
{
    syncFromMotorSide(_robot, _model);
    _model->getJointPosition(_q);
    _model->getJointVelocity(_dq);
}


bool QPPVMPlugin::close()
{
    _matlogger->flush();
}

void demo::QPPVMPlugin::syncFromMotorSide(XBot::RobotInterface::Ptr robot, XBot::ModelInterface::Ptr model)
{
    robot->getMotorPosition(_jidmap);
    model->setJointPosition(_jidmap);
    
    robot->getMotorVelocity(_jidmap);
    model->setJointVelocity(_jidmap);
    
    model->update();
}

