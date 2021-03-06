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

#define TRJ_TIME 3.0


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

bool QPPVMPlugin::init_control_plugin(  XBot::Handle::Ptr handle)
{
    _matlogger = XBot::MatLogger::getLogger("/tmp/qppvm_log");
    
    _set_ref = false;

    _robot = handle->getRobotInterface();
    //_model = XBot::ModelInterface::getModel(path_to_config_file);
    _model = XBot::ModelInterface::getModel(
        "/home/centauro/advr-superbuild/configs/ADVR_shared/centauro/configs/config_centauro_fixed_wrists.yaml");
    //_model->syncFrom(*_robot);
    
    _model->initLog(_matlogger, 30000);

    _model->getEffortLimits(_tau_max_const);
//     _tau_max_const.setConstant(_model->getJointNum(), 20.0); //PAY ATTENTION
    _tau_min_const = -_tau_max_const;
    std::cout<<"tau_max_const: "<<_tau_max_const<<std::endl;

    _tau_d.resize(_model->getJointNum());
    _tau_d.setZero(_tau_d.size());

    //sense();
    _model->computeNonlinearTerm(_h);
    _tau_max = _tau_max_const-_h;
    _tau_min = _tau_min_const-_h;

    _model->getRobotState("home", _q_home);
    _model->setJointPosition(_q_home);
    _model->setJointVelocity(Eigen::VectorXd(_q_home.size()).setConstant(0.0));
    _model->update();
    
    _q = _q_home;
    _q_ref = _q;

    _k.setZero(_robot->getJointNum());
    _d.setZero(_robot->getJointNum());

    Eigen::VectorXd k0, d0;
    _robot->getStiffness(k0);
    _robot->getDamping(d0);
    
    _k[_robot->getDofIndex("j_arm1_5")] = k0[_robot->getDofIndex("j_arm1_5")];
    _k[_robot->getDofIndex("j_arm1_6")] = k0[_robot->getDofIndex("j_arm1_6")];
    _k[_robot->getDofIndex("j_arm1_7")] = k0[_robot->getDofIndex("j_arm1_7")];
    _k[_robot->getDofIndex("j_arm2_5")] = k0[_robot->getDofIndex("j_arm2_5")];
    _k[_robot->getDofIndex("j_arm2_6")] = k0[_robot->getDofIndex("j_arm2_6")];
    _k[_robot->getDofIndex("j_arm2_7")] = k0[_robot->getDofIndex("j_arm2_7")];
    
    _d[_robot->getDofIndex("j_arm1_5")] = d0[_robot->getDofIndex("j_arm1_5")];
    _d[_robot->getDofIndex("j_arm1_6")] = d0[_robot->getDofIndex("j_arm1_6")];
    _d[_robot->getDofIndex("j_arm1_7")] = d0[_robot->getDofIndex("j_arm1_7")];
    _d[_robot->getDofIndex("j_arm2_5")] = d0[_robot->getDofIndex("j_arm2_5")];
    _d[_robot->getDofIndex("j_arm2_6")] = d0[_robot->getDofIndex("j_arm2_6")];
    _d[_robot->getDofIndex("j_arm2_7")] = d0[_robot->getDofIndex("j_arm2_7")];


    Eigen::VectorXd _k_matrix;
    Eigen::VectorXd _d_matrix;
    _model->getStiffness(_k_matrix);
    _model->getDamping(_d_matrix);
    
    
    _k_matrix = Eigen::VectorXd(_k_matrix.size()).setConstant(5.0);
    _d_matrix = Eigen::VectorXd(_d_matrix.size()).setConstant(2.0);

    
    std::cout<<"_d_matrix: \n"<<_d_matrix.transpose()<<std::endl;
    std::cout<<"_k_matrix: \n"<<_k_matrix.transpose()<<std::endl;

    _torque_limits.reset(new TorqueLimits(_tau_max, _tau_min));
    
    _joint_task.reset(new JointImpedanceCtrl(_q_home, *_model));
    _joint_task->setStiffness(_k_matrix.asDiagonal());
    _joint_task->setDamping(_d_matrix.asDiagonal());
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
                                                                       "arm1_7",
                                                                       "world",
                                                                       OpenSoT::Indices::range(0,2)) );
    
    Eigen::MatrixXd Kc(6,6); Kc.setIdentity(6,6); Kc = 700.*Kc; //Kc(1,1) = 50.0;
    Eigen::MatrixXd Dc(6,6); Dc.setIdentity(6,6); Dc = 70.0*Dc; //Dc(1,1) = 1.; 
    _ee_task_left->setStiffnessDamping(Kc, Dc);
    _ee_task_left->useInertiaMatrix(true);
    _ee_task_left->update(_q);
     
    _ee_task_right.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("RIGHT_ARM", 
                                                                       _q, 
                                                                       *_model, 
                                                                       "arm2_7",
                                                                       "world",
                                                                       OpenSoT::Indices::range(0,2)) );
    Kc.setIdentity(6,6); Kc = 700.*Kc; 
    Dc.setIdentity(6,6); Dc = 70.0*Dc;  
    _ee_task_right->setStiffnessDamping(Kc, Dc);
    _ee_task_right->useInertiaMatrix(true);
    _ee_task_right->update(_q);    
    
   _elbow_task_left.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("LEFT_ELBOW",
                                                                         _q,
                                                                          *_model,
                                                                          "arm1_4",
                                                                         "world",
                                                                        OpenSoT::Indices::range(0,2)) );
    
   _elbow_task_right.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("RIGHT_ELBOW",
                                                                         _q,
                                                                          *_model,
                                                                          "arm2_4",
                                                                         "world",
                                                                        OpenSoT::Indices::range(0,2)) );
    

//    _joint_limits.reset( new OpenSoT::constraints::torque::JointLimits(_q_home, _q_max, _q_min, *_model) );
//    _joint_limits->setGains(k0*10, d0*20);
//    _joint_limits->update(_q_home);

//     _joint_task->getConstraints().push_back(_torque_limits);
//     _joint_task->useInertiaMatrix(false);
     

    _autostack =  ((_ee_task_right + _ee_task_left)
                    // /(_elbow_task_left + _elbow_task_right))<< _torque_limits; 
                     /(_joint_task ))<< _torque_limits; 
                    
                    
    
//      QPOases_sot::Stack stack_of_tasks;
//      stack_of_tasks.push_back((_ee_task_left));
//      stack_of_tasks.push_back(_joint_task);
     

    _solver.reset(new QPOases_sot(_autostack->getStack(), _autostack->getBounds(), 1.0 ) );
    _solver->log(_matlogger);
    
//     qpOASES::Options opt;
//     opt.setToReliable();
//     opt.printLevel = qpOASES::PL_LOW;
//     
//     _solver->setOptions(0, opt);
//     _solver->setOptions(1, opt);

    return true;
}

void QPPVMPlugin::QPPVMControl(const double time)
{
    _tau_max = _tau_max_const-_h;
    _tau_min = _tau_min_const-_h;
    _torque_limits->setTorqueLimits(_tau_max, _tau_min);


//    _q_ref = _q_home;
//    _joint_task->setReference(_q_ref);


     //_torque_limits->update(_q);
     //_joint_task->update(_q);
//     _ee_task_left_pos->update(_q);
    
    
    if(_set_ref)
    {
        _ref = _start_pose;
        _ref.p.y(_start_pose.p.y()+0.15*std::sin(time-_start_time));
        _ref.p.z(_start_pose.p.z()+0.15*(1.0-std::cos(time-_start_time)));
        _ee_task_left->setReference(_ref);
    }
    
    
    _autostack->update(_q);
    _autostack->log(_matlogger);
    
//     if(time-_start_time >= 2.){
//       std::cout<<"_ee_task_left->getReference(): \n"<<_ee_task_left->getReference()<<std::endl;
//       std::cout<<"_ee_task_left->getActualPose(): \n"<<_ee_task_left->getActualPose()<<std::endl;
//       std::cout<<"_ee_task_left_pos->getb(): \n"<<_ee_task_left_pos->getb()<<std::endl;
//       std::cout<<"_ee_task_left->getb(): \n"<<_ee_task_left->getb()<<std::endl;
// //      std::cout<<"_ee_task_left->linearVelocityError: \n"<<_ee_task_left->linearVelocityError<<std::endl;
// //      std::cout<<"_ee_task_left->orientationVelocityError: \n"<<_ee_task_left->orientationVelocityError<<std::endl;
//       std::cout<<"_ee_task_left->positionError \n"<<_ee_task_left->positionError<<std::endl;
// //      std::cout<<"_ee_task_left->orientationError \n"<<_ee_task_left->orientationError<<std::endl;
//       std::cout<<"_ee_task_left->getSpringForce() \n"<<_ee_task_left->getSpringForce()<<std::endl;
//       std::cout<<"_ee_task_left->getDamperForce() \n"<<_ee_task_left->getDamperForce()<<std::endl;
// //     }
    
     
//      std::cout << "Left task error: \n" << _ee_task_left->getSpringForce() + _ee_task_left->getDamperForce() << std::endl;
     

    if(!_solver->solve(_tau_d)){
        _tau_d.setZero(_tau_d.size());
        std::cout<<"SOLVER ERROR!"<<std::endl;
    }
    _solver->log(_matlogger);
    
    
   
    _matlogger->add("tau_qp", _tau_d);
    
    _tau_d = _tau_d + _h;
    
     _matlogger->add("tau_desired", _tau_d);
}

void QPPVMPlugin::on_start(double time)
{
    sense();
    _model->computeNonlinearTerm(_h);
    
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
    
    
    _model->getPose(_ee_task_left->getDistalLink(), _start_pose);
    
    Eigen::VectorXd joint_spring, joint_damping,
                    left_spring, left_damping,
                    right_spring, right_damping;
    _joint_task->getSpringForce(joint_spring);
    _joint_task->getDampingForce(joint_damping);
    _ee_task_left->getSpringForce(left_spring);
    _ee_task_left->getDamperForce(left_damping);
    _ee_task_right->getSpringForce(right_spring);
    _ee_task_right->getDamperForce(right_damping);
    
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Right task error: \n" << right_spring + right_damping << std::endl;
    std::cout << "Left task error: \n" << left_spring + left_damping<< std::endl;
    std::cout << "Joint task error: \n" << joint_spring + joint_damping << std::endl;
    
    
}       


void QPPVMPlugin::control_loop(double time, double period)
{
    
    sense();
    _model->computeNonlinearTerm(_h);


    QPPVMControl(time);
    
    // set the joint effort on the model and then synchronize the effort on the robot
    _model->setJointEffort(_tau_d);

    _robot->setReferenceFrom(*_model, XBot::Sync::Effort);

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

