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

    _robot->model().getEffortLimits(_tau_max);
    _tau_min = -_tau_max;

    _tau_d.resize(_robot->model().getJointNum());
    _tau_d.setZero(_tau_d.size());

    sense();
    //_robot->model().computeNonlinearTerm(_h);
    _h.setZero(_tau_d.size());
    _tau_max = _tau_max-_h;
    _tau_min = _tau_min-_h;

    _q_ref = _q;
    _q_home = _q;
    _q0 = _q;

    _q_home[0] = 0.0;

    _q_home[1] = 0.0;
    _q_home[2] = -0.3;
    _q_home[3] = -0.8;
    _q_home[4] = -0.8;
    _q_home[5] = 0.0;
    _q_home[6] = -0.8;
    _q_home[7] = 0.0;

    _q_home[8] = 0.0;
    _q_home[9] = 0.3;
    _q_home[10] = 0.8;
    _q_home[11] = 0.8;
    _q_home[12] = 0.0;
    _q_home[13] = 0.8;
    _q_home[14] = 0.0;

    _homing_done = false;

    _homing_time = 4.;

    _k.setZero(_robot->getJointNum());
    _d.setZero(_robot->getJointNum());

    Eigen::VectorXd k0, d0;
    _robot->model().getStiffness(k0);
    _robot->model().getDamping(d0);

    Eigen::MatrixXd _k_matrix = k0.asDiagonal();
    Eigen::MatrixXd _d_matrix = d0.asDiagonal();
    
    std::cout<<"_d_matrix: "<<_d_matrix<<std::endl;

    _torque_limits.reset(new TorqueLimits(_tau_max, _tau_min));
    _joint_task.reset(new JointImpedanceCtrl(_q_home, _robot->model()));
    _joint_task->setStiffness(0.*_k_matrix);
    _joint_task->setDamping(.01*_d_matrix);
    _joint_task->useInertiaMatrix(true);
    _joint_task->update(_q);
    _robot->model().getJointLimits(_q_min, _q_max);
    Eigen::VectorXd q_range = _q_max - _q_min;
    _q_max -= q_range * 0.1;
    _q_min += q_range * 0.1;
    std::cout << "QMIN: " << _q_min.transpose() << std::endl;
    std::cout << "QHOME: " << _q_home.transpose() << std::endl;
    std::cout << "QMAX: " << _q_max.transpose() << std::endl;


    _ee_task_left.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("LEFT_ARM", 
                                                                       _q_home, 
                                                                       _robot->model(), 
                                                                       _robot->chain("left_arm").getTipLinkName(),
                                                                       "world"
                                                                       ) );
    Eigen::MatrixXd W(6,6); W.setIdentity(6,6); W = 1000*W; W(3,3) = 100; W(4,4) = 100; W(5,5) = 100;
    _ee_task_left->setStiffnessDamping(W, Eigen::VectorXd::Constant(6,.3).asDiagonal());
    _ee_task_left->update(_q_home);
    _ee_task_left->useInertiaMatrix(true);
    
    _ee_task_right.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("RIGHT_ARM", 
                                                                       _q_home, 
                                                                       _robot->model(), 
                                                                       _robot->chain("right_arm").getTipLinkName(),
                                                                       "world"
                                                                       ) );
    _ee_task_right->setStiffnessDamping(W, Eigen::VectorXd::Constant(6,.3).asDiagonal());;
    
    _ee_task_right->useInertiaMatrix(true);
    
    _elbow_task_left.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("LEFT_ELBOW", 
                                                                          _q, 
                                                                           _robot->model(), 
                                                                           _robot->chain("left_arm").getUrdfLinks()[4]->name,
                                                                          "world"
                                                                         ) );
    
    _elbow_task_right.reset( new OpenSoT::tasks::torque::CartesianImpedanceCtrl("RIGHT_ELBOW", 
                                                                          _q, 
                                                                           _robot->model(), 
                                                                           _robot->chain("right_arm").getUrdfLinks()[4]->name,
                                                                          "world"
                                                                         ) );
    

    _joint_limits.reset( new OpenSoT::constraints::torque::JointLimits(_q_home, _q_max, _q_min, _robot->model()) );
    _joint_limits->setGains(k0*10, d0*20);
    _joint_limits->update(_q_home);

//     _joint_task->getConstraints().push_back(_torque_limits);
//     _joint_task->setStiffnessDamping(k0.asDiagonal(), d0.asDiagonal());
//     _joint_task->useInertiaMatrix(false);
     

    _autostack = ( (_ee_task_left + _ee_task_right) 
                    /*/ (_elbow_task_left + _elbow_task_right)*/ 
                    / _joint_task ) << _torque_limits;
    
    
  //    QPOases_sot::Stack stack_of_tasks;
  //    stack_of_tasks.push_back(_ee_task_left);
  //    stack_of_tasks.push_back(_joint_task);
     

    _solver.reset(new QPOases_sot(_autostack->getStack(), _autostack->getBounds(), 1e5 ) );

    return true;
}

void QPPVMPlugin::QPPVMControl()
{
    _robot->model().getEffortLimits(_tau_max);
    _tau_min = -_tau_max;
    _tau_max = _tau_max-_h;
    _tau_min = _tau_min-_h;


    _q_ref = _q_home;

    _joint_task->setReference(_q_ref);
    _torque_limits->setTorqueLimits(_tau_max, _tau_min);

     //_torque_limits->update(_q);
     //_joint_task->update(_q);
     //_ee_task_left->update(_q);
    _autostack->update(_q);
     
     std::cout << "Left task error: \n" << _ee_task_left->getSpringForce() + _ee_task_left->getDamperForce() << std::endl;

    if(!_solver->solve(_tau_d)){
        _tau_d.setZero(_tau_d.size());
        std::cout<<"SOLVER ERROR!"<<std::endl;
        //printf("SOLVER ERROR! \n");
    }

    
    
    _tau_d = _tau_d + _h;
}

void QPPVMPlugin::on_start(double time)
{
    sense();
    
    _robot->model().getJointPosition(_q);
    
    
    
    _start_time = time;
    _robot->setStiffness(_k);
    _robot->setDamping(_d);
    _robot->move();

    Eigen::Affine3d left_ee_pose;
    _robot->model().getPose(_ee_task_left->getDistalLink(), left_ee_pose);
    
    Eigen::Affine3d right_ee_pose;
    _robot->model().getPose(_ee_task_right->getDistalLink(), right_ee_pose);
    
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
    _robot->model().computeNonlinearTerm(_h);
//     _h.setZero(_h.size());


    QPPVMControl();
    
    // set the joint effort on the model and then synchronize the effort on the robot
    _robot->model().setJointEffort(_tau_d);

    _robot->setReferenceFrom(_robot->model(), XBot::Sync::Effort);

    _matlogger->add("tau_desired", _tau_d);

   _robot->move();
}

void QPPVMPlugin::sense()
{
    _robot->sense();
    _robot->model().getJointPosition(_q);
    _robot->model().getJointVelocity(_dq);
}


bool QPPVMPlugin::close()
{

}

