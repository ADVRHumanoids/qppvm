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

SHLIBPP_DEFINE_SHARED_SUBCLASS(QPPVMPlugin_factory, demo::QPPVMPlugin, XBot::XBotControlPlugin);

using namespace demo;
using namespace OpenSoT::constraints::torque;
using namespace OpenSoT::tasks::torque;
using namespace OpenSoT::solvers;

bool QPPVMPlugin::init_control_plugin(  std::string path_to_config_file,
                                        XBot::SharedMemory::Ptr shared_memory,
                                        XBot::RobotInterface::Ptr robot)
{
    _robot = robot;
    Eigen::VectorXd tau_max;
    _robot->model().getEffortLimits(tau_max);
    _tau_d.resize(_robot->model().getJointNum());
    _tau_d.setZero(_tau_d.size());

    sense();
    
    _q_ref = _q;
    _k.setZero(_robot->getJointNum());
    _d.setZero(_robot->getJointNum());
    
    Eigen::MatrixXd _k_matrix = 1000 * Eigen::MatrixXd::Identity(_robot->model().getJointNum(), _robot->model().getJointNum());

    _torque_limits.reset(new TorqueLimits(tau_max));
    _joint_task.reset(new JointImpedanceCtrl(_q, _robot->model()));
    _joint_task->setStiffness(_k_matrix);
    _joint_task->setDamping(_k_matrix*0.001);
    _joint_task->useInertiaMatrix(true);
    _joint_task->update(_q);

    QPOases_sot::Stack stack_of_tasks;
    stack_of_tasks.push_back(_joint_task);

    _solver.reset(new QPOases_sot(stack_of_tasks, _torque_limits, 2e12));
}

void QPPVMPlugin::control_loop(double time, double period)
{
    sense();
    
    _q_ref(_robot->getDofIndex(_robot->chain("left_arm").getJointId(0))) = -90. * M_PI / 180. ;
    _joint_task->setReference(_q_ref);
    
    _torque_limits->update(_q);
    _joint_task->update(_q);

    if(!_solver->solve(_tau_d))
        _tau_d.setZero(_tau_d.size());

    _robot->setStiffness(_k);
    _robot->setDamping(_d);
    
    // set the joint effort on the model and then synchronize the effort on the robot
    _robot->model().setJointEffort(_tau_d);
    _robot->setReferenceFrom(_robot->model(), XBot::Sync::Effort);
    
//     _robot->printTracking();
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
 
