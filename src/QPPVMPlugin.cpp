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

using namespace demo;
using namespace OpenSoT::constraints::torque;
using namespace OpenSoT::tasks::torque;
using namespace OpenSoT::solvers;

bool QPPVMPlugin::init_control_plugin(std::string path_to_config_file,
                                 XBot::SharedMemory::Ptr shared_memory,
                                 XBot::RobotInterface::Ptr robot)
{
    _robot->sense();

    _robot = robot;
    Eigen::VectorXd tau_max;
    _robot->getEffortLimits(tau_max);

    _robot->getJointPosition(_q);

    _torque_limits.reset(new TorqueLimits(tau_max));
    _joint_task.reset(new JointImpedanceCtrl(_q, _robot));

    QPOases_sot::Stack stack_of_tasks;
    stack_of_tasks.push_back(_joint_task);

    _solver.reset(new QPOases_sot(stack_of_tasks, _torque_limits, 2e11));
}

void QPPVMPlugin::control_loop(double time, double period)
{

}

bool QPPVMPlugin::close()
{

}
 
