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

    Eigen::MatrixXd _k_matrix = 1000. * Eigen::MatrixXd::Identity(_robot->model().getJointNum(), _robot->model().getJointNum());

    _torque_limits.reset(new TorqueLimits(_tau_max, _tau_min));
    _joint_task.reset(new JointImpedanceCtrl(_q, _robot->model()));
    _joint_task->setStiffness(_k_matrix);
    _joint_task->setDamping(_k_matrix*0.01);
    _joint_task->useInertiaMatrix(false);
    _joint_task->update(_q);

    QPOases_sot::Stack stack_of_tasks;
    stack_of_tasks.push_back(_joint_task);

    _solver.reset(new QPOases_sot(stack_of_tasks, _torque_limits, 2e12));

    return true;
}

void QPPVMPlugin::QPPVMControl()
{
    _robot->model().getEffortLimits(_tau_max);
    _tau_min = -_tau_max;
    _tau_max = _tau_max-_h;
    _tau_min = _tau_min-_h;

    _joint_task->setReference(_q_ref);
    _torque_limits->setTorqueLimits(_tau_max, _tau_min);

    _torque_limits->update(_q);
    _joint_task->update(_q);

    if(!_solver->solve(_tau_d)){
        _tau_d.setZero(_tau_d.size());
        //std::cout<<"SOLVER ERROR!"<<std::endl;
        printf("SOLVER ERROR! \n");
    }

    _tau_d.setZero(_tau_d.size());
    _tau_d = _tau_d + _h;
}

void QPPVMPlugin::on_start(double time)
{
    _start_time = time;
}


void QPPVMPlugin::control_loop(double time, double period)
{
    sense();
    _robot->model().computeGravityCompensation(_h);
//     _h.setZero(_h.size());

    if(!_homing_done)
    {
        if( (time - _start_time) <= _homing_time ){
            _q_ref = _q0 + 0.5*(1-std::cos(3.1415*(time - _start_time)/_homing_time))*(_q_home-_q0);

            _robot->model().setJointPosition(_q_ref);
            _robot->setReferenceFrom(_robot->model(), XBot::Sync::Position);
        }
        else{
            //std::cout<<"Homing Done! Starting QPPVM Control!"<<std::endl;
            printf("Homing Done! Starting QPPVM Control! \n");
            _homing_done = true;}
    }

    QPPVMControl();

    if(_homing_done)
    {


//         _robot->setStiffness(_k);
//         _robot->setDamping(_d);
        // set the joint effort on the model and then synchronize the effort on the robot
        _robot->model().setJointEffort(_tau_d);

        _robot->setReferenceFrom(_robot->model(), XBot::Sync::Effort);
    }

    _matlogger->add("tau_desired", _tau_d);

//     _robot->printTracking();
//    _robot->move();
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

