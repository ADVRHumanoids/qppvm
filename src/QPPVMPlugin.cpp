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
#include <boost/make_shared.hpp>
#include <QPPVM_RT_plugin/ForceOptimization_ifopt.h>
#include <OpenSoT/constraints/force/WrenchLimits.h>


#define TRJ_TIME 3.0


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

    _sh_fb_pos = handle->getSharedMemory()->getSharedObject<Eigen::Vector3d>("/gazebo/floating_base_position");
    _sh_fb_rot = handle->getSharedMemory()->getSharedObject<Eigen::Quaterniond>("/gazebo/floating_base_orientation");
    _sh_fb_vel = handle->getSharedMemory()->getSharedObject<Eigen::Vector6d>("/gazebo/floating_base_velocity");
    
    _sh_fb_pos.set(Eigen::Vector3d::Zero());
    _sh_fb_rot.set(Eigen::Quaterniond::Identity());
    _sh_fb_vel.set(Eigen::Vector6d::Zero());

    _matlogger = XBot::MatLogger::getLogger("/tmp/qppvm_log");
    XBot::Logger::SetVerbosityLevel(XBot::Logger::Severity::HIGH);

    _set_ref = false;

    _robot = handle->getRobotInterface();
    _model = XBot::ModelInterface::getModel("configs/ADVR_shared/user_example/centauro_simple_example.yaml");

   
    _model->initLog(_matlogger, 1000000);

    _tau_d.setZero(_model->getJointNum());
    
    _model->computeNonlinearTerm(_h);

    _model->update();
  
    std::vector<std::string> links_in_contact;
    

    for(int i = 0; i < _robot->legs(); i++)
        links_in_contact.push_back(_robot->leg(i).getTipLinkName());

    _force_opt = boost::make_shared<ForceOptimization>(_model, links_in_contact);

    return true;
}

void QPPVMPlugin::QPPVMControl(const double time)
{
        
     _matlogger->add("tau_desired", _tau_d);
     
}

void QPPVMPlugin::on_start(double time)
{
    sense();

    _start_time = time;
     

}


void QPPVMPlugin::control_loop(double time, double period)
{

    sense();
    
    QPPVMControl(time);

     
    std::vector<Eigen::Vector6d> Fopt;
    Eigen::VectorXd tau_opt;
    
    std::vector<Eigen::Vector6d> Fref_ifopt; 
    
    //TODO: set Fref_ifopt from ifopt_contacts_node
    //TODO: rotation matrices for friction cones from ifopt_contacts_node
    //TODO: working with CARTESIO
     
    for(int i = 0; i < _robot->legs(); i++)
        Fref_ifopt.push_back(Eigen::Vector6d::Zero(6));
    
    _tau_d = _h;
    
   _force_opt->compute(_tau_d, Fref_ifopt, Fopt, tau_opt);


    _matlogger->add("tau_opt", tau_opt);
    
    for(int i = 0; i < 4; i++)
    {
        _matlogger->add("wrench_opt_" + std::to_string(i+1), Fopt[i]);
    }
    
    _force_opt->log(_matlogger);
    
    _robot->setEffortReference(tau_opt.tail(42));

    _matlogger->add("time_matlogger", time);

    _model->log(_matlogger, time);

    _robot->move();
}

void QPPVMPlugin::sense()
{
     Eigen::Affine3d w_T_fb;
     Eigen::Vector6d fb_twist;

     w_T_fb.translation() = _sh_fb_pos.get();
     w_T_fb.linear() = _sh_fb_rot.get().toRotationMatrix();
     fb_twist = _sh_fb_vel.get();


     _model->setFloatingBaseState(w_T_fb, fb_twist);
     _model->update();


     _model->syncFrom(*_robot, XBot::Sync::Position, XBot::Sync::Velocity, XBot::Sync::MotorSide);
     _model->getJointPosition(_q);
     _model->getJointVelocity(_dq);
     _model->computeNonlinearTerm(_h);

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


