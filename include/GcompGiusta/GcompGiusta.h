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

#ifndef __QPPVM_GCOMP_PLUGIN_H__
#define __QPPVM_GCOMP_PLUGIN_H__

#include <XCM/XBotControlPlugin.h>
#include <OpenSoT/solvers/iHQP.h>
#include <OpenSoT/solvers/BackEndFactory.h>

#include <OpenSoT/utils/AutoStack.h>
#include <OpenSoT/SubTask.h>
#include <geometry_msgs/Twist.h>
#include <OpenSoT/tasks/force/CoM.h>
#include <QPPVM_RT_plugin/ForceOptimization.h>

#include <std_msgs/Float64.h>

#include <XBotInterface/Logger.hpp>
#include <atomic>

#include <dynamic_reconfigure_advr/server.h>
#include <QPPVM_RT_plugin/QppvmConfig.h>

#include <cartesian_interface/CartesianPlugin/Utils.h>

#include <dynamic_reconfigure_advr/server.h>
#include <QPPVM_RT_plugin/QppvmConfig.h>

namespace demo {

    class GcompGiusta : public XBot::XBotControlPlugin {

    public:

        virtual bool init_control_plugin(XBot::Handle::Ptr handle);
        virtual void on_start(double time);
        virtual void control_loop(double time, double period);
        virtual bool close();

    protected:

    private:
        
        void cfg_callback(QPPVM_RT_plugin::QppvmConfig &config, uint32_t level);
        
        dynamic_reconfigure_advr::Server<QPPVM_RT_plugin::QppvmConfig> _server;
        
        XBot::RobotInterface::Ptr _robot;
        XBot::ModelInterface::Ptr _model;
        XBot::ImuSensor::ConstPtr _imu;
        
        
        demo::ForceOptimization::Ptr _forza_giusta;
        
        std::vector<std::string> _contact_links;
        std::vector<Eigen::Vector6d> _Fc;
        
        Eigen::VectorXd _gcomp, _tau_d, _k0, _k;
        Eigen::MatrixXd _JC;
        
        std::atomic_int _imp_gain;

        

    };

}


#endif

bool demo::GcompGiusta::init_control_plugin(XBot::Handle::Ptr handle)
{
    _contact_links = {"wheel_1", "wheel_2", "wheel_3", "wheel_4"};
    
    dynamic_reconfigure_advr::Server<QPPVM_RT_plugin::QppvmConfig>::CallbackType f;
    f = boost::bind(&GcompGiusta::cfg_callback, this, _1, _2);
    _server.setCallback(f);
    
    _imp_gain.store(1000);
    
    _robot = handle->getRobotInterface();
    _imu = _robot->getImu().begin()->second;
    
    _model = XBot::ModelInterface::getModel(handle->getPathToConfigFile());
    _model->computeGravityCompensation(_gcomp);
    _tau_d = _gcomp;
    
    _robot->getStiffness(_k0);
    _k = _k0;
    
    _forza_giusta = boost::make_shared<ForceOptimization>(_model, _contact_links, false);
    
}

void demo::GcompGiusta::on_start(double time)
{
    
}


void demo::GcompGiusta::control_loop(double time, double period)
{
    

    
    _model->syncFrom(*_robot, XBot::Sync::Position, XBot::Sync::MotorSide);
    _model->setFloatingBaseState(_imu);
    _model->update();
    
    _model->computeGravityCompensation(_gcomp);
    
    _forza_giusta->compute(_gcomp, _Fc, _tau_d);
    
    _model->setJointEffort(_tau_d);
    
    _robot->setReferenceFrom(*_model, XBot::Sync::Effort); // , XBot::Sync::Impedance);
    
    _k = _imp_gain.load() / 1000.0 * _k0;
    
    _robot->setStiffness(_k);
    
    _robot->move();
}

bool demo::GcompGiusta::close()
{

}


void demo::GcompGiusta::cfg_callback(QPPVM_RT_plugin::QppvmConfig& config, uint32_t level)
{
    Logger::info(Logger::Severity::HIGH, "Setting impedance gain to %f\n", config.impedance_gain);
    _imp_gain.store(config.impedance_gain * 1000);
}





