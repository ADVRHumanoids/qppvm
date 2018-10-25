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

#ifndef _INVDYN_PLUGIN_H_
#define _INVDYN_PLUGIN_H_

#include <std_msgs/Float64.h>
#include <XCM/XBotControlPlugin.h>
#include <OpenSoT/utils/InverseDynamics.h>
#include <OpenSoT/tasks/acceleration/Postural.h>
#include <OpenSoT/tasks/acceleration/Cartesian.h>
#include <OpenSoT/constraints/acceleration/DynamicFeasibility.h>
#include <OpenSoT/utils/AutoStack.h>
#include <OpenSoT/floating_base_estimation/qp_estimation.h>
#include <atomic>
#include <dynamic_reconfigure_advr/server.h>
#include <QPPVM_RT_plugin/InvDynConfig.h>
#include <cartesian_interface/CartesianPlugin/Utils.h>

namespace XBotPlugin {
    
    
class DynReconfigure {
    
public:
    
    DynReconfigure();

    dynamic_reconfigure_advr::Server<InvDyn::InvDynConfig> _server;
    
    void cfg_callback(InvDyn::InvDynConfig &config, uint32_t level);
    
    std::atomic<double> _impedance_gain;
    std::atomic<double> _joints_lambda, _waist_lambda, _feet_lambda; 
    std::atomic<double> _joints_lambda2, _waist_lambda2, _feet_lambda2; 

};

/**
 * @brief InvDynPlugin XBot RT Plugin
 *
 **/
class InvDynPlugin : public XBot::XBotControlPlugin
{

public:

    virtual bool init_control_plugin(XBot::Handle::Ptr handle);

    virtual bool close(){ _logger->flush(); return true; }

    virtual void on_start(double time);

    virtual void on_stop(double time){}

    virtual ~InvDynPlugin(){}

protected:

    virtual void control_loop(double time, double period);

private:
    
    void sync_model(double period);
    void set_gains();
    void solve();
    void integrate(double period);
    void sync_cartesian_ifc(double time, double period);

    XBot::SharedObject<Eigen::Vector3d> _sh_fb_pos;
    XBot::SharedObject<Eigen::Vector6d> _sh_fb_vel;
    OpenSoT::floating_base_estimation::qp_estimation::Ptr _fbest;

    OpenSoT::utils::InverseDynamics::Ptr _invdyn;

    XBot::Handle::Ptr _xbot_handle;
    
    XBot::RobotInterface::Ptr _robot;
    XBot::ModelInterface::Ptr _model;
    XBot::ImuSensor::ConstPtr _imu;

    double _start_time;

    Eigen::VectorXd _q0, _k0, _d0, _k, _d;

    XBot::MatLogger::Ptr _logger;
    std::vector<std::string> _contact_links;
    Eigen::VectorXd _x, _q, _qdot, _tau, _qddot;

    OpenSoT::tasks::acceleration::Cartesian::Ptr _waist_task;
    OpenSoT::tasks::acceleration::Postural::Ptr _postural_task;
    OpenSoT::constraints::acceleration::DynamicFeasibility::Ptr _dyn_feas;
    std::vector<OpenSoT::tasks::acceleration::Cartesian::Ptr> _feet_cartesian;

    OpenSoT::AutoStack::Ptr _autostack;
    OpenSoT::solvers::iHQP::Ptr _solver;
    
    DynReconfigure _dynreconfig;

    XBot::Cartesian::CartesianInterfaceImpl::Ptr _ci;
    XBot::Cartesian::Utils::SyncFromIO::Ptr _sync_from_nrt;
    bool _first_sync_done;
};

}






#endif // _GCOMP_PLUGIN_H_
