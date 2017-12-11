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

#ifndef _GCOMP_PLUGIN_H_
#define _GCOMP_PLUGIN_H_

#include <std_msgs/Float64.h>
#include <XCM/XBotControlPlugin.h>
#include <OpenSoT/tasks/acceleration/Postural.h>
#include <OpenSoT/tasks/acceleration/Cartesian.h>
#include <OpenSoT/constraints/acceleration/DynamicFeasibility.h>
#include <OpenSoT/tasks/force/CoM.h>
#include <OpenSoT/utils/AutoStack.h>
#include <atomic>

namespace XBotPlugin {

/**
 * @brief GcompPlugin XBot RT Plugin
 *
 **/
class GcompPlugin : public XBot::XBotControlPlugin
{

public:

    virtual bool init_control_plugin(XBot::Handle::Ptr handle);

    virtual bool close(){ _logger->flush(); return true; }

    virtual void on_start(double time);

    virtual void on_stop(double time){}

    virtual ~GcompPlugin(){}

protected:

    virtual void control_loop(double time, double period);

private:

    void sync_model();

    void callback_impedance(const std_msgs::Float64ConstPtr& msg)
    {
        int data = msg->data;
        data = data > 0 ? data : 0;
        data = data < 100 ? data : 100;

        _impedance_coeff.store(data);
    }
    
    void callback_postural(const std_msgs::Float64ConstPtr& msg)
    {
        int data = msg->data;
        data = data > 0 ? data : 0;

        _postural_gain.store(data);
    }

    XBot::Handle::Ptr _xbot_handle;
    
    XBot::RobotInterface::Ptr _robot;
    XBot::ModelInterface::Ptr _model;
    XBot::ImuSensor::ConstPtr _imu;

    XBot::SharedObject<Eigen::Vector3d> _sh_fb_pos;
    XBot::SharedObject<Eigen::Vector3d> _sh_fb_vel;

    double _start_time;

    Eigen::VectorXd _q0, _k0, _d0, _k, _d;
    Eigen::Vector3d _initial_com;

    XBot::MatLogger::Ptr _logger;
    std::vector<std::string> _contact_links;
    Eigen::VectorXd _x, _qddot_value, _q, _qdot, _tau, _tau_c;
    Eigen::MatrixXd _Jtmp;
    std::vector<Eigen::VectorXd> _wrench_value;
    OpenSoT::AffineHelper _qddot;
    std::vector<OpenSoT::AffineHelper> _wrenches;

    OpenSoT::tasks::acceleration::Cartesian::Ptr _waist_task;
    OpenSoT::tasks::force::CoM::Ptr _com_task;
    OpenSoT::tasks::acceleration::Postural::Ptr _postural_task;
    OpenSoT::constraints::acceleration::DynamicFeasibility::Ptr _dyn_feas;
    std::vector<OpenSoT::tasks::acceleration::Cartesian::Ptr> _feet_cartesian;

    OpenSoT::AutoStack::Ptr _autostack;
    OpenSoT::solvers::QPOases_sot::Ptr _solver;
    
    XBot::JointIdMap _q_ref_map, _qdot_ref_map;
    Eigen::VectorXd _q_ref, _qdot_ref, _qddot_ref;

    std::atomic_int _impedance_coeff, _postural_gain;
    XBot::RosUtils::SubscriberWrapper::Ptr _sub_impedance, _sub_postural;

};

}






#endif // _GCOMP_PLUGIN_H_
