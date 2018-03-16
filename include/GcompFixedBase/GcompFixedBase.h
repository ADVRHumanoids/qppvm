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

#ifndef _GCOMP_FIXED_BASE_PLUGIN_H_
#define _GCOMP_FIXED_BASE_PLUGIN_H_

#include <std_msgs/Float64.h>
#include <XCM/XBotControlPlugin.h>
#include <OpenSoT/tasks/acceleration/Postural.h>
#include <OpenSoT/tasks/acceleration/Cartesian.h>
#include <OpenSoT/constraints/acceleration/DynamicFeasibility.h>
#include <OpenSoT/tasks/force/CoM.h>
#include <OpenSoT/utils/AutoStack.h>
#include <atomic>
#include <functional>

namespace XBotPlugin {

/**
 * @brief GcompPlugin XBot RT Plugin
 *
 **/
class GcompFixedBase : public XBot::XBotControlPlugin
{

public:

    virtual bool init_control_plugin(XBot::Handle::Ptr handle);

    virtual bool close(){ _logger->flush(); return true; }

    virtual void on_start(double time);

    virtual void on_stop(double time);

    virtual ~GcompFixedBase(){}

protected:

    virtual void control_loop(double time, double period);

private:
    
    void callback(const std_msgs::Float64ConstPtr& msg, int leg_id);

    XBot::MatLogger::Ptr _logger;
    XBot::RobotInterface::Ptr _robot;
    XBot::ModelInterface::Ptr _model;
    XBot::ImuSensor::ConstPtr _imu;
    
    Eigen::VectorXd _gcomp, _tau_ref_leg, _leg_impedance;
    std::vector<Eigen::VectorXd> _klegs, _dlegs, _klegs_0, _dlegs_0, _tau_offset;
    
    std::vector<XBot::RosUtils::SubscriberWrapper::Ptr> _imp_sub;
    std::vector<std::unique_ptr<std::atomic<float>>> _impedance_gain;

};

}


#endif // _GCOMP_PLUGIN_H_


