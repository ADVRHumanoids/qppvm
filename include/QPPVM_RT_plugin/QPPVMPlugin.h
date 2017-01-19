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

#ifndef __QPPVM_PLUGIN_H__
#define __QPPVM_PLUGIN_H__

#include <XCM/XBotControlPlugin.h>
#include <OpenSoT/solvers/QPOases.h>
#include <OpenSoT/tasks/torque/CartesianImpedanceCtrl.h>
#include <OpenSoT/tasks/torque/JointImpedanceCtrl.h>
#include <OpenSoT/constraints/torque/TorqueLimits.h>


#include <XBotInterface/Logger.hpp>

namespace demo {

    class QPPVMPlugin : public XBot::XBotControlPlugin {
      
    public:
        
        QPPVMPlugin();
        
        virtual bool init_control_plugin(std::string path_to_config_file, 
                                         XBot::SharedMemory::Ptr shared_memory, 
                                         XBot::RobotInterface::Ptr robot);
        virtual void control_loop(double time, double period);
        virtual bool close();
            
    protected:
        
    private:
        
        XBot::RobotInterface::Ptr _robot;
        
        XBot::MatLogger::Ptr _matlogger;

        OpenSoT::solvers::QPOases_sot::Ptr _solver;

        OpenSoT::constraints::torque::TorqueLimits::Ptr _torque_limits;
        OpenSoT::tasks::torque::JointImpedanceCtrl::Ptr _joint_task;

        Eigen::VectorXd _q;
        Eigen::VectorXd _dq;
        Eigen::VectorXd _q_ref;
        
        Eigen::VectorXd _k;
        Eigen::VectorXd _d;

        Eigen::VectorXd _tau_d;
        Eigen::VectorXd _h;

        Eigen::VectorXd _tau_max;
        Eigen::VectorXd _tau_min;
        

        void sense();

    };
    
}
#endif
