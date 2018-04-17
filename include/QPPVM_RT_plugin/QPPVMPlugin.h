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
#include <OpenSoT/solvers/iHQP.h>
#include <OpenSoT/solvers/BackEndFactory.h>
#include <OpenSoT/tasks/torque/CartesianImpedanceCtrl.h>
#include <OpenSoT/tasks/torque/JointImpedanceCtrl.h>
#include <OpenSoT/constraints/torque/TorqueLimits.h>
#include <OpenSoT/constraints/torque/JointLimits.h>
#include <OpenSoT/utils/AutoStack.h>
#include <OpenSoT/SubTask.h>
#include <geometry_msgs/Twist.h>
#include <OpenSoT/tasks/force/CoM.h>
#include <QPPVM_RT_plugin/ForceOptimization.h>
#include <ForceAccPlugin/FloatingBaseEstimation.h>
#include <std_msgs/Float64.h>

#include <XBotInterface/Logger.hpp>
#include <atomic>

#include <dynamic_reconfigure/server.h>
#include <QPPVM_RT_plugin/QppvmConfig.h>

namespace demo {

    class QPPVMPlugin : public XBot::XBotControlPlugin {

    public:

        QPPVMPlugin();

        virtual bool init_control_plugin(XBot::Handle::Ptr handle);
        virtual void on_start(double time);
        virtual void control_loop(double time, double period);
        virtual bool close();

    protected:

    private:

        typedef OpenSoT::tasks::torque::CartesianImpedanceCtrl CartesianImpedanceTask;
        
        void set_gains();

        void cart_stiffness_callback(const geometry_msgs::TwistConstPtr& msg, int id);
        void cart_damping_callback(const geometry_msgs::TwistConstPtr& msg, int id);
        void impedance_gain_callback(const std_msgs::Float64ConstPtr& msg);
        void feedback_gain_callback(const std_msgs::Float64ConstPtr& msg);
        void cfg_callback(QPPVM_RT_plugin::QppvmConfig &config, uint32_t level);
        std::atomic<double> _impedance_gain, _stiffness_gain, _damping_gain, _joints_gain;
        
        dynamic_reconfigure::Server<QPPVM_RT_plugin::QppvmConfig> _server;

        std::vector<XBot::RosUtils::SubscriberWrapper::Ptr> _cart_stiffness_sub, _cart_damping_sub;
        XBot::RosUtils::PublisherWrapper::Ptr _fb_pub;
        XBot::RosUtils::SubscriberWrapper::Ptr _impedance_gain_sub, _feedback_gain_sub;
        
        std::vector<Eigen::Vector6d> _Fopt;
        Eigen::VectorXd _tau_opt;

        double _start_time;

        bool _set_ref;

        void syncFromMotorSide(XBot::RobotInterface::Ptr robot, XBot::ModelInterface::Ptr model);
        XBot::JointIdMap _jidmap;

        XBot::RobotInterface::Ptr _robot;
        XBot::ModelInterface::Ptr _model;

        XBot::MatLogger::Ptr _matlogger;

        OpenSoT::solvers::iHQP::Ptr _solver;

        OpenSoT::constraints::torque::TorqueLimits::Ptr _torque_limits;
        OpenSoT::tasks::torque::JointImpedanceCtrl::Ptr _joint_task;
        std::vector<CartesianImpedanceTask::Ptr> _leg_impedance_task;
        OpenSoT::constraints::torque::JointLimits::Ptr _joint_limits;
        CartesianImpedanceTask::Ptr _waist,_ee_task_left,_ee_task_right;

        OpenSoT::AutoStack::Ptr _autostack;


        Eigen::VectorXd _q;
        Eigen::VectorXd _dq;
        Eigen::VectorXd _q_ref, _q0;
        Eigen::VectorXd _q_home;

        Eigen::VectorXd _k, _k_ref;
        Eigen::VectorXd _d, _d_ref;

        Eigen::VectorXd _tau_d, _tau_offset;
        Eigen::VectorXd _h;

        Eigen::VectorXd _tau_max;
        Eigen::VectorXd _tau_min;
        Eigen::VectorXd _tau_max_const;
        Eigen::VectorXd _tau_min_const;

        Eigen::VectorXd _q_max;
        Eigen::VectorXd _q_min;

        Eigen::MatrixXd _Kc, _Dc;

        bool _homing_done;

        double _homing_time;

        KDL::Frame _start_pose;
        KDL::Frame _ref;


        XBot::SharedObject<Eigen::Vector3d> _sh_fb_pos;
        XBot::SharedObject<Eigen::Quaterniond> _sh_fb_rot;
        XBot::SharedObject<Eigen::Vector6d> _sh_fb_vel;


         Eigen::VectorXd tau_f1;
         Eigen::VectorXd tau_f2;
         Eigen::VectorXd tau_f3;
         Eigen::VectorXd tau_f4;
         
         Eigen::Vector3d com0;
         
         ForceOptimization::Ptr _force_opt;

        void sense(const double period);

        void QPPVMControl(const double time);
        
        estimation::FloatingBaseEstimator::Ptr _fb_estimator;
        XBot::ImuSensor::ConstPtr _imu;
        std::vector<std::string> _contact_links;

    };

}
#endif




