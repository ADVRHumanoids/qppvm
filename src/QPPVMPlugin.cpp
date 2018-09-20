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
#include <QPPVM_RT_plugin/ForceOptimization.h>
#include <OpenSoT/constraints/force/WrenchLimits.h>
#include <geometry_msgs/Pose.h>
#include <eigen_conversions/eigen_msg.h>

#define TRJ_TIME 3.0


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

bool QPPVMPlugin::init_control_plugin(  XBot::Handle::Ptr handle)
{
    
    
    dynamic_reconfigure_advr::Server<QPPVM_RT_plugin::QppvmConfig>::CallbackType f;

    f = boost::bind(&QPPVMPlugin::cfg_callback, this, _1, _2);
    _server.setCallback(f);

//     _sh_fb_pos = handle->getSharedMemory()->getSharedObject<Eigen::Vector3d>("/gazebo/floating_base_position");
//     _sh_fb_rot = handle->getSharedMemory()->getSharedObject<Eigen::Quaterniond>("/gazebo/floating_base_orientation");
//     _sh_fb_vel = handle->getSharedMemory()->getSharedObject<Eigen::Vector6d>("/gazebo/floating_base_velocity");
//     
//     _sh_fb_pos.set(Eigen::Vector3d::Zero());
//     _sh_fb_rot.set(Eigen::Quaterniond::Identity());
//     _sh_fb_vel.set(Eigen::Vector6d::Zero());

    _impedance_gain_sub = handle->getRosHandle()->subscribe("/qppvm/impedance_gain",
                                                     1,
                                                     &QPPVMPlugin::impedance_gain_callback, this);
    
    _feedback_gain_sub = handle->getRosHandle()->subscribe("/qppvm/feedback_gain",
                                                     1,
                                                     &QPPVMPlugin::feedback_gain_callback, this);
    
    _fb_pub = handle->getRosHandle()->advertise<geometry_msgs::Pose>("/qppvm/floating_base", 1);



    _matlogger = XBot::MatLogger::getLogger("/tmp/qppvm_log");

    _set_ref = false;

    _robot = handle->getRobotInterface();
    _imu = _robot->getImu().begin()->second;
//    _model = XBot::ModelInterface::getModel("external/qppvm/config/floating_base/centauro_simple_example.yaml");
//     _model = XBot::ModelInterface::getModel("configs/ADVR_shared/user_example/centauro_simple_example.yaml");
    _model = XBot::ModelInterface::getModel(handle->getPathToConfigFile());
    
    
    
    
    
    _contact_links = {"l_sole", "r_sole"};
    Eigen::MatrixXd contact_matrix(6,6); contact_matrix.setIdentity(6,6);
    _fb_estimator = std::make_shared<estimation::FloatingBaseEstimator>(_model, _imu, _contact_links, contact_matrix);
    

    for(int i = 0; i < _robot->legs(); i++)
    {

        auto callback_k = boost::bind(&QPPVMPlugin::cart_stiffness_callback,
                                     this, _1, i);

        auto pub_k = handle->getRosHandle()->subscribe<geometry_msgs::Twist>("/qppvm/leg_"+std::to_string(i+1)+"/stiffness",
                                                                           1,
                                                                           callback_k);


        auto callback_d = boost::bind(&QPPVMPlugin::cart_damping_callback,
                                     this, _1, i);

        auto pub_d = handle->getRosHandle()->subscribe<geometry_msgs::Twist>("/qppvm/leg_"+std::to_string(i+1)+"/damping",
                                                                           1,
                                                                           callback_d);

        _cart_stiffness_sub.push_back(pub_k);
        _cart_damping_sub.push_back(pub_d);
    }

    _model->initLog(_matlogger, 30000);

    _model->getEffortLimits(_tau_max_const);
    _tau_min_const = -_tau_max_const;

    for(unsigned int i = 0; i < 6; ++i)
    {
        _tau_max_const[i] = 100000.0;
        _tau_min_const[i] = -_tau_max_const[i];
    }

    _tau_d.setZero(_model->getJointNum());

    _model->computeNonlinearTerm(_h);

    _tau_max = _tau_max_const - _h;
    _tau_min = _tau_min_const - _h;

    _model->getRobotState("home", _q_home);
    _model->setJointPosition(_q_home);
    _model->setJointVelocity(Eigen::VectorXd::Zero(_model->getJointNum()));
    _model->setStiffness(Eigen::VectorXd::Zero(_model->getJointNum()));
    _model->setDamping(Eigen::VectorXd::Zero(_model->getJointNum()));
    _model->update();


//    _model->getRobotState("torque_offset", _tau_offset);

//    std::cout << _tau_offset << std::endl;

    _q = _q_home;
    _q_ref = _q;

    _k.setZero(_robot->getJointNum());
    _d.setZero(_robot->getJointNum());
    _robot->getStiffness(_k);
    _robot->getDamping(_d);

    Eigen::VectorXd _k_matrix(_model->getJointNum());
    Eigen::VectorXd _d_matrix(_model->getJointNum());
    _k_matrix.setZero(_k_matrix.size());
    _d_matrix.setZero(_d_matrix.size());
    _k_matrix.segment(6,_robot->getJointNum()) = 0.01*_k;//0.1*_k;
    _d_matrix.segment(6,_robot->getJointNum()) = 0.1*_d;//0.1*_d;


    std::cout<<"_d_matrix: \n"<<_d_matrix.transpose()<<std::endl;
    std::cout<<"_k_matrix: \n"<<_k_matrix.transpose()<<std::endl;

    /* Create torque limit constraint */
    _torque_limits = boost::make_shared<TorqueLimits>(_tau_max, _tau_min);
    std::cout<<"tau_max: "<<_tau_max.transpose()<<std::endl;
    std::cout<<"tau_min: "<<_tau_min.transpose()<<std::endl;


    /* Create torque joint impedance task */
    _joint_task = boost::make_shared<JointImpedanceCtrl>(_q_home, *_model);

    _joint_task->setStiffness(_k_matrix.asDiagonal());
    _joint_task->setDamping(_d_matrix.asDiagonal());
    _joint_task->useInertiaMatrix(true);
    _joint_task->update(_q);
    
    std::vector<std::string> links_in_contact;

    /* Create leg impedance tasks */
    for(int i = 0; i < _robot->legs(); i++)
    {
        auto imp_task = boost::make_shared<CartesianImpedanceTask>("LEG_"+std::to_string(i+1)+"_CART_IMP",
                                                                    _q,
                                                                   *_model,
                                                                   _robot->leg(i).getTipLinkName(),
                                                                   "world",
                                                                   OpenSoT::Indices::range(0,2)
                                                                   );

        _leg_impedance_task.push_back(imp_task);
        
        links_in_contact.push_back(_robot->leg(i).getTipLinkName());

        _Kc.setIdentity(6,6); _Kc = 500.*_Kc;
        _Dc.setIdentity(6,6); _Dc = 50.*_Dc;
//         _Kc.setIdentity(6,6); _Kc = 0.*_Kc;
//         _Dc.setIdentity(6,6); _Dc = 0.*_Dc;

        imp_task->setStiffnessDamping(_Kc, _Dc);
        imp_task->useInertiaMatrix(true);

    }

    _waist = boost::make_shared<CartesianImpedanceTask>("WAIST_CART_IMP",
                                                            _q,
                                                           *_model,
                                                           "Waist",
                                                           "world");
    Eigen::MatrixXd _Kw(6,6); _Kw.setIdentity(6,6);
    Eigen::MatrixXd _Dw(6,6); _Dw.setIdentity(6,6);

    double k_waist = 500.;
    double d_waist = 50.;
    
    _Kw(0,0) = k_waist;
    _Kw(1,1) = k_waist;
    _Kw(2,2) = k_waist;
    _Kw(3,3) = 0.5*k_waist;
    _Kw(4,4) = 0.5*k_waist;
    _Kw(5,5) = 0.5*k_waist;

    _Dw(0,0) = d_waist;
    _Dw(1,1) = d_waist;
    _Dw(2,2) = d_waist;
    _Dw(3,3) = d_waist;
    _Dw(4,4) = d_waist;
    _Dw(5,5) = d_waist;

    _waist->setStiffnessDamping(_Kw, _Dw);
      
    
    _ee_task_left = boost::make_shared<CartesianImpedanceTask>("LEFT_ARM",
                                                               _q,
                                                              *_model,
                                                              _robot->chain("left_arm").getTipLinkName(),
                                                              "world",
                                                               OpenSoT::Indices::range(0,2)
                                                              );
    
    _Kc.setIdentity(6,6); _Kc = 500.*_Kc;
    _Dc.setIdentity(6,6); _Dc = 50.*_Dc;
    
    _ee_task_left->setStiffnessDamping(_Kc, _Dc);
    _ee_task_left->useInertiaMatrix(true);
    
    
    _ee_task_right = boost::make_shared<CartesianImpedanceTask>("RIGHT_ARM",
                                                                _q,
                                                                *_model,
                                                                _robot->chain("right_arm").getTipLinkName(),
                                                                "world",
                                                                 OpenSoT::Indices::range(0,2)
                                                                );
    
    _Kc.setIdentity(6,6); _Kc = 500.*_Kc;
    _Dc.setIdentity(6,6); _Dc = 50.*_Dc;
    
    _ee_task_right->setStiffnessDamping(_Kc, _Dc);
    _ee_task_right->useInertiaMatrix(true);
    


    auto legs_impedance_aggr = _leg_impedance_task[0] + _leg_impedance_task[1];// + _leg_impedance_task[2] + _leg_impedance_task[3];
    auto ee_impedance_aggr = _ee_task_left +  _ee_task_right;

//     _autostack =  ( (legs_impedance_aggr + ee_impedance_aggr) / ( _joint_task) ) << _torque_limits;
    
//      _autostack =  ( legs_impedance_aggr / _waist /ee_impedance_aggr /  _joint_task ) << _torque_limits;
    
    _autostack =  ( legs_impedance_aggr / _waist /  _joint_task ) << _torque_limits;
        
//    _autostack =  ( legs_impedance_aggr /  _joint_task ) << _torque_limits;

//     _autostack.reset(new  OpenSoT::AutoStack(_joint_task));
//     _autostack<<_torque_limits;


    _solver = boost::make_shared<iHQP>(_autostack->getStack(), _autostack->getBounds(), 1.0);
    _solver->log(_matlogger);


    _force_opt = boost::make_shared<ForceOptimization>(_model, links_in_contact);
    
    _impedance_gain.store(1.0);
    _stiffness_Waist_gain.store(0.0);
    _damping_Waist_gain.store(0.0); 
    _stiffness_Feet_gain.store(0.0);
    _damping_Feet_gain.store(0.0);
    _joints_gain.store(0.0);


    /* Init cartesian ifc */
    YAML::Node yaml_file = YAML::LoadFile(handle->getPathToConfigFile());
    XBot::Cartesian::ProblemDescription ik_problem(yaml_file["CartesianInterface"]["problem_description"], _model);
    _ci = std::make_shared<XBot::Cartesian::CartesianInterfaceImpl>(_model, ik_problem);
    _sync_from_nrt = std::make_shared<XBot::Cartesian::Utils::SyncFromIO>("/xbotcore/cartesian_interface", handle->getSharedMemory());

    return true;
}

void QPPVMPlugin::QPPVMControl(const double time)
{
    
    set_gains();

    _tau_max = _tau_max_const - _h;
    _tau_min = _tau_min_const - _h;
    _torque_limits->setTorqueLimits(_tau_max, _tau_min);


    _autostack->update(_q);
    _autostack->log(_matlogger);

    if(!_solver->solve(_tau_d)){
        _tau_d.setZero(_tau_d.size());
        XBot::Logger::error("Unable to solve \n");
    }

    _solver->log(_matlogger);

    _matlogger->add("tau_qp", _tau_d);

    _matlogger->add("h", _h);

    _tau_d.noalias() = _tau_d + _h;

    _matlogger->add("tau_desired", _tau_d);
    

}

void demo::QPPVMPlugin::set_gains()
{
    double current_stiffness_Waist_gain = _stiffness_Waist_gain.load();
    double current_damping_Waist_gain = _damping_Waist_gain.load();
    double current_stiffness_Feet_gain = _stiffness_Feet_gain.load();
    double current_damping_Feet_gain = _damping_Feet_gain.load();
    double current_joint_gain = _joints_gain.load();
    _Kc.setIdentity(6,6); _Dc.setIdentity(6,6);
    _Dc *= current_damping_Waist_gain;
    _Kc *= current_stiffness_Waist_gain;
    
    _waist->setStiffnessDamping(_Kc, _Dc);
    
    _Kc.setIdentity(6,6); _Dc.setIdentity(6,6);
    _Dc *= current_damping_Feet_gain;
    _Kc *= current_stiffness_Feet_gain;
    
    for(auto imptask : _leg_impedance_task)
    {
        imptask->setStiffnessDamping(_Kc, _Dc);
    }
    
    _Kj.setIdentity(_model->getJointNum(),_model->getJointNum()); _Dj.setIdentity(_model->getJointNum(),_model->getJointNum());
    _Kj = _Kj * 100. * current_joint_gain;
    _Dj = _Dj * 10. * current_joint_gain;
    _joint_task->setStiffnessDamping(_Kj, _Dj);
}


void QPPVMPlugin::on_start(double time)
{
    _first_sync_done = false;
    
    _robot->setControlMode("torso", XBot::ControlMode::Idle());
    _robot->setControlMode("left_arm", XBot::ControlMode::Idle());
    _robot->setControlMode("right_arm", XBot::ControlMode::Idle());
    _robot->setControlMode("neck", XBot::ControlMode::Idle());
    
    
    
    sense(0); 

    _start_time = time;
    

     for(auto task : _leg_impedance_task)
     {
         Eigen::Affine3d task_pose;
         _model->getPose(task->getDistalLink(), task->getBaseLink(), task_pose);

         task->setReference(task_pose.matrix());
     }


     Eigen::Affine3d left_ee_pose;
     _model->getPose(_ee_task_left->getDistalLink(), _ee_task_left->getBaseLink(), left_ee_pose);

     Eigen::Affine3d right_ee_pose;
     _model->getPose(_ee_task_right->getDistalLink(),_ee_task_right->getBaseLink(), right_ee_pose);

     _ee_task_left->setReference(left_ee_pose.matrix());
     _ee_task_right->setReference(right_ee_pose.matrix());


     Eigen::Affine3d task_pose;
     _model->getPose("Waist", task_pose);
     _waist->setReference(task_pose.matrix());


     _joint_task->setReference(_q);

     _h.segment(0,6) = Eigen::Vector6d::Zero(); ///REMOVE LATER!
     _tau_max = _tau_max_const - _h;
     _tau_min = _tau_min_const - _h;
     _torque_limits->setTorqueLimits(_tau_max, _tau_min);


     _autostack->update(_q);
    

}


void QPPVMPlugin::control_loop(double time, double period)
{

    sync_cartesian_ifc(time, period);
    
    sense(period);
    
    QPPVMControl(time);
    
    _force_opt->compute(_tau_d, _Fopt, _tau_opt);

    _matlogger->add("tau_opt", _tau_opt);
    
    for(int i = 0; i < _robot->legs(); i++)
    {
        _matlogger->add("wrench_opt_" + std::to_string(i+1), _Fopt[i]);
    }
    
    _force_opt->log(_matlogger);
    _model->setJointEffort(_tau_opt);
    
    _k_ref = _impedance_gain.load() * _k;
    _d_ref = std::sqrt(2.0) * _impedance_gain.load() * _d;

    _robot->setReferenceFrom(*_model, XBot::Sync::Effort); // , XBot::Sync::Impedance);
    
    _robot->setStiffness(_k_ref);
    _robot->setDamping(_d_ref);


    _matlogger->add("time_matlogger", time);

    _model->log(_matlogger, time);

    _robot->move();
}

void QPPVMPlugin::sense(const double period)
{
//      Eigen::Affine3d w_T_fb;
//      Eigen::Vector6d fb_twist;
// 
//      w_T_fb.translation() = _sh_fb_pos.get();
//      w_T_fb.linear() = _sh_fb_rot.get().toRotationMatrix();
//      fb_twist = _sh_fb_vel.get();


//      _model->setFloatingBaseState(w_T_fb, fb_twist);
//      _model->update();
    

    _model->syncFrom(*_robot, XBot::Sync::Position, XBot::Sync::Velocity, XBot::Sync::MotorSide);
    
    _fb_estimator->update(period);
    _fb_estimator->log(_matlogger);
    
    _model->getJointPosition(_q);
    _model->getJointVelocity(_dq);
    _model->computeNonlinearTerm(_h);
    
    geometry_msgs::Pose pose_msg;
    Eigen::Affine3d T;
    _model->getFloatingBasePose(T);
    tf::poseEigenToMsg(T, pose_msg);
    _fb_pub->pushToQueue(pose_msg);
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

void demo::QPPVMPlugin::cart_stiffness_callback(const geometry_msgs::TwistConstPtr& msg, int id)
{
    _Kc.diagonal() << msg->linear.x, msg->linear.y, msg->linear.z,
                      msg->angular.x, msg->angular.y, msg->angular.z;

    if( (_Kc.array() < 0).any() )
    {
        XBot::Logger::warning("Invalid stiffness matrix received \n");
        return;
    }

    XBot::Logger::info(XBot::Logger::Severity::HIGH) << "New stiffness: " << _Kc.diagonal().transpose() << XBot::Logger::endl();

    _leg_impedance_task[id]->setStiffness(_Kc);
}

void demo::QPPVMPlugin::cart_damping_callback(const geometry_msgs::TwistConstPtr& msg, int id)
{
    _Dc.diagonal() << msg->linear.x, msg->linear.y, msg->linear.z,
                      msg->angular.x, msg->angular.y, msg->angular.z;

    if( (_Dc.diagonal().array() <= 0).any() )
    {
        XBot::Logger::warning("Invalid damping matrix received \n");
        return;
    }

    XBot::Logger::info(XBot::Logger::Severity::HIGH) << "New damping: " << _Dc.diagonal().transpose() << XBot::Logger::endl();

    _leg_impedance_task[id]->setDamping(_Dc);
}

void demo::QPPVMPlugin::sync_cartesian_ifc(double time, double period)
{
    /* Sync cartesian references from ROS */
    
    if(!_first_sync_done)
    {
        if(_sync_from_nrt->try_reset(_model, time))
        {
            _first_sync_done = true;
            XBot::Logger::info(Logger::Severity::HIGH, "Resetting NRT CI \n");
        }
    }
    
    if(_first_sync_done)
    {
        if(_sync_from_nrt->try_sync(time, period, _ci, _model))
        {
            _matlogger->add("sync_done", 1);
        }
        else
        {
            _matlogger->add("sync_done", 0);
        }
    }
    

    if(!_ci->update(time, period))
    {
        XBot::Logger::error("CartesianInterface: unable to solve \n");
        return;
    }
    
    
    /* Update qppvm references */
    Eigen::Affine3d T_ref;
    _ci->getPoseReference(_waist->getDistalLink(), T_ref);
    _waist->setReference(T_ref.matrix());
}

void demo::QPPVMPlugin::impedance_gain_callback(const std_msgs::Float64ConstPtr& msg)
{

    double impedance_gain = msg->data;
    impedance_gain = std::min(impedance_gain, 1.0);
    impedance_gain = std::max(impedance_gain, 0.0);

    Logger::info(Logger::Severity::HIGH, "Setting impedance gain to %f \n", impedance_gain);

    _impedance_gain.store(impedance_gain);
}

void demo::QPPVMPlugin::feedback_gain_callback(const std_msgs::Float64ConstPtr& msg)
{

    double feedback_gain = msg->data;
    feedback_gain = std::min(feedback_gain, 1.0);
    feedback_gain = std::max(feedback_gain, 0.0);

    Logger::info(Logger::Severity::HIGH, "Setting feedback gain to %f \n", feedback_gain);

    _stiffness_Waist_gain.store(feedback_gain);
    _damping_Waist_gain.store(std::sqrt(feedback_gain));
}


void demo::QPPVMPlugin::cfg_callback(QPPVM_RT_plugin::QppvmConfig& config, uint32_t level)
{
    _impedance_gain.store(config.impedance_gain);
    _stiffness_Waist_gain.store(config.stiffness_waist);
    _damping_Waist_gain.store(config.damping_waist);
    _stiffness_Feet_gain.store(config.stiffness_feet);
    _damping_Feet_gain.store(config.damping_feet);
    _joints_gain.store(config.joints_gain);
    
    Logger::info(Logger::Severity::HIGH, "\nSetting impedance gain to %f \n", config.impedance_gain);
    Logger::info(Logger::Severity::HIGH, "Setting joints gain to %f \n", config.joints_gain);
    Logger::info(Logger::Severity::HIGH, "Setting waist stiffness gain to %f and damping gain to %f\n", config.stiffness_waist, config.damping_waist);
    Logger::info(Logger::Severity::HIGH, "Setting feet stiffness gain to %f and damping gain to %f\n", config.stiffness_feet, config.damping_feet);
}
