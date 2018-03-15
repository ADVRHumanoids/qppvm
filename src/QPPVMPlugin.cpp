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
    _matlogger = XBot::MatLogger::getLogger("/tmp/qppvm_log");
    
    _set_ref = false;

    _robot = handle->getRobotInterface();
    _model = XBot::ModelInterface::getModel("external/qppvm/config/centauro_simple_example.yaml");
    
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
    
    _q = _q_home;
    _q_ref = _q;

    _k.setZero(_robot->getJointNum());
    _d.setZero(_robot->getJointNum());

    Eigen::VectorXd k0, d0;
    _robot->getStiffness(k0);
    _robot->getDamping(d0);

    Eigen::VectorXd _k_matrix;
    Eigen::VectorXd _d_matrix;
    
    _k_matrix = Eigen::VectorXd(_k_matrix.size()).setConstant(5.0);
    _d_matrix = Eigen::VectorXd(_d_matrix.size()).setConstant(2.0);

    
    std::cout<<"_d_matrix: \n"<<_d_matrix.transpose()<<std::endl;
    std::cout<<"_k_matrix: \n"<<_k_matrix.transpose()<<std::endl;

    /* Create torque limit constraint */
    _torque_limits = boost::make_shared<TorqueLimits>(_tau_max, _tau_min);
    
    /* Create torque joint impedance task */
    _joint_task = boost::make_shared<JointImpedanceCtrl>(_q_home, *_model);
    
    _joint_task->setStiffness(_k_matrix.asDiagonal());
    _joint_task->setDamping(_d_matrix.asDiagonal());
    _joint_task->useInertiaMatrix(true);
    _joint_task->update(_q);

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
        
        _Kc.setIdentity(6,6); _Kc = 700.*_Kc; 
        _Dc.setIdentity(6,6); _Dc = 70.0*_Dc; 
        
        imp_task->setStiffnessDamping(_Kc, _Dc);
        imp_task->useInertiaMatrix(true);
      
    }
    
    auto legs_impedance_aggr = _leg_impedance_task[0] + _leg_impedance_task[1] + _leg_impedance_task[2] + _leg_impedance_task[3];
    
    _autostack =  ( legs_impedance_aggr /_joint_task ) << _torque_limits; 

    _solver = boost::make_shared<QPOases_sot>(_autostack->getStack(), _autostack->getBounds(), 1.0);
    _solver->log(_matlogger);

    return true;
}

void QPPVMPlugin::QPPVMControl(const double time)
{
    _tau_max = _tau_max_const - _h;
    _tau_min = _tau_min_const - _h;
    _torque_limits->setTorqueLimits(_tau_max, _tau_min);

    _autostack->update(_q);
    
    if(!_solver->solve(_tau_d)){
        _tau_d.setZero(_tau_d.size());
        XBot::Logger::error("Unable to solve \n");
    }
    
    _solver->log(_matlogger);
   
    _matlogger->add("tau_qp", _tau_d);
    
    _tau_d = _tau_d + _h;
    
     _matlogger->add("tau_desired", _tau_d);
}

void QPPVMPlugin::on_start(double time)
{
    sense();
    
    _start_time = time;
    _robot->move();
    
    for(auto task : _leg_impedance_task)
    {
        Eigen::Affine3d task_pose;
        _model->getPose(task->getDistalLink(), task->getBaseLink(), task_pose);
        
        task->setReference(task_pose.matrix());
    }
    
    _joint_task->setReference(_q);
    
}       


void QPPVMPlugin::control_loop(double time, double period)
{
    
    sense();

    QPPVMControl(time);
    
    _model->setJointEffort(_tau_d);

    _robot->setReferenceFrom(*_model, XBot::Sync::Effort, XBot::Sync::Impedance);

    _matlogger->add("time_matlogger", time);
    
    _model->log(_matlogger, time);

    _robot->move();
}

void QPPVMPlugin::sense()
{
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
                      
    if( (_Dc.array() <= 0).any() )
    {
        XBot::Logger::warning("Invalid damping matrix received \n");
        return;
    }
    
    XBot::Logger::info(XBot::Logger::Severity::HIGH) << "New damping: " << _Dc.diagonal().transpose() << XBot::Logger::endl();
                  
    _leg_impedance_task[id]->setDamping(_Dc);
}

