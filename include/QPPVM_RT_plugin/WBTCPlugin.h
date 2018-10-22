#ifndef __WBTC_PLUGIN_H__
#define __WBTC_PLUGIN_H__

#include <XCM/XBotControlPlugin.h>
#include <OpenSoT/solvers/iHQP.h>
#include <OpenSoT/tasks/torque/JointImpedanceCtrl.h>
#include <OpenSoT/tasks/torque/CartesianImpedanceCtrl.h>
#include <OpenSoT/constraints/torque/TorqueLimits.h>
#include <OpenSoT/utils/AutoStack.h>

#include <XBotInterface/Logger.hpp>

#include <dynamic_reconfigure_advr/server.h>
#include <QPPVM_RT_plugin/QppvmConfig.h>

#include <cartesian_interface/CartesianPlugin/Utils.h>

#include<atomic>

namespace opensot{
class WBTCController{
public:
    typedef boost::shared_ptr<WBTCController> Ptr;

    WBTCController(XBot::ModelInterface& model);

    bool control(Eigen::VectorXd& tau);
    void setFilter(const double period, const double cut_off_freq);

    void log(XBot::MatLogger::Ptr logger);

    OpenSoT::tasks::torque::JointImpedanceCtrl::Ptr joint_impedance;
    OpenSoT::tasks::torque::CartesianImpedanceCtrl::Ptr LFoot, RFoot, Waist;
    OpenSoT::constraints::torque::TorqueLimits::Ptr torque_lims;


private:
    XBot::ModelInterface& _model;
    Eigen::VectorXd _h;
    Eigen::VectorXd _q, _qdot;
    Eigen::VectorXd _tau_opt;
    
    Eigen::VectorXd _qdot_filtered;
    double _alpha_filter;
    double _period;
    double _cut_off_freq;


    OpenSoT::AutoStack::Ptr _autostack;

    OpenSoT::solvers::iHQP::Ptr _solver;
    
    
    Eigen::MatrixXd J;
    Eigen::VectorXd spring, damper, force, Jtf;

};
}

namespace xbotcore {
class dynamic_reconf{
public:
    dynamic_reconf();

    dynamic_reconfigure_advr::Server<QPPVM_RT_plugin::QppvmConfig> _server;
    void cfg_callback(QPPVM_RT_plugin::QppvmConfig &config, uint32_t level);
    std::atomic<double> _impedance_gain; //Impedance GAINS in the DSPs
    std::atomic<double> _joints_gain; //Impedance GAINS in the joint impedance task
    std::atomic<double> _stiffness_Feet_gain, _damping_Feet_gain; //Impedance GAINS for the feet
    std::atomic<double> _stiffness_Waist_gain, _damping_Waist_gain;

};

class WBTCPlugin : public XBot::XBotControlPlugin {

public:
    WBTCPlugin();

    virtual bool init_control_plugin(XBot::Handle::Ptr handle);
    virtual void on_start(double time);
    virtual void control_loop(double time, double period);
    virtual bool close();

    void log();
    void set_dyn_reconfigure_gains();
    void sync_cartesian_ifc(double time, double period);

    opensot::WBTCController::Ptr controller;

private:
    XBot::MatLogger::Ptr _matlogger;
    XBot::RobotInterface::Ptr _robot;

    // To store impedance values in the dsps
    Eigen::VectorXd _k_dsp,_k_dsp_ref;
    Eigen::MatrixXd _Kj, _Dj, _Kj_ref, _Dj_ref;
    Eigen::VectorXd _Kj_vec, _Dj_vec;
    Eigen::VectorXd _d_dsp, _d_dsp_ref;
    Eigen::VectorXd _tau, _tau_ref, _tau_offset;

    Eigen::MatrixXd _K_Waist, _D_Waist;
    Eigen::MatrixXd _K_Foot, _D_Foot;
    Eigen::MatrixXd _K_Lfoot_ref, _D_Lfoot_ref, _K_Rfoot_ref, _D_Rfoot_ref, _K_Waist_ref, _D_Waist_ref;

    dynamic_reconf _dynamic_reconfigure;
    
    bool _use_offsets;
    
    /* Cartesian ifc variables */
    XBot::Cartesian::CartesianInterfaceImpl::Ptr _ci;
    XBot::Cartesian::Utils::SyncFromIO::Ptr _sync_from_nrt;
    bool _first_sync_done;
    XBot::ModelInterface::Ptr _model;


};

}

#endif
