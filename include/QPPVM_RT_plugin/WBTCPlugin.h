#ifndef __WBTC_PLUGIN_H__
#define __WBTC_PLUGIN_H__

#include <XCM/XBotControlPlugin.h>
#include <OpenSoT/solvers/iHQP.h>
#include <OpenSoT/tasks/torque/JointImpedanceCtrl.h>
#include <OpenSoT/constraints/torque/TorqueLimits.h>
#include <OpenSoT/utils/AutoStack.h>

#include <XBotInterface/Logger.hpp>

#include <dynamic_reconfigure_advr/server.h>
#include <QPPVM_RT_plugin/QppvmConfig.h>

#include<atomic>

namespace opensot{
class WBTCController{
public:
    typedef boost::shared_ptr<WBTCController> Ptr;

    WBTCController(XBot::ModelInterface& model);

    bool control(Eigen::VectorXd& tau);

    void log(XBot::MatLogger::Ptr logger);

    OpenSoT::tasks::torque::JointImpedanceCtrl::Ptr joint_impedance;
    OpenSoT::constraints::torque::TorqueLimits::Ptr torque_lims;


private:
    XBot::ModelInterface& _model;
    Eigen::VectorXd _h;
    Eigen::VectorXd _q;
    Eigen::VectorXd _tau_opt;


    OpenSoT::AutoStack::Ptr _autostack;

    OpenSoT::solvers::iHQP::Ptr _solver;

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

};

class WBTCPlugin : public XBot::XBotControlPlugin {

public:
    WBTCPlugin();

    virtual bool init_control_plugin(XBot::Handle::Ptr handle);
    virtual void on_start(double time);
    virtual void control_loop(double time, double period);
    virtual bool close();

    void log();

    opensot::WBTCController::Ptr controller;

private:
    XBot::MatLogger::Ptr _matlogger;
    XBot::RobotInterface::Ptr _robot;

    // To store impedance values in the dsps
    Eigen::VectorXd _k_dsp,_k_dsp_ref;
    Eigen::MatrixXd _Kj, _Dj, _Kj_ref, _Dj_ref;
    Eigen::VectorXd _d_dsp, _d_dsp_ref;
    Eigen::VectorXd _tau, _tau_ref, _tau_offset;

    dynamic_reconf _dynamic_reconfigure;


};

}

#endif
