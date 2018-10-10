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
    WBTCController(XBot::ModelInterface& model);

    void log(XBot::MatLogger::Ptr logger);

private:
    XBot::ModelInterface& _model;

};
}

namespace xbotcore {
class WBTCPlugin : public XBot::XBotControlPlugin {

public:
    WBTCPlugin();

    virtual bool init_control_plugin(XBot::Handle::Ptr handle);
    virtual void on_start(double time);
    virtual void control_loop(double time, double period);
    virtual bool close();

    void log();

    bool control();

private:
    XBot::MatLogger::Ptr _matlogger;
    XBot::RobotInterface::Ptr _robot;

    // To store impedance values in the dsps
    Eigen::VectorXd _k_dsp,_k_dsp_ref;
    Eigen::VectorXd _d_dsp, _d_dsp_ref;
    Eigen::VectorXd _h, _tau;

    //// DNAMIC RECONFIGURE TODO: PUT IT IN ANOTHER CLASS
    dynamic_reconfigure_advr::Server<QPPVM_RT_plugin::QppvmConfig> _server;
    void cfg_callback(QPPVM_RT_plugin::QppvmConfig &config, uint32_t level);
    std::atomic<double> _impedance_gain; //Impedance GAINS in the DSPs

};
}

#endif
