#ifndef __WBTC_PLUGIN_H__
#define __WBTC_PLUGIN_H__

#include <XCM/XBotControlPlugin.h>
#include <OpenSoT/solvers/iHQP.h>
#include <OpenSoT/tasks/torque/JointImpedanceCtrl.h>
#include <OpenSoT/constraints/torque/TorqueLimits.h>
#include <OpenSoT/utils/AutoStack.h>

#include <XBotInterface/Logger.hpp>

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

private:
    XBot::MatLogger::Ptr _matlogger;
    XBot::RobotInterface::Ptr _robot;

    Eigen::VectorXd _q, _qdot;

};
}

#endif
