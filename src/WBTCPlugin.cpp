#include <QPPVM_RT_plugin/WBTCPlugin.h>

REGISTER_XBOT_PLUGIN(WBTCPlugin, xbotcore::WBTCPlugin)

namespace opensot{

WBTCController::WBTCController(XBot::ModelInterface &model):
    _model(model)
{

}

void WBTCController::log(XBot::MatLogger::Ptr logger)
{

}

}

namespace xbotcore{

WBTCPlugin::WBTCPlugin()
{

}

bool WBTCPlugin::init_control_plugin(XBot::Handle::Ptr handle)
{
    dynamic_reconfigure_advr::Server<QPPVM_RT_plugin::QppvmConfig>::CallbackType f;
    f = boost::bind(&WBTCPlugin::cfg_callback, this, _1, _2);
    _server.setCallback(f);
    _impedance_gain.store(1.0);


    _matlogger = XBot::MatLogger::getLogger("/tmp/wbtc_log");

    _robot = handle->getRobotInterface();
    _robot->sense(true);

    _k_dsp.setZero(_robot->getJointNum());
    _d_dsp.setZero(_robot->getJointNum());
    _tau_ref.setZero(_robot->getJointNum());
    _tau.setZero(_robot->model().getJointNum());
    _tau_offset.setZero(_robot->getJointNum());
    _h.setZero(_robot->model().getJointNum());
    
    
    _robot->getRobotState("torque_offset", _tau_offset);

    _robot->getStiffness(_k_dsp);
    _robot->getDamping(_d_dsp);
    _k_dsp_ref = _k_dsp;
    _d_dsp_ref = _d_dsp;

    log();

    return true;
}

void WBTCPlugin::on_start(double time)
{
    _robot->sense(true);

    log();
}

void WBTCPlugin::control_loop(double time, double period)
{
    _robot->sense(true);

    if(control())
    {
        _k_dsp_ref = _impedance_gain.load() * _k_dsp;
        _d_dsp_ref = std::sqrt(2.0) * _impedance_gain.load() * _d_dsp;

        _robot->setStiffness(_k_dsp_ref);
        _robot->setDamping(_d_dsp_ref);

        _tau_ref = _tau.tail(_robot->model().getActuatedJointNum()) - _tau_offset;
        _robot->setEffortReference(_tau_ref);

        _robot->move();


    }

    log();
}

bool WBTCPlugin::control()
{
    _tau.setZero(_tau.size());

    _robot->model().computeNonlinearTerm(_h);

    _tau.noalias() = _tau +  _h;
    
    return true;
}

void WBTCPlugin::log()
{
    _robot->log(_matlogger, XBot::get_time_ns(),"wbtc");
    _matlogger->add("_h", _h);
    _matlogger->add("_tau_offset", _tau_offset);
}

bool WBTCPlugin::close()
{
    _matlogger->flush();
    return true;
}

void WBTCPlugin::cfg_callback(QPPVM_RT_plugin::QppvmConfig& config, uint32_t level)
{
    _impedance_gain.store(config.impedance_gain);
//    _stiffness_Waist_gain.store(config.stiffness_waist);
//    _damping_Waist_gain.store(config.damping_waist);
//    _stiffness_Feet_gain.store(config.stiffness_feet);
//    _damping_Feet_gain.store(config.damping_feet);
//    _joints_gain.store(config.joints_gain);

    Logger::info(Logger::Severity::HIGH, "\nSetting impedance gain to %f \n", config.impedance_gain);
//    Logger::info(Logger::Severity::HIGH, "Setting joints gain to %f \n", config.joints_gain);
//    Logger::info(Logger::Severity::HIGH, "Setting waist stiffness gain to %f and damping gain to %f\n", config.stiffness_waist, config.damping_waist);
//    Logger::info(Logger::Severity::HIGH, "Setting feet stiffness gain to %f and damping gain to %f\n", config.stiffness_feet, config.damping_feet);
}

}
