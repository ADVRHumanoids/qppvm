#include <QPPVM_RT_plugin/WBTCPlugin.h>

REGISTER_XBOT_PLUGIN(WBTCPlugin, xbotcore::WBTCPlugin)

namespace opensot{

WBTCController::WBTCController(XBot::ModelInterface &model):
    _model(model)
{
    _h.setZero(_model.getJointNum());
    _tau_opt.setZero(_model.getJointNum());

    _model.getJointPosition(_q);

    joint_impedance = boost::make_shared<OpenSoT::tasks::torque::JointImpedanceCtrl>(_q, _model);

    Eigen::VectorXd tau_lims;
    _model.getEffortLimits(tau_lims);
    torque_lims = boost::make_shared<OpenSoT::constraints::torque::TorqueLimits>(tau_lims, -tau_lims);

    _autostack = boost::make_shared<OpenSoT::AutoStack>(joint_impedance);
    _autostack<<torque_lims;

    _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 1.);

}

bool WBTCController::control(Eigen::VectorXd& tau)
{
    _model.computeNonlinearTerm(_h);
    _model.getJointPosition(_q);

    _autostack->update(_q);

    if(!_solver->solve(_tau_opt))
        _tau_opt.setZero(_tau_opt.size());

    tau = _tau_opt +  _h;

    return true;
}


void WBTCController::log(XBot::MatLogger::Ptr logger)
{
    logger->add("_h", _h);
    logger->add("_tau_opt", _tau_opt);
    _autostack->log(logger);
    _solver->log(logger);
}

}

namespace xbotcore{

WBTCPlugin::WBTCPlugin()
{

}

bool WBTCPlugin::init_control_plugin(XBot::Handle::Ptr handle)
{
    _matlogger = XBot::MatLogger::getLogger("/tmp/wbtc_log");

    _robot = handle->getRobotInterface();
    _robot->sense(true);

    _k_dsp.setZero(_robot->getJointNum());
    _d_dsp.setZero(_robot->getJointNum());
    _tau_ref.setZero(_robot->getJointNum());
    _tau.setZero(_robot->model().getJointNum());
    _tau_offset.setZero(_robot->getJointNum());    
    
    _robot->getRobotState("torque_offset", _tau_offset);

    _robot->getStiffness(_k_dsp);
    _robot->getDamping(_d_dsp);
    _k_dsp_ref = _k_dsp;
    _d_dsp_ref = _d_dsp;

    _Kj.setOnes(_k_dsp.size()+6, _k_dsp.size()+6);
    _Kj.diagonal().tail(_k_dsp.size()) = _k_dsp;
    _Kj.diagonal().head(6) = 1.*Eigen::Vector6d::Ones();
    _Dj.setOnes(_d_dsp.size()+6, _d_dsp.size()+6);
    _Dj.diagonal().tail(_d_dsp.size()) = _d_dsp;
    _Dj.diagonal().head(6) = 0.1*Eigen::Vector6d::Ones();

    log();

    return true;
}

void WBTCPlugin::on_start(double time)
{
    _robot->sense(true);

    controller = boost::make_shared<opensot::WBTCController>(_robot->model());

    controller->joint_impedance->setStiffnessDamping(_dynamic_reconfigure._joints_gain*_Kj,
                                                     _dynamic_reconfigure._joints_gain*_Dj);

    log();
}

void WBTCPlugin::control_loop(double time, double period)
{
    _robot->sense(true);

    controller->joint_impedance->setStiffnessDamping(_dynamic_reconfigure._joints_gain*_Kj,
                                                     _dynamic_reconfigure._joints_gain*_Dj);

    if(controller->control(_tau))
    {
        _k_dsp_ref = _dynamic_reconfigure._impedance_gain.load() * _k_dsp;
        _d_dsp_ref = std::sqrt(2.0) * _dynamic_reconfigure._impedance_gain.load() * _d_dsp;

        _robot->setStiffness(_k_dsp_ref);
        _robot->setDamping(_d_dsp_ref);

        _tau_ref = _tau.tail(_robot->model().getActuatedJointNum()) - _tau_offset;
        _robot->setEffortReference(_tau_ref);

        _robot->move();
    }

    log();
}

void WBTCPlugin::log()
{
    _robot->log(_matlogger, XBot::get_time_ns(),"wbtc");
    _matlogger->add("_Kj", _Kj);
    _matlogger->add("_Dj", _Dj);
    _matlogger->add("_tau_offset", _tau_offset);
}

bool WBTCPlugin::close()
{
    _matlogger->flush();
    return true;
}

dynamic_reconf::dynamic_reconf()
{
    dynamic_reconfigure_advr::Server<QPPVM_RT_plugin::QppvmConfig>::CallbackType f;
    f = boost::bind(&dynamic_reconf::cfg_callback, this, _1, _2);
    _server.setCallback(f);
    _impedance_gain.store(1.0);
    _joints_gain.store(0.0);
}

void dynamic_reconf::cfg_callback(QPPVM_RT_plugin::QppvmConfig& config, uint32_t level)
{
    _impedance_gain.store(config.impedance_gain);
//    _stiffness_Waist_gain.store(config.stiffness_waist);
//    _damping_Waist_gain.store(config.damping_waist);
//    _stiffness_Feet_gain.store(config.stiffness_feet);
//    _damping_Feet_gain.store(config.damping_feet);
    _joints_gain.store(config.joints_gain);

    Logger::info(Logger::Severity::HIGH, "\nSetting impedance gain to %f \n", config.impedance_gain);
    Logger::info(Logger::Severity::HIGH, "Setting joints gain to %f \n", config.joints_gain);
//    Logger::info(Logger::Severity::HIGH, "Setting waist stiffness gain to %f and damping gain to %f\n", config.stiffness_waist, config.damping_waist);
//    Logger::info(Logger::Severity::HIGH, "Setting feet stiffness gain to %f and damping gain to %f\n", config.stiffness_feet, config.damping_feet);
}

}
