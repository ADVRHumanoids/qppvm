#include <QPPVM_RT_plugin/WBTCPlugin.h>

REGISTER_XBOT_PLUGIN(WBTCPlugin, xbotcore::WBTCPlugin)

namespace opensot{

WBTCController::WBTCController(XBot::ModelInterface &model):
    _model(model)
{
    _h.setZero(_model.getJointNum());
    _tau_opt.setZero(_model.getJointNum());

    _model.getJointPosition(_q);
    _model.getJointVelocity(_qdot);
    _qdot_filtered = _qdot;
    _alpha_filter = 1.;

    joint_impedance = boost::make_shared<OpenSoT::tasks::torque::JointImpedanceCtrl>(_q, _model);
    std::list<unsigned int> id = {0, 1, 2};
    LFoot = boost::make_shared<OpenSoT::tasks::torque::CartesianImpedanceCtrl>("LFoot", _q, _model, "l_ankle", "Waist", id);
    RFoot = boost::make_shared<OpenSoT::tasks::torque::CartesianImpedanceCtrl>("RFoot", _q, _model, "r_ankle", "Waist", id);

    Eigen::VectorXd tau_lims;
    _model.getEffortLimits(tau_lims);

    if(_model.isFloatingBase())
    {
        for(unsigned int i = 0; i < 6; ++i)
            tau_lims[i] = 0.0;
    }
    torque_lims = boost::make_shared<OpenSoT::constraints::torque::TorqueLimits>(tau_lims, -tau_lims);
    

    _autostack = ((LFoot+RFoot)/joint_impedance)<<torque_lims;
//    _autostack = boost::make_shared<OpenSoT::AutoStack>(joint_impedance);
//    _autostack<<torque_lims;

    _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 1.);
    

    ///////
//     J.setZero(6, model.getJointNum());
//     Jtf.setZero(model.getJointNum());
//     spring.setZero(6);
//     damper.setZero(6);
//     force.setZero(6);
}

void WBTCController::setFilter(const double period, const double cut_off_freq)
{
    _period = period;
    _cut_off_freq = cut_off_freq;
    _alpha_filter = 1-std::exp(-2.*M_PI*cut_off_freq*period);
}

bool WBTCController::control(Eigen::VectorXd& tau)
{
    _model.getJointPosition(_q);
    _model.getJointVelocity(_qdot);
    
    //HERE WE FILTER and we apply the filtered velocities to the model
    _qdot_filtered.noalias() = _alpha_filter*_qdot + (1.-_alpha_filter)*_qdot_filtered;
    _model.setJointPosition(_q);
    _model.setJointVelocity(_qdot_filtered);
    _model.update();
    
    _model.computeNonlinearTerm(_h);

    _autostack->update(_q);

    if(!_solver->solve(_tau_opt))
        _tau_opt.setZero(_tau_opt.size());

    
    
//     _model.getRelativeJacobian("l_sole","Waist",J);
//     LFoot->getSpringForce(spring); LFoot->getDamperForce(damper);
//     force = spring+damper;
//     Jtf.noalias() = J.transpose()*force;
//     _tau_opt.segment(0,6) = Jtf.segment(0,6);
//     _model.getRelativeJacobian("r_sole","Waist",J);
//     RFoot->getSpringForce(spring); RFoot->getDamperForce(damper);
//     force = spring+damper;
//     Jtf.noalias() = J.transpose()*force;
//     _tau_opt.segment(6,6) = Jtf.segment(6,6);
    
    
    tau = _tau_opt +  _h;

    return true;
}


void WBTCController::log(XBot::MatLogger::Ptr logger)
{
    logger->add("h", _h);
    logger->add("tau_opt", _tau_opt);
    logger->add("q", _q);
    logger->add("qdot", _qdot);
    logger->add("qdot_filtered", _qdot_filtered);
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
    _model = XBot::ModelInterface::getModel(_robot->getPathToConfig());
    _robot->sense(true);

    _k_dsp.setZero(_robot->getJointNum());
    _d_dsp.setZero(_robot->getJointNum());
    _tau_ref.setZero(_robot->getJointNum());
    _tau.setZero(_model->getJointNum());
    _tau_offset.setZero(_robot->getJointNum());    
    
    if(!_robot->getRobotState("torque_offset", _tau_offset))
        _tau_offset.setZero(_robot->getJointNum());


    _robot->getStiffness(_k_dsp);
    _robot->getDamping(_d_dsp);
    _k_dsp_ref = _k_dsp;
    _d_dsp_ref = _d_dsp;

    _Kj_vec.setZero(_robot->getJointNum());
    _Dj_vec.setZero(_robot->getJointNum());

    if(!_robot->getRobotState("joint_stiffness", _Kj_vec))
    {
        XBot::Logger::error("joint_stiffness param missing\n");
        return false;
    }
    for(unsigned int i = 0; i < _Kj_vec.size(); ++i)
    {
        if(_Kj_vec[i] < 0.0)
        {
            XBot::Logger::error("joint_stiffness value can not be negative!\n");
            return false;
        }
    }

    if(!_robot->getRobotState("joint_damping", _Dj_vec))
    {
        XBot::Logger::error("joint_damping param missing\n");
        return false;
    }
    for(unsigned int i = 0; i < _Dj_vec.size(); ++i)
    {
        if(_Dj_vec[i] < 0.0)
        {
            XBot::Logger::error("joint_damping value can not be negative!\n");
            return false;
        }
    }


    if(_model->isFloatingBase())
    {
        _Kj.setZero(_k_dsp.size()+6, _k_dsp.size()+6);
        _Kj.block(6,6,_k_dsp.size(),_k_dsp.size()) = _Kj_vec.asDiagonal();
    
        _Dj.setZero(_d_dsp.size()+6, _d_dsp.size()+6);
        _Dj.block(6,6,_d_dsp.size(),_d_dsp.size()) = _Dj_vec.asDiagonal();
    }
    else
    {
        _Kj.setZero(_k_dsp.size(), _k_dsp.size());
        _Kj = _Kj_vec.asDiagonal();
    
        _Dj.setZero(_d_dsp.size(), _d_dsp.size());
        _Dj = _Dj_vec.asDiagonal();
    }
        
    
    _Kj_ref = _Kj;
    _Dj_ref = _Dj;

    _K_Foot.setIdentity(6,6); _D_Foot.setIdentity(6,6);
    _K_Lfoot_ref = _K_Rfoot_ref = _K_Foot;
    _D_Lfoot_ref = _D_Rfoot_ref = _D_Foot;
    
    _use_offsets = true;
    
    log();
    
    /* Init cartesian ifc */
    YAML::Node yaml_file = YAML::LoadFile(handle->getPathToConfigFile());
    XBot::Cartesian::ProblemDescription ik_problem(yaml_file["CartesianInterface"]["problem_description"], _model);
    _ci = std::make_shared<XBot::Cartesian::CartesianInterfaceImpl>(_model, ik_problem);
    _ci->enableOtg(0.002);
    _sync_from_nrt = std::make_shared<XBot::Cartesian::Utils::SyncFromIO>("/xbotcore/cartesian_interface", handle->getSharedMemory());

    return true;
}

void WBTCPlugin::on_start(double time)
{
    _model->syncFrom(*_robot, XBot::Sync::Position, XBot::Sync::Velocity, XBot::Sync::MotorSide);

    controller = boost::make_shared<opensot::WBTCController>(*_model);

    _Kj_ref = _dynamic_reconfigure._joints_gain*_Kj;
    _Dj_ref = _dynamic_reconfigure._joints_gain*_Dj;
    controller->joint_impedance->setStiffnessDamping(_Kj_ref, _Dj_ref);

    _K_Lfoot_ref = _K_Rfoot_ref = _dynamic_reconfigure._stiffness_Feet_gain*_K_Foot;
    _D_Lfoot_ref = _D_Rfoot_ref = _dynamic_reconfigure._damping_Feet_gain*_D_Foot;
    
    controller->LFoot->setStiffnessDamping(_K_Lfoot_ref, _D_Lfoot_ref);
    controller->RFoot->setStiffnessDamping(_K_Rfoot_ref, _D_Rfoot_ref);
                                                     

    log();
    controller->log(_matlogger);
    
    _first_sync_done = false;
    
    
}

void WBTCPlugin::set_dyn_reconfigure_gains()
{
    _Kj_ref = _dynamic_reconfigure._joints_gain*_Kj;
    _Dj_ref = _dynamic_reconfigure._joints_gain*_Dj;
    controller->joint_impedance->setStiffnessDamping(_Kj_ref, _Dj_ref);
    
    _K_Lfoot_ref = _K_Rfoot_ref = _dynamic_reconfigure._stiffness_Feet_gain*_K_Foot;
    _D_Lfoot_ref = _D_Rfoot_ref = _dynamic_reconfigure._damping_Feet_gain*_D_Foot;
    
    controller->LFoot->setStiffnessDamping(_K_Lfoot_ref, _D_Lfoot_ref);
    controller->RFoot->setStiffnessDamping(_K_Rfoot_ref, _D_Rfoot_ref);
}

void WBTCPlugin::control_loop(double time, double period)
{
    //controller->setFilter(period, 10.); //Matteo dice che 12 lo usano in COMAU
    
    _model->syncFrom(*_robot, XBot::Sync::Position, XBot::Sync::Velocity, XBot::Sync::MotorSide);

    set_dyn_reconfigure_gains();
    sync_cartesian_ifc(time, period);

    if(controller->control(_tau))
    {
        _k_dsp_ref = _dynamic_reconfigure._impedance_gain.load() * _k_dsp;
        _d_dsp_ref = std::sqrt(2.0) * _dynamic_reconfigure._impedance_gain.load() * _d_dsp;

        _robot->setStiffness(_k_dsp_ref);
        _robot->setDamping(_d_dsp_ref);

        if(_use_offsets)
            _tau_ref = _tau.tail(_model->getActuatedJointNum()) + _tau_offset;
        else
            _tau_ref = _tau.tail(_model->getActuatedJointNum());
        
        _robot->setEffortReference(_tau_ref);

        _robot->move();
    }

    log();
    controller->log(_matlogger);
}

void WBTCPlugin::log()
{
    _robot->log(_matlogger, XBot::get_time_ns(),"wbtc_");
    _matlogger->add("Kj", _Kj);
    _matlogger->add("Dj", _Dj);
    _matlogger->add("k_dsp", _k_dsp);
    _matlogger->add("d_dsp", _d_dsp);
    _matlogger->add("Kj_ref", _Kj_ref);
    _matlogger->add("Dj_ref", _Dj_ref);
    _matlogger->add("K_Lfoot_ref", _K_Lfoot_ref);
    _matlogger->add("D_Lfoot_ref", _D_Lfoot_ref);
    _matlogger->add("K_Rfoot_ref", _K_Rfoot_ref);
    _matlogger->add("D_Rfoot_ref", _D_Rfoot_ref);
    _matlogger->add("tau_offset", _tau_offset);
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
    _stiffness_Feet_gain.store(0.0);
    _damping_Feet_gain.store(0.0);
}

void dynamic_reconf::cfg_callback(QPPVM_RT_plugin::QppvmConfig& config, uint32_t level)
{
    _impedance_gain.store(config.impedance_gain);
//    _stiffness_Waist_gain.store(config.stiffness_waist);
//    _damping_Waist_gain.store(config.damping_waist);
    _stiffness_Feet_gain.store(config.stiffness_feet);
    _damping_Feet_gain.store(config.damping_feet);
    _joints_gain.store(config.joints_gain);

    Logger::info(Logger::Severity::HIGH, "\nSetting impedance gain to %f \n", config.impedance_gain);
    Logger::info(Logger::Severity::HIGH, "Setting joints gain to %f \n", config.joints_gain);
//    Logger::info(Logger::Severity::HIGH, "Setting waist stiffness gain to %f and damping gain to %f\n", config.stiffness_waist, config.damping_waist);
    Logger::info(Logger::Severity::HIGH, "Setting feet stiffness gain to %f and damping gain to %f\n", config.stiffness_feet, config.damping_feet);
}

void WBTCPlugin::sync_cartesian_ifc(double time, double period)
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
    static Eigen::MatrixXd T_ref_matrix(4,4);
    _ci->getPoseReference(controller->LFoot->getDistalLink(), T_ref);
    T_ref_matrix = T_ref.matrix();
    controller->LFoot->setReference(T_ref_matrix);
    _ci->getPoseReference(controller->RFoot->getDistalLink(), T_ref);
    T_ref_matrix = T_ref.matrix();
    controller->RFoot->setReference(T_ref_matrix);
}

}
