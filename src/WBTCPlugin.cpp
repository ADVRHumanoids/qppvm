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
    Waist = boost::make_shared<OpenSoT::tasks::torque::CartesianImpedanceCtrl>("Waist", _q, _model, "Waist", "world");//, id);
    LFoot = boost::make_shared<OpenSoT::tasks::torque::CartesianImpedanceCtrl>("LFoot", _q, _model, "l_ankle", "world");//, id);
    RFoot = boost::make_shared<OpenSoT::tasks::torque::CartesianImpedanceCtrl>("RFoot", _q, _model, "r_ankle", "world");//, id);
    

    Eigen::VectorXd tau_lims;
    _model.getEffortLimits(tau_lims);

    torque_lims = boost::make_shared<OpenSoT::constraints::torque::TorqueLimits>(tau_lims, -tau_lims);
    
    _autostack = ((LFoot+RFoot)/Waist/joint_impedance)<<torque_lims;


    _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 1.);

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
    _K_Foot.diagonal()<<1,1,1,0.25,0.25,0.25;
    _D_Foot.diagonal()<<1,1,1,0.25,0.25,0.25;
    _K_Lfoot_ref = _K_Rfoot_ref = _K_Foot;
    _D_Lfoot_ref = _D_Rfoot_ref = _D_Foot;
    _K_Waist_ref = _K_Waist = _K_Foot;
    _D_Waist_ref = _D_Waist = _D_Foot;
    
    _use_offsets = true;
    
    log();
    
    /* Init cartesian ifc */
    YAML::Node yaml_file = YAML::LoadFile(handle->getPathToConfigFile());
    XBot::Cartesian::ProblemDescription ik_problem(yaml_file["CartesianInterface"]["problem_description"], _model);
    _ci = std::make_shared<XBot::Cartesian::CartesianInterfaceImpl>(_model, ik_problem);
    _ci->enableOtg(0.002);
    _sync_from_nrt = std::make_shared<XBot::Cartesian::Utils::SyncFromIO>("/xbotcore/cartesian_interface", handle->getSharedMemory());


#ifdef USE_GAZEBO_GROUND_TRUTH
    _floating_base_ground_truth = std::make_shared<gazebo::FloatingBaseGroundTruth>(handle);
#endif

    return true;
}

void WBTCPlugin::on_start(double time)
{
    sense(0.0);

    controller = boost::make_shared<opensot::WBTCController>(*_model);

    std::vector<std::string> links_in_contact = {"l_ankle","r_ankle"};
    Eigen::Vector6d F; F.setZero(); _Fc.push_back(F); _Fc.push_back(F);
    forza_giusta = boost::make_shared<OpenSoT::utils::ForceOptimization>(_model, links_in_contact);

    _floating_base_differential_kineamtics = boost::make_shared<OpenSoT::floating_base_estimation::qp_estimation>(
                _model, (*(_robot->getImu().begin())).second, links_in_contact);

    _floating_base_forward_kinematics = boost::make_shared<OpenSoT::floating_base_estimation::kinematic_estimation>(
                _model, links_in_contact[0]);

    _Kj_ref = _dynamic_reconfigure._joints_gain*_Kj;
    _Dj_ref = _dynamic_reconfigure._joints_gain*_Dj;
    controller->joint_impedance->setStiffnessDamping(_Kj_ref, _Dj_ref);

    _K_Lfoot_ref = _K_Rfoot_ref = _dynamic_reconfigure._stiffness_Feet_gain*_K_Foot;
    _D_Lfoot_ref = _D_Rfoot_ref = _dynamic_reconfigure._damping_Feet_gain*_D_Foot;
    _K_Waist_ref = _dynamic_reconfigure._stiffness_Waist_gain*_K_Waist;
    _D_Waist_ref = _dynamic_reconfigure._damping_Waist_gain*_D_Waist;
    
    controller->LFoot->setStiffnessDamping(_K_Lfoot_ref, _D_Lfoot_ref);
    controller->RFoot->setStiffnessDamping(_K_Rfoot_ref, _D_Rfoot_ref);
    controller->Waist->setStiffnessDamping(_K_Waist_ref, _D_Waist_ref);
                                                     

    log();
    controller->log(_matlogger);
    
    _first_sync_done = false;
}

void WBTCPlugin::sense(const double dT)
{
    _model->syncFrom(*_robot, XBot::Sync::Position, XBot::Sync::Velocity, XBot::Sync::MotorSide);

    if(_floating_base_differential_kineamtics)
    {
        if(!_floating_base_differential_kineamtics->update(dT))
            XBot::Logger::error("_floating_base_differential_kineamtics->update() returned false!");
    }

    if(_floating_base_forward_kinematics)
    {
        _floating_base_forward_kinematics->update();
    }

#ifdef USE_GAZEBO_GROUND_TRUTH
       _floating_base_pose_gazebo =     _floating_base_ground_truth->getFloatingBasePose();
       _flaoting_base_velocity_gazebo = _floating_base_ground_truth->getFloatingBaseTwist();
//    _model->setFloatingBaseState(_floating_base_ground_truth->getFloatingBasePose(),
//                                _floating_base_ground_truth->getFloatingBaseTwist());
#endif



    _model->update();
}

void WBTCPlugin::set_dyn_reconfigure_gains()
{
    _Kj_ref = _dynamic_reconfigure._joints_gain*_Kj;
    _Dj_ref = std::sqrt(_dynamic_reconfigure._joints_gain)*_Dj;
    controller->joint_impedance->setStiffnessDamping(_Kj_ref, _Dj_ref);
    
    _K_Lfoot_ref = _K_Rfoot_ref = _dynamic_reconfigure._stiffness_Feet_gain*_K_Foot;
    _D_Lfoot_ref = _D_Rfoot_ref = _dynamic_reconfigure._damping_Feet_gain*_D_Foot;
    _K_Waist_ref = _dynamic_reconfigure._stiffness_Waist_gain*_K_Waist;
    _D_Waist_ref = _dynamic_reconfigure._damping_Waist_gain*_D_Waist;
    
    controller->LFoot->setStiffnessDamping(_K_Lfoot_ref, _D_Lfoot_ref);
    controller->RFoot->setStiffnessDamping(_K_Rfoot_ref, _D_Rfoot_ref);
    controller->Waist->setStiffnessDamping(_K_Waist_ref, _D_Waist_ref);
}

void WBTCPlugin::control_loop(double time, double period)
{
    controller->setFilter(period, 10.); //Matteo dice che 12 lo usano in COMAU
    
    sense(period);

    set_dyn_reconfigure_gains();
    sync_cartesian_ifc(time, period);

    if(controller->control(_tau)) //HERE WE COMPUTE TAU_BAR
    {
        if(forza_giusta->compute(_tau, _Fc, _tau_actuated)) //HERE WE COMPUTE TAU_ACTUATED
        {
            _k_dsp_ref = _dynamic_reconfigure._impedance_gain.load() * _k_dsp;
            _d_dsp_ref = std::sqrt(2.0) * _dynamic_reconfigure._impedance_gain.load() * _d_dsp;

            _robot->setStiffness(_k_dsp_ref);
            _robot->setDamping(_d_dsp_ref);

            if(_use_offsets)
                _tau_ref = _tau_actuated.tail(_model->getActuatedJointNum()) + _tau_offset;
            else
                _tau_ref = _tau_actuated.tail(_model->getActuatedJointNum());
        
            _robot->setEffortReference(_tau_ref);

            _robot->move();
        }
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
    _matlogger->add("K_Waist_ref", _K_Waist_ref);
    _matlogger->add("D_Waist_ref", _D_Waist_ref);
    _matlogger->add("tau_offset", _tau_offset);
    _matlogger->add("tau_actuated", _tau_actuated);

#ifdef USE_GAZEBO_GROUND_TRUTH
    if(_floating_base_ground_truth){
        _matlogger->add("floating_base_pose_gazebo", _floating_base_pose_gazebo.matrix());
    }   _matlogger->add("flaoting_base_velocity_gazebo",_flaoting_base_velocity_gazebo);
#endif
    if(_floating_base_differential_kineamtics)
        _floating_base_differential_kineamtics->log(_matlogger);
}

bool WBTCPlugin::close()
{
    _matlogger->flush();
    return true;
}

void WBTCPlugin::on_stop(double time)
{
    _matlogger->flush();
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

     _ci->getPoseReference(controller->Waist->getDistalLink(), T_ref);
     T_ref_matrix = T_ref.matrix();
     controller->Waist->setReference(T_ref_matrix);
}

}
