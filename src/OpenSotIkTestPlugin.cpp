#include <QPPVM_RT_plugin/OpenSotIkTestPlugin.h>

REGISTER_XBOT_PLUGIN(OpenSotIkTestPlugin, OpenSotIkTestPlugin)

bool OpenSotIkTestPlugin::init_control_plugin(std::string path_to_config_file,
                                              XBot::SharedMemory::Ptr shared_memory,
                                              XBot::RobotInterface::Ptr robot)
{
    _logger = XBot::MatLogger::getLogger("/tmp/OpenSotIkTestPlugin_logger");

    _robot = robot;
    _model = XBot::ModelInterface::getModel(path_to_config_file);

    _robot->sense();
    _robot->model().getJointPosition(_q0);

    _model->getRobotState("home", _qhome);

    _left_ref = shared_memory->get<Eigen::Affine3d>("w_T_left_ee");
    _right_ref = shared_memory->get<Eigen::Affine3d>("w_T_right_ee");

    _left_ref.reset(new Eigen::Affine3d);
    _right_ref.reset(new Eigen::Affine3d);
    
    std::vector<bool> active_joints(_model->getJointNum(), true);
    active_joints[_model->getDofIndex(_model->chain("torso").getJointId(0))] = false;

    /* Create cartesian tasks for both hands */
    _left_ee.reset( new OpenSoT::tasks::velocity::Cartesian("CARTESIAN_LEFT",
                                                            _qhome,
                                                            *_model,
                                                            _model->chain("left_arm").getTipLinkName(),
                                                            "world"
                                                            ) );
    _left_ee->setActiveJointsMask(active_joints);
//     _left_ee->setLambda(100);

    _right_ee.reset( new OpenSoT::tasks::velocity::Cartesian("CARTESIAN_RIGHT",
                                                             _qhome,
                                                             *_model,
                                                             _model->chain("right_arm").getTipLinkName(),
                                                             "world"
                                                             ) );
    _right_ee->setActiveJointsMask(active_joints);
//     _right_ee->setLambda(100);

    /* Create postural task */
    _postural.reset( new OpenSoT::tasks::velocity::Postural(_qhome) );
    Eigen::VectorXd weight;
    weight.setOnes((_model->getJointNum()));
    weight(0) = 100;
    _postural->setWeight(weight.asDiagonal());

    /* Create joint limits & velocity limits */
    Eigen::VectorXd qmin, qmax, qdotmax;
    _model->getJointLimits(qmin, qmax);
    _model->getVelocityLimits(qdotmax);
    double qdotmax_min = qdotmax.minCoeff();
    _final_qdot_lim = 2.0;

    _joint_lims.reset( new OpenSoT::constraints::velocity::JointLimits(_q0, qmax, qmin) );

    _joint_vel_lims.reset( new OpenSoT::constraints::velocity::VelocityLimits(qdotmax_min, 0.001, _model->getJointNum()) );

    /* Create autostack and set solver */
    _autostack = ( (/*_right_ee +*/ _left_ee) / _postural ) << _joint_lims << _joint_vel_lims;

    _solver.reset( new OpenSoT::solvers::QPOases_sot(_autostack->getStack(), _autostack->getBounds()) );
    
    // Logger
    
    Eigen::Affine3d left_pose, right_pose;
    _logger->add("left_ref_pos", _left_ref->translation());
    _logger->add("right_ref_pos", _right_ref->translation());
    _logger->add("left_actual_pos", left_pose.translation());
    _logger->add("right_actual_pos", right_pose.translation());
    
    _logger->add("left_ref_or", _left_ref->linear());
    _logger->add("right_ref_or", _right_ref->linear());
    _logger->add("left_actual_or", left_pose.linear());
    _logger->add("right_actual_or", right_pose.linear());
    _logger->add("computed_q", _q0);
    
    _logger->add("computed_qdot", _q0);
    
    _logger->add("time", 0.0);


    return true;
}

void OpenSotIkTestPlugin::on_start(double time)
{

    _model->syncFrom(*_robot);
    _model->getJointPosition(_q);
    

    _model->getPose(_left_ee->getDistalLink(), *_left_ref);
    _model->getPose(_right_ee->getDistalLink(), *_right_ref);
    
    _start_time = time;

    std::cout << "OpenSotIkTestPlugin STARTED!\nInitial q is " << _q.transpose() << std::endl;
    std::cout << "Home q is " << _qhome.transpose() << std::endl;
}

void OpenSotIkTestPlugin::control_loop(double time, double period)
{
    /* Model update */
    _model->setJointPosition(_q);
    _model->update();
    
    /* HACK: shape vel lim to avoid discontinuity */
    double alpha = 0;
    alpha = (time - _start_time)/10;
    alpha = alpha > 1 ? 1 : alpha;
    
    _joint_vel_lims->setVelocityLimits( (0 + alpha*(_final_qdot_lim - 0)) );

    /* Simple upward reference motion */
//     _right_ref->translation().y() += 0.05*period;

    /* Set cartesian tasks reference */
    _left_ee->setReference(_left_ref->matrix());
    _right_ee->setReference(_right_ref->matrix());

    /* Log data */
    Eigen::Affine3d left_pose, right_pose;
    _model->getPose(_left_ee->getDistalLink(), left_pose);
    _model->getPose(_right_ee->getDistalLink(), right_pose);
    
    _logger->add("left_ref_pos", _left_ref->translation());
    _logger->add("right_ref_pos", _right_ref->translation());
    _logger->add("left_actual_pos", left_pose.translation());
    _logger->add("right_actual_pos", right_pose.translation());
    
    _logger->add("left_ref_or", _left_ref->linear());
    _logger->add("right_ref_or", _right_ref->linear());
    _logger->add("left_actual_or", left_pose.linear());
    _logger->add("right_actual_or", right_pose.linear());
    _logger->add("computed_q", _q);
    
    _logger->add("time", time);


    /* Stack update and solve */
    _autostack->update(_q);

    _dq.setZero(_model->getJointNum());

    if( !_solver->solve(_dq) ){
        std::cerr << "UNABLE TO SOLVE" << std::endl;
        return;
    }
    
    _logger->add("computed_qdot", _dq/period);

    /* Update q */
    _q += _dq;
    

    /* Send command to motors */
    _robot->setReferenceFrom(*_model);
    _robot->move();

}

bool OpenSotIkTestPlugin::close()
{
    _logger->flush();
    return true;
}