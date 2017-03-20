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

    _left_ref = shared_memory->advertise<Eigen::Affine3d>("w_T_left_ee");
    _right_ref = shared_memory->advertise<Eigen::Affine3d>("w_T_right_ee");

    _left_ref.reset(new Eigen::Affine3d);
    _right_ref.reset(new Eigen::Affine3d);

    /* Create cartesian tasks for both hands */
    _left_ee.reset( new OpenSoT::tasks::velocity::Cartesian("CARTESIAN_LEFT",
                                                            _q0,
                                                            *_model,
                                                            _model->chain("left_arm").getTipLinkName(),
                                                            "world"
                                                            ) );
    _left_ee->setLambda(100);

    _right_ee.reset( new OpenSoT::tasks::velocity::Cartesian("CARTESIAN_RIGHT",
                                                             _q0,
                                                             *_model,
                                                             _model->chain("right_arm").getTipLinkName(),
                                                             "world"
                                                             ) );
    _right_ee->setLambda(100);

    /* Create postural task */
    _postural.reset( new OpenSoT::tasks::velocity::Postural(_qhome) );

    /* Create joint limits & velocity limits */
    Eigen::VectorXd qmin, qmax, qdotmax;
    _model->getJointLimits(qmin, qmax);
    _model->getVelocityLimits(qdotmax);
    double qdotmax_min = qdotmax.minCoeff();

    _joint_lims.reset( new OpenSoT::constraints::velocity::JointLimits(_q0, qmax, qmin) );

    _joint_vel_lims.reset( new OpenSoT::constraints::velocity::VelocityLimits(qdotmax_min, 0.001, _model->getJointNum()) );

    /* Create autostack and set solver */
    _autostack = ( (_right_ee + _left_ee) / _postural ) << _joint_lims << _joint_vel_lims;

    _solver.reset( new OpenSoT::solvers::QPOases_sot(_autostack->getStack(), _autostack->getBounds()) );

    return true;
}

void OpenSotIkTestPlugin::on_start(double time)
{

    _model->syncFrom(*_robot);
    _model->getJointPosition(_q);

    _model->getPose(_left_ee->getDistalLink(), *_left_ref);
    _model->getPose(_right_ee->getDistalLink(), *_right_ref);

    std::cout << "OpenSotIkTestPlugin STARTED!\nInitial q is " << _q.transpose() << std::endl;
}

void OpenSotIkTestPlugin::control_loop(double time, double period)
{
    /* Model update */
    _model->setJointPosition(_q);
    _model->update();

    /* Simple upward reference motion */
    _right_ref->translation().z() += 0.05*period;

    /* Set cartesian tasks reference */
    _left_ee->setReference(_left_ref->matrix());
    _right_ee->setReference(_right_ref->matrix());

    /* Log data */
    Eigen::Vector3d left_pos, right_pos;
    _model->getPointPosition(_left_ee->getDistalLink(), Eigen::Vector3d::Zero(), left_pos);
    _model->getPointPosition(_right_ee->getDistalLink(), Eigen::Vector3d::Zero(), right_pos);
    _logger->add("left_ref", _left_ref->translation());
    _logger->add("right_ref", _right_ref->translation());
    _logger->add("left_actual", left_pos);
    _logger->add("right_actual", right_pos);


    /* Stack update and solve */
    _autostack->update(_q);

    _dq.setZero(_model->getJointNum());

    if( !_solver->solve(_dq) ){
        std::cerr << "UNABLE TO SOLVE" << std::endl;
    }

    /* Update q */
    _q += period * _dq;

    /* Send command to motors */
    _robot->setReferenceFrom(*_model);
    _robot->move();

}

bool OpenSotIkTestPlugin::close()
{
    _logger->flush();
    return true;
}
