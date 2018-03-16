#include <GcompFixedBase/GcompFixedBase.h>

REGISTER_XBOT_PLUGIN_(XBotPlugin::GcompFixedBase)


bool XBotPlugin::GcompFixedBase::init_control_plugin(XBot::Handle::Ptr handle)
{
    _logger = XBot::MatLogger::getLogger("/tmp/gcomp_fixed_base_log");
    _robot = handle->getRobotInterface();
    _model = XBot::ModelInterface::getModel(handle->getPathToConfigFile());
    _imu = _robot->getImu().at("imu_link");

    _model->computeGravityCompensation(_gcomp);

    _tau_offset.assign(4, Eigen::VectorXd::Zero(6));
    _klegs.assign(4, Eigen::VectorXd::Zero(6));
    _dlegs.assign(4, Eigen::VectorXd::Zero(6));
    _klegs_0.assign(4, Eigen::VectorXd::Zero(6));
    _dlegs_0.assign(4, Eigen::VectorXd::Zero(6));

    for(int i = 0; i < 4; i++)
    {
        _robot->leg(i).getStiffness(_klegs_0[i]);
        _robot->leg(i).getDamping(_dlegs_0[i]);
        _impedance_gain.emplace_back(new std::atomic<float>(1.0));

        auto cb = std::bind(&GcompFixedBase::callback, this, std::placeholders::_1, i);

        auto sub = handle->getRosHandle()->subscribe<std_msgs::Float64>("/impedance_gain_leg_" + std::to_string(i),
                                                                        1, cb
                                                                        );

        _imp_sub.push_back(sub);

    }

    /*** HARD CODE HR LEG OFFSET ***/

    _tau_offset[0] << 0.0,  0.0, -7.0, 12.0, 0.0,  0.0;
    _tau_offset[1] << 0.0,  13.0, -2.0, 0.0, 0.0,  0.0;
    _tau_offset[2] << 0.0, -2.0,  2.0, -3.8, 0.0,  0.0;
    _tau_offset[3] << 0.0,  1.0,  -3.0, 0.0, 0.0,  0.0;

    _tau_ref_leg.setZero(6);
    _leg_impedance.setZero(6);

    return true;

}

void XBotPlugin::GcompFixedBase::on_start(double time)
{
    for(int i = 0; i < 4; i++)
    {
        _robot->leg(i).getStiffness(_klegs_0[i]);
        _robot->leg(i).getDamping(_dlegs_0[i]);
    }
}


void XBotPlugin::GcompFixedBase::control_loop(double time, double period)
{
    _model->syncFrom(*_robot, XBot::Sync::Position, XBot::Sync::Velocity);
    _model->setFloatingBaseState(_imu);
    _model->update();

    _model->computeGravityCompensation(_gcomp);
    _model->setJointEffort(_gcomp);


    for(int i = 0; i < 4; i++)
    {
        _model->leg(i).getJointEffort(_tau_ref_leg);
        _tau_ref_leg += _tau_offset[i];
        _model->leg(i).setJointEffort(_tau_ref_leg);

        _leg_impedance = _impedance_gain[i]->load() * _klegs_0[i];
        _leg_impedance(5) = _klegs_0[i](5);
        _robot->leg(i).setStiffness(_leg_impedance);

        _leg_impedance = std::sqrt(_impedance_gain[i]->load()) * _dlegs_0[i];
        _leg_impedance.array() += 5.0;
        _leg_impedance(5) = _dlegs_0[i](5);
        _robot->leg(i).setDamping(_leg_impedance);
    }

    for(int i = 0; i < 7; i++)
    {
        _model->chain("left_arm").setJointEffort(i, 0.0);
        _model->chain("right_arm").setJointEffort(i, 0.0);
    }

    _model->chain("torso").setJointEffort(0, 0.0);

    _robot->setReferenceFrom(*_model, XBot::Sync::Effort);
    _robot->move();
}

void XBotPlugin::GcompFixedBase::on_stop(double time)
{
}

void XBotPlugin::GcompFixedBase::callback(const std_msgs::Float64ConstPtr& msg, int leg_id)
{
    float gain = msg->data;
    if( gain < 0 || gain > 1){
        return;
    }

    _impedance_gain[leg_id]->store(gain);

}

