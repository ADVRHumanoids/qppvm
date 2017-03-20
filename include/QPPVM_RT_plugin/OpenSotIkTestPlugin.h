#include <XCM/XBotControlPlugin.h>
#include <OpenSoT/tasks/velocity/Cartesian.h>
#include <OpenSoT/tasks/velocity/Postural.h>
#include <OpenSoT/constraints/velocity/JointLimits.h>
#include <OpenSoT/constraints/velocity/VelocityLimits.h>
#include <OpenSoT/utils/AutoStack.h>


class OpenSotIkTestPlugin : public XBot::XBotControlPlugin {

public:

    virtual bool init_control_plugin(std::string path_to_config_file,
                                    XBot::SharedMemory::Ptr shared_memory,
                                    XBot::RobotInterface::Ptr robot);

    virtual void on_start(double time);

    virtual void control_loop(double time, double period);

    virtual bool close();


private:

    double _start_time;

    Eigen::VectorXd _q0, _q, _dq, _qhome;

    XBot::RobotInterface::Ptr _robot;
    XBot::ModelInterface::Ptr _model;

    XBot::SharedObject<Eigen::Affine3d> _left_ref, _right_ref;

    OpenSoT::tasks::velocity::Cartesian::Ptr _left_ee, _right_ee;
    OpenSoT::tasks::velocity::Postural::Ptr _postural;
    OpenSoT::constraints::velocity::JointLimits::Ptr _joint_lims;
    OpenSoT::constraints::velocity::VelocityLimits::Ptr _joint_vel_lims;

    OpenSoT::AutoStack::Ptr _autostack;
    OpenSoT::solvers::QPOases_sot::Ptr _solver;

    XBot::MatLogger::Ptr _logger;

};


