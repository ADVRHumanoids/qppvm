#include <ifopt_problem/ifopt_contacts_node.h>
#include <RobotInterfaceROS/ConfigFromParam.h>


ModelInterface::Ptr model;
boost::shared_ptr<cartesian_interface_handler> _ci;
Eigen::VectorXd q;


int main(int argc, char **argv)
{
    /* Init ROS node */
    ros::init(argc, argv, "ipopt_contcats_node");
    ros::NodeHandle nh;

    double rate;
    nh.param("rate", rate, 100.);

    ros::Rate loop_rate(rate);

    ConfigOptions config = XBot::ConfigOptionsFromParamServer();

    model = ModelInterface::getModel(config);
    q.setZero(model->getJointNum());

    std::map<std::string, double> map;
    if(nh.getParam("zeros", map))
    {
        ROS_INFO("Loading zeros param");
        setZeros(q, *model, map);
    }

    bool log;
    nh.param("log", log, false);

    XBot::MatLogger::Ptr logger;
    uint T = ros::Time::now().sec;
    std::stringstream ss;
    ss<<"/tmp/ipopt_contacts_node_"<<T;

    if(log)
        logger = XBot::MatLogger::getLogger(ss.str());

    _ci = boost::make_shared<cartesian_interface_handler>(model);

    double time = 0;
    double period = loop_rate.expectedCycleTime().toSec();

    while(ros::ok())
    {
        model->setJointPosition(q);
        model->update();

        _ci->ci_ros->run();

        _ci->ci->update(time, period);

        time += period;

        ros::spinOnce();
        loop_rate.sleep();
    }

    if(log)
        logger->flush();
    return 0;
}
