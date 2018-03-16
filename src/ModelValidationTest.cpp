#include <ros/ros.h>
#include <XBotInterface/RobotInterface.h>
#include <XBotInterface/Utils.h>

int current_target = 0;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "model_validation_test");
    ros::NodeHandle nh;
    
    auto robot = XBot::RobotInterface::getRobot(XBot::Utils::getXBotConfig());
    robot->sense();
    
    auto model_imu = XBot::ModelInterface::getModel(XBot::Utils::getXBotConfig());
    auto model_raw = XBot::ModelInterface::getModel(XBot::Utils::getXBotConfig());
    auto imu = robot->getImu().at("imu_link");
    
    auto logger = XBot::MatLogger::getLogger("/tmp/model_validation_log");
    
    Eigen::VectorXd tau_m, gcomp_raw, gcomp_imu, nl_raw, nl_imu;
    Eigen::VectorXd qm, qdotm, qref;
    
    Eigen::VectorXd fl_ref, 
                    fr_ref, 
                    hl_ref, 
                    hr_ref;
                    
    robot->chain("front_left_leg").getJointPosition(fl_ref);
    robot->chain("front_right_leg").getJointPosition(fr_ref);
    robot->chain("rear_left_leg").getJointPosition(hl_ref);
    robot->chain("rear_right_leg").getJointPosition(hr_ref);
    
    std::vector<Eigen::VectorXd> poses_fl, poses_fr, poses_hl, poses_hr;
    poses_fl.assign(3, Eigen::VectorXd(6));
    poses_fl[0] = fl_ref;
    poses_fl[1] << 0.0, -1.0, -0.4, 0.1, 0.0, 0.0;
    poses_fl[2] << 0.0, -1.5, 0.0, 0.0, 0.0, 0.0;
    
    poses_fr.assign(3, Eigen::VectorXd(6));
    poses_fr[0] = fr_ref;
    poses_fr[1] = -poses_fl[1];
    poses_fr[2] = -poses_fl[2];
    
    poses_hl.assign(3, Eigen::VectorXd(6));
    poses_hl[0] = hl_ref;
    poses_hl[1] = -poses_fl[1];
    poses_hl[2] = -poses_fl[2];
    
    poses_hr.assign(3, Eigen::VectorXd(6));
    poses_hr[0] = hr_ref;
    poses_hr[1] = poses_fl[1];
    poses_hr[2] = poses_fl[2];
    
    ros::Rate loop_rate(100);
    double time = 0.0;
    
    std::vector<double> time_v = {0.0, 20.0};
    
    while(ros::ok())
    {
        robot->sense();
        
        model_imu->syncFrom(*robot, XBot::Sync::All);
        model_imu->setFloatingBaseState(imu);
        
        model_raw->syncFrom(*robot, XBot::Sync::All);
        
        model_imu->computeGravityCompensation(gcomp_imu);
        model_imu->computeNonlinearTerm(nl_imu);
        
        model_raw->computeGravityCompensation(gcomp_raw);
        model_raw->computeNonlinearTerm(nl_raw);
        model_raw->getJointEffort(tau_m);
        
        model_imu->getJointPosition(qm);
        model_imu->getJointVelocity(qdotm);
        
        /* Manage references */
        double start_time, end_time;
        int segment_id = 0;
        for(int i = 0; i < time_v.size(); i++)
        {
            if(time >= time_v[i])
            {
                start_time = time_v[i];
                end_time = start_time + 10.0;
                segment_id = i;
            }
        }
        
        double alpha = 0, dalpha = 0, ddalpha = 0;
        XBot::Utils::FifthOrderPlanning(alpha, dalpha, ddalpha, 
                                        1.0, 
                                        start_time, end_time, time, 
                                        alpha, dalpha, ddalpha);
        
        /* Set reference to legs */
        fl_ref = (1-alpha) * poses_fl[segment_id] + alpha * poses_fl[segment_id+1];
        fr_ref = (1-alpha) * poses_fr[segment_id] + alpha * poses_fr[segment_id+1];
        hl_ref = (1-alpha) * poses_hl[segment_id] + alpha * poses_hl[segment_id+1];
        hr_ref = (1-alpha) * poses_hr[segment_id] + alpha * poses_hr[segment_id+1];
        
        robot->chain("front_left_leg").setPositionReference(fl_ref);
        robot->chain("front_right_leg").setPositionReference(fr_ref);
        robot->chain("rear_left_leg").setPositionReference(hl_ref);
        robot->chain("rear_right_leg").setPositionReference(hr_ref);
        
        
        robot->getPositionReference(qref);
        
        logger->add("gcomp_imu", gcomp_imu);
        logger->add("nl_imu", nl_imu);
        logger->add("gcomp_raw", gcomp_raw);
        logger->add("nl_raw", nl_raw);
        logger->add("tau_m", tau_m);
        logger->add("qm", qm);
        logger->add("qdotm", qdotm);
        logger->add("qref", qref);
        
        robot->move();
        
        time += loop_rate.expectedCycleTime().toSec();
        loop_rate.sleep();
        
    }
    
    
    logger->flush();
    
}