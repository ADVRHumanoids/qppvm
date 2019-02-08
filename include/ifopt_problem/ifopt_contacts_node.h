#ifndef _CONTACTS_NODE_H_
#define _CONTACTS_NODE_H_

#include <ros/ros.h>
#include <cartesian_interface/CartesianInterfaceImpl.h>
#include <cartesian_interface/ros/RosServerClass.h>
#include <cartesian_interface/problem/Com.h>

#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt_problem/ifopt_contacts.h>

using namespace XBot;
using namespace XBot::Cartesian;
using namespace ifopt;

static void setZeros(Eigen::VectorXd& q, ModelInterface& model, const std::map<std::string, double>& zeros)
{
    q.setZero(model.getJointNum());

    std::map<std::string, double>::const_iterator it = zeros.begin();
    for(it; it != zeros.end(); it++)
    {
        std::string joint_name = it->first;
        double value = it->second;

        int id = model.getDofIndex(joint_name);
        if(id == -1)
            ROS_WARN("Joint %s does not exists! Skip!", joint_name);
        else
        {
            ROS_INFO("Set %f to joint %s", value, joint_name.c_str());
            q[id] = value;
        }
    }

    model.setJointPosition(q);
    model.update();
}

class cartesian_interface_handler
{
public:
    cartesian_interface_handler(ModelInterface::Ptr model):
        _model(model)
    {
        create_cartesian_interface();
    }

    void create_cartesian_interface()
    {
        _lHand = std::make_shared<CartesianTask>("l_wrist");
        _rHand = std::make_shared<CartesianTask>("r_wrist");
        _lsole = std::make_shared<CartesianTask>("l_sole");
        _rsole = std::make_shared<CartesianTask>("r_sole");

        AggregatedTask aggr1;
        aggr1.push_back(_lsole);
        aggr1.push_back(_rsole);
        AggregatedTask aggr2;
        aggr2.push_back(_lHand);
        aggr2.push_back(_rHand);

        Stack stack;
        stack.push_back(aggr1);
        stack.push_back(aggr2);

        _pd = boost::make_shared<ProblemDescription>(stack);

        ci = std::make_shared<CartesianInterfaceImpl>(_model, *_pd);

        ci_ros = std::make_shared<RosServerClass>(ci, _model, _opt);
    }

    CartesianInterfaceImpl::Ptr ci;
    RosServerClass::Ptr ci_ros;
private:
    boost::shared_ptr<ProblemDescription> _pd;

    CartesianTask::Ptr _lHand;
    CartesianTask::Ptr _rHand;
    CartesianTask::Ptr _lsole;
    CartesianTask::Ptr _rsole;

    ModelInterface::Ptr _model;

    RosServerClass::Options _opt;
};

#endif
