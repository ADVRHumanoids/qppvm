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


class IFOPT_problem
{
public:
    IFOPT_problem(ModelInterface& model):
        _model(model)
    {
        _ext_w << 90, 0, 0, 0, 0, 0.0;
        _F_max.setOnes(); _F_max *= 100;
        _C << 1.5, 0.0, 1.0;
	_R << 2, 20, 1;
	_P << 8, 8, 4; 
	_mu = 0.2;
	_Wp = 10;
	
	_p_ref.setZero(12);
	_p_ref <<  0.5, -0.3, 0.0, 
		   0.5,  0.3, 0.0, 
		  -0.5, -0.3, 0.0, 
		  -0.5,  0.3, 0.0;
      
        if(!create_problem())
            throw std::runtime_error("Problems during creation of IPOPT problem!");
	
    }

    bool create_problem()
    {
        
	Eigen::VectorXd q;
        _model.getJointPosition(q);		
	
	auto p1 = std::make_shared<ExVariables>("p1");
	auto p2 = std::make_shared<ExVariables>("p2");
	auto p3 = std::make_shared<ExVariables>("p3");
	auto p4 = std::make_shared<ExVariables>("p4");
      
	auto F1 = std::make_shared<ExVariables>("F1");
	auto F2 = std::make_shared<ExVariables>("F2");
	auto F3 = std::make_shared<ExVariables>("F3");
	auto F4 = std::make_shared<ExVariables>("F4"); 
      
	auto n1 = std::make_shared<ExVariables>("n1");
	auto n2 = std::make_shared<ExVariables>("n2");
	auto n3 = std::make_shared<ExVariables>("n3");
	auto n4 = std::make_shared<ExVariables>("n4");
      
	auto com = std::make_shared<ExVariables>("com");
	
	auto static_constr = std::make_shared<StaticConstraint>();
	
	auto SE_p1 = std::make_shared<SuperEllipsoidConstraint>("p1");
	auto SE_p2 = std::make_shared<SuperEllipsoidConstraint>("p2");
	auto SE_p3 = std::make_shared<SuperEllipsoidConstraint>("p3");
	auto SE_p4 = std::make_shared<SuperEllipsoidConstraint>("p4");
	
	auto fr_F1 = std::make_shared<FrictionConstraint>("F1");
	auto fr_F2 = std::make_shared<FrictionConstraint>("F2");
	auto fr_F3 = std::make_shared<FrictionConstraint>("F3");
	auto fr_F4 = std::make_shared<FrictionConstraint>("F4");
	
	auto n_p1 = std::make_shared<NormalConstraint>("p1");
	auto n_p2 = std::make_shared<NormalConstraint>("p2");
	auto n_p3 = std::make_shared<NormalConstraint>("p3");
	auto n_p4 = std::make_shared<NormalConstraint>("p4");
      
	auto cost = std::make_shared<ExCost>();
	
	
	_nlp.AddVariableSet(p1); p1->SetBounds(Eigen::Vector3d( 0.1, -1.0, 0.0),Eigen::Vector3d( 2.0, -0.1, 0.4));
	_nlp.AddVariableSet(p2); p2->SetBounds(Eigen::Vector3d( 0.1,  0.1, 0.0),Eigen::Vector3d( 2.0,  1.0, 0.4));
	_nlp.AddVariableSet(p3); p3->SetBounds(Eigen::Vector3d(-2.0, -1.0, 0.0),Eigen::Vector3d(-0.1, -0.1, 0.4));
	_nlp.AddVariableSet(p4); p4->SetBounds(Eigen::Vector3d(-2.0,  0.1, 0.0),Eigen::Vector3d(-0.1,  1.0, 0.4));
  
	_nlp.AddVariableSet(F1); F1->SetBounds(-_F_max,_F_max);
	_nlp.AddVariableSet(F2); F2->SetBounds(-_F_max,_F_max);
	_nlp.AddVariableSet(F3); F3->SetBounds(-_F_max,_F_max);
	_nlp.AddVariableSet(F4); F4->SetBounds(-_F_max,_F_max);
  
	_nlp.AddVariableSet(n1);
	_nlp.AddVariableSet(n2);
	_nlp.AddVariableSet(n3);
	_nlp.AddVariableSet(n4); 
  
	_nlp.AddVariableSet(com); com->SetBounds(Eigen::Vector3d(-2.0, -1.0, 0.4),Eigen::Vector3d( 2.0, 1.0, 0.6));
	
	static_constr->SetExternalWrench(_ext_w); _nlp.AddConstraintSet(static_constr);
	
	SE_p1->SetParam(_C,_R,_P); _nlp.AddConstraintSet(SE_p1);
	SE_p2->SetParam(_C,_R,_P); _nlp.AddConstraintSet(SE_p2);
	SE_p3->SetParam(_C,_R,_P); _nlp.AddConstraintSet(SE_p3);
	SE_p4->SetParam(_C,_R,_P); _nlp.AddConstraintSet(SE_p4);
	
	fr_F1->set_mu(_mu); _nlp.AddConstraintSet(fr_F1);
	fr_F2->set_mu(_mu); _nlp.AddConstraintSet(fr_F2);
	fr_F3->set_mu(_mu); _nlp.AddConstraintSet(fr_F3);
	fr_F4->set_mu(_mu); _nlp.AddConstraintSet(fr_F4); 
	
	n_p1->SetParam(_C,_R,_P); _nlp.AddConstraintSet(n_p1);
	n_p2->SetParam(_C,_R,_P); _nlp.AddConstraintSet(n_p2);
	n_p3->SetParam(_C,_R,_P); _nlp.AddConstraintSet(n_p3);
	n_p4->SetParam(_C,_R,_P); _nlp.AddConstraintSet(n_p4); 
	
	cost->SetPosRef(_p_ref,_Wp); _nlp.AddCostSet(cost);
	
// 	solve(_x);
	
// 	_nlp.PrintCurrent();
		
        return true;
    }
    
    bool solve(Eigen::VectorXd& x)
    {       
	_ipopt.SetOption("derivative_test", "first-order");     
	_ipopt.Solve(_nlp);	
	
	x = _nlp.GetOptVariables()->GetValues();

//         std::cout<<"com: \n"<< x.tail(3).transpose() <<std::endl; 
// 	
// 	for(int i = 0; i < 4; i++)
// 	std::cout<<"p"<< i+1 << ": \n" << x.segment(i*3, 3).transpose() <<std::endl;
// 	
// 	for(int i = 0; i < 4; i++)
// 	std::cout<<"F"<< i+1 << ": \n" << x.segment(i*3 + 12, 3).transpose() <<std::endl;
	
	return true;
    }

    void log(XBot::MatLogger::Ptr logger)
    {
      for(int k = 0; k < _x.size(); k++)
      {
	logger->add("x_sol", _x[k]);
      }
    }
   
private:

    ModelInterface& _model;  
    
    Problem _nlp;
    IpoptSolver _ipopt;
  
    Eigen::VectorXd _x;
    Eigen::Vector3d _F_max; 
    Eigen::Vector6d _ext_w;
    Eigen::Vector3d _C, _R, _P; 
    Eigen::VectorXd _p_ref; 
    double _mu, _Wp;
    

};



#endif
