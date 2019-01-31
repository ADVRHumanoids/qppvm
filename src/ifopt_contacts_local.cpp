#include <iostream>

#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
// #include <ifopt_problem/ifopt_test.h>
#include <ifopt_problem/ifopt_contacts_local.h>

#include <XBotInterface/MatLogger.hpp>

using namespace ifopt;

int main()
{
    
   auto logger = XBot::MatLogger::getLogger("/tmp/ifopt_contacts_log");  
    
  // 1. define the problem
  Problem nlp;
  
  auto p1 = std::make_shared<ExVariables>("p1");
  auto p2 = std::make_shared<ExVariables>("p2");
  auto p3 = std::make_shared<ExVariables>("p3");
  auto p4 = std::make_shared<ExVariables>("p4");
  
  auto F1 = std::make_shared<ExVariables>("F1");
  auto F2 = std::make_shared<ExVariables>("F2");
  auto F3 = std::make_shared<ExVariables>("F3");
  auto F4 = std::make_shared<ExVariables>("F4");
  
//   auto n1 = std::make_shared<ExVariables>("n1");
//   auto n2 = std::make_shared<ExVariables>("n2");
//   auto n3 = std::make_shared<ExVariables>("n3");
//   auto n4 = std::make_shared<ExVariables>("n4");
  
  auto fr_F1 = std::make_shared<FrictionConstraint>("F1");
  auto fr_F2 = std::make_shared<FrictionConstraint>("F2");
  auto fr_F3 = std::make_shared<FrictionConstraint>("F3");
  auto fr_F4 = std::make_shared<FrictionConstraint>("F4");
  
  
  nlp.AddVariableSet(p1);
  nlp.AddVariableSet(p2);
  nlp.AddVariableSet(p3);
  nlp.AddVariableSet(p4);
  
  nlp.AddVariableSet(F1);
  nlp.AddVariableSet(F2);
  nlp.AddVariableSet(F3);
  nlp.AddVariableSet(F4);
  
//   nlp.AddVariableSet(n1);
//   nlp.AddVariableSet(n2);
//   nlp.AddVariableSet(n3);
//   nlp.AddVariableSet(n4);
  
   // BOUNDS 
  p1->SetBounds(Eigen::Vector3d( 0.1, -1.0, 0.0),Eigen::Vector3d( 2.0, -0.1, 0.4));
  p2->SetBounds(Eigen::Vector3d( 0.1,  0.1, 0.0),Eigen::Vector3d( 2.0,  1.0, 0.4));
  p3->SetBounds(Eigen::Vector3d(-2.0, -1.0, 0.0),Eigen::Vector3d(-0.1, -0.1, 0.4));
  p4->SetBounds(Eigen::Vector3d(-2.0,  0.1, 0.0),Eigen::Vector3d(-0.1,  1.0, 0.4));
 
  Eigen::Vector3d F_max; 
  F_max.setOnes();
  F_max *= 100;
  
  F1->SetBounds(-F_max,F_max);
  F2->SetBounds(-F_max,F_max);
  F3->SetBounds(-F_max,F_max);
  F4->SetBounds(-F_max,F_max);
  
  // SUPERELLIPSOID PARAMETERS
  Eigen::Vector3d C; C << 1.5, 0.0, 1.0;
  Eigen::Vector3d R; R << 2, 5, 1;
  Eigen::Vector3d P; P << 8, 8, 4; 
  
  
  // CENTROIDAL DYNAMICS
  auto static_constr = std::make_shared<StaticConstraint>();
  Eigen::Vector6d ext_w;
  ext_w << 80.0, 0, 0, 0, 0, 0.0;
  static_constr->SetExternalWrench(ext_w);
  
  Eigen::Vector3d com;
  com << -0.1, 0.0, 0.5;
  static_constr->SetCoM(com);
  
  static_constr->setParam_SE(C,R,P);
  
  nlp.AddConstraintSet(static_constr);
  
  
  // SUPERELLIPSOID ENVIRONMENT
  auto SE_p1 = std::make_shared<SuperEllipsoidConstraint>("p1");
  auto SE_p2 = std::make_shared<SuperEllipsoidConstraint>("p2");
  auto SE_p3 = std::make_shared<SuperEllipsoidConstraint>("p3");
  auto SE_p4 = std::make_shared<SuperEllipsoidConstraint>("p4");
  
  SE_p1->SetParam(C,R,P);
  SE_p2->SetParam(C,R,P);
  SE_p3->SetParam(C,R,P);
  SE_p4->SetParam(C,R,P);
  
  nlp.AddConstraintSet(SE_p1);
  nlp.AddConstraintSet(SE_p2);
  nlp.AddConstraintSet(SE_p3);
  nlp.AddConstraintSet(SE_p4);
   
  // FRICTION CONES
  double mu = .2;
  
  fr_F1->set_mu(mu);
  fr_F2->set_mu(mu);
  fr_F3->set_mu(mu);
  fr_F4->set_mu(mu);
 
  nlp.AddConstraintSet(fr_F1);
  nlp.AddConstraintSet(fr_F2);
  nlp.AddConstraintSet(fr_F3);
  nlp.AddConstraintSet(fr_F4); 
   
  // COST
  auto cost = std::make_shared<ExCost>();
  Eigen::VectorXd p_ref;  
  p_ref.setZero(12);
  
  p_ref <<  0.5, -0.3, 0.0, 
            0.5,  0.3, 0.0, 
           -0.5, -0.3, 0.0, 
           -0.5,  0.3, 0.0;
  
  double W_p = 10;
  
  cost->SetPosRef(p_ref,W_p);
  
  nlp.AddCostSet(cost);

  
  nlp.PrintCurrent();
  

  // 2. choose solver and options
  IpoptSolver ipopt;
//   ipopt.SetOption("linear_solver", "ma27");
//   ipopt.SetOption("jacobian_approximation", "exact");
//   ipopt.SetOption("jacobian_approximation", "finite-difference-values");
//   ipopt.SetOption("max_iter", 6000);
//   ipopt.SetOption("tol", 1e-3);
//   ipopt.SetOption("constr_viol_tol",1e-3);
//   ipopt.SetOption("mu_strategy", "adaptive");

  

  // 3 . solve
  ipopt.Solve(nlp);
  Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();
//   std::cout << x.transpose() << std::endl;
  for(int i = 0; i < 4; i++)
  std::cout<<"p"<< i+1 << ": \n" << x.segment(i*3, 3).transpose() <<std::endl;

  for(int i = 0; i < 4; i++)
  std::cout<<"F"<< i+1 << ": \n" << x.segment(i*3 + 12, 3).transpose() <<std::endl;

  // 4. test if solution correct
  double eps = 1e-5; //double precision
  assert(1.0-eps < x(0) && x(0) < 1.0+eps);
  assert(0.0-eps < x(1) && x(1) < 0.0+eps);
  
   for(int k = 0; k < x.size(); k++)
   {
        logger->add("x_sol", x[k]);
   }
  
    logger->flush();
  
}