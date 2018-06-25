#ifndef __QPPVM_FORCE_OPT_H__
#define __QPPVM_FORCE_OPT_H__

#include <OpenSoT/constraints/force/FrictionCone.h>
#include <OpenSoT/Task.h>
#include <OpenSoT/utils/Affine.h>
#include <OpenSoT/solvers/iHQP.h>
#include <OpenSoT/utils/AutoStack.h>
#include <OpenSoT/SubTask.h>
#include <OpenSoT/constraints/GenericConstraint.h>
#include <OpenSoT/constraints/TaskToConstraint.h>
#include <OpenSoT/tasks/MinimizeVariable.h>
#include <OpenSoT/Constraint.h>
#include <qpOASES/Options.hpp>
#include <OpenSoT/tasks/GenericTask.h>
#include <OpenSoT/utils/Piler.h>


namespace demo {
    
    class ForzaGiusta : public OpenSoT::Task<Eigen::MatrixXd, Eigen::VectorXd>
    {
        
    public:
        
        typedef boost::shared_ptr<ForzaGiusta> Ptr;
        
        ForzaGiusta(XBot::ModelInterface::Ptr model, 
                    const std::vector<OpenSoT::AffineHelper>& wrenches_dot,
                    std::vector<std::string> contact_links);
        
        void setFixedBaseTorque(const Eigen::VectorXd& fixed_base_torque, double& period);
        
        void setWrenches(const std::vector<Eigen::VectorXd>& wrenches_lasso);
        
    private:
        
        virtual void _update(const Eigen::VectorXd& x);

        virtual void _log(XBot::MatLogger::Ptr logger);
        
        XBot::ModelInterface::Ptr _model;
        std::vector<std::string> _contact_links;
        std::vector<bool> _enabled_contacts;
        std::vector<OpenSoT::AffineHelper> _wrenches_dot;
        OpenSoT::AffineHelper _task;
        Eigen::MatrixXd _J_i;
        Eigen::MatrixXd _Jfb_i;
        Eigen::VectorXd _x_ref_fb;
        double period_;
        std::vector<Eigen::VectorXd> _wrenches_lasso;
        
        
    };
    
ForzaGiusta::ForzaGiusta(XBot::ModelInterface::Ptr model, 
                         const std::vector< OpenSoT::AffineHelper >& wrenches_dot, 
                         std::vector< std::string > contact_links): 
    Task< Eigen::MatrixXd, Eigen::VectorXd >("FORZA_GIUSTA", wrenches_dot[0].getInputSize()),
    _model(model),
    _contact_links(contact_links),
    _wrenches_dot(wrenches_dot),
    _enabled_contacts(wrenches_dot.size(), true)
{
     _x_ref_fb.setZero(6);
    
    _wrenches_lasso.resize(_contact_links.size());
    for(int i = 0; i < _contact_links.size(); i++)
    _wrenches_lasso[i].setZero(6);

    
    period_=0.001;
    
    update(Eigen::VectorXd());
    
    _W.setIdentity(_A.rows(), _A.rows());
}

void ForzaGiusta::setFixedBaseTorque(const Eigen::VectorXd& fixed_base_torque, double& period)
{
    if(fixed_base_torque.size() != _model->getJointNum())
    {
        throw std::invalid_argument("fixed_base_torque != _model->getJointNum()");
    }
    
         _x_ref_fb = fixed_base_torque.head<6>();
     
//          _x_ref_fb << 0, 0, 1000, 0, 0, 0;

          
     period_=period;
    
    
}

void ForzaGiusta::setWrenches(const std::vector<Eigen::VectorXd>& wrenches_lasso)
{

    for(int i = 0; i < _contact_links.size(); i++)
     _wrenches_lasso[i] = wrenches_lasso[i];
     
        
    
}


void ForzaGiusta::_update(const Eigen::VectorXd& x)
{
     _task.setZero(_wrenches_dot[0].getInputSize(), 6);
    
    for(int i = 0; i < _enabled_contacts.size(); i++)
    {
        if(!_enabled_contacts[i]){
            continue;
        }
        else {
            
             _model->getJacobian(_contact_links[i], _J_i);
            
              _Jfb_i=_J_i.block<6,6>(0,0).transpose(); 
              
              _task = _task + (period_ * _Jfb_i *_wrenches_dot[i]) + (_Jfb_i * _wrenches_lasso[i]);
                                                 
//               _task = _task + (period_ * _wrenches_dot[i]) + ( _wrenches_lasso[i]);
//               _task = _task + (period_ * (i+1) * _wrenches_dot[i]) + ((i+1) * _wrenches_lasso[i]);
            
                     
        }
    }
    
  
    _A = _task.getM();
    _b = -_task.getq() + _x_ref_fb;
    
//      std::cout << "A:\n" << _A << std::endl;
//      std::cout << "b:\n" << _b << std::endl;
    
}

void ForzaGiusta::_log(XBot::MatLogger::Ptr logger)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


 class FrictionCone_wrenches_dot : public OpenSoT::Constraint<Eigen::MatrixXd, Eigen::VectorXd> 
    {
        
     public:
                typedef boost::shared_ptr<FrictionCone_wrenches_dot> Ptr;
                typedef std::pair<std::string, double> friction_cone;
                typedef std::vector<friction_cone> friction_cones;

        friction_cones _mu; //Friction Coefficient associated to each contact surface
        XBot::ModelInterface& _robot;

        Eigen::Matrix<double, 5, 3> _Ci;
        Eigen::Matrix<double, 5, 3> _Cdi;

        std::vector<Eigen::Affine3d> _wTl;

        int _n_of_contacts;

        Eigen::MatrixXd _A;
        Eigen::MatrixXd _A_d;
        Eigen::VectorXd _b;

        OpenSoT::AffineHelper _friction_cone;
        OpenSoT::AffineHelper _wrenches_dot;
        
        Eigen::VectorXd _wrenches_lasso;
        double period_;

    public:


        FrictionCone_wrenches_dot(const std::vector<OpenSoT::AffineHelper>& wrenches_dot,
                                  XBot::ModelInterface &robot,
                                  const friction_cones & mu);


        void update(const Eigen::VectorXd &x);

        void setMu(const friction_cones& mu){ _mu = mu;}
        
        void setWrenches(const std::vector<Eigen::VectorXd>& wrenches_lasso);

        int getNumberOfContacts(){return _n_of_contacts;}

    private:
        void computeAineq();
        void computeUpperBound();

            };
    
 

       FrictionCone_wrenches_dot::FrictionCone_wrenches_dot(const std::vector<OpenSoT::AffineHelper>& wrenches_dot,
                                                            XBot::ModelInterface &robot,
                                                            const friction_cones& mu):
           Constraint("friction_cone", wrenches_dot[0].getInputSize()),
           _robot(robot),
           _mu(mu),
           _Ci(5,3)
       {
           _n_of_contacts = _mu.size();

           _wrenches_dot.setZero(wrenches_dot[0].getInputSize(), 0);
           
           _wrenches_lasso.setZero(_n_of_contacts * 6);

           period_=0.001;

           for(auto& w : wrenches_dot){
               
               _wrenches_dot = _wrenches_dot / w;
           }
           

//            Eigen::Affine3d wTli;
//            for(unsigned int i = 0; i < _n_of_contacts; ++i)
//            {
//                _robot.getPose(_mu[i].first, wTli);
//                _wTl.push_back(wTli);
//            }

           update(Eigen::VectorXd::Zero(0));

       }

       void FrictionCone_wrenches_dot::computeAineq()
       {
           _Ci.setZero(5,3);

           _A.resize(5*_n_of_contacts, 6*_n_of_contacts);

           _A.setZero(_A.rows(), _A.cols());

           for(unsigned int i = 0; i < _n_of_contacts; ++i)
           {

                double __mu = _mu[i].second;
                
//                 __mu = std::sqrt(2.*__mu)/2.;

                _Ci(0,0) = 1.; _Ci(0,1) = 0.; _Ci(0,2) = -__mu;
                _Ci(1,0) = -1.; _Ci(1,1) = 0.; _Ci(1,2) = -__mu;
                _Ci(2,0) = 0.; _Ci(2,1) = 1.; _Ci(2,2) = -__mu;
                _Ci(3,0) = 0.; _Ci(3,1) = -1.; _Ci(3,2) = -__mu;
                _Ci(4,0) = 0.; _Ci(4,1) = 0.; _Ci(4,2) = -1.;

//                  _Ci = _Ci*_wTl[i].linear().transpose();

                _A.block<5,3>(5*i, 6*i) = _Ci;
           }
           
           _A_d.resize(5*_n_of_contacts, 6*_n_of_contacts);

           _A_d.setZero(_A_d.rows(), _A_d.cols());
           
           _Cdi.setZero(5,3);

           for(unsigned int i = 0; i < _n_of_contacts; ++i)
           {

                _Cdi(0,0) = 1.; _Cdi(0,1) = 0.; _Cdi(0,2) = 0;
                _Cdi(1,0) = -1.; _Cdi(1,1) = 0.; _Cdi(1,2) = 0;
                _Cdi(2,0) = 0.; _Cdi(2,1) = 1.; _Cdi(2,2) = 0;
                _Cdi(3,0) = 0.; _Cdi(3,1) = -1.; _Cdi(3,2) = 0;
                _Cdi(4,0) = 0.; _Cdi(4,1) = 0.; _Cdi(4,2) = -1.;

                _A_d.block<5,3>(5*i, 6*i) = _Cdi;
           }
           

       }

       void FrictionCone_wrenches_dot::computeUpperBound()
       {
           _b.resize(5*_n_of_contacts);
           _b.setZero(_b.rows());
           _bLowerBound.resize(5*_n_of_contacts);
           _bLowerBound = -1.0e20*_bLowerBound.setOnes(_bLowerBound.size());
       }

       void FrictionCone_wrenches_dot::setWrenches(const std::vector<Eigen::VectorXd>& wrenches_lasso)
       {
           
           _wrenches_lasso << wrenches_lasso[0], wrenches_lasso[1], wrenches_lasso[2],wrenches_lasso[3];
         
//            std::cout<<"_wrenches_lasso: \n"<<_wrenches_lasso<<std::endl;
 
       }

       void FrictionCone_wrenches_dot::update(const Eigen::VectorXd &x)
       {

               _n_of_contacts = _mu.size();

               computeAineq();
               computeUpperBound();

               _friction_cone =  (_A_d *  period_ * _wrenches_dot) + (_A *  _wrenches_lasso) - _b;   
               _Aineq = _friction_cone.getM();
               _bUpperBound = - _friction_cone.getq();
               
//                 std::cout<<"_Aineq: \n"<<_Aineq<<std::endl;

       }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


 
 
    class ForceOptimization 
    {
      
    public:
        
        typedef boost::shared_ptr<ForceOptimization> Ptr;
        
        ForceOptimization(XBot::ModelInterface::Ptr model,
                          std::vector<std::string> contact_links
                          );
        
        bool compute(const Eigen::VectorXd& fixed_base_torque, 
                     double& _start_time,
                     double&  time,
                     double& period,
                     std::vector<Eigen::VectorXd>& Fc,
                     std::vector<Eigen::VectorXd>& dFc,
                     Eigen::VectorXd& tau
                    );
        
        void log(XBot::MatLogger::Ptr logger);
        
        
    private:
        
        XBot::ModelInterface::Ptr _model;
        std::vector<std::string> _contact_links;
        std::vector< OpenSoT::AffineHelper > _aux_var;
        std::vector< OpenSoT::AffineHelper > _wrenches_dot;
        std::vector< OpenSoT::AffineHelper > _wrenches_min_var;
        std::vector< OpenSoT::AffineHelper > _slack_var;
        std::vector< OpenSoT::AffineHelper > _t_var;
        
        OpenSoT::constraints::force::FrictionCone::Ptr _friction_cone;
        ForzaGiusta::Ptr _forza_giusta;
        FrictionCone_wrenches_dot::Ptr _friction_constraint;
        OpenSoT::solvers::iHQP::Ptr _solver;
        OpenSoT::AutoStack::Ptr _autostack;
        OpenSoT::tasks::MinimizeVariable::Ptr min_variation;
        OpenSoT::tasks::MinimizeVariable::Ptr min_variation2;
        OpenSoT::tasks::MinimizeVariable::Ptr min_slack;
        OpenSoT::tasks::MinimizeVariable::Ptr task_fake;
        double start_time_, time_;
        
        
        Eigen::VectorXd _dx_value;
        Eigen::MatrixXd _JC;

        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_dot_bounds_ub;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_dot_bounds_lb;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> _wrench_lasso_ub;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> _wrench_lasso_lb;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_bounds;       
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> aux_var_bounds;
        
        std::vector<Eigen::VectorXd> Fc_old;
        std::vector<Eigen::VectorXd> dFc;
        
        Eigen::VectorXd wrench_ub, wrench_lb;
        Eigen::VectorXd inf_b;
        Eigen::VectorXd ub_slack, lb_slack, t_ub, t_lb;            
        Eigen::VectorXd wrench_dot_ub, wrench_dot_lb;
        Eigen::VectorXd opt_ub, opt_lb;
        
        OpenSoT::AffineHelper t_var;
                    

    };
    
}




demo::ForceOptimization::ForceOptimization(XBot::ModelInterface::Ptr model, 
                                           std::vector< std::string > contact_links):
    _model(model),
    _contact_links(contact_links)
{  
    /* Do we want to consider contact torques? */
    const bool optimize_contact_torque = false;
    
//////////////////////////////////////////////////////////////////////////////////////////// 

    /* Define optimization vector by stacking all contact wrenches */
    OpenSoT::OptvarHelper::VariableVector vars;

    for(auto cl : _contact_links)
    {
        vars.emplace_back(cl+"_wrenches_dot",6);
        vars.emplace_back(cl+"_t_ub",1);
        vars.emplace_back(cl+"_t_lb",1);
        vars.emplace_back(cl+"_slack_ub",3);
        vars.emplace_back(cl+"_slack_lb",3);
    }
    
    
    Fc_old.resize(_contact_links.size());
    for(int i = 0; i < _contact_links.size(); i++)
    Fc_old[i].setZero(6);

    OpenSoT::OptvarHelper opt(vars);
       
    /* Wrench bounds */
    wrench_ub.setZero(6); wrench_ub[0] =  1000; wrench_ub[1] =  1000; wrench_ub[2] =  1000;
    wrench_lb.setZero(6); wrench_lb[0] = -1000; wrench_lb[1] = -1000; wrench_lb[2] = -1000;
    
    inf_b.setOnes(6); inf_b*=1e30; 
               
    wrench_dot_ub.setOnes(6); wrench_dot_ub *=  1e3; 
    wrench_dot_lb.setOnes(6); wrench_dot_lb *= -1e3; 
    
    Eigen::VectorXd _aux_var_ub, _aux_var_lb;
    _aux_var_ub.setOnes(8);_aux_var_ub*=1e30;
    _aux_var_lb.setZero(8);
//     _aux_var_ub.setZero(2);
//     _aux_var_lb.setZero(2);

        
    
    /* Define affine mappings for all wrenches */
    for(auto cl : _contact_links){
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        /* Auxiliary variables */
  
        _aux_var.emplace_back(opt.getVariable(cl+"_t_ub")/opt.getVariable(cl+"_t_lb")/opt.getVariable(cl+"_slack_ub")/opt.getVariable(cl+"_slack_lb"));
//         _aux_var.emplace_back(opt.getVariable(cl+"_t_ub")/opt.getVariable(cl+"_t_lb"));///opt.getVariable(cl+"_slack_ub")/opt.getVariable(cl+"_slack_lb"));
        
               
        
        aux_var_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_bound",
                                                                                              _aux_var.back(),
                                                                                              _aux_var_ub,
                                                                                              _aux_var_lb,
                                                                                               OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
       
                              );    
        
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* forza giusta */        
      
        OpenSoT::AffineHelper wrench_dot = opt.getVariable(cl+"_wrenches_dot");
         _wrenches_dot.emplace_back(wrench_dot);
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        /* Wrench bounds */
 
        wrench_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_bound",
                                                                                             _wrenches_dot.back(),
                                                                                              wrench_ub,
                                                                                              wrench_lb,
                                                                                              OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
       
                              );   
        
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* Lasso wrench upper bound */ 
       
       Eigen::VectorXd zero,inf;
       zero.setZero(1); inf.setOnes(1);inf*=1e30;       
  
       OpenSoT::AffineHelper wrench_lasso_ub = 0.001 * wrench_dot.segment(2,1) + (-1.0) * opt.getVariable(cl+"_t_ub");
        
       _wrench_lasso_ub.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_dot_ub",
                                                                                               wrench_lasso_ub,
                                                                                               zero,
                                                                                               -inf,
                                                                                               OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
                               ); 
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
         /* Lasso wrench upper bound */   

      
        OpenSoT::AffineHelper wrench_lasso_lb =   0.001 * wrench_dot.segment(2,1) +  opt.getVariable(cl+"_t_lb");
        
         _wrench_lasso_lb.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_dot_lb",
                                                                                                 wrench_lasso_lb,
                                                                                                 inf,
                                                                                                 zero,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
         
                               );                                
        
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* wrench_dot upper bound with slack variables */ 
       
        Eigen::MatrixXd tmp;
        tmp.setZero(6,3);
        tmp.block<3,3>(0,0).setIdentity();
        
        
        OpenSoT::AffineHelper wrench_dot_slack_ub = wrench_dot + (-1.0) * tmp * opt.getVariable(cl+"_slack_ub");
        
         wrench_dot_bounds_ub.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_dot_ub",
                                                                                                 wrench_dot_slack_ub,
                                                                                                 wrench_dot_ub,
                                                                                                 -inf_b,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
                               ); 
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        /* wrench_dot lower bound with slack variables */  

      
        OpenSoT::AffineHelper wrench_dot_slack_lb =   wrench_dot + tmp * opt.getVariable(cl+"_slack_lb");
        
         wrench_dot_bounds_lb.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_dot_lb",
                                                                                                 wrench_dot_slack_lb,
                                                                                                 inf_b,
                                                                                                 wrench_dot_lb,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
         
                               );                                
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* min slack variables */ 
         
        _slack_var.emplace_back(opt.getVariable(cl+"_slack_ub")/opt.getVariable(cl+"_slack_lb"));
      
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* Lasso */ 

        _t_var.emplace_back(opt.getVariable(cl+"_t_ub") + opt.getVariable(cl+"_t_lb"));

 
    }  
    
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    /* Define friction cones */
    OpenSoT::constraints::force::FrictionCone::friction_cones friction_cones;

    for(auto cl : _contact_links)
    friction_cones.emplace_back(cl, .3);
         
    _friction_constraint = boost::make_shared<FrictionCone_wrenches_dot>(_wrenches_dot, *_model, friction_cones);
    
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    /* Construct forza giusta task */
    _forza_giusta = boost::make_shared<ForzaGiusta>(_model,_wrenches_dot, _contact_links);
    
     /* Construct min slack */          
     OpenSoT::AffineHelper  slack;
     slack = _slack_var[0]/_slack_var[1]/_slack_var[2]/_slack_var[3];
    
     min_slack = boost::make_shared<OpenSoT::tasks::MinimizeVariable>("MIN slack",slack);

    
     /* Construct Lasso */ 
//      t_var = _t_var[0] + _t_var[1] +  _t_var[2] +  _t_var[3];
     t_var = _t_var[0] + 20* _t_var[1] + 30 * _t_var[2] +  40* _t_var[3]; /* WEIGHTED LASSO */
         
     std::cout<<"g_Lasso: \n"<<t_var.getM().transpose()<<std::endl;
         

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
//    _autostack = boost::make_shared<OpenSoT::AutoStack>( 1e-2 * min_slack + _forza_giusta);
//               
// 
//    for(int i = 0; i < _contact_links.size(); i++)
//    {      
//     _autostack  << aux_var_bounds[i] << wrench_bounds[i] << _wrench_lasso_ub[i] << _wrench_lasso_lb[i] << wrench_dot_bounds_ub[i]<< wrench_dot_bounds_lb[i];
//    }
//    
//     _autostack << _friction_constraint;
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
    
  OpenSoT::tasks::Aggregated::Ptr aggr1 = (1e-2 * min_slack + _forza_giusta); 
 
  task_fake = boost::make_shared<OpenSoT::tasks::MinimizeVariable>("task_fake", t_var);
     

   for(int i = 0; i < _contact_links.size(); i++)
   {      
    task_fake  << _wrench_lasso_ub[i] << _wrench_lasso_lb[i];
   }
   
   _autostack = (aggr1/task_fake); 
//    _autostack = boost::make_shared<OpenSoT::AutoStack>( aggr1 );
   
   for(int i = 0; i < _contact_links.size(); i++)
   {      
    _autostack   << aux_var_bounds[i] <<  wrench_bounds[i] << wrench_dot_bounds_ub[i]<< wrench_dot_bounds_lb[i];
   }
   
    _autostack << _friction_constraint;
         
 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
    
    _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 0.0, OpenSoT::solvers::solver_back_ends::OSQP);
  

}


bool demo::ForceOptimization::compute(const Eigen::VectorXd& fixed_base_torque,
                                      double& _start_time,
                                      double& time,
                                      double& period,
                                      std::vector<Eigen::VectorXd>& Fc,
                                      std::vector<Eigen::VectorXd>& dFc,
                                      Eigen::VectorXd& tau)
{
    Fc.resize(_contact_links.size());
     
    dFc.resize(_contact_links.size());

    std::vector<Eigen::VectorXd> opt_c;
    opt_c.resize(_contact_links.size());
    
    start_time_= _start_time;
    time_= time;
    
    double period_= period;
    
    _forza_giusta->setFixedBaseTorque(fixed_base_torque, period_);  
    
    
    _autostack->update(Eigen::VectorXd());

/////////////////////////////////////////////////////////////////////////////////////     
/////////////////////////////////////////////////////////////////////////////////////  
    
//      Eigen::VectorXd g_Lasso =   1e-2 * t_var.getM().transpose();     
//      _solver->setUserG(g_Lasso);
// //      std::cout<<"g: \n" << _solver->getBackEnd(0)->getg() <<std::endl;
//       
//     if(!_solver->solve(_dx_value))
//     {
//         std::cout<<"FORCE OPTIMIZATION CANNOT SOLVE!"<<std::endl;
//         return false;
//     }
//     
/////////////////////////////////////////////////////////////////////////////////////     
/////////////////////////////////////////////////////////////////////////////////////     
 
    _solver->setUserSolveFlag(1, false);

    
    if(!_solver->solve(_dx_value))
    {
        std::cout<<"FORCE OPTIMIZATION CANNOT SOLVE!"<<std::endl;
        return false;
    }
    
    Eigen::MatrixXd tmpH;
    tmpH.setZero(_solver->getBackEnd(0)->getNumVariables(),_solver->getBackEnd(0)->getNumVariables());
   
    _solver->getBackEnd(1)->updateTask(tmpH, t_var.getM().transpose());
    
//      std::cout<<"H: \n" << _solver->getBackEnd(1)->getH() <<std::endl;
//      std::cout<<"g: \n" << _solver->getBackEnd(1)->getg() <<std::endl;
    
    _solver->getBackEnd(1)->solve();
    _dx_value = _solver->getBackEnd(1)->getSolution();
    
/////////////////////////////////////////////////////////////////////////////////////     
/////////////////////////////////////////////////////////////////////////////////////  
    
    tau = fixed_base_torque;  
    
    for(int i = 0; i < _contact_links.size(); i++)
    {
         _wrenches_dot[i].getValue(_dx_value, opt_c[i]);
         
        dFc[i]=opt_c[i].head(6);
                        
        Fc[i] =   Fc_old[i] + dFc[i] * period;
        Fc_old[i] = Fc[i];
        
         std::cout<<"Fc"<< i+1 << ": \n" << Fc[i].head(3) <<std::endl;
               
        
        _model->getJacobian(_contact_links[i], _JC);
        
         tau.noalias() -= _JC.transpose()*Fc[i];                          
         
         /* Wrench bounds */
         wrench_bounds[i]->setBounds((1.0/period_)*(wrench_ub-Fc[i]),(1.0/period_)*(wrench_lb-Fc[i]));
 
         /* Lasso bounds */       
         Eigen::VectorXd zero(1),inf(1);
         zero.setZero(); inf.setOnes();inf*=1e30;   
         
         _wrench_lasso_ub[i]->setBounds(zero - Fc[i].segment(2,1), -inf - Fc[i].segment(2,1));
         _wrench_lasso_lb[i]->setBounds( inf - Fc[i].segment(2,1), zero - Fc[i].segment(2,1));
                     
    }
  
    
    _forza_giusta->setWrenches(Fc);
    _friction_constraint->setWrenches(Fc);
     



    
    return true;
    
}

void demo::ForceOptimization::log(XBot::MatLogger::Ptr logger)
{
    _autostack->log(logger);
    _solver->log(logger);
}


#endif
