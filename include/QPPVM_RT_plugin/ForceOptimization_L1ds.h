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
     
//          _x_ref_fb << 1100, 1100, 1100, 0, 0, 0;

          
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
                                                 
//              _task = _task + (period_ * _wrenches_dot[i]) + ( _wrenches_lasso[i]);
            
                     
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
        std::vector< OpenSoT::AffineHelper > _wrenches_d;
        std::vector< OpenSoT::AffineHelper > _wrenches_dot;
        std::vector< OpenSoT::AffineHelper > _wrenches_min_var;
        std::vector< OpenSoT::AffineHelper > _slack_var;
        
        OpenSoT::constraints::force::FrictionCone::Ptr _friction_cone;
        ForzaGiusta::Ptr _forza_giusta;
        OpenSoT::solvers::iHQP::Ptr _solver;
        OpenSoT::AutoStack::Ptr _autostack;
        OpenSoT::tasks::MinimizeVariable::Ptr min_variation;
        OpenSoT::tasks::GenericTask::Ptr min_variation_GT;
        OpenSoT::tasks::MinimizeVariable::Ptr min_slack, min_slack1, min_slack2, min_slack3, min_slack4;
        double start_time_, time_;
        
        
        Eigen::VectorXd _dx_value;
        Eigen::MatrixXd _JC;
        
        Eigen::VectorXd b; 
        bool _change_ref1, _change_ref2;
        

        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_dot_bounds_ub;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_dot_bounds_lb;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_d_bounds;
        
        std::vector<Eigen::VectorXd> Fc_old;
        std::vector<Eigen::VectorXd> dFc;
        std::vector<Eigen::VectorXd> wrenches_lasso; 
        Eigen::VectorXd err_lb, err_ub;  
            

    };
    
}




demo::ForceOptimization::ForceOptimization(XBot::ModelInterface::Ptr model, 
                                           std::vector< std::string > contact_links):
    _model(model),
    _contact_links(contact_links)
{  
    /* Do we want to consider contact torques? */
    const bool optimize_contact_torque = false;
    _change_ref1 = false;
    _change_ref2 = false;
    
    /* Define optimization vector by stacking all contact wrenches */
    OpenSoT::OptvarHelper::VariableVector vars;

    for(auto cl : _contact_links)
    {
        vars.emplace_back(cl+"_wrenches_p",6);
        vars.emplace_back(cl+"_wrenches_m",6);
        vars.emplace_back(cl+"_slack_ub",3);
        vars.emplace_back(cl+"_slack_lb",3);
    }
    
    
    Fc_old.resize(_contact_links.size());
    for(int i = 0; i < _contact_links.size(); i++)
    Fc_old[i].setZero(12);

    OpenSoT::OptvarHelper opt(vars);
       
    /* Wrench bounds */
    Eigen::VectorXd wrench_p_ub(6), wrench_p_lb(6), wrench_m_ub(6), wrench_m_lb(6);
    wrench_p_lb.setZero();wrench_m_lb.setZero();
    wrench_p_ub.setZero();wrench_p_ub[0] = 25; wrench_p_ub[1] = 25; wrench_p_ub[2] = 500;
    wrench_m_ub.setZero();wrench_m_ub[0] = 25; wrench_m_ub[1] = 25; 
    
    Eigen::VectorXd inf_b(6);
    inf_b.setOnes();inf_b*=1e30; 
    
    Eigen::VectorXd ub_slack(3), lb_slack(3);
    ub_slack.setOnes();ub_slack*=1e20;
    lb_slack.setZero();
               
    Eigen::VectorXd wrench_dot_ub(6), wrench_dot_lb(6);
    wrench_dot_ub.setOnes();wrench_dot_ub *=  1e6; 
    wrench_dot_lb.setOnes();wrench_dot_lb *= -1e6; 
    
    Eigen::VectorXd opt_ub(18), opt_lb(18);
    opt_ub << wrench_p_ub, wrench_m_ub, ub_slack, ub_slack;
    opt_lb << wrench_p_lb, wrench_m_lb, lb_slack, lb_slack;
        
    
    /* Define affine mappings for all wrenches */
    for(auto cl : _contact_links){
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  
        _wrenches_d.emplace_back(opt.getVariable(cl+"_wrenches_p")/opt.getVariable(cl+"_wrenches_m")/opt.getVariable(cl+"_slack_ub")/opt.getVariable(cl+"_slack_lb"));
        
               
        
        wrench_d_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_bound",
                                                                                             _wrenches_d.back(),
                                                                                              opt_ub,
                                                                                              opt_lb,
                                                                                              OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
       
                              );    
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* forza giusta */        
      
        OpenSoT::AffineHelper wrench_var = opt.getVariable(cl+"_wrenches_p") + (-1.0) * opt.getVariable(cl+"_wrenches_m");
         _wrenches_dot.emplace_back(wrench_var);
        
        
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* wrench_dot upper bound with slack variables */ 
       
        Eigen::MatrixXd tmp(6,3);
        tmp.setZero();
        tmp.block<3,3>(0,0).setIdentity();
        
        
        OpenSoT::AffineHelper wrench_var_slack_ub = wrench_var + (-1.0) * tmp * opt.getVariable(cl+"_slack_ub");
        
         wrench_dot_bounds_ub.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_dot_ub",
                                                                                                 wrench_var_slack_ub,
                                                                                                 wrench_dot_ub,
                                                                                                 -inf_b,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
                               ); 
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        /* wrench_dot lower bound with slack variables */  

      
        OpenSoT::AffineHelper wrench_var_slack_lb =   wrench_var + tmp * opt.getVariable(cl+"_slack_lb");
        
         wrench_dot_bounds_lb.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_dot_lb",
                                                                                                 wrench_var_slack_lb,
                                                                                                 inf_b,
                                                                                                 wrench_dot_lb,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
         
                               );                 
       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
       /* min variation task */        
         
        _wrenches_min_var.emplace_back(opt.getVariable(cl+"_wrenches_p")/opt.getVariable(cl+"_wrenches_m"));
        
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* min slack variables */ 
         
        _slack_var.emplace_back(opt.getVariable(cl+"_slack_ub")/opt.getVariable(cl+"_slack_lb"));
                
    }  
    
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    /* Construct forza giusta task */
    _forza_giusta = boost::make_shared<ForzaGiusta>(_model,_wrenches_dot, _contact_links);

    
     /* Construct min slack */          
     std::vector<OpenSoT::AffineHelper>  slack;
     slack.emplace_back(_slack_var[0]/_slack_var[1]/_slack_var[2]/_slack_var[3]);
    
     min_slack = boost::make_shared<OpenSoT::tasks::MinimizeVariable>("MIN slack",slack[0]);
    
   
           
    /* Construct min variation task */    
   min_variation = boost::make_shared<OpenSoT::tasks::MinimizeVariable>("MIN_VARIATION", _wrenches_min_var[2]/_wrenches_min_var[3]);
     std::cout << "min_variation: \n" << _wrenches_min_var[2]/_wrenches_min_var[3] << std::endl;
   
   
//     min_variation = boost::make_shared<OpenSoT::tasks::MinimizeVariable>("MIN_VARIATION", _wrenches_d[3]);
//    min_variation = boost::make_shared<OpenSoT::tasks::MinimizeVariable>("MIN_VARIATION", _wrenches_dot[1]);
        
     
     min_variation_GT = boost::make_shared<OpenSoT::tasks::GenericTask>("MIN_VARIATION GT",_wrenches_min_var[2].getM(),_wrenches_min_var[2].getq());
//       std::cout << "min_variation_GT M: \n" << _wrenches_min_var[2].getM() << std::endl;
//       std::cout << "min_variation_GT q: \n" << _wrenches_min_var[2].getq() << std::endl;
     
    
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    

           _autostack = boost::make_shared<OpenSoT::AutoStack>( 1e-2 * min_slack + _forza_giusta +  1e-4 * min_variation);
//         _autostack = boost::make_shared<OpenSoT::AutoStack>( 1e-2 * min_slack + _forza_giusta + 1e-4 * min_variation_GT );
          
//          _autostack/min_variation;
          
          
   
   for(int i = 0; i < _contact_links.size(); i++)
   {      
    _autostack  <<  wrench_d_bounds[i] << wrench_dot_bounds_ub[i]<< wrench_dot_bounds_lb[i]; 
   }

    
    
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
    
//       _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 1.0);             
     _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 1e-10, OpenSoT::solvers::solver_back_ends::OSQP);
   
   
/*for(unsigned int i = 0; i < _autostack->getStack().size(); ++i)
{
   boost::any any_opt;
   _solver->getOptions(i, any_opt);
   qpOASES::Options opt_;
   opt_ = boost::any_cast<qpOASES::Options>(any_opt);
   opt_.setToDefault();
   opt_.printLevel = qpOASES::PL_NONE;
//     opt_.enableRegularisation = qpOASES::BT_TRUE;
   _solver->setOptions(i, opt_);
   
   //_solver->getOptions(i, any_opt);
   //opt_ = boost::any_cast<qpOASES::Options>(any_opt);
   //opt_.print();
}   
  */ 
   
    
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
    
    wrenches_lasso.resize(_contact_links.size());

    std::vector<Eigen::VectorXd> opt_c;
    opt_c.resize(_contact_links.size());
    
    start_time_= _start_time;
    time_= time;
    
    double period_= period;
    
    _forza_giusta->setFixedBaseTorque(fixed_base_torque, period_);  
    
    
    _autostack->update(Eigen::VectorXd());
    
    
    
    if(!_solver->solve(_dx_value))
    {
        std::cout<<"FORCE OPTIMIZATION CANNOT SOLVE!"<<std::endl;
        return false;
    }
    
      
    
    tau = fixed_base_torque;  
    
    for(int i = 0; i < _contact_links.size(); i++)
    {
         _wrenches_d[i].getValue(_dx_value, opt_c[i]);
         
        dFc[i]=opt_c[i].head(12);
                        
        Fc[i] =   Fc_old[i] + dFc[i] * period;
        Fc_old[i] = Fc[i];
        
//          std::cout<<"opt_c"<<opt_c[i].head(12)<<std::endl;
//          std::cout<<"slack: \n"<<opt_c[i].tail(6)<<std::endl;
         std::cout<<"Fc"<< i+1 << ": \n" << (Fc[i].head(6)-Fc[i].tail<6>()).head(3) <<std::endl;
               
        
        _model->getJacobian(_contact_links[i], _JC);
        
         tau.noalias() -= _JC.transpose()*(Fc[i].head(6)-Fc[i].tail(6));
                
         wrenches_lasso[i] = Fc[i].head(6)-Fc[i].tail<6>();
         
         
          /* Wrench bounds */
         Eigen::VectorXd wrench_p_ub(6), wrench_p_lb(6), wrench_m_ub(6), wrench_m_lb(6);
         wrench_p_lb.setZero();wrench_m_lb.setZero();
         wrench_p_ub.setZero(); wrench_p_ub[0] = 25; wrench_p_ub[1] = 25; wrench_p_ub[2] = 500;
         wrench_m_ub.setZero(); wrench_m_ub[0] = 25; wrench_m_ub[1] = 25; 
        
         Eigen::VectorXd ub_slack(3), lb_slack(3);
         ub_slack.setOnes();ub_slack*=1e20*period_;
         lb_slack.setZero();
                      
         Eigen::VectorXd opt_ub(18), opt_lb(18);
         opt_ub << wrench_p_ub, wrench_m_ub, ub_slack, ub_slack;
         opt_lb << wrench_p_lb, wrench_m_lb, lb_slack, lb_slack;
                       
         Eigen::VectorXd tmp_F(18);
         tmp_F.setZero();
         tmp_F.head(12)=Fc[i];
          
          wrench_d_bounds[i]->setBounds((1.0/period_)*(opt_ub-tmp_F),(1.0/period_)*(opt_lb-tmp_F));

                     
    }
  
    
    _forza_giusta->setWrenches(wrenches_lasso);
     
    Eigen::VectorXd tmpt(24);
    tmpt.setZero();
    tmpt.head(12)=Fc[2];
    min_variation->setReference(- (1.0/period_) * tmpt);
    
    
//      min_variation->setReference(- (1.0/period_) * Fc[2]);
//      min_variation_GT->setb( (1.0/period_) * Fc[2]);
//      min_variation->setReference(- (1.0/period_)*Fc[3]);
     

     
     if(time_-start_time_ >= 0.5)
     {
   
     std::cout<<"WRENCH_DOT BOUNDS UPDATE"<<std::endl;    
         
     Eigen::VectorXd wrench_dot_ub(6), wrench_dot_lb(6);
     wrench_dot_ub.setOnes();wrench_dot_ub *=  500; 
     wrench_dot_lb.setOnes();wrench_dot_lb *= -500;  
     Eigen::VectorXd inf_b(6);
     inf_b.setOnes();inf_b*=1e30; 
     
     for(int i = 0; i < _contact_links.size(); i++){             
     wrench_dot_bounds_ub[i]->setBounds(wrench_dot_ub, -inf_b);
     wrench_dot_bounds_lb[i]->setBounds(inf_b, wrench_dot_lb);
     }
     
     Eigen::VectorXd tmpt(24);
     tmpt.setZero();
     tmpt.tail(12)=Fc[3];
     min_variation->setReference(- (1.0/period_) * tmpt);
     
//      min_variation->setReference(- (1.0/period_)*Fc[3]);

          
     }


    
    return true;
    
}

void demo::ForceOptimization::log(XBot::MatLogger::Ptr logger)
{
    _autostack->log(logger);
    _solver->log(logger);
}


#endif
