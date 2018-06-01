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
//     _x_ref_fb.setZero(1);
    
    _wrenches_lasso.resize(_contact_links.size());
    for(int i = 0; i < _contact_links.size(); i++)
    _wrenches_lasso[i].setZero(6);
//     _wrenches_lasso[i].setZero(1);

    
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
     
//       _x_ref_fb << fixed_base_torque[0], fixed_base_torque[1], fixed_base_torque[2], 0, 0, 0;  
//       _x_ref_fb << 0, 0, 100, 0, 0, 0;
//       _x_ref_fb <<100;
          
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
//      _task.setZero(_wrenches_dot[0].getInputSize(), 1);
    
    for(int i = 0; i < _enabled_contacts.size(); i++)
    {
        if(!_enabled_contacts[i]){
            continue;
        }
        else {
            
             _model->getJacobian(_contact_links[i], _J_i);
            
              _Jfb_i=_J_i.block<6,6>(0,0).transpose(); 
              
              _task = _task + (period_ * _Jfb_i *_wrenches_dot[i]) + (_Jfb_i * _wrenches_lasso[i]);
                                                 
//            _task = _task + (period_ *_wrenches_dot[i]) + ( _wrenches_lasso[i]);
              
         
              
             /* Weight matrix */                                  
//              Eigen::MatrixXd P_i(6,6); 
//              P_i.setIdentity();
//              P_i.block<3,3>(3,3).setZero();
//              
//               _task = _task + (_Jfb_i + (i+1)*P_i) * (period_ *  _wrenches_dot[i]) + _Jfb_i * _wrenches_lasso[i];
                          
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
        
        OpenSoT::constraints::force::FrictionCone::Ptr _friction_cone;
        ForzaGiusta::Ptr _forza_giusta;
        OpenSoT::solvers::iHQP::Ptr _solver;
        OpenSoT::AutoStack::Ptr _autostack;
        OpenSoT::tasks::GenericTask::Ptr min_variation;
        OpenSoT::tasks::GenericTask::Ptr wrenches_dot_task;
        double start_time_, time_;
        
        
        Eigen::VectorXd _dx_value;
        Eigen::MatrixXd _JC;
        
        Eigen::VectorXd b; 
        bool _change_ref1, _change_ref2;
        

        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_dot_bounds;
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

    for(auto cl : _contact_links){
        vars.emplace_back(cl,12);
//         vars.emplace_back(cl,2);
    }
    
    
    Fc_old.resize(_contact_links.size());
    for(int i = 0; i < _contact_links.size(); i++)
    Fc_old[i].setZero(12);
    
    
//     for(int i = 0; i < _contact_links.size(); i++)
//     Fc_old[i].setZero(2);

    OpenSoT::OptvarHelper opt(vars);
       
    
    /* Wrench bounds */
    Eigen::VectorXd wrench_d_ub(12), wrench_d_lb(12);
    wrench_d_ub << 25, 25, 500, 0,  0, 0,  25, 25, 0.0, 0,  0, 0;
    wrench_d_lb <<  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    
    
    Eigen::VectorXd wrench_dot_ub(6), wrench_dot_lb(6);
    wrench_dot_ub <<  1e6, 1e6,  1e6, 1e6,  1e6,  1e6;
    wrench_dot_lb << -1e6, -1e6,  -1e6, -1e6,  -1e6,  -1e6;
    
       

//      /* wrench_d bounds */
//     Eigen::VectorXd wrench_d_ub(2), wrench_d_lb(2);
//     wrench_d_ub << 60, 0;
//     wrench_d_lb <<  0, 0;

//      /* wrench_dot bounds */
//     Eigen::VectorXd wrench_dot_ub(1), wrench_dot_lb(1);    
//     wrench_dot_ub <<  50;
//     wrench_dot_lb << -50;
   

    /* Define affine mappings for all wrenches */
    for(auto cl : _contact_links){

        
        _wrenches_d.emplace_back(opt.getVariable(cl));
        
               
        
        wrench_d_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_bound",
                                                                                             _wrenches_d.back(),
                                                                                             (1.0/.001)*wrench_d_ub,
                                                                                             (1.0/.001)*wrench_d_lb,
                                                                                             OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
        
                               );
        
        
        
        Eigen::MatrixXd tmpM(6,12);
        tmpM<<Eigen::MatrixXd::Identity(6,6),-Eigen::MatrixXd::Identity(6,6);

        
//         Eigen::MatrixXd tmpM(1,2);
//         tmpM<<Eigen::MatrixXd::Identity(1,1),-Eigen::MatrixXd::Identity(1,1);
        
        OpenSoT::AffineHelper wrench_var = tmpM*opt.getVariable(cl);
        _wrenches_dot.emplace_back(wrench_var);
        
//         std::cout << "opt var " << cl << ": \n" << opt.getVariable(cl) << std::endl;
//         std::cout << "wrench var " << cl << ": \n" << wrench_var << std::endl;
        
        
        
         wrench_dot_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_dot_bound",
                                                                                                 wrench_var,
                                                                                                 wrench_dot_ub,
                                                                                                 wrench_dot_lb,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
                               );
        
    
       
        
    }
    

    
    
    
    /* Construct forza giusta task */
    _forza_giusta = boost::make_shared<ForzaGiusta>(_model,_wrenches_dot, _contact_links);

           
    /* Construct min variation task */
    Eigen::MatrixXd A(12,48); 
    A.setZero();
    A.block<12,12>(0,32).setIdentity();
    A*=.001;
    
    b.setZero(12);

           
    min_variation = boost::make_shared<OpenSoT::tasks::GenericTask>("MIN_VARIATION",A,b);
    
//        _autostack = boost::make_shared<OpenSoT::AutoStack>(_forza_giusta);

//      _autostack = boost::make_shared<OpenSoT::AutoStack>(_forza_giusta + min_variation);  
     
     _autostack = (_forza_giusta/min_variation);
         
                   
   
   for(int i = 0; i < _contact_links.size(); i++)
    {    
         
         _autostack  <<  wrench_d_bounds[i] << wrench_dot_bounds[i];
       
    }
    
    
//     /* Define friction cones */
//     OpenSoT::constraints::force::FrictionCone::friction_cones friction_cones;
// 
//     for(auto cl : _contact_links)
//     {
//         friction_cones.emplace_back(cl, 0.3);
//     }
//     
//     auto friction_constraint = boost::make_shared<OpenSoT::constraints::force::FrictionCone>(_wrenches_dot, *_model, friction_cones);
//     
//     _autostack << friction_constraint;
        
       
    
    
        _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 0.0);             
//         _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 1e-10, OpenSoT::solvers::solver_back_ends::OSQP);
   
   
for(unsigned int i = 0; i < _autostack->getStack().size(); ++i)
{
   boost::any any_opt;
   _solver->getOptions(i, any_opt);
   qpOASES::Options opt_;
   opt_ = boost::any_cast<qpOASES::Options>(any_opt);
   opt_.setToDefault();
   opt_.printLevel = qpOASES::PL_NONE;
//    opt_.enableRegularisation = qpOASES::BT_TRUE;
   _solver->setOptions(i, opt_);
   
   //_solver->getOptions(i, any_opt);
   //opt_ = boost::any_cast<qpOASES::Options>(any_opt);
   //opt_.print();
}   
   
   
    
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

    
    
    start_time_= _start_time;
    time_= time;
    
    double period_= period;
    
    _forza_giusta->setFixedBaseTorque(fixed_base_torque, period_);  
    
    
    _autostack->update(Eigen::VectorXd());
    
    
    
    if(!_solver->solve(_dx_value))
    {
        std::cout<<"FORCE OPTIMIZATION COULD NOT SOLVE!"<<std::endl;
        return false;
    }
    
      
    
    tau = fixed_base_torque;  
    
    for(int i = 0; i < _contact_links.size(); i++)
    {
         _wrenches_d[i].getValue(_dx_value, dFc[i]);
                        
        Fc[i] =   Fc_old[i] + dFc[i] * period;
        Fc_old[i] = Fc[i];
        
//         std::cout<<"Fc: \n"<<Fc[i]<<std::endl;
               
        
        _model->getJacobian(_contact_links[i], _JC);
        
         tau.noalias() -= _JC.transpose()*(Fc[i].head(6)-Fc[i].tail(6));
        
        
        
           wrenches_lasso[i] = Fc[i].head(6)-Fc[i].tail<6>();
//            wrenches_lasso[i] = Fc[i].head(1)-Fc[i].tail<1>();
         

                              
          /* wrench_dot bounds */ 
         Eigen::VectorXd wrench_dot_ub(6), wrench_dot_lb(6);
         wrench_dot_ub <<  1e6,  1e6,  1e6,  1e6,  1e6,  1e6;
         wrench_dot_lb << -1e6, -1e6, -1e6, -1e6, -1e6, -1e6;
         
            
//       Eigen::VectorXd wrench_dot_ub(1), wrench_dot_lb(1);        
//       wrench_dot_ub <<  50;
//       wrench_dot_lb << -50;
        
         wrench_dot_bounds[i]->setBounds(wrench_dot_ub, wrench_dot_lb);
         
         /* wrench_d bounds */  
        Eigen::VectorXd wrench_d_ub(12), wrench_d_lb(12);
        wrench_d_ub << 25, 25, 500, 0,  0, 0,  25, 25, 0.0, 0,  0, 0;
        wrench_d_lb <<  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        

//     Eigen::VectorXd wrench_d_ub(2), wrench_d_lb(2);
//     wrench_d_ub << 60, 0;
//     wrench_d_lb <<  0, 0; 
         

                
        wrench_d_bounds[i]->setBounds((1.0/period_)*(wrench_d_ub-Fc[i]),(1.0/period_)*(wrench_d_lb-Fc[i]));

                     
    }
  
    
    _forza_giusta->setWrenches(wrenches_lasso);
       
    
    b.setZero(12);

    min_variation->setb(b - Fc[3]);
      

    
    if(time_-start_time_ >= 1.0 && !_change_ref1)
    {       
    std::cout<<"change leg"<<std::endl;
     _change_ref1 = true;
    }
   
    if(time_-start_time_ >= 1.0)
    {
    Eigen::MatrixXd A(12,48); 
    A.setZero();
    A.block<12,12>(0,24).setIdentity();
//     A.block<12,12>(0,0).setIdentity();
    A*=.001;    
              
    min_variation->setb(b - Fc[2]);
//     min_variation->setb(b - Fc[0]);
    min_variation->setA(A);    
    }
    
    
   if(time_-start_time_ >= 1.0 && time_-start_time_ <= 2.0)
   {  
    Eigen::VectorXd wrench_dot_ub(6), wrench_dot_lb(6);
//     wrench_dot_ub <<   4e3,  4e3,  4e3,  4e3,  4e3,  4e3;
//     wrench_dot_lb <<  -4e3, -4e3, -4e3, -4e3, -4e3, -4e3;
    
    wrench_dot_ub <<   500,  500,  500,  4e3,  4e3,  4e3;
    wrench_dot_lb <<  -500, -500, -500, -4e3, -4e3, -4e3;
    
    for(int i = 0; i < _contact_links.size(); i++)             
    wrench_dot_bounds[i]->setBounds(wrench_dot_ub, wrench_dot_lb); 
   }
    
    

    
    return true;
    
}

void demo::ForceOptimization::log(XBot::MatLogger::Ptr logger)
{
    _autostack->log(logger);
    _solver->log(logger);
}


#endif
