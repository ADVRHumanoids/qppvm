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
                    const std::vector<OpenSoT::AffineHelper>& wrenches,
                    std::vector<std::string> contact_links);
        
        void setFixedBaseTorque(const Eigen::VectorXd& fixed_base_torque);
        
    private:
        
        virtual void _update(const Eigen::VectorXd& x);

        virtual void _log(XBot::MatLogger::Ptr logger);
        
        XBot::ModelInterface::Ptr _model;
        std::vector<std::string> _contact_links;
        std::vector<bool> _enabled_contacts;
        std::vector<OpenSoT::AffineHelper> _wrenches;
        OpenSoT::AffineHelper _task;
        Eigen::MatrixXd _J_i;
        Eigen::MatrixXd _Jfb_i;
        Eigen::Vector6d _x_ref_fb;
        
        
    };
    
ForzaGiusta::ForzaGiusta(XBot::ModelInterface::Ptr model, 
                         const std::vector< OpenSoT::AffineHelper >& wrenches, 
                         std::vector< std::string > contact_links): 
    Task< Eigen::MatrixXd, Eigen::VectorXd >("FORZA_GIUSTA", wrenches[0].getInputSize()),
    _model(model),
    _contact_links(contact_links),
    _wrenches(wrenches),
    _enabled_contacts(wrenches.size(), true)
{
    _x_ref_fb.setZero();
    
    update(Eigen::VectorXd());
}

void ForzaGiusta::setFixedBaseTorque(const Eigen::VectorXd& fixed_base_torque)
{
    if(fixed_base_torque.size() != _model->getJointNum())
    {
        throw std::invalid_argument("fixed_base_torque != _model->getJointNum()");
    }
    
     _x_ref_fb = fixed_base_torque.head<6>();
    
    
}


void ForzaGiusta::_update(const Eigen::VectorXd& x)
{
     _task.setZero(_wrenches[0].getInputSize(), 6);
    
    for(int i = 0; i < _enabled_contacts.size(); i++)
    {
        if(!_enabled_contacts[i]){
            continue;
        }
        else {
            _model->getJacobian(_contact_links[i], _J_i);
             _Jfb_i=_J_i.block<6,6>(0,0).transpose(); 
//              _Jfb_i.setIdentity(6,6); 
                  
            Eigen::MatrixXd _Jfb_L1_i(6,12);
            _Jfb_L1_i <<  _Jfb_i, - _Jfb_i;
            
            _task = _task + _Jfb_L1_i * _wrenches[i];
        }
    }
    
    _A = _task.getM();
    _b = -_task.getq() + _x_ref_fb;
}

void ForzaGiusta::_log(XBot::MatLogger::Ptr logger)
{
    OpenSoT::Task< Eigen::MatrixXd, Eigen::VectorXd >::_log(logger);
}



class ForzaGiusta_Lasso : public OpenSoT::Task<Eigen::MatrixXd, Eigen::VectorXd>
    {
        
    public:
        
        typedef boost::shared_ptr<ForzaGiusta_Lasso> Ptr;
        
        ForzaGiusta_Lasso(XBot::ModelInterface::Ptr model, 
                    const std::vector<OpenSoT::AffineHelper>& wrenches,
                    std::vector<std::string> contact_links);
        
        void setFixedBaseTorque(const Eigen::VectorXd& fixed_base_torque);
        
    private:
        
        virtual void _update(const Eigen::VectorXd& x);

        virtual void _log(XBot::MatLogger::Ptr logger);
        
        XBot::ModelInterface::Ptr _model;
        std::vector<std::string> _contact_links;
        std::vector<bool> _enabled_contacts;
        std::vector<OpenSoT::AffineHelper> _wrenches;
        OpenSoT::AffineHelper _task;
        Eigen::MatrixXd _J_i;
        Eigen::MatrixXd _Jfb_i;
        Eigen::Vector6d _x_ref_fb;
        
        
    };
    
ForzaGiusta_Lasso::ForzaGiusta_Lasso(XBot::ModelInterface::Ptr model, 
                         const std::vector< OpenSoT::AffineHelper >& wrenches, 
                         std::vector< std::string > contact_links): 
    Task< Eigen::MatrixXd, Eigen::VectorXd >("FORZA_GIUSTA_Lasso", wrenches[0].getInputSize()),
    _model(model),
    _contact_links(contact_links),
    _wrenches(wrenches),
    _enabled_contacts(wrenches.size(), true)
{
   
    update(Eigen::VectorXd());
}



void ForzaGiusta_Lasso::_update(const Eigen::VectorXd& x)
{
  
     _task.setZero(_wrenches[0].getInputSize(), 1);
    
    for(int i = 0; i < _enabled_contacts.size(); i++)
    {
        if(!_enabled_contacts[i]){
            continue;
        }
        else {
                  
           Eigen::MatrixXd A(1,12);
           A.setZero();
           A.topRows(1)<< 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
                                              
            _task = _task + A* _wrenches[i];
        }
    }
    
    
    
    _A = _task.getM();
    _b = -_task.getq();
}

void ForzaGiusta_Lasso::_log(XBot::MatLogger::Ptr logger)
{
    OpenSoT::Task< Eigen::MatrixXd, Eigen::VectorXd >::_log(logger);
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
        std::vector< OpenSoT::AffineHelper > _wrenches;
        std::vector< OpenSoT::AffineHelper > _wrenches_dot;
        
        OpenSoT::constraints::force::FrictionCone::Ptr _friction_cone;
        ForzaGiusta::Ptr _forza_giusta;
        ForzaGiusta_Lasso::Ptr _forza_giusta_lasso;
        OpenSoT::solvers::iHQP::Ptr _solver;
        OpenSoT::AutoStack::Ptr _autostack;
        OpenSoT::tasks::GenericTask::Ptr min_variation;
        double start_time_, time_;
        
        
        Eigen::VectorXd _x_value;
        Eigen::MatrixXd _JC;
        
        Eigen::VectorXd b; 
        bool _change_ref;
        

        double period_;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_dot_bounds;

        
    };
    
}




demo::ForceOptimization::ForceOptimization(XBot::ModelInterface::Ptr model, 
                                           std::vector< std::string > contact_links):
    _model(model),
    _contact_links(contact_links)
{ 
    /* Do we want to consider contact torques? */
    const bool optimize_contact_torque = false;
    _change_ref = false;
    
    /* Define optimization vector by stacking all contact wrenches */
    OpenSoT::OptvarHelper::VariableVector vars;

    for(auto cl : _contact_links){
        vars.emplace_back(cl,12);
    }
    
    period_=.001;

    OpenSoT::OptvarHelper opt(vars);
       
    
    /* Wrench bounds */
    Eigen::VectorXd wrench_ub(12), wrench_lb(12);
    wrench_ub << 1000,  1000, 1000, 0,  0, 0,  1000,  1000, 0, 0,  0, 0;
    wrench_lb <<  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    
    Eigen::VectorXd wrench_dot_ub(6), wrench_dot_lb(6);
    wrench_dot_ub << 1e6,  1e6, 1e6, 1e6,  1e6, 1e6;
    wrench_dot_lb << -1e6,  -1e6, -1e6, -1e6,  -1e6, -1e6; 
       

    std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_bounds;

     

    /* Define affine mappings for all wrenches */
    for(auto cl : _contact_links){

        
        _wrenches.emplace_back(opt.getVariable(cl));
        
               
        
        wrench_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_bound",
                                                                                             _wrenches.back(),
                                                                                             wrench_ub,
                                                                                             wrench_lb,
                                                                                             OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
        
                               );
        
        
       //_wrenches_dot.emplace_back(opt.getVariable(cl).head(6) + (-1) * opt.getVariable(cl).tail(6));
        
        Eigen::MatrixXd tmpM(6,12);
        tmpM<<Eigen::MatrixXd::Identity(6,6),-Eigen::MatrixXd::Identity(6,6);
        OpenSoT::AffineHelper wrench_var = tmpM*opt.getVariable(cl);
        _wrenches_dot.emplace_back(wrench_var);
        
//         _wrenches_dot.emplace_back(tmpM*opt.getVariable(cl));        

        
        
        wrench_dot_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_dot_bound",
                                                                                                _wrenches_dot.back(),
                                                                                                 wrench_dot_ub,
                                                                                                 wrench_dot_lb,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
                               );
        
    
       
        
    }
    
    
    for(int i = 0; i < _contact_links.size(); i++)
    {      
         std::cout<<"UPPER BOUND: \n"<<wrench_dot_bounds[i]->getbUpperBound()<<std::endl;
         std::cout<<"LOWER BOUND: \n"<<wrench_dot_bounds[i]->getbLowerBound()<<std::endl;
         std::cout<<"A: \n"<<wrench_dot_bounds[i]->getAineq()<<std::endl;
         
         
    }
         
    
    /* Define friction cones */
//     OpenSoT::constraints::force::FrictionCone::friction_cones friction_cones;
// 
//     for(auto cl : _contact_links)
//     {
//         friction_cones.emplace_back(cl, 0.3);
//     }
//     
//     auto friction_constraint = boost::make_shared<OpenSoT::constraints::force::FrictionCone>(_wrenches, *_model, friction_cones);
    
    
    
    /* Construct forza giusta task */
    _forza_giusta = boost::make_shared<ForzaGiusta>(_model, _wrenches, _contact_links);
//      for(unsigned int i = 0; i < 4; ++i)
//         _forza_giusta<< wrench_dot_bounds[i];
    
    
    
    /* Construct min variation task */
    int size = 48;
    Eigen::MatrixXd A(size,size); A.setIdentity();
    //Eigen::VectorXd b(size); 
    
    b.setZero(48);
    
    b[2] = 1000;//450.59;
    b[14] = 1000.;//74.85;
    b[26] = 1000.;//384.47;
    b[38] = 1000.;
    
     min_variation = boost::make_shared<OpenSoT::tasks::GenericTask>("MIN_VARIATION",A,b);
//      for(unsigned int i = 0; i < 4; ++i)
//         min_variation<< wrench_bounds[i];
        
      
        _autostack = (_forza_giusta/min_variation);
     
//          _autostack = boost::make_shared<OpenSoT::AutoStack>(_forza_giusta);      
       

    /* Construct forza giusta Lasso task - !!!UNNECESSARY!!! */
//       _forza_giusta_lasso = boost::make_shared<ForzaGiusta_Lasso>(_model, _wrenches, _contact_links);
// 
//       _autostack = boost::make_shared<OpenSoT::AutoStack>(_forza_giusta);       
//       Eigen::VectorXd err_lb(1), err_ub(1);                  
//       err_lb << 0;
//       err_ub << 2000;
//       _autostack << boost::make_shared<OpenSoT::constraints::TaskToConstraint>(_forza_giusta_lasso,err_lb,err_ub);
      
                    
   
   for(int i = 0; i < _contact_links.size(); i++)
    {    
         
         _autostack  <<  wrench_bounds[i]<<wrench_dot_bounds[i];

       
    }
    
    
       
    _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 1.0);             
//     _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 0.0, OpenSoT::solvers::solver_back_ends::OSQP);
   
   
for(unsigned int i = 0; i < _autostack->getStack().size(); ++i)
{
   boost::any any_opt;
   _solver->getOptions(i, any_opt);
   qpOASES::Options opt_;
   opt_ = boost::any_cast<qpOASES::Options>(any_opt);
   opt_.setToDefault();
   opt_.printLevel = qpOASES::PL_NONE;
   //opt_.enableRegularisation = qpOASES::BT_TRUE;
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
    
    start_time_= _start_time;
    time_= time;
    period_ = period; 
    
    _forza_giusta->setFixedBaseTorque(fixed_base_torque);
       
    
    _autostack->update(Eigen::VectorXd());
    
    
    
    if(!_solver->solve(_x_value))
    {
        std::cout<<"FORCE OPTIMIZATION COULD NOT SOLVE!"<<std::endl;
        return false;
    }
    
      
    
    tau = fixed_base_torque;  
    
    for(int i = 0; i < _contact_links.size(); i++)
    {
        _wrenches[i].getValue(_x_value, Fc[i]);
                       
        _model->getJacobian(_contact_links[i], _JC);
        
         tau.noalias() -= _JC.transpose()*(Fc[i].head(6)-Fc[i].tail(6));
         
         Eigen::VectorXd wrench_dot_ub(6), wrench_dot_lb(6);
         wrench_dot_ub << 1e6,  1e6, 1e6, 1e6,  1e6, 1e6;
         wrench_dot_lb << -1e6,  -1e6, -1e6, -1e6,  -1e6, -1e6;
         
//         std::cout<<"Fc: \n"<<Fc[i]<<std::endl;
         
         
//         std::cout<<"UPPER BOUND: \n"<<wrench_dot_bounds[i]->getbUpperBound()<<std::endl;
//         std::cout<<"LOWER BOUND: \n"<<wrench_dot_bounds[i]->getbLowerBound()<<std::endl;
//          std::cout<<"A: \n"<<wrench_dot_bounds[i]->getAineq()<<std::endl;
         
         
          wrench_dot_bounds[i]->setBounds(wrench_dot_ub * period_ + (Fc[i].head(6)-Fc[i].tail(6)), wrench_dot_lb * period_  + (Fc[i].head(6)-Fc[i].tail(6)));   
         
        


    }
    
    if(time_-start_time_ >= 2.0)
    {
    Eigen::VectorXd wrench_dot_ub(6), wrench_dot_lb(6);
    wrench_dot_ub << 1e6,  1e6, 100, 1e6,  1e6, 1e6;
    wrench_dot_lb << -1e6,  -1e6, -100, -1e6,  -1e6, -1e6;
    for(int i = 0; i < _contact_links.size(); i++)
    wrench_dot_bounds[i]->setBounds(wrench_dot_ub * period_ + (Fc[i].head(6)-Fc[i].tail(6)), wrench_dot_lb * period_  + (Fc[i].head(6)-Fc[i].tail(6))); 
    }
    
    
    
    min_variation->setb(_x_value);
   
    if(time_-start_time_ >= 2.0 && !_change_ref)
    {
    b.setZero(48);
    b[2] = 1000;//450.59;
    b[14] = 1000.;//74.85;
    b[26] = 1000.;//384.47;
    b[38] = 0.;
    
    Eigen::VectorXd wrench_dot_ub(6), wrench_dot_lb(6);
    wrench_dot_ub << 1e6,  1e6, 100, 1e6,  1e6, 1e6;
    wrench_dot_lb << -1e6,  -1e6, -100, -1e6,  -1e6, -1e6;
    for(int i = 0; i < _contact_links.size(); i++)
    wrench_dot_bounds[i]->setBounds(wrench_dot_ub * period_ + (Fc[i].head(6)-Fc[i].tail(6)), wrench_dot_lb * period_  + (Fc[i].head(6)-Fc[i].tail(6))); 
    
    min_variation->setb(b);
    
    _change_ref = true;
    }
    
    
    
    return true;
    
}

void demo::ForceOptimization::log(XBot::MatLogger::Ptr logger)
{
//     _autostack->log(logger);
//     _solver->log(logger);
}


#endif
