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
                    const std::vector<OpenSoT::AffineHelper>& forces_dot,
                    std::vector<std::string> contact_links);
        
        void setFixedBaseTorque(const Eigen::VectorXd& fixed_base_torque, double& period);
        
        void setforces(const std::vector<Eigen::VectorXd>& forces);
        
    private:
        
        virtual void _update(const Eigen::VectorXd& x);

        virtual void _log(XBot::MatLogger::Ptr logger);
        
        XBot::ModelInterface::Ptr _model;
        std::vector<std::string> _contact_links;
        std::vector<bool> _enabled_contacts;
        std::vector<OpenSoT::AffineHelper> _forces_dot;
        OpenSoT::AffineHelper _task;
        Eigen::MatrixXd _J_i;
        Eigen::MatrixXd _Jfb_i;
        Eigen::VectorXd _x_ref_fb;
        double period_;
        std::vector<Eigen::VectorXd> _forces;
        
        
    };
    
ForzaGiusta::ForzaGiusta(XBot::ModelInterface::Ptr model, 
                         const std::vector< OpenSoT::AffineHelper >& forces_dot, 
                         std::vector< std::string > contact_links): 
    Task< Eigen::MatrixXd, Eigen::VectorXd >("FORZA_GIUSTA", forces_dot[0].getInputSize()),
    _model(model),
    _contact_links(contact_links),
    _forces_dot(forces_dot),
    _enabled_contacts(forces_dot.size(), true)
{
     _x_ref_fb.setZero(6);
    
    _forces.resize(_contact_links.size());
    for(int i = 0; i < _contact_links.size(); i++)
    _forces[i].setZero(3);

    
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
     
          
     period_=period;
    
    
}

void ForzaGiusta::setforces(const std::vector<Eigen::VectorXd>& forces)
{

    for(int i = 0; i < _contact_links.size(); i++)
     _forces[i] = forces[i];
     
        
    
}


void ForzaGiusta::_update(const Eigen::VectorXd& x)
{
     _task.setZero(_forces_dot[0].getInputSize(), 6);
    
    for(int i = 0; i < _enabled_contacts.size(); i++)
    {
        if(!_enabled_contacts[i]){
            continue;
        }
        else {
            
             _model->getJacobian(_contact_links[i], _J_i);
            
              _Jfb_i=_J_i.block<6,6>(0,0).transpose().block<6,3>(0,0);
                           
              _task = _task + (period_ * _Jfb_i *_forces_dot[i]) + (_Jfb_i * _forces[i]);
                                                 
            
                     
        }
    }
    
  
    _A = _task.getM();
    _b = -_task.getq() + _x_ref_fb;
    

    
}

void ForzaGiusta::_log(XBot::MatLogger::Ptr logger)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


 class FrictionCone_forces_dot : public OpenSoT::Constraint<Eigen::MatrixXd, Eigen::VectorXd> 
    {
        
     public:
                typedef boost::shared_ptr<FrictionCone_forces_dot> Ptr;
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
        OpenSoT::AffineHelper _forces_dot;
        
        Eigen::VectorXd _forces;
        double period_;

    public:


        FrictionCone_forces_dot(const std::vector<OpenSoT::AffineHelper>& forces_dot,
                                  XBot::ModelInterface &robot,
                                  const friction_cones & mu);


        void update(const Eigen::VectorXd &x);

        void setMu(const friction_cones& mu){ _mu = mu;}
        
        void setforces(const std::vector<Eigen::VectorXd>& forces);

        int getNumberOfContacts(){return _n_of_contacts;}

    private:
        void computeAineq();
        void computeUpperBound();

            };
    
 

       FrictionCone_forces_dot::FrictionCone_forces_dot(const std::vector<OpenSoT::AffineHelper>& forces_dot,
                                                            XBot::ModelInterface &robot,
                                                            const friction_cones& mu):
           Constraint("friction_cone", forces_dot[0].getInputSize()),
           _robot(robot),
           _mu(mu),
           _Ci(5,3)
       {
           _n_of_contacts = _mu.size();

           _forces_dot.setZero(forces_dot[0].getInputSize(), 0);
           
           _forces.setZero(_n_of_contacts * 3);

           period_=0.001;

           for(auto& w : forces_dot){
               
               _forces_dot = _forces_dot / w;
           }
           

//            Eigen::Affine3d wTli;
//            for(unsigned int i = 0; i < _n_of_contacts; ++i)
//            {
//                _robot.getPose(_mu[i].first, wTli);
//                _wTl.push_back(wTli);
//            }

           update(Eigen::VectorXd::Zero(0));

       }

       void FrictionCone_forces_dot::computeAineq()
       {
           _Ci.setZero(5,3);

           _A.resize(5*_n_of_contacts, 3*_n_of_contacts);

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

                _A.block<5,3>(5*i, 3*i) = _Ci;
           }
           
           _A_d.resize(5*_n_of_contacts, 3*_n_of_contacts);

           _A_d.setZero(_A_d.rows(), _A_d.cols());
           
           _Cdi.setZero(5,3);

           for(unsigned int i = 0; i < _n_of_contacts; ++i)
           {

                _Cdi(0,0) = 1.; _Cdi(0,1) = 0.; _Cdi(0,2) = 0;
                _Cdi(1,0) = -1.; _Cdi(1,1) = 0.; _Cdi(1,2) = 0;
                _Cdi(2,0) = 0.; _Cdi(2,1) = 1.; _Cdi(2,2) = 0;
                _Cdi(3,0) = 0.; _Cdi(3,1) = -1.; _Cdi(3,2) = 0;
                _Cdi(4,0) = 0.; _Cdi(4,1) = 0.; _Cdi(4,2) = -1.;

                _A_d.block<5,3>(5*i, 3*i) = _Cdi;
           }
           

       }

       void FrictionCone_forces_dot::computeUpperBound()
       {
           _b.resize(5*_n_of_contacts);
           _b.setZero(_b.rows());
           _bLowerBound.resize(5*_n_of_contacts);
           _bLowerBound = -1.0e20*_bLowerBound.setOnes(_bLowerBound.size());
       }

       void FrictionCone_forces_dot::setforces(const std::vector<Eigen::VectorXd>& forces)
       {
           
           _forces << forces[0],forces[1],forces[2],forces[3];
         
 
       }

       void FrictionCone_forces_dot::update(const Eigen::VectorXd &x)
       {

               _n_of_contacts = _mu.size();

               computeAineq();
               computeUpperBound();

               _friction_cone =  (_A_d *  period_ * _forces_dot) + (_A *  _forces) - _b;   
               _Aineq = _friction_cone.getM();
               _bUpperBound = - _friction_cone.getq();
               
                std::cout<<"forces: \n"<<_forces<<std::endl;

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
        std::vector< OpenSoT::AffineHelper > _force_dot;
        std::vector< OpenSoT::AffineHelper > _wrenches_min_var;
        std::vector< OpenSoT::AffineHelper > _slack_var;
        std::vector< OpenSoT::AffineHelper > _t_var;
        
        OpenSoT::constraints::force::FrictionCone::Ptr _friction_cone;
        ForzaGiusta::Ptr _forza_giusta;
        FrictionCone_forces_dot::Ptr _friction_constraint;
        OpenSoT::solvers::iHQP::Ptr _solver;
        OpenSoT::AutoStack::Ptr _autostack;
        OpenSoT::tasks::MinimizeVariable::Ptr min_variation;
        OpenSoT::tasks::MinimizeVariable::Ptr min_variation2;
        OpenSoT::tasks::MinimizeVariable::Ptr min_slack;
        OpenSoT::tasks::GenericTask::Ptr task_Lasso;
        double start_time_, time_;
        
        
        Eigen::VectorXd _dx_value;
        Eigen::MatrixXd _JC;

        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> force_dot_bounds_ub;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> force_dot_bounds_lb;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> _force_lasso_ub;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> _force_lasso_lb;
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> force_bounds;       
        std::vector<OpenSoT::constraints::GenericConstraint::Ptr> aux_var_bounds;
        
        std::vector<Eigen::VectorXd> Fc_old;
        std::vector<Eigen::VectorXd> dFc;
        std::vector<Eigen::VectorXd> Wc;
        
        Eigen::VectorXd force_ub, force_lb;
        Eigen::VectorXd inf_b;
        Eigen::VectorXd ub_slack, lb_slack, t_ub, t_lb;            
        Eigen::VectorXd force_dot_ub, force_dot_lb;
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
        vars.emplace_back(cl+"_force_dot",3);
        vars.emplace_back(cl+"_t",1);
        vars.emplace_back(cl+"_slack_ub",3);
        vars.emplace_back(cl+"_slack_lb",3);
    }
    
    
    Fc_old.resize(_contact_links.size());
    Wc.resize(_contact_links.size());
    for(int i = 0; i < _contact_links.size(); i++)
    {
    Fc_old[i].setZero(3);
    Wc[i].setZero(6);
    }

    OpenSoT::OptvarHelper opt(vars);
       
    /* Wrench bounds */
    force_ub.setOnes(3); force_ub *= 1e3;
    force_lb.setOnes(3); force_lb *=-1e3; 
    
    inf_b.setOnes(3); inf_b*=1e30; 
               
    force_dot_ub.setOnes(3); force_dot_ub *=  1e3; 
    force_dot_lb.setOnes(3); force_dot_lb *= -1e3; 
    
    Eigen::VectorXd _aux_var_ub, _aux_var_lb;
    _aux_var_ub.setOnes(7);_aux_var_ub*=1e30;
    _aux_var_lb.setZero(7);

      
    
    /* Define affine mappings for all forces */
    for(auto cl : _contact_links){
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        /* Auxiliary variables */
  
        _aux_var.emplace_back(opt.getVariable(cl+"_t")/opt.getVariable(cl+"_slack_ub")/opt.getVariable(cl+"_slack_lb"));
             
        
        aux_var_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"force_bound",
                                                                                              _aux_var.back(),
                                                                                              _aux_var_ub,
                                                                                              _aux_var_lb,
                                                                                               OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
       
                              );    
        
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* forza giusta */        
       
        OpenSoT::AffineHelper force_dot = opt.getVariable(cl+"_force_dot");
        _force_dot.emplace_back(force_dot);
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        /* Wrench bounds */
 
        force_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"force_bound",
                                                                                             _force_dot.back(),
                                                                                              force_ub,
                                                                                              force_lb,
                                                                                              OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
       
                              );   
        
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* Lasso wrench upper bound */ 
       
       Eigen::VectorXd zero,inf;
       zero.setZero(1); inf.setOnes(1);inf*=1e30;       
  
       OpenSoT::AffineHelper force_lasso_ub = 0.001 * force_dot.tail(1) - opt.getVariable(cl+"_t");
        
       _force_lasso_ub.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"force_dot_ub",
                                                                                               force_lasso_ub,
                                                                                               zero,
                                                                                               -inf,
                                                                                               OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
                               ); 
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
         /* Lasso wrench upper bound */   

      
        OpenSoT::AffineHelper force_lasso_lb =   0.001 * force_dot.tail(1) +  opt.getVariable(cl+"_t");
        
         _force_lasso_lb.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"force_dot_lb",
                                                                                                 force_lasso_lb,
                                                                                                 inf,
                                                                                                 zero,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
         
                               );                                
        
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* wrench_dot upper bound with slack variables */ 
       

        OpenSoT::AffineHelper force_dot_slack_ub = force_dot - opt.getVariable(cl+"_slack_ub");
        
         force_dot_bounds_ub.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"wrench_dot_ub",
                                                                                                 force_dot_slack_ub,
                                                                                                 force_dot_ub,
                                                                                                 -inf_b,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
                               ); 
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        /* wrench_dot lower bound with slack variables */  

      
        OpenSoT::AffineHelper force_dot_slack_lb = force_dot + opt.getVariable(cl+"_slack_lb");
        
         force_dot_bounds_lb.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"force_dot_lb",
                                                                                                 force_dot_slack_lb,
                                                                                                 inf_b,
                                                                                                 force_dot_lb,
                                                                                                 OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
         
                               );                                
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* min slack variables */ 
         
        _slack_var.emplace_back(opt.getVariable(cl+"_slack_ub")/opt.getVariable(cl+"_slack_lb"));
      
       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
       /* Lasso */ 

        _t_var.emplace_back(opt.getVariable(cl+"_t"));

 
    }  
    
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    /* Define friction cones */
    OpenSoT::constraints::force::FrictionCone::friction_cones friction_cones;

    for(auto cl : _contact_links)
    friction_cones.emplace_back(cl, .3);
         
    _friction_constraint = boost::make_shared<FrictionCone_forces_dot>(_force_dot, *_model, friction_cones);
    
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

     /* Construct forza giusta task */
     _forza_giusta = boost::make_shared<ForzaGiusta>(_model,_force_dot, _contact_links);
    
     /* Construct min slack */          
     OpenSoT::AffineHelper  slack;
     slack = _slack_var[0]/_slack_var[1]/_slack_var[2]/_slack_var[3];
    
     min_slack = boost::make_shared<OpenSoT::tasks::MinimizeVariable>("MIN slack",slack);

  
     /* Construct Lasso task */
     t_var = _t_var[0]/_t_var[1]/_t_var[2]/_t_var[3]; 

     int size = t_var.getOutputSize();
     Eigen::MatrixXd tmp_A(size,size); tmp_A.setZero(size,size);
     Eigen::VectorXd tmp_b(size); tmp_b.setZero(size);
  
     task_Lasso = boost::make_shared<OpenSoT::tasks::GenericTask>("Lasso task",tmp_A,tmp_b,t_var);
     Eigen::VectorXd c_Lasso(size);
     c_Lasso.setOnes(size);
//      c_Lasso << 1.0, 20.0, 30.0, 40.0; /* WEIGHTED LASSO */
     task_Lasso->setc(c_Lasso);
     task_Lasso->setHessianType( OpenSoT::HST_ZERO);
     
     std::cout << "A_task_lasso:\n" << task_Lasso->getA() << std::endl;
     std::cout << "b_task_lasso:\n" << task_Lasso->getb() << std::endl;
     std::cout << "c_task_lasso:\n" << task_Lasso->getc() << std::endl;
     std::cout << "HessianAtype:\n" << task_Lasso->getHessianAtype() << std::endl;
 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
//    /* LASSO cost function (1 level)*/
//     OpenSoT::tasks::Aggregated::Ptr aggr1 = (1e-2 * min_slack + _forza_giusta + task_Lasso); 
//     
//     _autostack = boost::make_shared<OpenSoT::AutoStack>( aggr1 );
//    
//     for(int i = 0; i < _contact_links.size(); i++)
//     {      
//      _autostack   << aux_var_bounds[i] <<  force_bounds[i] << force_dot_bounds_ub[i] << force_dot_bounds_lb[i] << _force_lasso_ub[i] << _force_lasso_lb[i];
//     }
//    
//     _autostack << _friction_constraint;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
   /* LASSO null space (2 levels)*/
     
   for(int i = 0; i < _contact_links.size(); i++)
   {      
   task_Lasso  << _force_lasso_ub[i] << _force_lasso_lb[i];
   } 
     
   OpenSoT::tasks::Aggregated::Ptr aggr2 = (1e-2 * min_slack + _forza_giusta); 
   _autostack = boost::make_shared<OpenSoT::AutoStack>( aggr2 );
//    _autostack = (aggr2/task_Lasso); 

   
   for(int i = 0; i < _contact_links.size(); i++)
   {      
    _autostack   << aux_var_bounds[i] <<  force_bounds[i] << force_dot_bounds_ub[i]<< force_dot_bounds_lb[i];
   }
   
    _autostack << _friction_constraint;
    
    
     /* FrontEnd (iHQP) with multiple solvers */
     OpenSoT::solvers::solver_back_ends solver_1 = OpenSoT::solvers::solver_back_ends::qpOASES;
     OpenSoT::solvers::solver_back_ends solver_2 = OpenSoT::solvers::solver_back_ends::OSQP;
     OpenSoT::solvers::solver_back_ends solver_3 = OpenSoT::solvers::solver_back_ends::CBC;
     
     std::vector<OpenSoT::solvers::solver_back_ends> solver_vector(2);
     solver_vector[0]=solver_1; solver_vector[1]=solver_2;
   
    _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 1.0, solver_vector);
//      _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 0.0, OpenSoT::solvers::solver_back_ends::OSQP);


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
 
    
    if(!_solver->solve(_dx_value))
    {
        std::cout<<"FORCE OPTIMIZATION CANNOT SOLVE!"<<std::endl;
        return false;
    }
    
/////////////////////////////////////////////////////////////////////////////////////     
/////////////////////////////////////////////////////////////////////////////////////  
    
    tau = fixed_base_torque;  
    
    for(int i = 0; i < _contact_links.size(); i++)
    {
         _force_dot[i].getValue(_dx_value, opt_c[i]);
         
        dFc[i]=opt_c[i];
                        
        Fc[i] =   Fc_old[i] + dFc[i] * period;
        Fc_old[i] = Fc[i];
        
         std::cout<<"Fc"<< i+1 << ": \n" << Fc[i] <<std::endl;
         
         Wc[i].setZero(6);
         Wc[i].head(3)=Fc[i];              
        
        _model->getJacobian(_contact_links[i], _JC);
        
         tau.noalias() -= _JC.transpose()*Wc[i];        
         
         /* Wrench bounds */
         force_bounds[i]->setBounds((1.0/period_)*(force_ub-Fc[i]),(1.0/period_)*(force_lb-Fc[i]));
 
         /* Lasso bounds */       
         Eigen::VectorXd zero(1),inf(1);
         zero.setZero(); inf.setOnes();inf*=1e30;   
         
         _force_lasso_ub[i]->setBounds(zero - Fc[i].tail(1), -inf - Fc[i].tail(1));
         _force_lasso_lb[i]->setBounds( inf - Fc[i].tail(1), zero - Fc[i].tail(1));
                     
    }
  
    
    _forza_giusta->setforces(Fc);
    _friction_constraint->setforces(Fc);
     

    
    return true;
    
}

void demo::ForceOptimization::log(XBot::MatLogger::Ptr logger)
{
    _autostack->log(logger);
    _solver->log(logger);
}


#endif
