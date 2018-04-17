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
#include <OpenSoT/tasks/force/FloatingBase.h>
#include <OpenSoT/constraints/force/CoP.h>

namespace demo {
    
    class ForceOptimization 
    {
      
    public:
        
        typedef boost::shared_ptr<ForceOptimization> Ptr;
        
        ForceOptimization(XBot::ModelInterface::Ptr model, 
                          std::vector<std::string> contact_links
                          );
        
        bool compute(const Eigen::VectorXd& fixed_base_torque, 
                     std::vector<Eigen::Vector6d>& Fc,
                     Eigen::VectorXd& tau
                    );
        
        void log(XBot::MatLogger::Ptr logger);
        
        
    private:
        
        XBot::ModelInterface::Ptr _model;
        std::vector<std::string> _contact_links;
        std::vector< OpenSoT::AffineHelper > _wrenches;
        
        OpenSoT::constraints::force::FrictionCone::Ptr _friction_cone;
        OpenSoT::tasks::force::FloatingBase::Ptr _forza_giusta;
        OpenSoT::solvers::iHQP::Ptr _solver;
        OpenSoT::AutoStack::Ptr _autostack;
        
        Eigen::VectorXd _x_value;
        Eigen::MatrixXd _JC;
        Eigen::VectorXd _fc_i;
        Eigen::VectorXd _Fci;
        
    };
    
}




demo::ForceOptimization::ForceOptimization(XBot::ModelInterface::Ptr model, 
                                           std::vector< std::string > contact_links):
    _model(model),
    _contact_links(contact_links)
{
    /* Do we want to consider contact torques? */
    const bool optimize_contact_torque = true;
    
    /* Define optimization vector by stacking all contact wrenches */
    OpenSoT::OptvarHelper::VariableVector vars;

    for(auto cl : _contact_links){
        vars.emplace_back(cl, optimize_contact_torque ? 6 : 3); // put 6 for full wrench
    }

    OpenSoT::OptvarHelper opt(vars);
    
    /* Wrench bounds */
    Eigen::VectorXd wrench_ub(6), wrench_lb(6);
    wrench_ub << 1000, 1000, 1000, 50, 50, 50;
    wrench_lb << -1000, -1000, 0, -50, -50, -50;

    std::vector<OpenSoT::constraints::GenericConstraint::Ptr> wrench_bounds;
    std::list<OpenSoT::solvers::iHQP::TaskPtr> min_wrench_tasks;

    /* Define affine mappings for all wrenches */
    for(auto cl : _contact_links){

        _wrenches.emplace_back(opt.getVariable(cl) /
                               OpenSoT::AffineHelper::Zero(opt.getSize(), optimize_contact_torque ? 0 : 3)
                              );
        
        
        wrench_bounds.push_back( boost::make_shared<OpenSoT::constraints::GenericConstraint>(cl+"_bound",
                                                                                             _wrenches.back(),
                                                                                             wrench_ub,
                                                                                             wrench_lb,
                                                                                             OpenSoT::constraints::GenericConstraint::Type::CONSTRAINT)
                               );
        
        auto min_wrench = boost::make_shared<OpenSoT::tasks::MinimizeVariable>("MIN_" + cl + "_WRENCH", 
            _wrenches.back()
        );
        
        min_wrench_tasks.push_back(min_wrench);
        
    }
    
    /* Define friction cones */
    OpenSoT::constraints::force::FrictionCone::friction_cones friction_cones;

    for(auto cl : _contact_links)
    {
        friction_cones.emplace_back(cl, 0.3);
    }
    
    auto friction_constraint = boost::make_shared<OpenSoT::constraints::force::FrictionCone>(_wrenches, *_model, friction_cones);
    
    
    
    /* Construct forza giusta task */
    _forza_giusta = boost::make_shared<OpenSoT::tasks::force::FloatingBase>(*_model, _wrenches, _contact_links);
    
    /* Min wrench aggregated */
    auto min_force_aggr = boost::make_shared<OpenSoT::tasks::Aggregated>(min_wrench_tasks, opt.getSize());

    Eigen::Vector2d X_LIMS; X_LIMS.setZero(); X_LIMS[0] = -0.08; X_LIMS[1] = 0.08;
    Eigen::Vector2d Y_LIMS; Y_LIMS.setZero(); Y_LIMS[0] = -0.05; Y_LIMS[1] = 0.05;
    auto CoP = boost::make_shared<OpenSoT::constraints::force::CoP>(*_model, _wrenches, _contact_links,
                                                                    X_LIMS, Y_LIMS);
    
    /* Define optimization problem */
//    _autostack = boost::make_shared<OpenSoT::AutoStack>(min_force_aggr);
//    _autostack << boost::make_shared<OpenSoT::constraints::TaskToConstraint>(_forza_giusta);
    _autostack = boost::make_shared<OpenSoT::AutoStack>(_forza_giusta);
    //_autostack<<CoP;
    _autostack<<wrench_bounds[0]<<wrench_bounds[1];
    
    _solver = boost::make_shared<OpenSoT::solvers::iHQP>(_autostack->getStack(), _autostack->getBounds(), 1.0);
    
    /* Initialize solution */
    _x_value.setZero(opt.getSize());
    
}


bool demo::ForceOptimization::compute(const Eigen::VectorXd& fixed_base_torque, 
                                      std::vector<Eigen::Vector6d>& Fc,
                                      Eigen::VectorXd& tau)
{
    Fc.resize(_contact_links.size());
    
    _forza_giusta->setFloatingBaseTorque(fixed_base_torque.head<6>());
    _autostack->update(_x_value);
    
    if(!_solver->solve(_x_value))
    {
        return false;
    }
    
    tau = fixed_base_torque;
    
    for(int i = 0; i < _contact_links.size(); i++)
    {
        
        _wrenches[i].getValue(_x_value, _Fci);
        Fc[i] = _Fci;
        
        _model->getJacobian(_contact_links[i], _JC);
        _fc_i.noalias() = _JC.transpose()*_Fci;
        tau -= _fc_i;
    }
    
    
    return true;
    
}

void demo::ForceOptimization::log(XBot::MatLogger::Ptr logger)
{
    _autostack->log(logger);
//     _solver->log(logger);
}


#endif
