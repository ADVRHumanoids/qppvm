#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <Eigen/Geometry> 

namespace Eigen {
    typedef Matrix<double, 6, 1> Vector6d;
}

namespace ifopt {
    
using Eigen::Vector2d;


class ExVariables : public VariableSet {
public:
  // Every variable set has a name, here "var_set1". this allows the constraints
  // and costs to define values and Jacobians specifically w.r.t this variable set.
    
  ExVariables(const std::string& name) : VariableSet(3, name)
  {
    // the initial values where the NLP starts iterating from
    x0_ = 0.0;
    x1_ = 0.0;
    x2_ = 0.0;;
    
    lb_.setConstant(-1000.0);
    ub_.setConstant(1000.0);
    
  }

  // Here is where you can transform the Eigen::Vector into whatever
  // internal representation of your variables you have (here two doubles, but
  // can also be complex classes such as splines, etc..
  void SetVariables(const VectorXd& x) override
  {
    x0_ = x(0);
    x1_ = x(1);
    x2_ = x(2);
  };

  // Here is the reverse transformation from the internal representation to
  // to the Eigen::Vector
  VectorXd GetValues() const override
  {
      Eigen::VectorXd temp; 
      temp.setZero(3);
      temp << x0_, x1_, x2_;
      
      return temp;

  };
  
  void SetBounds(const Eigen::Vector3d& lower, const Eigen::Vector3d& upper)
  {
   
      lb_ = lower;
      ub_ = upper;
      
      if( ((ub_ - lb_).array() < 0).any() )
      {
          throw std::invalid_argument("Inconsistent bounds");
      }
      
  }

  // Each variable has an upper and lower bound set here
  VecBound GetBounds() const override
  {
    VecBound bounds(GetRows());
    bounds.at(0) = Bounds(lb_(0), ub_(0));
    bounds.at(1) = Bounds(lb_(1), ub_(1));
    bounds.at(2) = Bounds(lb_(2), ub_(2));
    return bounds;
  }

private:
  double x0_, x1_, x2_;
  Eigen::Vector3d lb_, ub_;
  
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// STATIC RELATION: cost term 

class StaticCost: public CostTerm {
public:
    
    StaticCost() : StaticCost("static_cost_term") 
    {
        _wrench_ext.setZero();
        _mg << 0.0, 0.0, -100;
        _com.setZero();
    }
    StaticCost(const std::string& name) : CostTerm(name) {}

    
    double GetCost() const override
    {
        
        double value = 0;
        
        Eigen::Vector6d cost;
        cost.setZero();
        
        for(int i : {1, 2, 3, 4})
        {
            Eigen::Vector3d Fi = GetVariables()->GetComponent("F" + std::to_string(i))->GetValues();
            Eigen::Vector3d pi = GetVariables()->GetComponent("p" + std::to_string(i))->GetValues();
            cost.head<3>() += Fi;
            cost.tail<3>() += (pi-_com).cross(Fi);
        }
        
        
        cost -= _wrench_ext;
        cost.head<3>() += _mg;
        
//         std::cout <<"cost: " << cost.transpose() << std::endl;     
        
        value = cost.dot(cost);
        
        return value;
        
    };
    
    void SetExternalWrench(const Eigen::Vector6d& w)
    {
        _wrench_ext = w;
    }
    
    void SetCoM(const Eigen::Vector3d& com)
    {
        _com = com;
    }

    void FillJacobianBlock (std::string var_set, Jacobian& jac) const override
    {
        jac.setZero();
    }

private:

    Eigen::Vector6d _wrench_ext;
    Eigen::Vector3d _mg, _com;

    
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ExCost: public CostTerm {
public:
    
    ExCost() : ExCost("cost_term") 
    {
        _W_p = 0;
        _p_ref.setZero(12);
    }
    ExCost(const std::string& name) : CostTerm(name) {}

    
    void SetPosRef(const Eigen::VectorXd& p_ref, const double& W_p)
    {
        _W_p = W_p; 
        _p_ref = p_ref;
    }
    
    
    double GetCost() const override
    {
        
            double value = 0;
                
            for(int i : {1, 2, 3, 4})
            {
                Eigen::Vector3d Fi = GetVariables()->GetComponent("F" + std::to_string(i))->GetValues();
                Eigen::Vector3d pi = GetVariables()->GetComponent("p" + std::to_string(i))->GetValues();
                
                value += _W_p*(pi -_p_ref.segment(3*(i-1),3)).norm()/2.0 + Fi.norm()/2.0;
            } 
            
            return value;
    };

    void FillJacobianBlock (std::string var_set, Jacobian& jac) const override
    {
        jac.setZero();
        for(int i = 0; i < 4; i++)
        {
            if(var_set == ("F" + std::to_string(i+1)))
            {
                Eigen::Vector3d Fi = GetVariables()->GetComponent("F" + std::to_string(i+1))->GetValues();
                
                jac.coeffRef(0, 0) = Fi.x();
                jac.coeffRef(0, 1) = Fi.y();
                jac.coeffRef(0, 2) = Fi.z();
                
            }
            
            if(var_set == ("p" + std::to_string(i+1)))
            {
                Eigen::Vector3d pi = GetVariables()->GetComponent("p" + std::to_string(i+1))->GetValues();
                
                jac.coeffRef(0, 0) = _W_p*(pi.x()-_p_ref(3*i));
                jac.coeffRef(0, 1) = _W_p*(pi.y()-_p_ref(3*i+1));
                jac.coeffRef(0, 2) = _W_p*(pi.z()-_p_ref(3*i+2));
                
            }
        }
    }

private:

   double _W_p;
   Eigen::VectorXd _p_ref;

    
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class StaticConstraint : public ConstraintSet {
    
public:

    StaticConstraint() : ConstraintSet(6, "StaticConstraint")
    {
        _wrench_ext.setZero();
        _mg << 0.0, 0.0, -100;
        _com.setZero();
    }

    // The constraint value minus the constant value "1", moved to bounds.
    VectorXd GetValues() const override
    {
        Eigen::Vector6d value;
        value.setZero();
        
        for(int i : {1, 2, 3, 4})
        {
            Eigen::Vector3d Fi = GetVariables()->GetComponent("F" + std::to_string(i))->GetValues();
            Eigen::Vector3d pi = GetVariables()->GetComponent("p" + std::to_string(i))->GetValues();
            value.head<3>() += Fi;
            value.tail<3>() += (pi-_com).cross(Fi);
        }
        
        
        value -= _wrench_ext;
        value.head<3>() += _mg;
        
        return value;
        
    };
    
    void SetExternalWrench(const Eigen::Vector6d& w)
    {
        _wrench_ext = w;
    }
    
    void SetCoM(const Eigen::Vector3d& com)
    {
        _com = com;
    }

    // The only constraint in this set is an equality constraint to 1.
    // Constant values should always be put into GetBounds(), not GetValues().
    // For inequality constraints (<,>), use Bounds(x, inf) or Bounds(-inf, x).
    VecBound GetBounds() const override
    {
        VecBound b(GetRows());
        for(int i = 0; i < 6; i++)
        {            
            b.at(i) = Bounds(.0, .0);  
        }
           
        return b;
    }

    // This function provides the first derivative of the constraints.
    // In case this is too difficult to write, you can also tell the solvers to
    // approximate the derivatives by finite differences and not overwrite this
    // function, e.g. in ipopt.cc::use_jacobian_approximation_ = true
    void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override
    {
        jac_block.setZero();
        for(int i = 0; i < 4; i++)
        {
            if(var_set == ("F" + std::to_string(i+1)))
            {
                jac_block.coeffRef(0, 0) = 1.0;
                jac_block.coeffRef(1, 1) = 1.0;
                jac_block.coeffRef(2, 2) = 1.0;
                
                Eigen::Vector3d pi = GetVariables()->GetComponent("p" + std::to_string(i+1))->GetValues();
                jac_block.coeffRef(3, 1) = -(pi.z()-_com.z());
                jac_block.coeffRef(3, 2) =   pi.y()-_com.y();
                jac_block.coeffRef(4, 0) =   pi.z()-_com.z();
                jac_block.coeffRef(4, 2) = -(pi.x()-_com.x());
                jac_block.coeffRef(5, 0) = -(pi.y()-_com.y());
                jac_block.coeffRef(5, 1) =   pi.x()-_com.x();
            }
            
            if(var_set == ("p" + std::to_string(i+1)))
            {
                Eigen::Vector3d Fi = GetVariables()->GetComponent("F" + std::to_string(i+1))->GetValues();
                jac_block.coeffRef(3, 1) =  Fi.z();
                jac_block.coeffRef(3, 2) = -Fi.y();
                jac_block.coeffRef(4, 0) = -Fi.z();
                jac_block.coeffRef(4, 2) =  Fi.x();
                jac_block.coeffRef(5, 0) =  Fi.y();
                jac_block.coeffRef(5, 1) = -Fi.x();
            }
        }
    }
    
private:
    
    Eigen::Vector6d _wrench_ext;
    Eigen::Vector3d _mg, _com;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class SuperEllipsoidConstraint : public ConstraintSet {
    
public:

    SuperEllipsoidConstraint(const std::string& position_name) :
        ConstraintSet(1, "SuperEllipsoidConstraint " + position_name),
        pname_(position_name)
    {
        _C.setZero();
        _R.setOnes();
        _P << 8, 8, 4;
    }

    // The constraint value minus the constant value "1", moved to bounds.
    VectorXd GetValues() const override
    {
       
        Eigen::VectorXd value;
        value.setZero(1);
        
        Eigen::Vector3d p = GetVariables()->GetComponent(pname_)->GetValues(); 
             
        for(int i = 0; i < 3; i++)
        {            
           value(0) += pow((p(i)-_C(i))/_R(i),_P(i));
        }
        
        value(0) -= 1;
                
        return value;
        
    };
    
    void SetParam(const Eigen::Vector3d& C, const Eigen::Vector3d& R, const Eigen::Vector3d& P)
    {
        _C = C; // superellipsoid center
        _R = R; // radial axis
        _P = P; // r s t parameters (r=s super-ellipsoid)
    }

    // The only constraint in this set is an equality constraint to 1.
    // Constant values should always be put into GetBounds(), not GetValues().
    // For inequality constraints (<,>), use Bounds(x, inf) or Bounds(-inf, x).
    VecBound GetBounds() const override
    {
        VecBound b(GetRows());
         
        b.at(0) = Bounds(0.0, 0.0);  
          
        return b;
    }

    // This function provides the first derivative of the constraints.
    // In case this is too difficult to write, you can also tell the solvers to
    // approximate the derivatives by finite differences and not overwrite this
    // function, e.g. in ipopt.cc::use_jacobian_approximation_ = true
    void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override
    {
        jac_block.setZero();
        
        if(var_set == pname_)
        {
            Eigen::Vector3d p = GetVariables()->GetComponent(pname_)->GetValues();
                          
            jac_block.coeffRef(0, 0) = _P.x()/pow(_R.x(),_P.x()) * pow(p.x()-_C.x(),_P.x()-1);
            jac_block.coeffRef(0, 1) = _P.y()/pow(_R.y(),_P.y()) * pow(p.y()-_C.y(),_P.y()-1);
            jac_block.coeffRef(0, 2) = _P.z()/pow(_R.z(),_P.z()) * pow(p.z()-_C.z(),_P.z()-1);
                               
        }
    }
    
private:
    
    Eigen::Vector3d _C, _R, _P;
    std::string pname_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class FrictionConstraint : public ConstraintSet {
    
public:

    FrictionConstraint(const std::string& force_name) :
        ConstraintSet(2, "FrictionConstraint " + force_name),
        fname_(force_name)
    {
        mu_= 1; 
              
        _C.setZero();
        _R.setOnes();
        _P << 8, 8, 4;

 
    }
    
    
    void set_mu(const double& mu)
    {
//         mu_ = std::sqrt(2.*mu)/2.;
        mu_ = mu;

    }
    
//     void set_rot(const Eigen::Vector3d& norm_vec)
//     {
// 
//       if (norm_vec == Eigen::Vector3d::UnitX())
//       {
//         R = Eigen::AngleAxisd(-0.5*M_PI, Eigen::Vector3d::UnitY());          
//       }
//       else if (norm_vec == - Eigen::Vector3d::UnitX())
//       {
//         R = Eigen::AngleAxisd(0.5*M_PI, Eigen::Vector3d::UnitY());
//       }      
//       else if (norm_vec == Eigen::Vector3d::UnitY())
//       {
//         R = Eigen::AngleAxisd(-0.5*M_PI, Eigen::Vector3d::UnitX());
//       }
//       else if (norm_vec == - Eigen::Vector3d::UnitY())
//       {
//         R = Eigen::AngleAxisd(0.5*M_PI, Eigen::Vector3d::UnitX());
//       }
//       else if (norm_vec == - Eigen::Vector3d::UnitZ())
//       {
//         R = Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitX());  
//       }
//       
//     }
     
    void setParam_SE(const Eigen::Vector3d& C, const Eigen::Vector3d& R, const Eigen::Vector3d& P)
    {
        
        _C = C; // superellipsoid center
        _R = R; // radial axis
        _P = P; // r s t parameters (r=s super-ellipsoid)
              
    }

    VectorXd GetValues() const override
    {
            Eigen::VectorXd value;
            
            value.setZero(2);
        
//             value.setZero(5);
            
            Eigen::Vector3d F = GetVariables()->GetComponent(fname_)->GetValues();  
  
            int k; 
            for(int i = 0; i < 4; i++)
            {
                if(fname_ == ("F" + std::to_string(i+1)))
                {
                    k = i + 1;
                }
            }                               
                       
            Eigen::Vector3d p = GetVariables()->GetComponent("p" + std::to_string(k))->GetValues();
            Eigen::Vector3d  n_SE = GetVariables()->GetComponent("n" + std::to_string(k))->GetValues();

//             std::cout <<"k: " << k << std::endl;              
//             std::cout <<"p: " << p.transpose() << std::endl; 

//             F << -0.3690, 0.0000, 25.4988;
//             p << 0.1000, -0.3000, 0.0147;
            
                                       
//             Eigen::Vector3d n_SE; 
// 
//             n_SE.x() = - _P.x() / pow(_R.x(),_P.x()) * pow(p.x()-_C.x(),_P.x()-1);
//             n_SE.y() = - _P.y() / pow(_R.y(),_P.y()) * pow(p.y()-_C.y(),_P.y()-1);
//             n_SE.z() = - _P.z() / pow(_R.z(),_P.z()) * pow(p.z()-_C.z(),_P.z()-1);
            
      
//             n_SE.normalize();   
            
//             std::cout <<"n_SE" << n_SE.transpose() << std::endl;
                       
//             Eigen::Matrix3d R;
//             R.setIdentity();
//                      
//             R.coeffRef(0, 0) =  n_SE.y()/((n_SE.head(2)).norm()); 
//             R.coeffRef(0, 1) = -n_SE.x()/((n_SE.head(2)).norm());  
//             
//             R.coeffRef(1, 0) =  (n_SE.x()*n_SE.z())/((n_SE.head(2)).norm());  
//             R.coeffRef(1, 1) =  (n_SE.y()*n_SE.z())/((n_SE.head(2)).norm());  
//             R.coeffRef(1, 2) = -(n_SE.head(2)).norm();  
//             
//             R.coeffRef(2, 0) = n_SE.x();  
//             R.coeffRef(2, 1) = n_SE.y();  
//             R.coeffRef(2, 2) = n_SE.z();  
                        
//             std::cout <<"R" << R << std::endl;
                        
//             Eigen::MatrixXd A; A.setZero(5,3);
//             A.coeffRef(0, 0) =  1.0;
//             A.coeffRef(0, 2) = -mu_;
//         
//             A.coeffRef(1, 0) = -1.0;
//             A.coeffRef(1, 2) = -mu_;
//         
//             A.coeffRef(2, 1) =  1.0;
//             A.coeffRef(2, 2) = -mu_;
//         
//             A.coeffRef(3, 1) = -1.0;
//             A.coeffRef(3, 2) = -mu_;
//         
//             A.coeffRef(4, 2) = -1.0;     
                    
//             value = A*R*F;
            
//             std::cout <<"F_local " << (R*F).transpose() << std::endl;
//             std::cout <<"value " << value.transpose() << std::endl;
            
            value(0) = -F.dot(n_SE); 
            value(1) = (F-(n_SE.dot(F))*n_SE).norm() - mu_*(F.dot(n_SE));
                  
//             std::cout <<"value " << value.transpose() << std::endl;        
  
        return value;
        
    };

    // Constant values should always be put into GetBounds(), not GetValues().
    // For inequality constraints (<,>), use Bounds(x, inf) or Bounds(-inf, x).
    VecBound GetBounds() const override
    {
        VecBound b(GetRows());
        for(int i = 0; i < 2; i++)
        {            
             b.at(i) = BoundSmallerZero;               
        }
                            
        return b;
    }

    // This function provides the first derivative of the constraints.
    // In case this is too difficult to write, you can also tell the solvers to
    // approximate the derivatives by finite differences and not overwrite this
    // function, e.g. in ipopt.cc::use_jacobian_approximation_ = true
    void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override
    {
        
       if(var_set == fname_)
       {     
            
           Eigen::Vector3d F = GetVariables()->GetComponent(fname_)->GetValues();  
  
            int k; 
            for(int i = 0; i < 4; i++)
            {
                if(fname_ == ("F" + std::to_string(i+1)))
                {
                    k = i + 1;
                }
            }                               
                                   
            Eigen::Vector3d  n_SE = GetVariables()->GetComponent("n" + std::to_string(k))->GetValues();
            
            Eigen::VectorXd fr;
            
            fr.setZero(3);
        
            jac_block.setZero();
        
            jac_block.coeffRef(0, 3) = -n_SE.x();
            jac_block.coeffRef(0, 4) = -n_SE.y();
            jac_block.coeffRef(0, 5) = -n_SE.z();
            jac_block.coeffRef(0, 6) = -F.x();
            jac_block.coeffRef(0, 7) = -F.y();
            jac_block.coeffRef(0, 8) = -F.z();
            
            jac_block.coeffRef(1, 3) = - mu_*n_SE.x();
            jac_block.coeffRef(1, 4) = - mu_*n_SE.y();
            jac_block.coeffRef(1, 5) = - mu_*n_SE.z();
            jac_block.coeffRef(1, 6) = - mu_*F.x();
            jac_block.coeffRef(1, 7) = - mu_*F.y();
            jac_block.coeffRef(1, 8) = - mu_*F.z();
            
       }
        
//         jac_block.coeffRef(0, 0) =  1.0;
//         jac_block.coeffRef(0, 2) = -mu_;
//         
//         jac_block.coeffRef(1, 0) = -1.0;
//         jac_block.coeffRef(1, 2) = -mu_;
//         
//         jac_block.coeffRef(2, 1) =  1.0;
//         jac_block.coeffRef(2, 2) = -mu_;
//         
//         jac_block.coeffRef(3, 1) = -1.0;
//         jac_block.coeffRef(3, 2) = -mu_;
//         
//         jac_block.coeffRef(4, 2) = -1.0;
//                             
//         Eigen::MatrixXd J = (jac_block * R);
//         jac_block = J.sparseView();

      
    }
    
private:
    
  double mu_;
  std::string fname_;
  Eigen::Vector3d _C, _R, _P;
  bool SE_flag;
    
};



} // namespace opt


