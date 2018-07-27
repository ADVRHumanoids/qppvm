#include <OpenSoT/solvers/BackEndFactory.h>
#include <OpenSoT/solvers/GLPKBackEnd.h>

int main()
{
    auto logger = XBot::MatLogger::getLogger("/tmp/glpk_test_log");
    
    auto be = OpenSoT::solvers::BackEndFactory(OpenSoT::solvers::solver_back_ends::GLPK, 
                                               20,
                                               26,
                                               OpenSoT::HessianType::HST_ZERO,
                                               1e-6
                                               );
    
    Eigen::MatrixXd Aeq;
    Aeq.setZero(6, 20);
    Aeq.leftCols(10).setRandom();
    Eigen::VectorXd beq;
    beq.setRandom(6);
    
    Eigen::VectorXd xmin, xmax;
    xmax.setConstant(20, 10.);
    xmin = -xmax;
    
    xmin.tail(10).setZero();
    xmax.tail(10).setOnes();
    
    const double M = 15;
    Eigen::MatrixXd Alog, eyeN;
    eyeN.setIdentity(10,10);
    Alog.setZero(20,20);
    Alog << eyeN, -M*eyeN,
           -eyeN, -M*eyeN;
           
    Eigen::VectorXd lAlog, uAlog;
    lAlog.setConstant(20, -1e30);
    uAlog.setConstant(20, 0);
    
    Eigen::MatrixXd A(26,20);
    Eigen::VectorXd lA(26), uA(26);
    
    A.topRows(6) = Aeq;
    A.bottomRows(20) = Alog;
    lA.head(6) = uA.head(6) = beq;
    lA.tail(20) = lAlog;
    uA.tail(20) = uAlog;
    
    
    Eigen::VectorXd c;
    c.setZero(20);
    c.tail(10).setOnes();
    
    std::cout << "\n\nA:\n" << A << "\n";
    std::cout << "\n\nlA:\n" << lA << "\n";
    std::cout << "\n\nuA:\n" << uA << "\n";
    
    logger->log("A", A);
    logger->log("lA", lA);
    logger->log("uA", uA);
    logger->log("xmin", xmin);
    logger->log("xmax", xmax);
    
    
    OpenSoT::solvers::GLPKBackEnd::GLPKBackEndOptions opt_GLPK;
    for(int i = 10; i < 20; i++)
    {
       opt_GLPK.var_id_kind_.push_back(std::pair<int,int>(i,GLP_BV)); 
    }
    
    be->setOptions(opt_GLPK);
    
    
    
    be->initProblem(Eigen::MatrixXd::Zero(20,20), c, A, lA, uA, xmin, xmax);
    std::cout << "\n\n\n" << be->getSolution() << std::endl;
    
    
    
    for(int i = 0; i < 5; i++)
    {
        A.block(0,0,6,10).setRandom();
        beq.setRandom(6);
        lA.head(6) = uA.head(6) = beq;
        
        be->updateConstraints(A, lA, uA);
        be->solve();
        std::cout << "\n\nSolution:\n" << be->getSolution() << std::endl;
        std::cout << "\n\nResidual:\n" << A.topRows(6)*be->getSolution() - beq << std::endl;
    }
    
    logger->flush();
    
    return 0;

    
}
