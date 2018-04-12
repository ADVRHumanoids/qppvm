#ifndef __QPPVM_FLOATING_BASE_EST_H__
#define __QPPVM_FLOATING_BASE_EST_H__

#include <XBotInterface/ModelInterface.h>
#include <OpenSoT/utils/Piler.h>

namespace estimation {

    class FloatingBaseEstimator {

    public:

        typedef std::shared_ptr<FloatingBaseEstimator> Ptr;

        FloatingBaseEstimator(XBot::ModelInterface::Ptr model,
                              XBot::ImuSensor::ConstPtr imu,
                              std::vector<std::string> contact_links,
                              const Eigen::MatrixXd& contact_matrix
                              );



        void update(double dt);

        void log(XBot::MatLogger::Ptr logger);

        Eigen::Vector3d getVelocity() const;


    private:

        typedef Eigen::Matrix<double, -1, -1, 0, 30, 3> LimitedSizeMatrix;

        Eigen::Vector3d _qdot_est;
        Eigen::VectorXd _qdot, _q;
        Eigen::MatrixXd _Jtmp, _KJtmp;
        OpenSoT::utils::MatrixPiler _Jc;

        LimitedSizeMatrix _A;
        Eigen::VectorXd _b;
        Eigen::ColPivHouseholderQR<LimitedSizeMatrix> _qr_solver;

        XBot::ModelInterface::Ptr _model;
        XBot::ImuSensor::ConstPtr _imu;
        std::vector<std::string> _contact_links;
        Eigen::MatrixXd _contact_matrix;



    };



}

Eigen::Vector3d estimation::FloatingBaseEstimator::getVelocity() const
{
    return _qdot_est;
}

void estimation::FloatingBaseEstimator::update(double dt)
{
    /* Update angular part from IMU */
    _model->getJointPosition(_q);


    _model->setFloatingBaseState(_imu);
    _model->update();
    _model->getJointPosition(_q);


    /* Compute contact jacobian */
    _Jc.reset(_model->getJointNum());

    for(const auto& cl : _contact_links)
    {
          _model->getJacobian(cl, cl, _Jtmp);
          _KJtmp.noalias() = _contact_matrix*_Jtmp;
          _Jc.pile(_KJtmp);
    }

    /* Solve for linear velocity */
    int na = _model->getActuatedJointNum();
    _model->getJointVelocity(_qdot);
    _model->getJointPosition(_q);

    _A = _Jc.generate_and_get().leftCols<3>();
    _b.noalias() = -_Jc.generate_and_get().rightCols(na+3)*_qdot.tail(na+3);
    _qr_solver.compute(_A);

    _qdot_est = _qr_solver.solve(_b);

    _qdot.head<3>() = _qdot_est;
    _q.head<3>() += _qdot_est * dt;

    _model->setJointPosition(_q);
    _model->setJointVelocity(_qdot);
    _model->update();


}

estimation::FloatingBaseEstimator::FloatingBaseEstimator(XBot::ModelInterface::Ptr model,
                                                         XBot::ImuSensor::ConstPtr imu,
                                                         std::vector< std::string > contact_links,
                                                         const Eigen::MatrixXd& contact_matrix):
    _model(model),
    _imu(imu),
    _contact_links(contact_links),
    _contact_matrix(contact_matrix)
{
    update(0);
}


void estimation::FloatingBaseEstimator::log(XBot::MatLogger::Ptr logger)
{
    logger->add("fb_est_qdot", _qdot_est);
    logger->add("fb_est_A", _A);
    logger->add("fb_est_b", _b);
    logger->add("contact_speed", _Jc.generate_and_get()*_qdot);
}





#endif