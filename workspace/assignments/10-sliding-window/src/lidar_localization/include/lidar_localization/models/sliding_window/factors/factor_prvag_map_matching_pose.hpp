/*
 * @Description: ceres residual block for map matching pose measurement
 * @Author: Ge Yao
 * @Date: 2020-11-29 15:47:49
 */
#ifndef LIDAR_LOCALIZATION_MODELS_SLIDING_WINDOW_FACTOR_PRVAG_MAP_MATCHING_POSE_HPP_
#define LIDAR_LOCALIZATION_MODELS_SLIDING_WINDOW_FACTOR_PRVAG_MAP_MATCHING_POSE_HPP_

#include <ceres/ceres.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>

#include <sophus/so3.hpp>

#include "glog/logging.h"

namespace sliding_window {

class FactorPRVAGMapMatchingPose : public ceres::SizedCostFunction<6, 15> {
public:
  static const int INDEX_P = 0;
  static const int INDEX_R = 3;

  FactorPRVAGMapMatchingPose(void){};

  void SetMeasurement(const Eigen::VectorXd &m) { m_ = m; }

  void SetInformation(const Eigen::MatrixXd &I) { I_ = I; }

  virtual bool Evaluate(double const *const *parameters, double *residuals,
                        double **jacobians) const {
    //
    // parse parameters:
    //
    // pose
    Eigen::Map<const Eigen::Vector3d> pos(&parameters[0][INDEX_P]);
    Eigen::Map<const Eigen::Vector3d> log_ori(&parameters[0][INDEX_R]);
    const Sophus::SO3d ori = Sophus::SO3d::exp(log_ori);

    //
    // parse measurement:
    //
    const Eigen::Vector3d &pos_prior = m_.block<3, 1>(INDEX_P, 0);
    const Eigen::Vector3d &log_ori_prior = m_.block<3, 1>(INDEX_R, 0);
    const Sophus::SO3d ori_prior = Sophus::SO3d::exp(log_ori_prior);

    //
    // get square root of information matrix:
    //
    Eigen::LLT<Eigen::Matrix<double, 6, 6>> LowerI(I_);
    Eigen::Matrix<double, 6, 6> sqrt_info = LowerI.matrixL().transpose();
    //
    // compute residual:
    //
    Eigen::Map<Eigen::Matrix<double, 6, 1>> res(residuals);
    res.segment<3>(INDEX_P) = pos - pos_prior;
    res.segment<3>(INDEX_R) = (ori_prior.inverse() * ori).log();
    //
    // compute jacobians:
    //
    if (jacobians) {
      const Eigen::Matrix3d J_r_inv = JacobianRInv(res.segment<3>(INDEX_R));
      if (jacobians[0]) {
        // implement jacobian computing:
        Eigen::Map<Eigen::Matrix<double, 6, 15, Eigen::RowMajor>> jacobian(
            jacobians[0]);
        jacobian.setZero();
        jacobian.block<3, 3>(INDEX_P, INDEX_P).setIdentity();
        jacobian.block<3, 3>(INDEX_R, INDEX_R) = J_r_inv;
        jacobian = sqrt_info * jacobian;
      }
    }

    //
    // correct residual by square root of information matrix:
    //
    res = sqrt_info * res;

    return true;
  }

private:
  static Eigen::Matrix3d JacobianRInv(const Eigen::Vector3d &w) {
    Eigen::Matrix3d J_r_inv = Eigen::Matrix3d::Identity();

    double theta = w.norm();
    double half_theta = 0.5 * theta;

    if (theta > 1e-5) {
      Eigen::Vector3d k = w.normalized();
      Eigen::Matrix3d K = Sophus::SO3d::hat(k);
      const double half_theta_cot_half_theta =
          half_theta * std::cos(half_theta) / std::sin(half_theta);
      Eigen::Matrix3d I33 = Eigen::Matrix3d::Identity();
      J_r_inv = half_theta_cot_half_theta * I33 +
                (1 - half_theta_cot_half_theta) * k * k.transpose() +
                half_theta * K;
      // Eigen::Matrix3d temp =
      //     I33 + 0.5 * K +
      //     (1.0 - (1.0 + std::cos(theta)) * theta / (2.0 * std::sin(theta))) *
      //         K * K;
      // LOG(INFO) << "Jr_inv: \n";
      // LOG(INFO) << J_r_inv;
      // LOG(INFO) << "temp: \n";
      // LOG(INFO) << temp;
      // J_r_inv =
      //     J_r_inv + 0.5 * K +
      //     (1.0 - (1.0 + std::cos(theta)) * theta / (2.0 * std::sin(theta))) *
      //         K * K;
    }

    return J_r_inv;
  }

  Eigen::VectorXd m_;
  Eigen::MatrixXd I_;
};

} // namespace sliding_window

#endif // LIDAR_LOCALIZATION_MODELS_SLIDING_WINDOW_FACTOR_PRVAG_MAP_MATCHING_POSE_HPP_
