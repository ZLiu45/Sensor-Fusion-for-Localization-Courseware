/*
 * @Description: IMU pre-integrator for LIO mapping, implementation
 * @Author: Ge Yao
 * @Date: 2020-11-29 15:47:49
 */

#include "lidar_localization/models/pre_integrator/imu_pre_integrator.hpp"

#include "lidar_localization/global_defination/global_defination.h"

#include "glog/logging.h"

namespace lidar_localization {

namespace {
Eigen::Matrix3d skewMat(const Eigen::Vector3d &rec) {
  Eigen::Matrix3d res = Eigen::Matrix3d::Zero();
  res << 0, -rec.z(), rec.y(), rec.z(), 0.0, -rec.x(), -rec.y(), rec.x(), 0.0;
  return res;
}
} // namespace

IMUPreIntegrator::IMUPreIntegrator(const YAML::Node &node) {
  //
  // parse config:
  //
  // a. earth constants:
  EARTH.GRAVITY_MAGNITUDE = node["earth"]["gravity_magnitude"].as<double>();
  // b. process noise:
  COV.MEASUREMENT.ACCEL =
      node["covariance"]["measurement"]["accel"].as<double>();
  COV.MEASUREMENT.GYRO = node["covariance"]["measurement"]["gyro"].as<double>();
  COV.RANDOM_WALK.ACCEL =
      node["covariance"]["random_walk"]["accel"].as<double>();
  COV.RANDOM_WALK.GYRO = node["covariance"]["random_walk"]["gyro"].as<double>();

  // prompt:
  LOG(INFO) << std::endl
            << "IMU Pre-Integration params:" << std::endl
            << "\tgravity magnitude: " << EARTH.GRAVITY_MAGNITUDE << std::endl
            << std::endl
            << "\tprocess noise:" << std::endl
            << "\t\tmeasurement:" << std::endl
            << "\t\t\taccel.: " << COV.MEASUREMENT.ACCEL << std::endl
            << "\t\t\tgyro.: " << COV.MEASUREMENT.GYRO << std::endl
            << "\t\trandom_walk:" << std::endl
            << "\t\t\taccel.: " << COV.RANDOM_WALK.ACCEL << std::endl
            << "\t\t\tgyro.: " << COV.RANDOM_WALK.GYRO << std::endl
            << std::endl;

  // a. gravity constant:
  state.g_ = Eigen::Vector3d(0.0, 0.0, EARTH.GRAVITY_MAGNITUDE);

  // b. process noise:
  Q_.block<3, 3>(INDEX_M_ACC_PREV, INDEX_M_ACC_PREV) =
      Q_.block<3, 3>(INDEX_M_ACC_CURR, INDEX_M_ACC_CURR) =
          COV.MEASUREMENT.ACCEL * Eigen::Matrix3d::Identity();
  Q_.block<3, 3>(INDEX_M_GYR_PREV, INDEX_M_GYR_PREV) =
      Q_.block<3, 3>(INDEX_M_GYR_CURR, INDEX_M_GYR_CURR) =
          COV.MEASUREMENT.GYRO * Eigen::Matrix3d::Identity();
  Q_.block<3, 3>(INDEX_R_ACC_PREV, INDEX_R_ACC_PREV) =
      COV.RANDOM_WALK.ACCEL * Eigen::Matrix3d::Identity();
  Q_.block<3, 3>(INDEX_R_GYR_PREV, INDEX_R_GYR_PREV) =
      COV.RANDOM_WALK.GYRO * Eigen::Matrix3d::Identity();

  // c. process equation, state propagation:
  F_.block<3, 3>(INDEX_ALPHA, INDEX_BETA) = Eigen::Matrix3d::Identity();
  F_.block<3, 3>(INDEX_THETA, INDEX_B_G) = -Eigen::Matrix3d::Identity();

  // d. process equation, noise input:
  B_.block<3, 3>(INDEX_THETA, INDEX_M_GYR_PREV) =
      B_.block<3, 3>(INDEX_THETA, INDEX_M_GYR_CURR) =
          0.50 * Eigen::Matrix3d::Identity();
  B_.block<3, 3>(INDEX_B_A, INDEX_R_ACC_PREV) =
      B_.block<3, 3>(INDEX_B_G, INDEX_R_GYR_PREV) = Eigen::Matrix3d::Identity();
}

/**
 * @brief  reset IMU pre-integrator
 * @param  init_imu_data, init IMU measurements
 * @return true if success false otherwise
 */
bool IMUPreIntegrator::Init(const IMUData &init_imu_data) {
  // reset pre-integrator state:
  ResetState(init_imu_data);

  // mark as inited:
  is_inited_ = true;

  return true;
}

/**
 * @brief  update IMU pre-integrator
 * @param  imu_data, current IMU measurements
 * @return true if success false otherwise
 */
bool IMUPreIntegrator::Update(const IMUData &imu_data) {
  if (imu_data_buff_.front().time < imu_data.time) {
    // set buffer:
    imu_data_buff_.push_back(imu_data);

    // update state mean, covariance and Jacobian:
    UpdateState();

    // move forward:
    imu_data_buff_.pop_front();
  }

  return true;
}

/**
 * @brief  reset IMU pre-integrator using new init IMU measurement
 * @param  init_imu_data, new init IMU measurements
 * @param  output pre-integration result for constraint building as
 * IMUPreIntegration
 * @return true if success false otherwise
 */
bool IMUPreIntegrator::Reset(const IMUData &init_imu_data,
                             IMUPreIntegration &imu_pre_integration) {
  // one last update:
  Update(init_imu_data);

  // set output IMU pre-integration:
  imu_pre_integration.T_ = init_imu_data.time - time_;

  // set gravity constant:
  imu_pre_integration.g_ = state.g_;

  // set measurement:
  imu_pre_integration.alpha_ij_ = state.alpha_ij_;
  imu_pre_integration.theta_ij_ = state.theta_ij_;
  imu_pre_integration.beta_ij_ = state.beta_ij_;
  imu_pre_integration.b_a_i_ = state.b_a_i_;
  imu_pre_integration.b_g_i_ = state.b_g_i_;
  // set information:
  imu_pre_integration.P_ = P_;
  // set Jacobian:
  imu_pre_integration.J_ = J_;

  // reset:
  ResetState(init_imu_data);

  return true;
}

/**
 * @brief  reset pre-integrator state using IMU measurements
 * @param  void
 * @return void
 */
void IMUPreIntegrator::ResetState(const IMUData &init_imu_data) {
  // reset time:
  time_ = init_imu_data.time;

  // a. reset relative translation:
  state.alpha_ij_ = Eigen::Vector3d::Zero();
  // b. reset relative orientation:
  state.theta_ij_ = Sophus::SO3d();
  // c. reset relative velocity:
  state.beta_ij_ = Eigen::Vector3d::Zero();
  // d. set init bias, acceleometer:
  state.b_a_i_ =
      Eigen::Vector3d(init_imu_data.accel_bias.x, init_imu_data.accel_bias.y,
                      init_imu_data.accel_bias.z);
  // d. set init bias, gyroscope:
  state.b_g_i_ =
      Eigen::Vector3d(init_imu_data.gyro_bias.x, init_imu_data.gyro_bias.y,
                      init_imu_data.gyro_bias.z);

  // reset state covariance:
  P_ = MatrixP::Zero();

  // reset Jacobian:
  J_ = MatrixJ::Identity();

  // reset buffer:
  imu_data_buff_.clear();
  imu_data_buff_.push_back(init_imu_data);
}

/**
 * @brief  update pre-integrator state: mean, covariance and Jacobian
 * @param  void
 * @return void
 */
void IMUPreIntegrator::UpdateState(void) {
  static double dt = 0.0;

  static Eigen::Vector3d w_mid = Eigen::Vector3d::Zero();
  static Eigen::Vector3d a_mid = Eigen::Vector3d::Zero();

  static Sophus::SO3d prev_theta_ij = Sophus::SO3d();
  static Sophus::SO3d curr_theta_ij = Sophus::SO3d();
  static Sophus::SO3d d_theta_ij = Sophus::SO3d();

  static Eigen::Matrix3d dR_inv = Eigen::Matrix3d::Identity();
  static Eigen::Matrix3d prev_R = Eigen::Matrix3d::Identity();
  static Eigen::Matrix3d curr_R = Eigen::Matrix3d::Identity();
  static Eigen::Matrix3d prev_R_a_hat = Eigen::Matrix3d::Identity();
  static Eigen::Matrix3d curr_R_a_hat = Eigen::Matrix3d::Identity();

  //
  // parse measurements:
  //
  // get measurement handlers:
  const IMUData &prev_imu_data = imu_data_buff_.at(0);
  const IMUData &curr_imu_data = imu_data_buff_.at(1);

  // get time delta:
  dt = curr_imu_data.time - prev_imu_data.time;

  // get measurements:
  const Eigen::Vector3d prev_w(
      prev_imu_data.angular_velocity.x - state.b_g_i_.x(),
      prev_imu_data.angular_velocity.y - state.b_g_i_.y(),
      prev_imu_data.angular_velocity.z - state.b_g_i_.z());
  const Eigen::Vector3d curr_w(
      curr_imu_data.angular_velocity.x - state.b_g_i_.x(),
      curr_imu_data.angular_velocity.y - state.b_g_i_.y(),
      curr_imu_data.angular_velocity.z - state.b_g_i_.z());

  const Eigen::Vector3d prev_a(
      prev_imu_data.linear_acceleration.x - state.b_a_i_.x(),
      prev_imu_data.linear_acceleration.y - state.b_a_i_.y(),
      prev_imu_data.linear_acceleration.z - state.b_a_i_.z());
  const Eigen::Vector3d curr_a(
      curr_imu_data.linear_acceleration.x - state.b_a_i_.x(),
      curr_imu_data.linear_acceleration.y - state.b_a_i_.y(),
      curr_imu_data.linear_acceleration.z - state.b_a_i_.z());

  //
  //
  // 1. get w_mid:
  w_mid = (prev_w + curr_w) * 0.5;
  // // 2. update relative orientation, so3:
  prev_R = state.theta_ij_.matrix();
  d_theta_ij = Sophus::SO3d::exp(w_mid * dt);
  curr_theta_ij = state.theta_ij_ * d_theta_ij;
  curr_R = curr_theta_ij.matrix();
  state.theta_ij_ = curr_theta_ij;
  // // 3. get a_mid:
  a_mid = (prev_R * prev_a + curr_R * curr_a) * 0.5;
  // 4. update relative translation:
  state.alpha_ij_ += state.beta_ij_ * dt + 0.5 * a_mid * dt * dt;
  // 5. update relative velocity:
  state.beta_ij_ += a_mid * dt;
  //
  //
  // 1. intermediate results:
  F_.setIdentity();
  dR_inv = d_theta_ij.inverse().matrix();
  Eigen::Matrix3d skew_w_mid_t = skewMat(w_mid) * dt;
  prev_R_a_hat = prev_R * skewMat(prev_a) * dt;
  curr_R_a_hat = curr_R * skewMat(curr_a) * dt;
  const Eigen::Matrix3d kId33 = Eigen::Matrix3d::Identity();
  //
  //
  // F12 & F32:
  F_.setIdentity();
  F_.block<3, 3>(INDEX_ALPHA, INDEX_THETA) =
      -0.25 * prev_R_a_hat * dt - 0.25 * curr_R_a_hat * dR_inv * dt;
  F_.block<3, 3>(INDEX_BETA, INDEX_THETA) =
      -0.5 * prev_R_a_hat - 0.5 * curr_R_a_hat * dR_inv;
  F_.block<3, 3>(INDEX_ALPHA, INDEX_BETA) = kId33 * dt;
  // F14 & F34:
  F_.block<3, 3>(INDEX_ALPHA, INDEX_B_A) = -0.25 * (curr_R + prev_R) * dt * dt;
  F_.block<3, 3>(INDEX_BETA, INDEX_B_A) = -0.5 * (curr_R + prev_R) * dt;
  // F15 & F35:
  F_.block<3, 3>(INDEX_ALPHA, INDEX_B_G) = 0.25 * curr_R_a_hat * dt * dt;
  F_.block<3, 3>(INDEX_BETA, INDEX_B_G) = 0.5 * curr_R_a_hat * dt;
  // F22:
  F_.block<3, 3>(INDEX_BETA, INDEX_BETA) = -Sophus::SO3d::hat(w_mid);
  //
  //
  B_.setZero();
  // G11 & G31:
  B_.block<3, 3>(INDEX_ALPHA, INDEX_M_ACC_PREV) = 0.25 * prev_R * dt * dt;
  B_.block<3, 3>(INDEX_BETA, INDEX_M_ACC_PREV) = 0.5 * prev_R * dt;
  // G12 & G32:
  B_.block<3, 3>(INDEX_ALPHA, INDEX_M_GYR_PREV) =
      -0.25 * curr_R_a_hat * dt * 0.5 * dt;
  B_.block<3, 3>(INDEX_BETA, INDEX_M_GYR_PREV) = -0.5 * curr_R_a_hat * 0.5 * dt;
  // G13 & G33:
  B_.block<3, 3>(INDEX_ALPHA, INDEX_M_ACC_CURR) = 0.25 * curr_R * dt * dt;
  B_.block<3, 3>(INDEX_BETA, INDEX_M_ACC_CURR) = 0.5 * curr_R * dt;
  // G14 & G34:
  B_.block<3, 3>(INDEX_ALPHA, INDEX_M_GYR_CURR) =
      B_.block<3, 3>(INDEX_ALPHA, INDEX_M_GYR_PREV);
  B_.block<3, 3>(INDEX_BETA, INDEX_M_GYR_CURR) =
      B_.block<3, 3>(INDEX_BETA, INDEX_M_ACC_PREV);

  // G22 & G24
  B_.block<3, 3>(INDEX_THETA, INDEX_M_GYR_PREV) = 0.5 * kId33 * dt;
  B_.block<3, 3>(INDEX_THETA, INDEX_M_GYR_CURR) = 0.5 * kId33 * dt;

  // G45 & G56
  B_.block<3, 3>(INDEX_B_A, INDEX_R_ACC_PREV) = kId33 * dt;
  B_.block<3, 3>(INDEX_B_G, INDEX_R_GYR_PREV) = kId33 * dt;
  P_ = F_ * P_ * F_.transpose() + B_ * Q_ * B_.transpose();

  J_ = F_ * J_;
}

} // namespace lidar_localization