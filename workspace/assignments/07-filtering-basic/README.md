# Multi-Sensor Fusion for Localization & Mapping: Filtering Basic -- 多传感器融合定位与建图: 基于滤波的融合方法I

深蓝学院, 多传感器融合定位与建图, 第7章Filtering Basic代码框架.

---

## Overview

本作业旨在加深对**基于滤波的融合方法**的理解.

### 代码实现
1. GetVelocityDelta():
```
linear_acc_mid = 0.5 * (linear_acc_curr + linear_acc_prev) - accl_bias_;
velocity_delta = T * 0.5 *
               (GetUnbiasedLinearAcc(linear_acc_curr, R_curr) +
                GetUnbiasedLinearAcc(linear_acc_prev, R_prev));
```
2. GetOdomEstimation(): 
```
// get deltas:
Eigen::Vector3d angular_delta = Eigen::Vector3d::Zero();
if (!GetAngularDelta(1, 0, angular_delta, angular_vel_mid)) {
LOG(FATAL) << "update to get angular delta";
return;
}
// update orientation:
Eigen::Matrix3d R_prev, R_curr;
UpdateOrientation(angular_delta, R_curr, R_prev);

// get velocity delta:
Eigen::Vector3d velocity_delta = Eigen::Vector3d::Zero();
double delta_t = 0.0;
if (!GetVelocityDelta(1, 0, R_curr, R_prev, delta_t, velocity_delta,
                    linear_acc_mid)) {
LOG(FATAL) << "update to get velocity delta";
return;
}
// save mid-value unbiased linear acc for error-state update:
// update position:
UpdatePosition(delta_t, velocity_delta);
```
3. SetProcessEquation(): 
```
void ErrorStateKalmanFilter::SetProcessEquation(const Eigen::Matrix3d &C_nb,
                                                const Eigen::Vector3d &f_n,
                                                const Eigen::Vector3d &w_b) {
  // TODO: set process / system equation:
  // a. set process equation for delta vel:

  F_.block<3, 3>(kIndexErrorVel, kIndexErrorOri) = -C_nb * skewMatrix(f_n);
  F_.block<3, 3>(kIndexErrorVel, kIndexErrorAccel) = -C_nb;

  B_.block<3, 3>(kIndexErrorVel, kIndexNoiseAccel) = C_nb;

  // b. set process equation for delta ori:
  F_.block<3, 3>(kIndexErrorOri, kIndexErrorOri) = -skewMatrix(w_b);
}
```
4. UpdateErrorEstimation(): 
```
void ErrorStateKalmanFilter::UpdateErrorEstimation(
    const double &T, const Eigen::Vector3d &linear_acc_mid,
    const Eigen::Vector3d &angular_vel_mid) {
  MatrixF F_1st;
  MatrixF F_2nd;
  // TODO: update process equation:
  UpdateProcessEquation(linear_acc_mid, angular_vel_mid);
  // TODO: get discretized process equations:
  F_1st = F_ * T;
  F_2nd = 0.5 * F_ * F_ * T * T;  // 2nd order taylor expansion
  MatrixF F_k = MatrixF::Identity() + F_1st + F_2nd;
  MatrixB B_k = B_;
  double sqrt_t = std::sqrt(T);
  B_k.block<3, 3>(kIndexErrorVel, kIndexNoiseAccel) *= T;
  B_k.block<3, 3>(kIndexErrorOri, kIndexNoiseGyro) *= T;
  B_k.block<3, 3>(kIndexErrorAccel, kIndexNoiseBiasAccel) *= sqrt_t;
  B_k.block<3, 3>(kIndexErrorGyro, kIndexNoiseBiasGyro) *= sqrt_t;
  // TODO: perform Kalman prediction
  X_ = F_k * X_;
  P_ = F_k * P_ * F_k.transpose() + B_k * Q_ * B_k.transpose(); }
```
5. CorrectErrorEstimation && CorrectErrorEstimationPose: 
```
void ErrorStateKalmanFilter::CorrectErrorEstimationPose(
    const Eigen::Matrix4d &T_nb, Eigen::VectorXd &Y, Eigen::MatrixXd &G,
    Eigen::MatrixXd &K) {
  //
  // TODO: set measurement:
  YPose_.head(3) = pose_.block<3, 1>(0, 3) - T_nb.block<3, 1>(0, 3);
  Eigen::Matrix3d delta_R =
      T_nb.block<3, 3>(0, 0).transpose() * pose_.block<3, 3>(0, 0);
  YPose_.tail(3) = Sophus::SO3d::vee(delta_R - Eigen::Matrix3d::Identity());
  // TODO: set measurement equation:
  G = GPose_;
  Y = G * X_;
  // TODO: set Kalman gain:
  MatrixRPose S = G * P_ * G.transpose() + RPose_;
  K = P_ * G.transpose() * S.inverse();
}

/**
 * @brief  correct error estimation
 * @param  measurement_type, measurement type
 * @param  measurement, input measurement
 * @return void
 */
void ErrorStateKalmanFilter::CorrectErrorEstimation(
    const MeasurementType &measurement_type, const Measurement &measurement) {
  //
  // TODO: understand ESKF correct workflow
  //
  Eigen::VectorXd Y;
  Eigen::MatrixXd G, K;
  switch (measurement_type) {
  case MeasurementType::POSE:
    CorrectErrorEstimationPose(measurement.T_nb, Y, G, K);
    break;
  default:
    break;
  }

  // TODO: perform Kalman correct:
  // LOG(INFO) << "K: \n" << K;
  P_ = ((MatrixP::Identity() - K * G) * P_).eval();
  const VectorYPose deltaY = YPose_ - Y;
  X_ = X_ + K * deltaY;
```
6. EliminateError(): 
```
void ErrorStateKalmanFilter::EliminateError(void) {
  //
  // TODO: correct state estimation using the state of ESKF
  //
  // a. position:
  pose_.block<3, 1>(0, 3) = pose_.block<3, 1>(0, 3) - X_.head(3);
  // b. velocity:
  vel_ = vel_ - X_.segment<3>(kIndexErrorVel);
  // c. orientation:
  Eigen::Quaterniond dq;
  dq.vec() << -0.5 * X_.segment<3>(kIndexErrorOri);
  dq.w() = 1.0;
  Eigen::Quaterniond curr_q(pose_.block<3, 3>(0, 0));
  curr_q = curr_q * dq;
  curr_q.normalize();

  pose_.block<3, 3>(0, 0) = curr_q.toRotationMatrix();

  // d. gyro bias:
  // if (IsCovStable(kIndexErrorGyro)) {
  gyro_bias_ -= X_.block<3, 1>(kIndexErrorGyro, 0);

  // e. accel bias:
  // if (IsCovStable(kIndexErrorAccel)) {
  accl_bias_ -= X_.block<3, 1>(kIndexErrorAccel, 0);
  // }
```
Result: 
![image](https://user-images.githubusercontent.com/11698181/153745043-b06038d5-9c8c-447c-aaa9-4dc60df493b7.png)

### 调试参数
#### Before tuning parameters: 
comparing ground_truth vs laser.txt vs Comparing ground_truth vs fused: 
##### using laser: 
![image](https://user-images.githubusercontent.com/11698181/153746111-3ec04de5-a8fa-4592-83c8-eb71fb557424.png)
##### using fused: 
![image](https://user-images.githubusercontent.com/11698181/153746122-a900d334-4c26-4858-9084-66862255667e.png)

#### fused data is very tight vs laser data. We should increase the noise of measurements and decrease the sigma of imu process noise.
With the help from IMU, the trajectory is much smoother and has better performance in terms of max error and min errors. 
![image](https://user-images.githubusercontent.com/11698181/153746535-b4e1d0f6-b2cc-4cc9-9a4f-1433a65dbbec.png)
##### using laser: 
![image](https://user-images.githubusercontent.com/11698181/153746575-63c4fc17-c2a6-4dd4-8783-7f3928c27dd7.png)
##### using fused: 
![image](https://user-images.githubusercontent.com/11698181/153746596-76fa6495-1099-468a-bb43-cfd7ad2a3828.png)


### 不考虑随机游走模型时，工程实现。
```
b_a_{k+1} = b_a_{k}
b_g_{k+1} = b_g_{k}
```
The bias doesn't have to update. 
```
void ErrorStateKalmanFilter::UpdateErrorEstimation(
    const double &T, const Eigen::Vector3d &linear_acc_mid,
    const Eigen::Vector3d &angular_vel_mid) {
  MatrixF F_1st;
  MatrixF F_2nd;
  // TODO: update process equation:
  UpdateProcessEquation(linear_acc_mid, angular_vel_mid);
  // TODO: get discretized process equations:
  F_1st = F_ * T;
  F_2nd = 0.5 * F_ * F_ * T * T;  // 2nd order taylor expansion
  MatrixF F_k = MatrixF::Identity() + F_1st + F_2nd;
  MatrixB B_k = B_;
  double sqrt_t = std::sqrt(T);
  B_k.block<3, 3>(kIndexErrorVel, kIndexNoiseAccel) *= T;
  B_k.block<3, 3>(kIndexErrorOri, kIndexNoiseGyro) *= T;
  B_k.block<3, 3>(kIndexErrorAccel, kIndexNoiseBiasAccel) *= 0;
  B_k.block<3, 3>(kIndexErrorGyro, kIndexNoiseBiasGyro) *= 0;
  // TODO: perform Kalman prediction
  X_ = F_k * X_;
  P_ = F_k * P_ * F_k.transpose() + B_k * Q_ * B_k.transpose(); }
```
Without IMU bias update, the fused result is worse than the previous one in terms of maximum error and std. But the overall performance was not impact too much in terms of mean and rmse. This is because the imu propagation was only used for a short segment so most of the cases, laser alignment can cover the errors from the IMU biases. However, under some hard cases, the bias will have a large impact, which will cause a large error. 
##### laser: 
![image](https://user-images.githubusercontent.com/11698181/153747026-96ee76fb-7503-4062-91e6-61e8ded5e29e.png)
##### fusion 
![image](https://user-images.githubusercontent.com/11698181/153747038-89e92f3c-1a58-47db-94ba-227fcba88f6c.png)

### 不同噪声设置情况下的结果对比

