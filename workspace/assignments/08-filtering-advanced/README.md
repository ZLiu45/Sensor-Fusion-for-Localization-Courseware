# Multi-Sensor Fusion for Localization & Mapping: Filtering Advanced -- 多传感器融合定位与建图: 基于滤波的融合方法II

深蓝学院, 多传感器融合定位与建图, 第8章Filtering Advanced代码框架.

---

## Overview

本作业旨在加深对**基于滤波的融合方法**的理解.

本章作业要求如下: 在上一讲作业里实现的滤波方案的基础上

1. 实现**融合运动模型**的滤波方法
2. 对比加入运动模型约束前后,滤波精度的变化. 由于运动模型约束更多的是改善速度的波动,而且是y向和z向的波动,因此要求展示结果时,提供b系y向和z向速度误差的曲线与指标.

注:同样由于kitti数据集质量的问题,效果的改善不一定在所有路段都能体现,可以挑选效果好的路段重点展示。

---

### 新模型代码实现：
1. correctErrorStateEstimation()
```
  switch (measurement_type) {
  case MeasurementType::POSE: {
    CorrectErrorEstimationPose(measurement.T_nb, Y, G, K);
    P_ = ((MatrixP::Identity() - K * G) * P_).eval();
    const VectorYPose deltaY = YPose_ - Y;
    X_ = X_ + K * deltaY;
    break;
  }
  case MeasurementType::POSE_VEL: {
    // TODO: register new correction logic here:
    CorrectErrorEstimationPoseVel(measurement.T_nb, measurement.v_b,
                                  measurement.w_b, Y, G, K);
    P_ = ((MatrixP::Identity() - K * G) * P_).eval();
    const VectorYPoseVel deltaY = YPoseVel_ - Y;
    X_ = X_ + K * deltaY;
    break;
  }
```
2. CorrectErrorEstimationPoseVel()
```
void ErrorStateKalmanFilter::CorrectErrorEstimationPoseVel(
    const Eigen::Matrix4d &T_nb, const Eigen::Vector3d &v_b,
    const Eigen::Vector3d &w_b, Eigen::VectorXd &Y, Eigen::MatrixXd &G,
    Eigen::MatrixXd &K) {
  //
  Eigen::Matrix3d w_R_b = pose_.block<3, 3>(0, 0);
  YPoseVel_.setZero();
  YPoseVel_.head(3) = pose_.block<3, 1>(0, 3) - T_nb.block<3, 1>(0, 3);
  YPoseVel_.segment<3>(3) = w_R_b.transpose() * vel_ - Eigen::Vector3d(v_b.x(), 0.0, 0.0);

  Eigen::Matrix3d delta_R =
      T_nb.block<3, 3>(0, 0).transpose() * pose_.block<3, 3>(0, 0);
  YPoseVel_.tail<3>(3) =
      Sophus::SO3d::vee(delta_R - Eigen::Matrix3d::Identity());

  // YPoseVel_(kDimMeasurementPose) = 0.0;
  // TODO: set measurement equation:
  G = GPoseVel_;
  G.block<3, 3>(3, 3) = w_R_b.transpose();
  G.block<3, 3>(3, 6) = skewMatrix(v_b);
  Y = G * X_;
  // TODO: set Kalman gain:
  MatrixRPoseVel S = G * P_ * G.transpose() + RPoseVel_;
  K = P_ * G.transpose() * S.inverse();
}
```

### 结果评价：
#### 整体结果： 
##### Not using velocity measurement:![evo_laser_v0_0](https://user-images.githubusercontent.com/11698181/154829537-52ed68e3-4347-406f-ad4c-f91d20d62299.png)

1. laser: 
![evo_![evo_laser2](https://user-images.githubusercontent.com/11698181/154829455-62f3521d-f80a-4513-92e1-ec12e3d636d0.png)
laser](https://user-images.githubusercontent.com/11698181/154829452-70a6a24b-794e-4b99-9fa9-de1d0d5cc6e1.png)
```
max: 13.22
mean: 8.65
median: 8.63
min: 6.26 
rmse: 8.67 
```
2. fused: 
![evo_fused](https://user-images.githubusercontent.com/11698181/154829505-34ea58f2-a10b-4121-b2c6-d07fab1ff300.png)
![evo_fused2](https://user-images.githubusercontent.com/11698181/154829507-59bdd035-ea03-4b9a-99e0-2164f5b976eb.png)
```
max: 13.03
mean: 8.66
median: 8.63
min: 6.51
rmse: 8.68 
```
##### using velocity measurement: 
1. laser 
![evo_laser_v0_0](https://user-images.githubusercontent.com/11698181/154829544-14080eee-5f6d-433a-a07a-8461f0e1f78c.png)
![evo_laser_v0](https://user-images.githubusercontent.com/11698181/154829548-0437b41d-da57-42ec-92c9-e3129a1240e5.png)
```
max: 12.36
mean: 7.81 
median: 7.77
min: 5.43 
rmse: 7.83 
```
2. fused: 
![evo_fused_v0_0](https://user-images.githubusercontent.com/11698181/154829582-adf2be15-3c4a-4a5b-8aad-4aabee667d16.png)
![evo_fused_v0](https://user-images.githubusercontent.com/11698181/154829584-60054da5-de84-4a3b-bc35-0e659e7c28a0.png)
```
max: 12.18
mean: 7.81
median: 7.77
min: 5.28
rmse: 7.82
```
### 仿真实现：
1. CorrectErrorEstimationPosiVel():
```
void ErrorStateKalmanFilter::CorrectErrorEstimationPosiVel(
    const Eigen::Matrix4d &T_nb, const Eigen::Vector3d &v_b,
    const Eigen::Vector3d &w_b, Eigen::VectorXd &Y, Eigen::MatrixXd &G,
    Eigen::MatrixXd &K) {
  // parse measurement:
  Eigen::Matrix3d w_R_b = pose_.block<3, 3>(0, 0);
  YPosiVel_.head(3) = pose_.block<3, 1>(0, 3) - T_nb.block<3, 1>(0, 3);
  YPosiVel_.tail(3) = w_R_b.transpose() * vel_ - v_b;
  // set measurement equation:
  G = GPosiVel_;
  G.block<3, 3>(3, 3) = w_R_b.transpose();
  G.block<3, 3>(3, 6) = skewMatrix(v_b);
  Y = G * X_;
  // set Kalman gain:
  MatrixRPosiVel S = G * P_ * G.transpose() + RPosiVel_;
  K = P_ * G.transpose() * S.inverse();
}
```
2. correctErrorStateEstimation():
