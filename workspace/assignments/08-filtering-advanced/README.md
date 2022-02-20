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
  YPoseVel_.segment<3>(3) = w_R_b.transpose() * vel_ - v_b;

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
#### 路段截取： 
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
