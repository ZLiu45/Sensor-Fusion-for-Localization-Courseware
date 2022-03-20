# Multi-Sensor Fusion for Localization & Mapping: Sliding Window -- 多传感器融合定位与建图: 基于图优化的定位方法

深蓝学院, 多传感器融合定位与建图, 第10章Graph Optimization for Localization through Sliding Window代码框架.

---

## Overview

本作业旨在加深对**基于图优化的定位, 滑动窗口, 方法**的理解.

补全基于滑动窗口的融合定位方法的实现, 并分别与:

* 不加融合
* EKF融合

的效果做对比.V

---

### 补全代码，且功能正常

* **IMU Pre-Integration**

* **Ceres Factors**
    * **Map Matching / GNSS Position Prior** [here](https://github.com/ZLiu45/Sensor-Fusion-for-Localization-Courseware/blob/zliu/sliding_window/workspace/assignments/10-sliding-window/src/lidar_localization/include/lidar_localization/models/sliding_window/factors/factor_prvag_map_matching_pose.hpp)
    ```
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
    ```
    * **IMU Pre-Integration** 
    * **Lidar Odometry** 
    * **Sliding Window Marginalization** 

* **Module Hyper Params.**
    * **Sliding Window Config** [here](src/lidar_localization/config/matching/sliding_window.yaml)

### 实现功能的基础上，性能在部分路段比EKF有改善。


### 不同窗口长度下的融合结果，并对效果及原因做对比分析

