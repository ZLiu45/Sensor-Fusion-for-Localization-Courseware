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
    * **IMU Pre-Integration** [here](https://github.com/ZLiu45/Sensor-Fusion-for-Localization-Courseware/blob/zliu/sliding_window/workspace/assignments/10-sliding-window/src/lidar_localization/include/lidar_localization/models/sliding_window/factors/factor_prvag_imu_pre_integration.hpp)
    ```
    const Eigen::Vector3d &alpha_ij = m_.block<3, 1>(INDEX_P, 0);
    const Eigen::Vector3d &theta_ij = m_.block<3, 1>(INDEX_R, 0);
    const Eigen::Vector3d &beta_ij = m_.block<3, 1>(INDEX_V, 0);
    const Sophus::SO3d ori_ij_est = Sophus::SO3d::exp(theta_ij);
    //
    // get square root of information matrix:
    Matrix15d sqrt_info =
        Eigen::LLT<Eigen::Matrix<double, 15, 15>>(I_).matrixL().transpose();

    //
    // compute residual:
    //

    const Sophus::SO3d i_ori_world = ori_i.inverse();
    const Eigen::Vector3d world_p_ij =
        (pos_j - pos_i - vel_i * T_ + 0.5 * g_ * T_ * T_);
    const Eigen::Vector3d world_v_ij = (vel_j - vel_i + g_ * T_);
    const Sophus::SO3d ori_ij = (ori_i.inverse() * ori_j);
    static const Eigen::Matrix3d kI33 = Eigen::Matrix3d::Identity();

    Eigen::Map<Eigen::Matrix<double, 15, 1>> res(residuals);
    res.segment<3>(INDEX_P) = i_ori_world * world_p_ij - alpha_ij;
    res.segment<3>(INDEX_R) = (ori_ij_est.inverse() * ori_ij).log();
    res.segment<3>(INDEX_V) = i_ori_world * world_v_ij - beta_ij;
    res.segment<3>(INDEX_A) = b_a_j - b_a_i;
    res.segment<3>(INDEX_G) = b_g_j - b_g_i;
    //
    // compute jacobians:
    //
    if (jacobians) {
      // compute shared intermediate results:
      Eigen::Matrix3d J_r_inv = JacobianRInv(res.segment<3>(INDEX_R));
      if (jacobians[0]) {
        Eigen::Map<Eigen::Matrix<double, 15, 15, Eigen::RowMajor>> jacobian_i(
            jacobians[0]);
        jacobian_i.setZero();
        // a. residual, position:
        jacobian_i.block<3, 3>(INDEX_P, INDEX_P) = -i_ori_world.matrix();
        jacobian_i.block<3, 3>(INDEX_P, INDEX_R) =
            Sophus::SO3d::hat(i_ori_world * world_p_ij).matrix();
        jacobian_i.block<3, 3>(INDEX_P, INDEX_V) = -i_ori_world.matrix() * T_;
        jacobian_i.block<3, 3>(INDEX_P, INDEX_A) =
            -J_.block<3, 3>(INDEX_P, INDEX_A);
        jacobian_i.block<3, 3>(INDEX_P, INDEX_G) =
            -J_.block<3, 3>(INDEX_P, INDEX_G);

        // b. residual, orientation:
        jacobian_i.block<3, 3>(INDEX_R, INDEX_R) =
            -J_r_inv * ori_ij.matrix().transpose();
        jacobian_i.block<3, 3>(INDEX_R, INDEX_G) =
            -J_r_inv *
            (Sophus::SO3d::exp(res.block<3, 1>(INDEX_R, 0)))
                .matrix()
                .inverse() *
            J_.block<3, 3>(INDEX_R, INDEX_G);

        // c. residual, velocity:
        jacobian_i.block<3, 3>(INDEX_V, INDEX_R) =
            Sophus::SO3d::hat(i_ori_world * world_v_ij);
        jacobian_i.block<3, 3>(INDEX_V, INDEX_V) = -i_ori_world.matrix();
        jacobian_i.block<3, 3>(INDEX_V, INDEX_A) =
            -J_.block<3, 3>(INDEX_V, INDEX_A);
        jacobian_i.block<3, 3>(INDEX_V, INDEX_G) =
            -J_.block<3, 3>(INDEX_V, INDEX_G);
        // d. residual, bias accel:
        jacobian_i.block<3, 3>(INDEX_A, INDEX_A) = -kI33;
        // d. residual, bias accel:
        jacobian_i.block<3, 3>(INDEX_G, INDEX_G) = -kI33;

        jacobian_i = sqrt_info * jacobian_i;
      }

      if (jacobians[1]) {
        Eigen::Map<Eigen::Matrix<double, 15, 15, Eigen::RowMajor>> jacobian_j(
            jacobians[1]);
        jacobian_j.setZero();
        // a. residual, position:
        jacobian_j.block<3, 3>(INDEX_P, INDEX_P) = i_ori_world.matrix();
        // b. residual, orientation:
        jacobian_j.block<3, 3>(INDEX_R, INDEX_R) = J_r_inv;
        // c. residual, velocity:
        jacobian_j.block<3, 3>(INDEX_V, INDEX_V) = i_ori_world.matrix();
        // d. residual, bias accel:
        jacobian_j.block<3, 3>(INDEX_A, INDEX_A) = kI33;
        // d. residual, bias accel:
        jacobian_j.block<3, 3>(INDEX_G, INDEX_G) = kI33;

        jacobian_j = sqrt_info * jacobian_j;
      }
    }

    //
    // orrect residual by square root of information matrix:
    //
    res = sqrt_info * res;
    ```
    * **Lidar Odometry** [here](https://github.com/ZLiu45/Sensor-Fusion-for-Localization-Courseware/blob/zliu/sliding_window/workspace/assignments/10-sliding-window/src/lidar_localization/include/lidar_localization/models/sliding_window/factors/factor_prvag_relative_pose.hpp)
    ```
        //
    Eigen::LLT<Eigen::Matrix<double, 6, 6>> LowerI(I_);
    Eigen::Matrix<double, 6, 6> sqrt_info = LowerI.matrixL().transpose();

    // compute residual:
    //
    Eigen::Map<Eigen::Matrix<double, 6, 1>> res(residuals);

    const Sophus::SO3d i_ori_world = ori_i.inverse();
    const Sophus::SO3d ori_ij_est = i_ori_world * ori_j;
    const Eigen::Vector3d pos_ij_est = i_ori_world * (pos_j - pos_i);

    res.segment<3>(INDEX_P) = pos_ij_est - pos_ij;
    res.segment<3>(INDEX_R) = (ori_ij.inverse() * ori_ij_est).log();

    //
    // compute jacobians:
    //
    if (jacobians) {
      // compute shared intermediate results:
      const Eigen::Matrix3d J_r_inv = JacobianRInv(res.segment<3>(INDEX_R));

      if (jacobians[0]) {
        // implement computing:
        Eigen::Map<Eigen::Matrix<double, 6, 15, Eigen::RowMajor>> jac_i(
            jacobians[0]);
        jac_i.setZero();
        jac_i.block<3, 3>(INDEX_P, INDEX_P) = -i_ori_world.matrix();
        jac_i.block<3, 3>(INDEX_P, INDEX_R) =
            Sophus::SO3d::hat(pos_ij_est).matrix();
        jac_i.block<3, 3>(INDEX_R, INDEX_R) =
            -J_r_inv * ori_ij_est.matrix().transpose();

        jac_i = sqrt_info * jac_i;
      }

      if (jacobians[1]) {
        // implement computing:
        Eigen::Map<Eigen::Matrix<double, 6, 15, Eigen::RowMajor>> jac_j(
            jacobians[1]);
        jac_j.setZero();
        jac_j.block<3, 3>(INDEX_P, INDEX_P) = i_ori_world.matrix();
        jac_j.block<3, 3>(INDEX_R, INDEX_R) = J_r_inv;

        jac_j = sqrt_info * jac_j;
      }
    }

    //
    // correct residual by square root of information matrix:
    //
    res = sqrt_info * res;
    ```
    * **Sliding Window Marginalization** [here](https://github.com/ZLiu45/Sensor-Fusion-for-Localization-Courseware/blob/zliu/sliding_window/workspace/assignments/10-sliding-window/src/lidar_localization/include/lidar_localization/models/sliding_window/factors/factor_prvag_marginalization.hpp)
    #### SetResMapMatchingPose
    ```
    void SetResMapMatchingPose(const ceres::CostFunction *residual,
                             const std::vector<double *> &parameter_blocks) {
    // init:
    ResidualBlockInfo res_map_matching_pose(residual, parameter_blocks);
    Eigen::VectorXd residuals;
    std::vector<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        jacobians;

    // compute:
    Evaluate(res_map_matching_pose, residuals, jacobians);
    const Eigen::MatrixXd &J_m = jacobians.at(0);
    //
    // Update H:
    //
    // a. H_mm:
    H_.block<15, 15>(INDEX_M, INDEX_M) += J_m.transpose() * J_m;

    //
    // Update b:
    //
    // a. b_m:
    b_.block<15, 1>(INDEX_M, 0) += J_m.transpose() * residuals;
  }
    ```
    #### SetResRelativePose
    ```
    void SetResRelativePose(const ceres::CostFunction *residual,
                          const std::vector<double *> &parameter_blocks) {
    // init:
    ResidualBlockInfo res_relative_pose(residual, parameter_blocks);
    Eigen::VectorXd residuals;
    std::vector<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        jacobians;

    // compute:
    Evaluate(res_relative_pose, residuals, jacobians);
    const Eigen::MatrixXd &J_m = jacobians.at(0);
    const Eigen::MatrixXd &J_r = jacobians.at(1);

    //
    // Update H:
    const Eigen::MatrixXd H_mm = J_m.transpose() * J_m;
    H_.block<15, 15>(INDEX_M, INDEX_M) += H_mm;
    // b. H_mr:
    const Eigen::MatrixXd H_mr = J_m.transpose() * J_r;
    H_.block<15, 15>(INDEX_M, INDEX_R) += H_mr;
    // c. H_rm:
    H_.block<15, 15>(INDEX_R, INDEX_M) += H_mr.transpose();
    // d. H_rr:
    const Eigen::MatrixXd H_rr = J_r.transpose() * J_r;
    H_.block<15, 15>(INDEX_R, INDEX_R) += H_rr;

    //
    // Update b:
    //
    // a. b_m:
    b_.block<15, 1>(INDEX_M, 0) += J_m.transpose() * residuals;
    // a. b_r:
    b_.block<15, 1>(INDEX_R, 0) += J_r.transpose() * residuals;
  }
    ```
   #### SetResImuPrevIntegration
   ```
   void SetResIMUPreIntegration(const ceres::CostFunction *residual,
                               const std::vector<double *> &parameter_blocks) {
    // init:
    ResidualBlockInfo res_imu_pre_integration(residual, parameter_blocks);
    Eigen::VectorXd residuals;
    std::vector<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        jacobians;

    // compute:
    Evaluate(res_imu_pre_integration, residuals, jacobians);
    const Eigen::MatrixXd &J_m = jacobians.at(0);
    const Eigen::MatrixXd &J_r = jacobians.at(1);

    //
    // Update H:
    //
    // a. H_mm:
    const Eigen::MatrixXd H_mm = J_m.transpose() * J_m;
    H_.block<15, 15>(INDEX_M, INDEX_M) += H_mm;
    // b. H_mr:
    const Eigen::MatrixXd H_mr = J_m.transpose() * J_r;
    H_.block<15, 15>(INDEX_M, INDEX_R) += H_mr;
    // c. H_rm:
    H_.block<15, 15>(INDEX_R, INDEX_M) += H_mr.transpose();
    // d. H_rr:
    const Eigen::MatrixXd H_rr = J_r.transpose() * J_r;
    H_.block<15, 15>(INDEX_R, INDEX_R) += H_rr;

    //
    // Update b:
    //
    // a. b_m:
    b_.block<15, 1>(INDEX_M, 0) += J_m.transpose() * residuals;
    // a. b_r:
    b_.block<15, 1>(INDEX_R, 0) += J_r.transpose() * residuals;
  }
   ```
   #### marginalize: 
   ```
   void Marginalize(const double *raw_param_r_0) {
    x_0_ = Eigen::Map<const Eigen::Matrix<double, 15, 1>>(raw_param_r_0);

    // implement marginalization logic ???
    const Eigen::MatrixXd &H_mm =
        0.5 * (H_.block<15, 15>(INDEX_M, INDEX_M) +
               H_.block<15, 15>(INDEX_M, INDEX_M).transpose());
    const Eigen::MatrixXd &H_rr = H_.block<15, 15>(INDEX_R, INDEX_R);
    const Eigen::MatrixXd &H_rm = H_.block<15, 15>(INDEX_R, INDEX_M);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(H_mm);
    Eigen::MatrixXd H_mm_inv =
        Eigen::VectorXd((saes.eigenvalues().array() > 1e-4)
                            .select(saes.eigenvalues().array().inverse(), 0))
            .asDiagonal() *
        saes.eigenvectors().transpose();
    Eigen::MatrixXd H_rm_mm_inv = H_rm * H_mm_inv;

    Eigen::MatrixXd H_rr_prime = H_rr - H_rm_mm_inv * H_rm.transpose();
    Eigen::VectorXd b_r_prime =
        b_.block<15, 1>(INDEX_R, 0) - H_rm_mm_inv * b_.block<15, 1>(INDEX_M, 0);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(H_rr_prime);
    Eigen::VectorXd S =
        Eigen::VectorXd((saes2.eigenvalues().array() > 1e-5)
                            .select(saes2.eigenvalues().array(), 0));
    Eigen::VectorXd S_inv =
        Eigen::VectorXd((saes2.eigenvalues().array() > 1e-5)
                            .select(saes2.eigenvalues().array().inverse(), 0));
    Eigen::VectorXd S_sqrt = S.cwiseSqrt();
    Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();

    J_ = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
    e_ = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose() * b_r_prime;
  }
   ```
   #### ceres sliding window: [here](https://github.com/ZLiu45/Sensor-Fusion-for-Localization-Courseware/blob/zliu/sliding_window/workspace/assignments/10-sliding-window/src/lidar_localization/src/models/sliding_window/ceres_sliding_window.cpp)
   ```
   if (!residual_blocks_.map_matching_pose.empty() &&
        !residual_blocks_.relative_pose.empty() &&
        !residual_blocks_.imu_pre_integration.empty()) {
      auto &key_frame_m = optimized_key_frames_.at(N - kWindowSize - 1);
      auto &key_frame_r = optimized_key_frames_.at(N - kWindowSize - 0);
      const ceres::CostFunction *factor_map_matching_pose =
          GetResMapMatchingPose(residual_blocks_.map_matching_pose.front());
      const ceres::CostFunction *factor_relative_pose =
          GetResRelativePose(residual_blocks_.relative_pose.front());
      const ceres::CostFunction *factor_imu_pre_integration =
          GetResIMUPreIntegration(residual_blocks_.imu_pre_integration.front());

      sliding_window::FactorPRVAGMarginalization *factor_marginalization =
          new sliding_window::FactorPRVAGMarginalization();
      factor_marginalization->SetResMapMatchingPose(
          factor_map_matching_pose, std::vector<double *>{key_frame_m.prvag});

      factor_marginalization->SetResRelativePose(
          factor_relative_pose,
          std::vector<double *>{key_frame_m.prvag, key_frame_r.prvag});

      factor_marginalization->SetResIMUPreIntegration(
          factor_imu_pre_integration,
          std::vector<double *>{key_frame_m.prvag, key_frame_r.prvag});

      factor_marginalization->Marginalize(key_frame_r.prvag);
      // add marginalization factor into sliding window
      problem.AddResidualBlock(factor_marginalization, NULL, key_frame_r.prvag);

      residual_blocks_.map_matching_pose.pop_front();
      residual_blocks_.relative_pose.pop_front();
      residual_blocks_.imu_pre_integration.pop_front();
    }

    LOG(INFO) << "set marginalization factor";
    // b.2. map matching pose constraint:
    if (!residual_blocks_.map_matching_pose.empty()) {
      for (const auto &residual_map_matching_pose :
           residual_blocks_.map_matching_pose) {
        auto &key_frame =
            optimized_key_frames_.at(residual_map_matching_pose.param_index);

        sliding_window::FactorPRVAGMapMatchingPose *factor_map_matching_pose =
            GetResMapMatchingPose(residual_map_matching_pose);

        // add map matching factor into sliding window
        problem.AddResidualBlock(factor_map_matching_pose,
                                 config_.loss_function_ptr.get(),
                                 key_frame.prvag);
      }
    }
    LOG(INFO) << "set map matching factor";

    // b.3. relative pose constraint:
    if (!residual_blocks_.relative_pose.empty()) {
      for (const auto &residual_relative_pose :
           residual_blocks_.relative_pose) {
        auto &key_frame_i =
            optimized_key_frames_.at(residual_relative_pose.param_index_i);
        auto &key_frame_j =
            optimized_key_frames_.at(residual_relative_pose.param_index_j);

        sliding_window::FactorPRVAGRelativePose *factor_relative_pose =
            GetResRelativePose(residual_relative_pose);

        // add relative pose factor into sliding window
        problem.AddResidualBlock(factor_relative_pose,
                                 config_.loss_function_ptr.get(),
                                 key_frame_i.prvag, key_frame_j.prvag);
      }
    }
    LOG(INFO) << "set relative pose factor";

    // TODO: b.4. IMU pre-integration constraint
    if (!residual_blocks_.imu_pre_integration.empty()) {
      for (const auto &residual_imu_pre_integration :
           residual_blocks_.imu_pre_integration) {
        auto &key_frame_i = optimized_key_frames_.at(
            residual_imu_pre_integration.param_index_i);
        auto &key_frame_j = optimized_key_frames_.at(
            residual_imu_pre_integration.param_index_j);

        sliding_window::FactorPRVAGIMUPreIntegration
            *factor_imu_pre_integration =
                GetResIMUPreIntegration(residual_imu_pre_integration);

        // TODO: add IMU factor into sliding window
        problem.AddResidualBlock(factor_imu_pre_integration,
                                 nullptr,
                                 key_frame_i.prvag, key_frame_j.prvag);
      }
   ```
   
   #### sliding window [here](https://github.com/ZLiu45/Sensor-Fusion-for-Localization-Courseware/blob/zliu/sliding_window/workspace/assignments/10-sliding-window/src/lidar_localization/src/matching/back_end/sliding_window.cpp)
   ```
   bool SlidingWindow::Update(void) {
  static KeyFrame last_key_frame_ = current_key_frame_;

  // add node for new key frame pose:
  //
  // fix the pose of the first key frame for lidar only mapping:
  if (sliding_window_ptr_->GetNumParamBlocks() == 0) {
    // TODO: add init key frame
    sliding_window_ptr_->AddPRVAGParam(current_key_frame_, true);
  } else {
    // TODO: add current key frame
    sliding_window_ptr_->AddPRVAGParam(current_key_frame_, false);
  }

  // get num. of vertices:
  const int N = sliding_window_ptr_->GetNumParamBlocks();

  // get param block ID, current:
  const int param_index_j = N - 1;

  //
  // add unary constraints:
  //
  //
  // a. map matching / GNSS position:
  //
  if (N > 0 && measurement_config_.source.map_matching) {
    // get prior position measurement:
    Eigen::Matrix4d prior_pose = current_map_matching_pose_.pose.cast<double>();

    // TODO: add constraint, GNSS position:
    sliding_window_ptr_->AddPRVAGMapMatchingPoseFactor(
        param_index_j, prior_pose, measurement_config_.noise.map_matching);
  }

  //
  // add binary constraints:
  //
  if (N > 1) {
    // get param block ID, previous:
    const int param_index_i = N - 2;

    //
    // a. lidar frontend:
    //
    // get relative pose measurement:
    Eigen::Matrix4d relative_pose =
        (last_key_frame_.pose.inverse() * current_key_frame_.pose)
            .cast<double>();
    // TODO: add constraint, lidar frontend / loop closure detection:
    sliding_window_ptr_->AddPRVAGRelativePoseFactor(
        param_index_i, param_index_j, relative_pose,
        measurement_config_.noise.lidar_odometry);
    //
    // b. IMU pre-integration:
    //
    if (measurement_config_.source.imu_pre_integration) {
      // TODO: add constraint, IMU pre-integraion:
      sliding_window_ptr_->AddPRVAGIMUPreIntegrationFactor(
          param_index_i, param_index_j, imu_pre_integration_);
    }
  }

  // move forward:
  last_key_frame_ = current_key_frame_;

   ```
* **Module Hyper Params.**
    * **Sliding Window Config** [here](src/lidar_localization/config/matching/sliding_window.yaml)

