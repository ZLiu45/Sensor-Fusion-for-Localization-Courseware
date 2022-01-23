# Sensor Fusion: Lidar Odometry -- 多传感器融合定位与建图: 惯性导航原理

## jacobian derivation 
[Assignment_5.pdf](https://github.com/ZLiu45/Sensor-Fusion-for-Localization-Courseware/files/7919972/Assignment_5.pdf)

## calibration with autodiff 
![task2](https://user-images.githubusercontent.com/11698181/150626499-53dd45c4-493e-40bc-a988-5e4784588135.png)

### code in calibration.cpp: 
```
for (int th_mult = 2; th_mult <= 10; th_mult++) {
    std::vector<imu_tk::DataInterval> static_intervals;
    std::vector<imu_tk::TriadData_<_T>> static_samples;
    std::vector<double> acc_calib_params(9);

    //
    // TODO: implement lower triad model here
    //
    acc_calib_params[0] = init_acc_calib_.misXZ();
    acc_calib_params[1] = -init_acc_calib_.misXY();
    acc_calib_params[2] = init_acc_calib_.misYX();

    acc_calib_params[3] = init_acc_calib_.scaleX();
    acc_calib_params[4] = init_acc_calib_.scaleY();
    acc_calib_params[5] = init_acc_calib_.scaleZ();

    acc_calib_params[6] = init_acc_calib_.biasX();
    acc_calib_params[7] = init_acc_calib_.biasY();
    acc_calib_params[8] = init_acc_calib_.biasZ();
   ....
   }
```

```
acc_calib_ = CalibratedTriad_<_T>(
      //
      // TODO: implement lower triad model here
      //
      0,  0,  0, 
      min_cost_calib_params[0], min_cost_calib_params[1], min_cost_calib_params[2], 
      min_cost_calib_params[3], min_cost_calib_params[4], min_cost_calib_params[5],
      min_cost_calib_params[6], min_cost_calib_params[7], min_cost_calib_params[8]);
```

```
 CalibratedTriad_<_T2> calib_triad(
        // TODO: implement lower triad model here
        // mis_yz, mis_zy, mis_zx:
        _T2(0), _T2(0),_T2(0),
        // mis_xz, mis_xy, mis_yx:
        params[0], params[1], params[2],
        //    s_x,    s_y,    s_z:
        params[3], params[4], params[5],
        //    b_x,    b_y,    b_z:
        params[6], params[7], params[8]);
```
## Calibration with Analytical Jacobians
```
virtual bool Evaluate(_T_a const *const *params, _T_a *residuals,
                        _T_a **jacobians) const {
    // compute the residuals
    CalibratedTriad_<_T_a> calib_triad(
        // mis_yz, mis_zy, mis_zx:
        _T_a(0), _T_a(0), _T_a(0),
        // mis_xz, mis_xy, mis_yx:
        (*params)[0], (*params)[1], (*params)[2],
        //    s_x,    s_y,    s_z:
        (*params)[3], (*params)[4], (*params)[5],
        //    b_x,    b_y,    b_z:
        (*params)[6], (*params)[7], (*params)[8]);
    Eigen::Matrix<_T_a, 3, 1> raw_samp(_T_a(sample_(0)), _T_a(sample_(1)),
                                       _T_a(sample_(2)));

    Eigen::Matrix<_T_a, 3, 3> S = Eigen::Matrix<_T_a, 3, 3>::Zero();
    S(1, 0) = (*params)[0];
    S(2, 0) = (*params)[1];
    S(2, 1) = (*params)[2];
    Eigen::Matrix<_T_a, 3, 3> I = Eigen::Matrix<_T_a, 3, 3>::Identity();
    Eigen::Matrix<_T_a, 3, 3> I_minus_S = I - S;

    Eigen::Matrix<_T_a, 3, 1> K_vec((*params)[3], (*params)[4], (*params)[5]);

    Eigen::Matrix<_T_a, 3, 1> ba((*params)[6], (*params)[7], (*params)[8]);
    Eigen::Matrix<_T_a, 3, 1> unbiased_samp = raw_samp - ba;

    Eigen::Matrix<_T_a, 3, 3> K_prime;
    K_prime.setZero();
    K_prime.diagonal() << 1.0 / K_vec.x(), 1.0 / K_vec.y(), 1.0 / K_vec.z();
    Eigen::Matrix<_T_a, 3, 1> true_samp = I_minus_S * K_prime * unbiased_samp;
    residuals[0] = 0.5 * (_T_a(g_mag_) - true_samp.norm());

    // compute jacobians
    if (jacobians) {
      if (jacobians[0]) {

        Eigen::Map<Eigen::Matrix<_T_a, 1, 9>> jac(jacobians[0]);
        jac.setZero();
        // jacobian over mis_alignment
        Eigen::Matrix<_T_a, 3, 3> jac_a_S;
        jac_a_S.setZero();
        jac_a_S(1, 0) = -unbiased_samp.x() / K_vec.x();
        jac_a_S(2, 1) = -unbiased_samp.x() / K_vec.x();
        jac_a_S(2, 2) = -unbiased_samp.y() / K_vec.y();
        jac.leftCols(3) = -true_samp.transpose() * jac_a_S;

        // jacobian over Scalar
        Eigen::Matrix<_T_a, 3, 3> jac_a_K;
        jac_a_K.setZero();
        jac_a_K(0, 0) = -unbiased_samp.x()/(K_vec.x() * K_vec.x());
        jac_a_K(1, 0) = -S(1, 0) * jac_a_K(0, 0); 
        jac_a_K(1, 1) = -unbiased_samp.y()/(K_vec.y() * K_vec.y());
        jac_a_K(2, 0) = -S(2, 0) * jac_a_K(0, 0);  
        jac_a_K(2, 1) = -S(2, 1) * jac_a_K(1, 1); 
        jac_a_K(2, 2) = -unbiased_samp.z()/(K_vec.z() * K_vec.z());
        jac.block(0, 3, 1, 3) = -true_samp.transpose() * jac_a_K; 

        // jacobian over bias
        Eigen::Matrix<_T_a, 3, 3> jac_a_ba;
        jac_a_ba.setZero(); 
        jac_a_ba.diagonal() << -1/K_vec.x(), -1/K_vec.y(), -1/K_vec.z(); 
        jac_a_ba(1, 0) = S(1, 0)/K_vec.x(); 
        jac_a_ba(2, 0) = S(2, 0)/K_vec.x(); 
        jac_a_ba(2, 1) = S(2, 1)/K_vec.y();
        jac.rightCols(3) = -true_samp.transpose() * jac_a_ba;
      }
    }
    return true;
  }
};
```
![image](https://user-images.githubusercontent.com/11698181/150672367-e073e8b9-6763-40ca-ad8e-00e3f3657fc2.png)
