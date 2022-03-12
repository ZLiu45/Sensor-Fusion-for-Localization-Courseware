# Multi-Sensor Fusion for Localization & Mapping: Graph Optimization

## code: 
### imu-pre-integrator: 
```
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
```

```
// update F and B 
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
```
### edge_prvag_imu_pre_integration.hpp 
```
    // update pre-integration measurement caused by bias change:
    if (v0->isUpdated()) {
      Eigen::Vector3d d_b_a_i, d_b_g_i;
      v0->getDeltaBiases(d_b_a_i, d_b_g_i);
      updateMeasurement(d_b_a_i, d_b_g_i);
    }

    Eigen::Vector3d rel_p = _measurement.segment<3>(INDEX_P);
    Eigen::Vector3d rel_v = _measurement.segment<3>(INDEX_V);
    Eigen::Vector3d rel_theta = _measurement.segment<3>(INDEX_R);
    const Sophus::SO3d rel_ori_ij = Sophus::SO3d::exp(rel_theta);

    // //
    // // compute error:
    // //
    _error.block<3, 1>(INDEX_P, 0) =
        ori_i.inverse() * (0.5 * g_ * T_ * T_ + pos_j - pos_i - vel_i * T_) -
        rel_p;
	LOG(INFO) << "position error: " << _error.block<3, 1>(INDEX_P, 0).transpose();
    _error.block<3, 1>(INDEX_V, 0) =
        ori_i.inverse() * (g_ * T_ + vel_j - vel_i) - rel_v;
    _error.block<3, 1>(INDEX_R, 0) =
        (Sophus::SO3d::exp(rel_theta).inverse() * ori_i.inverse() * ori_j)
            .log();
    _error.block<3, 1>(INDEX_A, 0) = b_a_j - b_a_i;
    _error.block<3, 1>(INDEX_G, 0) = b_g_j - b_g_i;
```

### vertex_prvag.hpp 
```
    _estimate.pos +=
        Eigen::Vector3d(update[PRVAG::INDEX_POS], update[PRVAG::INDEX_POS + 1],
                        update[PRVAG::INDEX_POS + 2]);
    _estimate.ori = _estimate.ori *
                    Sophus::SO3d::exp(Eigen::Vector3d(
                        update[PRVAG::INDEX_ORI], update[PRVAG::INDEX_ORI + 1],
                        update[PRVAG::INDEX_ORI + 2]));
    _estimate.vel +=
        Eigen::Vector3d(update[PRVAG::INDEX_VEL], update[PRVAG::INDEX_VEL + 1],
                        update[PRVAG::INDEX_VEL + 2]);
    Eigen::Vector3d d_b_a_i(Eigen::Vector3d(update[PRVAG::INDEX_B_A],
                                            update[PRVAG::INDEX_B_A + 1],
                                            update[PRVAG::INDEX_B_A + 2]));
    Eigen::Vector3d d_b_g_i(Eigen::Vector3d(update[PRVAG::INDEX_B_G],
                                            update[PRVAG::INDEX_B_G + 1],
                                            update[PRVAG::INDEX_B_G + 2]));
    _estimate.b_a += d_b_a_i;
    _estimate.b_g += d_b_g_i;
    updateDeltaBiases(d_b_a_i, d_b_g_i);
```
## evaluation:
### using IMU pre-integration 

#### plots From RViz: 
<img src="https://user-images.githubusercontent.com/11698181/158010558-6d9b7b04-c3ab-400c-b48a-742bb0faa225.png" width="640"/>

#### evo analaysis
<img src="https://user-images.githubusercontent.com/11698181/158030449-f33f583a-ec85-4b53-8a41-054c393d4a14.png" width="640"/>
<img src="https://user-images.githubusercontent.com/11698181/158030500-103f9545-6fe7-4a1e-b1f0-32a20bcee173.png" width="640"/>

#### evo rpe: 
| SE(3) | Umeyama alignment |
| --- | ----------- |
| max | 0.377172 |
| mean | 0.1216 |
| median | 0.105 |
| min | 0.0058 |
| rmse | 0.1452 |
| std | 0.079281 |
<img src="https://user-images.githubusercontent.com/11698181/158030656-c297613a-deb6-4664-8895-13ab7cf05978.png" width="640"/>
<img src="https://user-images.githubusercontent.com/11698181/158030690-ea865be3-e4de-4680-83f8-3b214efd0766.png" width="640"/>

#### evo ape: 
| SE(3) | Umeyama alignment |
| --- | ----------- |
| max | 12.361094 |
| mean | 2.86433 |
| median | 2.258796 |
| min | 0.096506 |
| rmse | 3.46 |
| std | 1.9409 |
<img src="https://user-images.githubusercontent.com/11698181/158030894-822b6434-24cf-4011-98eb-49616aa4bcbd.png" width="640"/>
<img src="https://user-images.githubusercontent.com/11698181/158030909-d9aa4207-6a1e-4e3a-8cc5-3411ed87c47e.png" width="640"/>

### with IMU pre-integration

#### plots From RViz: 
<img src="https://user-images.githubusercontent.com/11698181/158031441-ebdcf906-6ff0-40a7-b20c-d3ad4a9af925.png" width="640"/>

#### evo analaysis
#### evo rpe: 
| SE(3) | Umeyama alignment |
| --- | ----------- |
| max | 0.37133 |
| mean | 0.114145 |
| median | 0.107 |
| min | 0.00575 |
| rmse | 0.1286 |
| std | 0.059 |

<img src="https://user-images.githubusercontent.com/11698181/158031687-2d1f7c1e-3e8b-43c0-a34c-e01fa57ecff3.png" width="640"/>
<img src="https://user-images.githubusercontent.com/11698181/158031664-78884130-32c5-4f3b-8c4a-c3ac17121c81.png" width="640"/>

#### evo ape: 
| SE(3) | Umeyama alignment |
| --- | ----------- |
| max | 14.119667 |
| mean | 3.773 |
| median | 2.556 |
| min | 0.0404 |
| rmse | 5.0596 |
| std | 3.3669 |

<img src="https://user-images.githubusercontent.com/11698181/158031572-db567046-2d76-455a-af3b-ccaf0f062347.png" width="640"/>
<img src="https://user-images.githubusercontent.com/11698181/158031598-d7b784d7-2fbd-49ea-baef-4a844f90ec31.png" width="640"/>

## wheel encoder jacobian update:
