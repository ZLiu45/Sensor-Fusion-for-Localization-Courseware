# Sensor Fusion: Lidar Odometry -- 多传感器融合定位与建图: 惯性导航解算

## 中值法
```
while (imu_data_buff_.size() > 1) {
  size_t prev_index, curr_index;
  if (integration_method_ == "euler") {
      prev_index = 0; 
      curr_index = prev_index; 
  } else if (integration_method_ == "mid_point") {
      prev_index = 0; 
      curr_index = 1; 
  }
  Eigen::Vector3d angular_delta;
  if (!GetAngularDelta(curr_index, prev_index, angular_delta)) {
    return false;
  }
  Eigen::Matrix3d R_prev = pose_.block<3, 3>(0, 0);
  Eigen::Matrix3d R_curr = Eigen::Matrix3d::Identity();

  UpdateOrientation(angular_delta, R_curr, R_prev);
  Eigen::Quaterniond curr_q(R_curr);
  Eigen::Vector3d vel_delta;
  double delta_t;
  if (!GetVelocityDelta(curr_index, prev_index, R_curr, R_prev, delta_t,
                        vel_delta)) {
    return false;
  }

  UpdatePosition(delta_t, vel_delta);
  imu_data_buff_.pop_front();
}
```
## using EVO for evaluation
![image](https://user-images.githubusercontent.com/11698181/152662207-391f2d49-9a93-4c96-84c9-c2d5c9f09f4a.png)
![image](https://user-images.githubusercontent.com/11698181/152662224-f93ae0d8-c563-4622-97b8-7149cfe7ad04.png)


## 欧拉法

