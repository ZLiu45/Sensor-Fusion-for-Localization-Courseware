# Sensor Fusion: Lidar Odometry -- 多传感器融合定位与建图: 惯性导航解算

## 中值法
```
while (imu_data_buff_.size() > 1) {
  size_t prev_index = 0; 
  size_t curr_index = 1; 
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
In `GetAngularDelta()`: 
```
if (integration_method_ == "euler") {
    angular_delta = delta_t * angular_vel_prev;  
}
```
In `GetVelocityDelta()`: 
```
if (integration_method_ == "euler") {
      velocity_delta = delta_t * linear_acc_prev;
}
```
![image](https://user-images.githubusercontent.com/11698181/152662590-9976d41f-792f-4fe2-a4d0-b630768d2c6d.png)
![image](https://user-images.githubusercontent.com/11698181/152662594-982aae8b-38b4-40a5-818f-4b70bbad0825.png)

## GNSS_INS_SIM
### 使用高精度的IMU
#### motion file 
```
1,0,0,0,0,0,0,20,1
5,0,2,0,3,0,0,25,1
3,-180,0,0,0,0,0,25,1
1,0,0,0,0,0,0,25,1
3,180,0,0,0,0,0,25,1
1,0,0,0,0,0,0,25,1
3,-180,0,0,0,0,0,25,1
1,0,0,0,0,0,0,25,1
5,0,0,0,0,0,0,25,1
1,0,0,0,0,0,0,10,1
3,90,-2,0,0,0,0,25,1
1,0,0,0,0,0,0,25,1
3,180,0,0,0,0,0,25,1
1,0,0,0,0,0,0,25,1
3,-180,0,0,0,0,0,25,1
1,0,0,0,0,0,0,25,1
3,180,0,0,0,0,0,25,1
1,0,0,0,0,0,0,25,1
```
#### trajectory plots 
![image](https://user-images.githubusercontent.com/11698181/152671127-407aedb7-754e-458d-ad37-9bf70ccedb4d.png)
![image](https://user-images.githubusercontent.com/11698181/152671143-5059dbfc-40e0-4b59-a129-2df5bdc5f1dc.png)

#### 中值法
![image](https://user-images.githubusercontent.com/11698181/152671337-d617b1c6-099d-418d-99cc-6549d69ac5d6.png)
![image](https://user-images.githubusercontent.com/11698181/152671353-c655b4c9-61ad-49f9-9aa6-42840a2d65cd.png)

#### 欧拉法 
![image](https://user-images.githubusercontent.com/11698181/152671634-04275dfd-336d-40d6-bbcd-506e0f690255.png)
![image](https://user-images.githubusercontent.com/11698181/152671645-535457e7-d67f-4ce4-b624-97b2dffb9c3a.png)

### 使用低精度的IMU
#### 中值法


#### 欧拉法 
