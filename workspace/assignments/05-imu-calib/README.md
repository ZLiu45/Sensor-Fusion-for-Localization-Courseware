# Sensor Fusion: Lidar Odometry -- 多传感器融合定位与建图: 惯性导航原理

## jacobian derivation 
[Assignment_5.pdf](https://github.com/ZLiu45/Sensor-Fusion-for-Localization-Courseware/files/7917712/Assignment_5.pdf)

## calibration with autodiff 
![task2](https://user-images.githubusercontent.com/11698181/150626499-53dd45c4-493e-40bc-a988-5e4784588135.png)

### code: 
```
for (int th_mult = 2; th_mult <= 10; th_mult++) {
    std::vector<imu_tk::DataInterval> static_intervals;
    std::vector<imu_tk::TriadData_<_T>> static_samples;
    std::vector<double> acc_calib_params(9);

    //
    // TODO: implement lower triad model here
    //
    acc_calib_params[0] = init_acc_calib_.misXZ();
    acc_calib_params[1] = init_acc_calib_.misXY();
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
