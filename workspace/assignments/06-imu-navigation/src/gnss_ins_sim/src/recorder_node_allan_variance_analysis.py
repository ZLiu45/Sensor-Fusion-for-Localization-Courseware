#!/usr/bin/python

import os

import rospkg
import rospy
import rosbag

import math
import numpy as np

from pprint import pprint

from gnss_ins_sim.sim import imu_model
from gnss_ins_sim.sim import ins_sim

from std_msgs.msg import String
from sensor_msgs.msg import Imu
from sensor_msgs.msg import NavSatFix
from nav_msgs.msg import Odometry


def get_gnss_ins_sim(motion_def_file, fs_imu, fs_gps):

    # set IMU model:
    D2R = math.pi/180.0
    imu_err = 'low-accuracy'
    # imu_err = {
    #     # 1. gyro:
    #     # a. random noise:
    #     # gyro angle random walk, deg/rt-hr
    #     'gyro_arw': np.array([0.75, 0.75, 0.75]),
    #     # gyro bias instability, deg/hr
    #     'gyro_b_stability': np.array([10.0, 10.0, 10.0]),
    #     # gyro bias isntability correlation time, sec
    #     'gyro_b_corr': np.array([100.0, 100.0, 100.0]),
    #     # b. deterministic error:
    #     'gyro_b': np.array([0.0, 0.0, 0.0]),
    #     'gyro_k': np.array([1.0, 1.0, 1.0]),
    #     'gyro_s': np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    #     # 2. accel:
    #     # a. random noise:
    #     # accel velocity random walk, m/s/rt-hr
    #     'accel_vrw': np.array([0.05, 0.05, 0.05]),
    #     # accel bias instability, m/s2
    #     'accel_b_stability': np.array([2.0e-4, 2.0e-4, 2.0e-4]),
    #     # accel bias isntability correlation time, sec
    #     'accel_b_corr': np.array([100.0, 100.0, 100.0]),
    #     # b. deterministic error:
    #     'accel_b': np.array([0.0e-3, 0.0e-3, 0.0e-3]),
    #     'accel_k': np.array([1.0, 1.0, 1.0]),
    #     'accel_s': np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    #     # 3. mag:
    #     'mag_si': np.eye(3) + np.random.randn(3, 3)*0.0,
    #     'mag_hi': np.array([10.0, 10.0, 10.0])*0.0,
    #     'mag_std': np.array([0.1, 0.1, 0.1])
    # }
    # # generate GPS and magnetometer data:
    imu = imu_model.IMU(accuracy=imu_err, axis=9, gps=True)
    print("create imu model")
    # # init simulation:
    sim = ins_sim.Sim(
        [fs_imu, fs_gps, fs_imu],
        motion_def_file,
        ref_frame=1,
        imu=imu,
        mode=None,
        env=None,
        algorithm=None
    )

    # run:
    sim.run(1)
    imu_step_size = 1.0 / fs_imu
    gps_step_size = 1.0/fs_gps

    return {
        "ref_pos": sim.dmgr.get_data_all("ref_pos").data,
        "gyro": sim.dmgr.get_data_all("gyro").data[0],
        "accel": sim.dmgr.get_data_all("accel").data[0],
        "imu_stamp": imu_step_size * np.asarray(range(sim.dmgr.get_data_all("gyro").data[0].shape[0])),
        "gps_stamp": gps_step_size * np.asarray(range(sim.dmgr.get_data_all("ref_pos").data.shape[0])),
    }

def gnss_ins_sim_recorder():
    """
    Record simulated GNSS/IMU data as ROS bag
    """
    # ensure gnss_ins_sim_node is unique:
    rospy.init_node('gnss_ins_sim_recorder_node')

    # parse params:
    motion_def_name = rospy.get_param(
        '/gnss_ins_sim_recorder_node/motion_file')
    sample_freq_imu = rospy.get_param(
        '/gnss_ins_sim_recorder_node/sample_frequency/imu')
    sample_freq_gps = rospy.get_param(
        '/gnss_ins_sim_recorder_node/sample_frequency/gps')
    topic_name_imu = rospy.get_param('/gnss_ins_sim_recorder_node/topic_name')
    rosbag_output_path = rospy.get_param(
        '/gnss_ins_sim_recorder_node/output_path')
    rosbag_output_name = rospy.get_param(
        '/gnss_ins_sim_recorder_node/output_name')

    # # generate simulated data:
    motion_def_path = os.path.join(
        "/workspace/assignments/06-imu-navigation/src/gnss_ins_sim/config/motion_def", "allan_variance_analysis.csv")
    imu_simulator_data = get_gnss_ins_sim(motion_def_path, sample_freq_imu, sample_freq_gps)
    
    timestamp_start = rospy.Time.now()

    with rosbag.Bag(
        os.path.join(rosbag_output_path, rosbag_output_name), 'w'
    ) as bag:
        # get timestamp base:

        for (gyro, accel, stamp) in zip(imu_simulator_data["gyro"], imu_simulator_data["accel"], imu_simulator_data["imu_stamp"]):
            msg = Imu()
            # a. set header:
            msg.header.frame_id = 'ORI'
            msg.header.stamp = timestamp_start + \
                rospy.Duration.from_sec(stamp)
            # b. set orientation estimation:
            msg.orientation.x = 0.0
            msg.orientation.y = 0.0
            msg.orientation.z = 0.0
            msg.orientation.w = 1.0
            # c. gyro:
            msg.angular_velocity.x = gyro[0]
            msg.angular_velocity.y = gyro[1]
            msg.angular_velocity.z = gyro[2]
            msg.linear_acceleration.x = accel[0]
            msg.linear_acceleration.y = accel[1]
            msg.linear_acceleration.z = accel[2]
            # write:
            bag.write(topic_name_imu, msg, msg.header.stamp)

    # write rel_pos to txt file 
    gt_lla_path = os.path.join("/workspace/data/gnss_ins_sim/allan_variance_analysis", "gt_sim.txt")
    stamps = imu_simulator_data["gps_stamp"].reshape(imu_simulator_data["gps_stamp"].shape[0], 1)
    stamps = stamps + timestamp_start.to_sec()
    gt_data = np.concatenate((stamps, imu_simulator_data["ref_pos"]), axis=1)
    np.savetxt(gt_lla_path, gt_data, fmt ='%10.10f', delimiter=" ")
    

if __name__ == '__main__':
    gnss_ins_sim_recorder()
