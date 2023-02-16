import sys
sys.path.append('..')
import numpy as np
import design_projects.chap05.model_coef as TF
import parameters.aerosonde_parameters as MAV


#### TODO #####
gravity = MAV.gravity  # gravity constant
Va0 = TF.Va_trim
rho = MAV.rho # density of air
sigma = 0  # low pass filter gain for derivative

#----------roll loop-------------
# get transfer function data for delta_a to phi
wn_roll = 14.5
zeta_roll = 0.707
wn_roll = 14
zeta_roll = 0.707
roll_kp = 0
roll_kd = 0

#----------course loop-------------
wn_course = 1.45
zeta_course = 0.707
course_kp = 0
course_ki = 0

#----------yaw damper-------------
yaw_damper_p_wo = 0
yaw_damper_kr = 0

#----------pitch loop-------------
wn_pitch = 5.8
zeta_pitch = 0.707
pitch_kp = 0
pitch_kd = 0

k_p_theta = wn_pitch**2-TF.a_theta2/TF.a_theta3
K_theta_DC = k_p_theta*TF.a_theta3/(TF.a_theta2 + k_p_theta*TF.a_theta3)


#----------altitude loop-------------
wn_altitude = 0.058
zeta_altitude = 0.707
altitude_kp = 0
altitude_ki = 0
altitude_zone = 0

#---------airspeed hold using throttle---------------
wn_airspeed_throttle = 0.58
zeta_airspeed_throttle = 0.707
airspeed_throttle_kp = 1
airspeed_throttle_ki = .01
