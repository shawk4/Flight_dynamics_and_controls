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
roll_kp = wn_roll**2/TF.a_phi2
roll_kd = (2*zeta_roll*wn_roll-TF.a_phi1)/TF.a_phi2

#----------course loop-------------
wn_course = 1.45
zeta_course = 0.707
course_kp = 2*zeta_course*wn_course*Va0/gravity
course_ki = wn_course**2*Va0/gravity

#----------yaw damper-------------
yaw_damper_p_wo = 0.45
yaw_damper_kr = 0.2

#----------pitch loop-------------
wn_pitch = 30
zeta_pitch = 0.707
pitch_kp = (wn_pitch**2-TF.a_theta2)/TF.a_theta3
pitch_kd = (2*zeta_pitch*wn_pitch-TF.a_theta1)/TF.a_theta3

K_theta_DC = pitch_kp*TF.a_theta3/(TF.a_theta2 + pitch_kp*TF.a_theta3)


#----------altitude loop-------------
wn_altitude = 0.58
zeta_altitude = 0.707
altitude_kp = (2*zeta_altitude*wn_altitude)/(K_theta_DC*Va0)
altitude_ki = wn_altitude**2/(K_theta_DC*Va0)
altitude_zone = 10

#---------airspeed hold using throttle---------------
wn_airspeed_throttle = 5 #0.58
zeta_airspeed_throttle = 0.707
airspeed_throttle_kp = (2*zeta_airspeed_throttle*wn_airspeed_throttle-TF.a_V1)/TF.a_V2
airspeed_throttle_ki = wn_airspeed_throttle**2/TF.a_V2
pass
