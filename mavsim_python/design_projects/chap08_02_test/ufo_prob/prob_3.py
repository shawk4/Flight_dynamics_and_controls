import sys
sys.path.append('../../..')
import numpy as np
from control.pi_control import PIControl
from control.pd_control_with_rate import PDControlWithRate
import matplotlib.pyplot as plt
from models.mav_dynamics_control import MavDynamics
import parameters.simulation_parameters as SIM
from scipy import signal
# import control

Boeing = MavDynamics(SIM.ts_simulation)

#----------pitch loop-------------
a_theta1 = .668 
a_theta2 = 1.27
a_theta3 = -2.08

wn_pitch = 3.0
zeta_pitch = 0.6
pitch_kp = (wn_pitch**2-a_theta2)/a_theta3
pitch_kd = (2*zeta_pitch*wn_pitch-a_theta1)/a_theta3
print("Pitch_kp:" + str(pitch_kp))
print("Pitch_kd:" + str(pitch_kd))

K_theta_DC = pitch_kp*a_theta3/(a_theta2 + pitch_kp*a_theta3)

#----------altitude loop-------------
Va0 = 830

wn_altitude = 0.2 # 15 times slower than pitch
zeta_altitude = 0.9
altitude_kp = (2*zeta_altitude*wn_altitude)/(K_theta_DC*Va0)
altitude_ki = wn_altitude**2/(K_theta_DC*Va0)
print("altitude_kp:" + str(altitude_kp))
print("altitude_ki:" + str(altitude_ki))


lti1 = signal.lti([-2.08],[1.0,0.668,1.27])

de_q = signal.lti([a_theta3,0.0],[1.0,a_theta1,a_theta2])

Pc = signal.lti([pitch_kd, altitude_kp],[1.0])
Pp = signal.lti([],[])

hc = signal.lti([altitude_kp*altitude_ki,altitude_kp],[altitude_ki])
hp = signal.lti([],[])

tf = signal.sos2tf([de_q,Pc,Pp,hc,hp])
time, y = signal.step(tf)

plt.plot(time, y)
plt.grid
plt.show()


# instantiate longitudinal controllers
# pitch_from_elevator = PDControlWithRate(
#                 kp=pitch_kp,
#                 kd=pitch_kd,
#                 limit=np.radians(45))
# altitude_from_pitch = PIControl(
#                 kp=altitude_kp,
#                 ki=altitude_ki,
#                 limit=np.radians(30))

# h_c = 100
# altitude = 20000
# theta = 0
# q = 0
# altitude_hist = [20000]
# e_hist = [theta]
# timeStep = [0]
# for i in range(100000):
#     # longitudinal autopilot
#     # -------autopilot-------------
#     # might normaly saturate altitude command here
#     delta_h = altitude_from_pitch.update(h_c, altitude)
#     delta_e = pitch_from_elevator.update(delta_h, theta, q) # theta_c = comanded theta q = pitch_rate
#     # might saturate throttle command here
#     altitude += delta_h
#     theta += delta_e

#     # -------physical system-------------
#     # Boeing.update(delta_e)  # propagate the MAV dynamics

#     altitude_hist.append(altitude)
#     e_hist.append(delta_e)
#     timeStep.append(i)

# plt.plot(timeStep, altitude_hist, 'r--')
# # plt.plot(timeStep, e_hist, 'g.')
# plt.show()




    