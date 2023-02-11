"""
mavsimPy
    Homework check for chapter 3
"""
import sys
sys.path.append('../..')
import numpy as np
from models.mav_dynamics_forces import MavDynamics
import parameters.simulation_parameters as SIM
import parameters.aerosonde_parameters as MAV

state = np.array([[5], [2], [-20], [5],
                        [0], [0], [1], [0], [0],
                        [0], [1], [0.5], [0], [0], [0]])
forces_moments = np.array([[10, 5, 0, 0, 14, 0]]).T
mav = MavDynamics(SIM.ts_simulation)
x_dot = mav._derivatives(state, forces_moments)

print("State Derivatives: Case 1")
print("north_dot: ", x_dot[0])
print("east_dot: ", x_dot[1])
print("down_dot: ", x_dot[2])
print("u_dot: ", x_dot[3])
print("v_dot: " , x_dot[4])
print("w_dot: " , x_dot[5])
print("e0_dot: " , x_dot[6])
print("e1_dot: " , x_dot[7])
print("e2_dot: " , x_dot[8])
print("e3_dot: " , x_dot[9])
print("p_dot: " , x_dot[10])
print("q_dot: " , x_dot[11])
print("r_dot: " , x_dot[12])
print(" ")


state = np.array([[5], [2], [-20], [0],
                        [3], [6], [1], [.6], [0],
                        [.2], [0], [0], [3], [0], [0]])
forces_moments = np.array([[10, 5, 0, 0, 14, 0]]).T
mav = MavDynamics(SIM.ts_simulation)
x_dot = mav._derivatives(state, forces_moments)

print("State Derivatives: Case 2")
print("north_dot: ", x_dot[0])
print("east_dot: ", x_dot[1])
print("down_dot: ", x_dot[2])
print("u_dot: ", x_dot[3])
print("v_dot: " , x_dot[4])
print("w_dot: " , x_dot[5])
print("e0_dot: " , x_dot[6])
print("e1_dot: " , x_dot[7])
print("e2_dot: " , x_dot[8])
print("e3_dot: " , x_dot[9])
print("p_dot: " , x_dot[10])
print("q_dot: " , x_dot[11])
print("r_dot: " , x_dot[12])

# State Derivatives: Case 1
# north_dot:  [5.]
# east_dot:  [0.]
# down_dot:  [0.]
# u_dot:  [0.90909091]
# v_dot:  [0.45454545]
# w_dot:  [2.5]
# e0_dot:  [-0.]
# e1_dot:  [0.5]
# e2_dot:  [0.25]
# e3_dot:  [0.]
# p_dot:  [0.06073576]
# q_dot:  [12.22872247]
# r_dot:  [-0.08413156]

# State Derivatives: Case 2
# north_dot:  [0.08746356]
# east_dot:  [-1.96793003]
# down_dot:  [2.79883382]
# u_dot:  [9.90909091]
# v_dot:  [0.45454545]
# w_dot:  [0.]
# e0_dot:  [-0.3]
# e1_dot:  [0.]
# e2_dot:  [-0.9]
# e3_dot:  [1.5]
# p_dot:  [0.]
# q_dot:  [13.28951542]
# r_dot:  [0.]


