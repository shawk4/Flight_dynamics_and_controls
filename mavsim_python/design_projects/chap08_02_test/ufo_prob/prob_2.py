import sys
sys.path.append('../../..')
import numpy as np
from tools.rotations import Quaternion2Rotation, Quaternion2Euler, Euler2Rotation

# Givens:
vel = np.array([25,1,-3]).T # inertial velocties or Vg expressed in the body frame [u,v,w]
att = np.array([-10,5,105]).T *np.pi/180# inertial attitude expressed in the body frame [phi, theta, psi]
wind_steady_state = np.array([2,-5,-1]).T # wind in inertial or world frame

# convert wind vector from world to body frame 
R_bv = Euler2Rotation(att.item(0),att.item(1),att.item(2)).T
wind = R_bv @ wind_steady_state
wn = wind.item(0)
we = wind.item(1)
wd = wind.item(2)

# velocity vector relative to the airmass ([ur , vr, wr]= ?)
# V_ba = np.array([[u-uw, v-vw, w-ww]])
ur = vel[0]-wn
vr = vel[1]-we
wr = vel[2]-wd

# compute airspeed Va
Va = np.sqrt(ur**2 + vr**2 + wr**2)
print("Va:" + str(Va))

# compute angle of attack alpha
alpha = np.arctan(wr/ur)
print("alpha:" + str(alpha*180/np.pi))

# compute sideslip angle beta
beta = np.arcsin(vr/(Va))
print("Beta:" + str(beta*180/np.pi))

# compute air mass refrenced flight path angle gamma_a
gamma_a = att[1] - alpha # theta - alpha
print("gamma_a:" + str(gamma_a*180/np.pi))

# compute flight path angle gamma 
Vg = np.sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2)
p_d_dot = (R_bv.T @ vel)[2]
gamma = np.arcsin(-p_d_dot/Vg)
print("gamma:" + str(gamma*180/np.pi))




