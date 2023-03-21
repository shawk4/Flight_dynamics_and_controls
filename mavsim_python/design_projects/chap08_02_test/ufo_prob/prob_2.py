import sys
sys.path.append('../../..')
import numpy as np
from tools.rotations import Quaternion2Rotation, Quaternion2Euler, Euler2Rotation

# Givens:
vel = np.array([25,1,0]).T # inertial velocties or Vg expressed in the body frame [u,v,w]
att = np.array([-10,5,105]).T # inertial attitude expressed in the body frame [phi, theta, psi]
wn = 2
we = -5
wd = -1

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
print("alpha:" + str(alpha))

# compute sideslip angle beta
beta = np.arcsin(vr/(Va))
print("Beta:" + str(beta))

# compute air mass refrenced flight path angle gamma_a
gamma_a = att[1] - alpha # theta - alpha
print("gamma_a:" + str(gamma_a))

# compute flight path angle gamma !!! come back
gamma = alpha
print("gamma:" + str(gamma))




