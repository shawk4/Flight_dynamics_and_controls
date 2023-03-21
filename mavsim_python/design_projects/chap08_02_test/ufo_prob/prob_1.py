import sys
sys.path.append('../../..')
import numpy as np
from tools.rotations import Quaternion2Rotation, Quaternion2Euler, Euler2Rotation

# convert Euler to Quaternion
quaternion = Euler2Rotation(np.deg2rad(-25.0), np.deg2rad(-2.0), np.deg2rad(110.0))
print("Quaternion")
print(quaternion)

# part 2 rotate velocitys from body to inertial frame
vel = np.array([35,-2,-5]).T # inertial velocties or Vg expressed in the body frame
R = Quaternion2Rotation(quaternion)
print("North, East, Down")
print(R @ vel)


