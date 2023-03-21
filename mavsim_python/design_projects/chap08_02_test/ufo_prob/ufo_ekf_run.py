# ufo_ekf_run.py

import numpy as np
import matplotlib.pyplot as plt
from ufo_observer import EkfStateObserver

# load the data
ufo_data = np.loadtxt(open("ufo_data.txt", "rb"), delimiter=",")
time = ufo_data[:,0]
force = ufo_data[:,1]
z = ufo_data[:,2]
M = time.shape[0]

xhat = np.zeros((2,M))

# create ufo_ekf object of EkfStateObserver class 
ufo_ekf = EkfStateObserver()

for i in range(0, M):
    tmp = ufo_ekf.update(np.array([ force[i], z[i], time[i] ]))
    xhat[np.ix_([0, 1], [i])] = tmp

vel = xhat[0,:]
pos = xhat[1,:]

pos_error = z - pos
mean_pos_error = np.mean(pos_error)
print('Mean position error = ', mean_pos_error, 'm')

max_velocity = np.max(vel)
print('Maximum velocity = ', max_velocity, 'm/s')

plt.figure(1)
plt.subplot(311)
plt.plot(time, vel)
plt.ylabel('velocity (m/s)')
plt.grid()
plt.subplot(312)
plt.plot(time, z, time, pos)
plt.xlabel('time (s)')
plt.ylabel('position (m)')
plt.grid()
plt.subplot(313)
plt.plot(time, z-pos)
plt.xlabel('time (s)')
plt.ylabel('position error (m)')
plt.grid()
plt.show()

