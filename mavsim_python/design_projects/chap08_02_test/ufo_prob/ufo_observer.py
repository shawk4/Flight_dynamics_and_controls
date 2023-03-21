# ufo_ekf.py
import numpy as np
from math import fmod

class EkfStateObserver:
    # implement continuous-discrete EFK for UFO
    def __init__(self):
        # define model parameters
        self.m = 10
        self.b0 = 0.6
        self.b1 = .2
        self.g = 9.81
        self.Ts = 0.05
        # process and measurement noise
        self.Q = np.array([[1.0**2, 0.0 ],
                           [0.0, 0.1**2 ]])
        self.R = np.array([[0.3**2, 0.0 ],
                           [0.0, 0.3**2 ]])
        # initialize state and covariance
        self.vel0 = 0
        self.pos0 = 0
        self.xhat = np.array([ [self.vel0], [self.pos0] ])
        self.Px = self.P = np.diag([0, 0])

    
    def f(self, x, measurement):
        # system dynamics for propagation model: xdot = f(x, u)
        ##### TODO #####
        m  =  self.m
        b0 = self.b0
        b1 = self.b1
        g  =  self.g
        Ts = self.Ts

        v = x.item(0)
        pos = x.item(1)
        F = measurement.item(1)

        f_ = np.zeros((2,1))
        f_[0] = -b0/m*v - b1/m*v**3 -g + 1/m*F
        f_[1] =  v
        return f_

    def h(self, x, measurement):
        # measurement model y
        v = x.item(0)
        z = x.item(1)
        h_ = np.array([ [v], # x-accel 
                        [z]])# z-accel
        return h_

    def update(self,inp):
        # ekf algorithm for ufo
        F = inp[0]
        z_m = inp[1]
        t = inp[2]

        # prediction step
        N = 10
        Tp = self.Ts/10
        for j in range(0, N):
            self.xhat = self.xhat + Tp*self.f(self.xhat,inp)
            A = np.array([[-self.b0/self.m,0],
                          [1,0]])
            Ad = np.identity(2) + A*Tp + A@A*Tp**2
            self.Px = Ad @ self.Px @ Ad.T + Tp**2*self.Q 

        # correction step
                # measurement updates
        h = self.h(self.xhat, inp)
        Ci = np.array([[1,0],
                       [0,1]])
        y = np.array([[F, z_m]]).T
        S_inv = np.linalg.inv(self.R + Ci @ self.Px @ Ci.T)
        if (np.fmod(t,20*self.Ts) == 0): # np.fmod(t,20*self.Ts) == 0  # True
            Li = self.Px@Ci.T@S_inv
            self.Px = (np.identity(2) - Li@Ci) @ self.Px @ (np.identity(2) - Li @ Ci).T + Li @ self.R @ Li.T
            self.xhat = self.xhat + Li@(y - h)
        

        # return state estimate
        return self.xhat