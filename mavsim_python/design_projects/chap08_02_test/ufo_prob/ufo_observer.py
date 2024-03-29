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
        self.R = 0.3**2
        # self.R = np.array([[0.3**2, 0.0 ],
        #                    [0.0, 0.3**2 ]])

        # initialize state and covariance
        self.vel0 = 0
        self.pos0 = 20
        self.xhat = np.array([ [self.vel0], [self.pos0] ])
        self.Px = self.P = np.diag([1.0, 1.0])

    
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
            # self.xhat = self.xhat + Tp*self.f(self.xhat,inp) # !!! here vhat is used insted of x hat I bet this causes my zig zags
            vhat = self.xhat.item(0)
            f1 = -(self.b0/self.m)*vhat - (self.b1/self.m)*vhat**3 - self.g + F/self.m
            f2 = vhat
            f_ufo = np.array([[f1],[f2]])
            self.xhat = self.xhat + Tp*f_ufo
            A = np.array([[-self.b0/self.m - 3*self.b1/self.m*vhat**2,0.0],
                          [1,0.0]])
            # Ad = np.identity(2) + A*Tp + A@A*Tp**2                       # why don't we need to make this discrete? 
            self.Px = self.Px + Tp* (A @ self.Px + self.Px @ A.T + self.Q) # why is this handeled differently?

        # correction step
        # measurement updates
        # h = self.h(self.xhat, inp)
        Ci = np.array([[0,1]])
        # y = np.array([[F, z_m]]).T
        # S_inv = np.linalg.inv(self.R + Ci @ self.Px @ Ci.T)
        if True: # np.fmod(t,20*self.Ts) == 0  # True
            zhat = self.xhat.item(1)
            Li = self.Px@Ci.T/(self.R + Ci@self.Px@Ci.T)

            # this is also handled differently and I don't really know why the second part is added or left off?
            # @ (np.identity(2) - Li @ Ci).T + Li @ self.R @ Li.T 
            self.Px = (np.identity(2) - Li@Ci) @ self.Px 
            self.xhat = self.xhat + Li*(z_m-zhat)#Li@(y - h)
        

        # return state estimate
        return self.xhat