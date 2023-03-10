"""
observer
    - Beard & McLain, PUP, 2012
    - Last Update:
        3/2/2019 - RWB
"""
import sys
import numpy as np
from scipy import stats
sys.path.append('..')
import parameters.control_parameters as CTRL
import parameters.simulation_parameters as SIM
import parameters.sensor_parameters as SENSOR
from tools.wrap import wrap
from message_types.msg_state import MsgState
from message_types.msg_sensors import MsgSensors

import parameters.aerosonde_parameters as MAV # !!! should have already been here or what?

class Observer:
    def __init__(self, ts_control, initial_measurements = MsgSensors()):
        # initialized estimated state message
        self.estimated_state = MsgState()
        # use alpha filters to low pass filter gyros and accels
        # alpha = Ts/(Ts + tau) where tau is the LPF time constant

        ##### TODO #####
        self.lpf_gyro_x = AlphaFilter(alpha=0.9, y0=initial_measurements.gyro_x)
        self.lpf_gyro_y = AlphaFilter(alpha=.9, y0=initial_measurements.gyro_y)
        self.lpf_gyro_z = AlphaFilter(alpha=.9, y0=initial_measurements.gyro_z)
        self.lpf_accel_x = AlphaFilter(alpha=.5, y0=initial_measurements.accel_x)
        self.lpf_accel_y = AlphaFilter(alpha=.5, y0=initial_measurements.accel_y)
        self.lpf_accel_z = AlphaFilter(alpha=.5, y0=initial_measurements.accel_z)
        # use alpha filters to low pass filter absolute and differential pressure
        self.lpf_abs = AlphaFilter(alpha=0.9, y0=initial_measurements.abs_pressure)
        self.lpf_diff = AlphaFilter(alpha=0.7, y0=initial_measurements.diff_pressure)
        # ekf for phi and theta
        self.attitude_ekf = EkfAttitude()
        # ekf for pn, pe, Vg, chi, wn, we, psi
        self.position_ekf = EkfPosition()

    def update(self, measurement):
        ##### TODO #####
        # estimates for p, q, r are low pass filter of gyro minus bias estimate
        self.estimated_state.p = self.lpf_gyro_x.update(measurement.gyro_x)
        self.estimated_state.q = self.lpf_gyro_y.update(measurement.gyro_y)
        self.estimated_state.r = self.lpf_gyro_z.update(measurement.gyro_z)

        # invert sensor model to get altitude and airspeed
        self.estimated_state.altitude = self.lpf_abs.update(measurement.abs_pressure)/MAV.rho/MAV.gravity
        self.estimated_state.Va = np.sqrt(2/MAV.rho*self.lpf_diff.update(measurement.diff_pressure))

        # estimate phi and theta, (roll and pitch) with simple ekf
        self.attitude_ekf.update(measurement, self.estimated_state)

        # estimate pn, pe, Vg, chi, wn, we, psi
        self.position_ekf.update(measurement, self.estimated_state)

        # not estimating these
        self.estimated_state.alpha = self.estimated_state.theta#0.0
        self.estimated_state.beta = 0.0
        self.estimated_state.bx = 0.0
        self.estimated_state.by = 0.0
        self.estimated_state.bz = 0.0
        return self.estimated_state


class AlphaFilter:
    # alpha filter implements a simple low pass filter
    # y[k] = alpha * y[k-1] + (1-alpha) * u[k]
    def __init__(self, alpha=0.7, y0=0.0):
        self.alpha = alpha  # filter parameter
        self.y = y0  # initial condition

    def update(self, u):
        ##### TODO #####
        self.y = self.alpha * self.y + (1 - self.alpha) * u
        return self.y


class EkfAttitude:
    # implement continous-discrete EKF to estimate roll and pitch angles
    def __init__(self):        
        ##### TODO #####
        self.Q = np.diag([0.1, 0.1])
        self.Q_gyro = np.diag([SENSOR.gyro_sigma, SENSOR.gyro_sigma, SENSOR.gyro_sigma])
        self.R_accel = np.diag([SENSOR.accel_sigma, SENSOR.accel_sigma, SENSOR.accel_sigma])
        # self.Q_gyro = np.diag([0.1, 0.1, 0.1])
        # self.R_accel = np.diag([0.000001, 0.000001, 0.000001])
        self.N = 10  # number of prediction step per sample
        self.xhat = np.array([[0.0], [0.0]]) # initial state: phi, theta
        self.P = np.diag([0, 0])
        self.Ts = SIM.ts_control/self.N
        self.gate_threshold = 999 #stats.chi2.isf(q=?, df=?) !!! turned off for now

    def update(self, measurement, state):
        self.propagate_model(measurement, state)
        self.measurement_update(measurement, state)
        state.phi = self.xhat.item(0)
        state.theta = self.xhat.item(1)


    def f(self, x, measurement, state):
        # system dynamics for propagation model: xdot = f(x, u)
        ##### TODO #####
        phi = x.item(0)
        theta = x.item(1)
        f_ = np.zeros((2,1))
        p = measurement.gyro_x
        q = measurement.gyro_y
        r = measurement.gyro_z #!!! What in the world am I doing with state ??? it is my currnent unpropagated state
        f_[0] =  p + q*np.sin(phi)*np.tan(theta) + r*np.cos(phi)*np.tan(theta)
        f_[1] =  q*np.cos(phi)-r*np.sin(phi)
        return f_

    def h(self, x, measurement, state):
        # measurement model y
        ##### TODO #####
        phi = x.item(0)
        theta = x.item(1)
        g = MAV.gravity
        Va = np.sqrt(2/MAV.rho*measurement.diff_pressure)
        p = measurement.gyro_x
        q = measurement.gyro_y
        r = measurement.gyro_z
        #predicting acceleration measurements
        h_ = np.array([ [q*Va*np.sin(theta)+g*   np.sin(theta) ], # x-accel 
                        [r*Va*np.cos(theta)-p*Va*np.sin(theta)-g*np.cos(theta)*np.sin(phi) ], # y-accel
                        [-q*Va*np.cos(theta)-  g*np.cos(theta)*np.cos(phi) ]])# z-accel
        return h_

    def propagate_model(self, measurement, state):
        # model propagation
        ##### TODO #####
        Tp = self.Ts #/self.N
        for i in range(0, self.N):
            # Tp = Tout/self.N
            self.xhat[:2] = self.xhat[:2] + Tp*self.f(self.xhat,measurement,state)
            A = jacobian(self.f,self.xhat,measurement,state)
            Ad = np.identity(2) + A*Tp + A@A*Tp**2
            self.P = Ad@self.P@Ad.T + Tp**2*self.Q 
            # self.P = np.zeros((2,2))

    def measurement_update(self, measurement, state):
        # measurement updates
        h = self.h(self.xhat, measurement, state)
        Ci = jacobian(self.h, self.xhat, measurement, state)
        y = np.array([[measurement.accel_x, measurement.accel_y, measurement.accel_z]]).T
        ##### TODO #####
        S_inv = np.zeros((3,3))
        if (y-h).T @ S_inv @ (y-h) < self.gate_threshold:
            Li = self.P@Ci.T@np.linalg.inv(self.R_accel+ Ci@self.P@Ci.T)
            self.P = (np.identity(2) - Li@Ci) @ self.P @ (np.identity(2) - Li@Ci).T + Li@self.R_accel@Li.T
            self.xhat = self.xhat + Li@(y - h)
            # self.P = np.zeros((2,2))
            # self.xhat = np.zeros((2,1))


class EkfPosition:
    # implement continous-discrete EKF to estimate pn, pe, Vg, chi, wn, we, psi
    def __init__(self):
        self.Q = np.diag([
                    0.1,  # pn
                    0.1,  # pe
                    0.1,  # Vg
                    0.1, # chi
                    0.1, # wn
                    0.1, # we
                    0.0001, #0.0001, # psi
                    ])
        self.R_gps = np.diag([
                    SENSOR.gps_n_sigma,  # y_gps_n
                    SENSOR.gps_e_sigma,  # y_gps_e
                    SENSOR.gps_Vg_sigma,  # y_gps_Vg
                    SENSOR.gps_course_sigma,  # y_gps_course
                    ])
        self.R_pseudo = np.diag([
                    0.01, # pseudo measurement #1
                    0.01 # pseudo measurement #2
                    ])
        self.N = 10  # number of prediction step per sample
        self.Ts = (SIM.ts_control / self.N)
        self.xhat = np.array([[0.0], [0.0], [25.0], [0.0], [0.0], [0.0], [0.0]])
        self.P = np.diag([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.gps_n_old = 0
        self.gps_e_old = 0
        self.gps_Vg_old = 0
        self.gps_course_old = 0
        self.pseudo_threshold = 0 #stats.chi2.isf(q=?, df=?)
        self.gps_threshold = 100000 # don't gate GPS

    def update(self, measurement, state):
        self.propagate_model(measurement, state)
        self.measurement_update(measurement, state)
        state.north = self.xhat.item(0)
        state.east = self.xhat.item(1)
        state.Vg = self.xhat.item(2)
        state.chi = self.xhat.item(3)
        state.wn = self.xhat.item(4)
        state.we = self.xhat.item(5)
        state.psi = self.xhat.item(6)

    def f(self, x, measurement, state):
        # system dynamics for propagation model: xdot = f(x, u)
        # xhat estimated states
        Vg  = x.item(2)
        chi = x.item(3)
        wn  = x.item(4)
        we  = x.item(5)
        psi = x.item(6)
        # u commanded inputs (state?)
        Va = np.sqrt(2/MAV.rho*measurement.diff_pressure)
        phi = state.phi
        theta = state.theta
        q = measurement.gyro_y
        r = measurement.gyro_z
        # other
        g = MAV.gravity
        phi_dot = q*np.sin(phi)/np.cos(theta)+r*np.cos(phi)/np.cos(theta)

        f_ = np.array([[Vg*np.cos(chi)],
                       [Vg*np.sin(chi)],
                       [Va/Vg*phi_dot*(-wn*np.sin(phi)+we*np.cos(phi))],
                       [g/Vg*np.tan(phi)*np.cos(chi-psi)],
                       [0.0],
                       [0.0],
                       [phi_dot]])
        return f_

    def h_gps(self, x, measurement, state):
        # measurement model for gps measurements
        pn = x.item(0)
        pe = x.item(1)
        Vg = x.item(2)
        chi = x.item(3)

        h_ = np.array([
            [pn], #pn
            [pe], #pe
            [Vg], #Vg
            [chi], #chi
        ])
        return h_

    def h_pseudo(self, x, measurement, state):
        # measurement model for wind triangale pseudo measurement
        # x hat estimated states
        Vg  = x.item(2)
        chi = x.item(3)
        wn  = x.item(4)
        we  = x.item(5)
        
        # u commanded inputs
        Va = np.sqrt(2/MAV.rho*measurement.diff_pressure)
        phi = state.phi

        h_ = np.array([
            [Va*np.cos(phi)+wn-Vg*np.cos(chi)],  # wind triangle x
            [Va*np.sin(phi)+we-Vg*np.sin(chi)],  # wind triangle y
        ])
        return h_

    def propagate_model(self, measurement, state):
        # model propagation
        Tp = self.Ts #/self.N !!! what is this???
        for i in range(0, self.N):
            # propagate model
            # self.xhat = np.zeros((7,1))
            self.xhat = self.xhat + Tp*self.f(self.xhat,measurement,state)
            # compute Jacobian
            A = jacobian(self.f,self.xhat,measurement,state)
            # convert to discrete time models
            Ad = np.identity(7) + A*Tp + A@A*Tp**2
            # update P with discrete time model
            self.P = Ad@self.P@Ad.T + Tp**2*self.Q 
            # self.P = np.zeros((7,7))

    def measurement_update(self, measurement, state):
        # always update based on wind triangle pseudo measurement
        h = self.h_pseudo(self.xhat, measurement, state)
        Ci = jacobian(self.h_pseudo, self.xhat, measurement, state)
        y = np.array([[0, 0]]).T
        S_inv = np.zeros((2,2))
        if (y-h).T @ S_inv @ (y-h) < self.pseudo_threshold:
            Li = self.P@Ci.T@np.linalg.inv(self.R_pseudo+Ci@self.P@Ci.T)
            self.P = (np.identity(2) - Li@Ci) @ self.P @ (np.identity(2)-Li@Ci).T + Li@self.R_pseudo@Li.T
            self.xhat = self.xhat + Li @ (y - h) 
            # self.P = np.zeros((7,7))
            # self.xhat = np.zeros((7,1))

        # only update GPS when one of the signals changes
        if (measurement.gps_n != self.gps_n_old) \
            or (measurement.gps_e != self.gps_e_old) \
            or (measurement.gps_Vg != self.gps_Vg_old) \
            or (measurement.gps_course != self.gps_course_old):

            h = self.h_gps(self.xhat, measurement, state)
            C = jacobian(self.h_gps, self.xhat, measurement, state)
            y_chi = wrap(measurement.gps_course, h[3, 0])
            y = np.array([[measurement.gps_n,
                           measurement.gps_e,
                           measurement.gps_Vg,
                           y_chi]]).T
            S_inv = np.zeros((4,4))
            if (y-h).T @ S_inv @ (y-h) < self.gps_threshold:
                Li = self.P@C.T@np.linalg.inv(self.R_gps + C@self.P@C.T)
                self.P = (np.identity(7) - Li@C) @ self.P
                self.xhat = self.xhat + Li @ (y - h)
                # self.P = np.zeros((7,7))
                # self.xhat = np.zeros((7,1))

            # update stored GPS signals
            self.gps_n_old = measurement.gps_n
            self.gps_e_old = measurement.gps_e
            self.gps_Vg_old = measurement.gps_Vg
            self.gps_course_old = measurement.gps_course


def jacobian(fun, x, measurement, state):
    # compute jacobian of fun with respect to x
    f = fun(x, measurement, state)
    m = f.shape[0]
    n = x.shape[0]
    eps = 0.0001  # deviation
    J = np.zeros((m, n))
    for i in range(0, n):
        x_eps = np.copy(x)
        x_eps[i][0] += eps
        f_eps = fun(x_eps, measurement, state)
        df = (f_eps - f) / eps
        J[:, i] = df[:, 0]
    return J