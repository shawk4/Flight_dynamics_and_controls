"""
mavDynamics 
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
mavsim_python
    - Beard & McLain, PUP, 2012
    - Update history:  
        2/24/2020 - RWB
"""
import sys
sys.path.append('..')
import numpy as np


# load message types
from message_types.msg_state import MsgState
from message_types.msg_delta import MsgDelta
import parameters.aerosonde_parameters as MAV
from tools.rotations import Quaternion2Rotation, Quaternion2Euler, Euler2Rotation

class MavDynamics:
    def __init__(self, Ts):
        self._ts_simulation = Ts
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        # We will also need a variety of other elements that are functions of the _state and the wind.
        # self.true_state is a 19x1 vector that is estimated and used by the autopilot to control the aircraft:
        # true_state = [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        self._state = np.array([[MAV.north0],  # (0)
                               [MAV.east0],   # (1)
                               [MAV.down0],   # (2)
                               [MAV.u0],    # (3)
                               [MAV.v0],    # (4)
                               [MAV.w0],    # (5)
                               [MAV.e0],    # (6)
                               [MAV.e1],    # (7)
                               [MAV.e2],    # (8)
                               [MAV.e3],    # (9)
                               [MAV.p0],    # (10)
                               [MAV.q0],    # (11)
                               [MAV.r0],    # (12)
                               [0],   # (13)
                               [0],   # (14)
                               ])
        # initialize true_state message
        self.true_state = MsgState()


    ###################################
    # public functions
    def update(self, forces_moments):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''

        # Integrate ODE using Runge-Kutta RK4 algorithm
        time_step = self._ts_simulation
        k1 = self._derivatives(self._state[0:13], forces_moments)
        k2 = self._derivatives(self._state[0:13] + time_step/2.*k1, forces_moments)
        k3 = self._derivatives(self._state[0:13] + time_step/2.*k2, forces_moments)
        k4 = self._derivatives(self._state[0:13] + time_step*k3, forces_moments)
        self._state[0:13] += time_step/6 * (k1 + 2*k2 + 2*k3 + k4)

        # normalize the quaternion
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        normE = np.sqrt(e0**2+e1**2+e2**2+e3**2)
        self._state[6][0] = self._state.item(6)/normE
        self._state[7][0] = self._state.item(7)/normE
        self._state[8][0] = self._state.item(8)/normE
        self._state[9][0] = self._state.item(9)/normE

        # update the message class for the true state
        self._update_true_state()


    def external_set_state(self, new_state):
        self._state = new_state

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):
        """
        for the dynamics xdot = f(x, u), returns f(x, u)
        """

        ##### TODO #####
        # Airframe parameters from Appendix E
        # physical parameters
        m = 11 #kg
        Jx = 0.824 # kg-m^2
        Jy = 1.135 # kg-m^2
        Jz = 1.759 # kg-m^2
        Jxz = 0.12 # kg-m^2
        S = .55 # m^2
        b = .29 # m
        c = 0.19 # m
        ro = 1.268 # kg/m^3
        e = 0.9
        # Motor parameters
        Vmax = 44.4 # V
        Dprop = 0.5 # m
        # ...

        # Changeable parameters for each state include:
        moments = 0.0
        products_of_inertia = 0
        initial_conditions = 0

        # Extract the States
        north = state.item(0)[0]
        east = state.item(1)[1]
        down = state.item(2)[2]
        u = state.item(3)
        v = state.item(4)
        w = state.item(5)
        pos = np.array([north,east,down]).T # inertial positions
        Q = np.array([state.item(6),state.item(7),state.item(8),state.item(9)]).T
        p = state.item(10)
        q = state.item(11)
        r = state.item(12)
        # rotation_rate = np.array([p,q,r]).T


        # Extract Forces/Moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)
        f_xyz = np.array([fx,fy,fz]).T #forces along primary axis in the body frame
        l = forces_moments.item(3) # external moment applied to the airfram about body x axis
        My = forces_moments.item(4) # external moment applied to the airfram about body y axis
        n = forces_moments.item(5) # external moment applied to the airfram about body z axis

        theta = 0
        psi = 0
        phi = 0
        #rotation from the vehicle frame to the body frame. Is the vehicle frame the same as the inertial frame???
        R_v_b = np.array([
            [np.cos(theta)*np.cos(psi), np.sin(phi)*np.sin(theta)*np.cos(psi)-np.cos(theta)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.cos(psi)+np.sin(theta)*np.sin(psi)],
            [np.cos(theta)*np.cos(psi), np.sin(phi)*np.sin(theta)*np.cos(psi)+np.cos(theta)*np.sin(psi), np.cos(phi)*np.sin(theta)*np.cos(psi)-np.sin(theta)*np.sin(psi)],
            [-np.sin(theta)           , np.sin(phi)*np.cos(theta)                                      , np.cos(phi)*np.cos(theta)]])
        # R_yaw_pitch_roll = np.array([
        #     [1, np.sin(phi)*np.tan(theta)   , np.cos(phi)*np.tan(theta)],
        #     [0, np.cos(theta)               , -np.sin(phi)],
        #     [0, np.sin(phi)*1/np.cos(theta), np.cos(phi)*1/np.cos(theta)]])
        Q_yaw_pitch_roll = np.array([
            [0,-p,-q,-r],
            [p, 0, r,-q],
            [q,-r, 0, p],
            [r, q,-p, 0]])
        G = Jx*Jz-Jxz**2
        G1 = Jxz*(Jx-Jy+Jz)/G
        G2 = (Jz*(Jz-Jy)+Jxz**2)/G
        G3 = Jz/G
        G4 = Jxz/G
        G5 = Jz-Jx/Jy
        G6 = Jxz/Jy
        G7 = ((Jx-Jy)*Jx+Jxz**2)/G
        G8 = Jx/G

        U = np.array([u,v,w]).T # inertial velocties or Vg expressed in the body frame
        # Position Kinematics
        pos_dot = R_v_b * U

        # Position Dynamics
        U_dot = np.cross(pos,U) + 1/m*f_xyz


        # rotational kinematics
        # euler_dot = R_yaw_pitch_roll @ rotation_rate
        # E_dot = euler_dot # convert euler angles to quaternion
        Q_dot = .5* Q_yaw_pitch_roll @ Q


        # rotatonal dynamics
        p_dot = G1*p*q-G2*q*r         + G3*l+G4*n
        q_dot = G5*p*r-G6*(p**2-r**2) + 1/Jy*My
        r_dot =  G7*p*q-G1*q*r        + G4*l+G8*n
        

        # collect the derivative of the states
        x_dot = np.array([[
            pos_dot[0],pos_dot[1],pos_dot[2],
            U_dot[0],U_dot[1],U_dot[2],
            Q_dot[0],Q_dot[1],Q_dot[2],Q_dot[3],
            p_dot,q_dot,r_dot ]]).T
        # x_dot = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0]]).T
        return x_dot

    def _update_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = 0
        self.true_state.alpha = 0
        self.true_state.beta = 0
        self.true_state.phi = phi
        self.true_state.theta = theta
        self.true_state.psi = psi
        self.true_state.Vg = 0
        self.true_state.gamma = 0
        self.true_state.chi = 0
        self.true_state.p = self._state.item(10)
        self.true_state.q = self._state.item(11)
        self.true_state.r = self._state.item(12)
        self.true_state.wn = 0
        self.true_state.we = 0
        self.true_state.bx = 0
        self.true_state.by = 0
        self.true_state.bz = 0
        self.true_state.camera_az = 0
        self.true_state.camera_el = 0
