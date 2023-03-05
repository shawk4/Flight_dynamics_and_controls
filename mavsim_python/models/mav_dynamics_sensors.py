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
from message_types.msg_sensors import MsgSensors
from message_types.msg_delta import MsgDelta

import parameters.aerosonde_parameters as MAV
import parameters.sensor_parameters as SENSOR
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
        # store wind data for fast recall since it is used at various points in simulation
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.u0
        self._alpha = 0 #!!!hmmm... 
        self._beta = 0  #!!!hmmm... 
        # initialize true_state message
        self.true_state = MsgState()
        # initialize the sensors message
        self._sensors = MsgSensors()
        # random walk parameters for GPS
        self._gps_eta_n = 0.
        self._gps_eta_e = 0.
        self._gps_eta_h = 0.
        # timer so that gps only updates every ts_gps seconds
        self._t_gps = 999.  # large value ensures gps updates at initial time.
        # update velocity data and forces and moments
        self._update_velocity_data()
        self._forces_moments(delta=MsgDelta())

        self.vn_prev = 0
        self.ve_prev = 0
        self.vd_prev = 0


    ###################################
    # public functions
    def update(self, delta, wind):
        '''
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        '''
        # get forces and moments acting on rigid bod
        forces_moments = self._forces_moments(delta)

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

        # update the airspeed, angle of attack, and side slip angles using new state
        self._update_velocity_data(wind)
        # update the message class for the true state
        self._update_true_state()


    def sensors(self):
        "Return value of sensors on MAV: gyros, accels, absolute_pressure, dynamic_pressure, GPS"
        ##### TODO #####
        # simulate rate gyros(units are rad / sec)
        self._sensors.gyro_x = self._state.item(10) + SENSOR.gyro_x_bias + np.random.randn()*SENSOR.gyro_sigma # !!! suspect this is the wrong p q and r self._state.item(10) 
        self._sensors.gyro_y = self._state.item(11) + SENSOR.gyro_y_bias + np.random.randn()*SENSOR.gyro_sigma
        self._sensors.gyro_z = self._state.item(12) + SENSOR.gyro_z_bias + np.random.randn()*SENSOR.gyro_sigma

        # simulate accelerometers(units of g)
        mass = MAV.mass
        g = MAV.gravity
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        self._sensors.accel_x = self._forces.item(0)/mass + g*np.sin(theta)             + np.random.randn()*SENSOR.accel_sigma
        self._sensors.accel_y = self._forces.item(1)/mass - g*np.cos(theta)*np.sin(phi) + np.random.randn()*SENSOR.accel_sigma
        self._sensors.accel_z = self._forces.item(2)/mass - g*np.cos(theta)*np.cos(phi) + np.random.randn()*SENSOR.accel_sigma

        # simulate magnetometers
        # magnetic field in provo has magnetic declination of 12.5 degrees
        # and magnetic inclination of 66 degrees not required.
        delta = 12.5
        phi = phi + SENSOR.mag_beta + np.random.randn()*SENSOR.mag_sigma 
        m_0 = phi - delta

        # phi_m = -np.arctan2(self._state.item(1), self._state.item(0))
        # heading = delta + phi_m
        # R = Euler2Rotation(0,theta,psi)
        rx = np.array([np.cos(theta), np.sin(theta)*np.sin(phi),np.sin(theta)*np.cos(phi)])
        ry = np.array([0,             np.cos(phi),            -np.sin(phi)])
        rz = np.array([-np.sin(theta),np.cos(theta)*np.sin(phi),np.cos(theta)*np.cos(phi)])
        self._sensors.mag_x = (rx*m_0).item(0)
        self._sensors.mag_y = (ry*m_0).item(1)
        self._sensors.mag_z = (rz*m_0).item(2) # !!!hmmm...

        # simulate pressure sensors
        Va = self._Va
        # h_AGL = 1387 h of aircraft
        beta_abs = 0 #!!!hmmm... 
        beta = 0     #!!!hmmm... 
        self._sensors.abs_pressure = MAV.rho*g*self._state.item(2) + beta_abs + np.random.randn()*SENSOR.abs_pres_sigma
        self._sensors.diff_pressure = (MAV.rho*Va**2)/2 + beta + SENSOR.diff_pres_sigma*np.random.randn() 

        
        # simulate GPS sensor
        if self._t_gps >= SENSOR.ts_gps:
            wind_n = self._wind.item(0)
            wind_e = self._wind.item(1)
            
            
            self._gps_eta_n = np.random.randn()*SENSOR.gps_n_sigma
            self._gps_eta_e = np.random.randn()*SENSOR.gps_e_sigma
            self._gps_eta_h = np.random.randn()*SENSOR.gps_h_sigma
            vn = np.exp(-SENSOR.gps_k*SENSOR.ts_gps)*self.vn_prev + self._gps_eta_h # !!!!!!!!!!!!!!!!! vn_prev used before declared
            ve = np.exp(-SENSOR.gps_k*SENSOR.ts_gps)*self.ve_prev + self._gps_eta_h
            vd = np.exp(-SENSOR.gps_k*SENSOR.ts_gps)*self.vd_prev + self._gps_eta_h
            self._sensors.gps_n =  self._state.item(1) + self.vn_prev
            self._sensors.gps_e =  self._state.item(2) + self.ve_prev
            self._sensors.gps_h = -self._state.item(3) + self.vd_prev
            self._sensors.gps_Vg = np.sqrt((Va*np.cos(psi)+wind_n)**2 + (Va*np.sin(psi)+wind_e)**2) + np.random.randn()*SENSOR.gps_Vg_sigma**2 #!!! why is this one squared?
            self._sensors.gps_course = np.arctan2((Va*np.sin(psi)+wind_e),(Va*np.cos(psi)+wind_n)) + np.random.randn()*SENSOR.gps_course_sigma**2 #!!! what is Omega p q and r?
            self._t_gps = 0.
            self.vn_prev = vn
            self.ve_prev = ve
            self.vd_prev = vd
        else:
            self._t_gps += self._ts_simulation
        return self._sensors

    def external_set_state(self, new_state):
        self._state = new_state

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):
        """
        for the dynamics xdot = f(x, u), returns f(x, u)
        """
        ##### TODO #####
        #  -- Copy From mav_dynamic_control.py --
        ##### TODO #####
        #  -- Copy From mav_dynamic_forces.py --
        # Airframe parameters from Appendix E
        # physical parameters
        mass = MAV.mass # 11 kg
        # Jx = MAV.Jx # 0.824 # kg-m^2
        Jy = MAV.Jy # 1.135 # kg-m^2
        # Jz = MAV.Jz # 1.759 # kg-m^2
        # Jxz = MAV.Jxz # 0.12 # kg-m^2

        # Changeable parameters for each state include:
        # moments = 0.0
        # products_of_inertia = 0
        # initial_conditions = 0

        # Extract the States
        north = state.item(0)
        east = state.item(1)
        down = state.item(2)
        u = state.item(3)
        v = state.item(4)
        w = state.item(5)
        # pos = np.array([north,east,down]).T # inertial positions
        Q = np.array([state.item(6),state.item(7),state.item(8),state.item(9)]).T
        e0 = state.item(6)
        e1 = state.item(7)
        e2 = state.item(8)
        e3 = state.item(9)
        p = state.item(10)
        q = state.item(11)
        r = state.item(12)
        rotation_rate = np.array([p,q,r]).T


        # Extract Forces/Moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)
        f_xyz = np.array([fx,fy,fz]).T #forces along primary axis in the body frame
        l = forces_moments.item(3) # external moment applied to the airfram about body x axis
        My = forces_moments.item(4) # external moment applied to the airfram about body y axis
        n = forces_moments.item(5) # external moment applied to the airfram about body z axis

        R_bi = Quaternion2Rotation(np.array([e0,e1,e2,e3]))

        Q_yaw_pitch_roll = np.array([
            [0,-p,-q,-r],
            [p, 0, r,-q],
            [q,-r, 0, p],
            [r, q,-p, 0]])
        G = MAV.gamma   # Jx*Jz-Jxz**2
        G1 = MAV.gamma1 # Jxz*(Jx-Jy+Jz)/G
        G2 = MAV.gamma2 # (Jz*(Jz-Jy)+Jxz**2)/G
        G3 = MAV.gamma3 # Jz/G
        G4 = MAV.gamma4 # Jxz/G
        G5 = MAV.gamma5 # Jz-Jx/Jy
        G6 = MAV.gamma6 # Jxz/Jy
        G7 = MAV.gamma7 # ((Jx-Jy)*Jx+Jxz**2)/G
        G8 = MAV.gamma8 # Jx/G

        U = np.array([u,v,w]).T # inertial velocties or Vg expressed in the body frame
        # Position Kinematics
        pos_dot = R_bi @ U

        # Position Dynamics
        U_dot = np.cross(-rotation_rate,U) + 1/mass*f_xyz


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

    def _update_velocity_data(self, wind=np.zeros((6,1))):
        steady_state = wind[0:3]
        gust = wind[3:6]

        ##### TODO #####
        #  -- Copy From mav_dynamic_control.py --
        # convert wind vector from world to body frame (self._wind = ?)
        R_bv = Quaternion2Rotation(np.array([self._state.item(6),self._state.item(7),self._state.item(8),self._state.item(9)]))
        # velocity vector relative to the airmass ([ur , vr, wr]= ?)
        self._wind = R_bv @ steady_state + gust

        # velocity vector relative to the airmass ([ur , vr, wr]= ?)
        # V_ba = np.array([[u-uw, v-vw, w-ww]])
        ur = self._state[3]-self._wind[0]
        vr = self._state[4]-self._wind[1]
        wr = self._state[5]-self._wind[2]
        # compute airspeed (self._Va= ?)
        self._Va = np.sqrt(ur**2 + vr**2 + wr**2).item()

        # compute angle of attack (self._alpha = ?)
        self._alpha = np.arctan(wr/ur).item()
        
        # compute sideslip angle (self._beta = ?)
        self._beta = np.arcsin(vr/(self._Va)).item()

    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """
        ##### TODO ###### 
        # #  -- Copy From mav_dynamic_control.py --
        # extract states (phi, theta, psi, p, q, r)
        e0 = self._state[6]
        e1 = self._state[7]
        e2 = self._state[8]
        e3 = self._state[9]
        p = self._state[10]
        q = self._state[11]
        r = self._state[12]
        mass = MAV.mass
        g = MAV.gravity
        alpha = self._alpha# set in update velocity 
        beta = self._beta #set in update velocity
        Va = self._Va
        c = MAV.c
        b = MAV.b
        alpha0 = MAV.alpha0
        rho = MAV.rho
        S = MAV.S_wing
        delta_a, delta_e, delta_r, delta_t = delta.aileron, delta.elevator, delta.rudder, delta.throttle 
        # compute gravitaional forces ([fg_x, fg_y, fg_z])
        fg_x = mass*g*2*(e1*e3-e2*e0)
        fg_y = mass*g*2*(e2*e3+e1*e0)
        fg_z = mass*g*(e3**2+e0**2-e1**2-e2**2)
        # compute Lift and Drag coefficients (CL, CD)
        sigma_alpha = (1+np.exp(-MAV.M*(alpha-alpha0))+ np.exp(MAV.M*(alpha+alpha0))) / ((1+np.exp(-MAV.M*(alpha-alpha0)))*(1+np.exp(MAV.M*(alpha+alpha0))))
        CL = (1-sigma_alpha)*(MAV.C_L_0+MAV.C_L_alpha*alpha) + sigma_alpha*(2*np.sign(alpha)*(np.sin(alpha))**2*alpha*np.cos(alpha))
        CD = MAV.C_D_p + (((MAV.C_L_0 + MAV.C_L_alpha * alpha))**2) / (np.pi * MAV.e * MAV.AR)
        # compute Lift and Drag Forces (F_lift, F_drag)
        F_lift = 0.5 * rho * Va**2 * S * (CL + MAV.C_L_q * c/(2*Va)*q + MAV.C_L_delta_e*delta_e)
        F_drag = 0.5 * rho * Va**2 * S * (CD + MAV.C_D_q * c/(2*Va)*q + MAV.C_D_delta_e*delta_e)
        # propeller thrust and torque
        # thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, delta.throttle)
        thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, delta_t)

        # compute longitudinal forces in body frame (fx, fz)
        fx = np.cos(alpha) * -F_drag + - np.sin(alpha) * - F_lift
        fz = np.sin(alpha) * -F_drag + np.cos(alpha) * - F_lift
        C1 = 0.5*rho*Va**2*S
        # compute lateral forces in body frame (fy)
        fy = C1*(MAV.C_Y_0 + MAV.C_Y_beta*beta + MAV.C_Y_p*(b/2/Va)*p + MAV.C_Y_r*(b/2/Va)*r + MAV.C_Y_delta_a*delta_a + MAV.C_Y_delta_r*delta_r)

        # compute lateral torques in body frame (Mx, Mz)
        # compute logitudinal torque in body frame (My)
        # roll
        Mx = C1*b*(MAV.C_ell_0 + MAV.C_ell_beta*beta + MAV.C_ell_p*(b/2/Va)*p + MAV.C_ell_r*(b/2/Va)*r + MAV.C_ell_delta_a*delta_a + MAV.C_ell_delta_r*delta_r)
        # pitch
        My = C1*c*(MAV.C_m_0 + MAV.C_m_alpha*alpha + MAV.C_m_q*(c/2/Va)*q + MAV.C_m_delta_e*delta_e)
        # yaw
        Mz = C1*b*(MAV.C_n_0 + MAV.C_n_beta*beta + MAV.C_n_p*(b/2/Va)*p + MAV.C_n_r*(b/2/Va)*r + MAV.C_n_delta_a*delta_a + MAV.C_n_delta_r*delta_r)
        
        # self._forces[0] = fx #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! don't seem right
        # self._forces[1] = fy
        # self._forces[2] = fz

        Fx = fx + thrust_prop + fg_x
        Fy = fy + fg_y
        Fz = fz + fg_z
        self._forces[0] = Fx
        self._forces[1] = Fy
        self._forces[2] = Fz
        # forces_moments = np.array([[0, 0, 0, 0, 0, 0]]).T
        forces_moments = np.array([[Fx, Fy, Fz, Mx, My, Mz]]).T
        return forces_moments

    def _motor_thrust_torque(self, Va, delta_t):
        # compute thrust and torque due to propeller
        ##### TODO #####
        #  -- Copy From mav_dynamic_control.py --
        # map delta_t throttle command(0 to 1) into motor input voltage
        v_in = MAV.V_max*delta_t
        a = MAV.C_Q0 * MAV.rho * np.power(MAV.D_prop,5)/((2.0*np.pi)**2)
        b = MAV.C_Q1 * MAV.rho * np.power(MAV.D_prop,4)/ (2.0*np.pi) * Va + MAV.KQ*MAV.KV/MAV.R_motor
        c = MAV.C_Q2 * MAV.rho * np.power(MAV.D_prop,3) * Va**2 -(MAV.KQ/MAV.R_motor) * v_in + MAV.KQ * MAV.i0
        # use quadratic formula to solve for motor speed, then only consider positve roots 
        Omega_op = (-b+np.sqrt(b**2 - 4*a*c))/(2*a)
        # find the advance ratio propeller speed to airspeed
        J_op = 2*np.pi*Va/(Omega_op*MAV.D_prop)
        # find thrust and torque coefficients
        C_T = MAV.C_T2*J_op**2 + MAV.C_T1*J_op+MAV.C_T0
        C_Q = MAV.C_Q2*J_op**2 + MAV.C_Q1*J_op+MAV.C_Q0
        # add thrust and torque caused by the propeller
        n = Omega_op / (2*np.pi)

        # Angular speed of propeller (omega_p = ?) ??? in radians ???
        # Omega_p = 2*np.pi*n

        # thrust and torque due to propeller
        thrust_prop =  MAV.rho * n**2 * np.power(MAV.D_prop,4) * C_T
        torque_prop =  MAV.rho * n**2 * np.power(MAV.D_prop,5) * C_Q

        return thrust_prop.item(), torque_prop.item()

    def _update_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        pdot = Quaternion2Rotation(self._state[6:10]) @ self._state[3:6]
        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = self._Va
        self.true_state.alpha = self._alpha
        self.true_state.beta = self._beta
        self.true_state.phi = phi
        self.true_state.theta = theta
        self.true_state.psi = psi
        self.true_state.Vg = np.linalg.norm(pdot)
        self.true_state.gamma = np.arcsin(pdot.item(2) / self.true_state.Vg)
        self.true_state.chi = np.arctan2(pdot.item(1), pdot.item(0))
        self.true_state.p = self._state.item(10)
        self.true_state.q = self._state.item(11)
        self.true_state.r = self._state.item(12)
        self.true_state.wn = self._wind.item(0)
        self.true_state.we = self._wind.item(1)
        self.true_state.bx = SENSOR.gyro_x_bias
        self.true_state.by = SENSOR.gyro_y_bias
        self.true_state.bz = SENSOR.gyro_z_bias
        self.true_state.camera_az = self._state.item(13)
        self.true_state.camera_el = self._state.item(14)
