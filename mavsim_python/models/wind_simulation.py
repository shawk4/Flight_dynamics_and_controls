"""
Class to determine wind velocity at any given moment,
calculates a steady wind speed and uses a stochastic
process to represent wind gusts. (Follows section 4.4 in uav book)
"""
import sys
sys.path.append('..')
from tools.transfer_function import TransferFunction
import numpy as np
import parameters.aerosonde_parameters as MAV


class WindSimulation:
    def __init__(self, Ts, gust_flag = True, steady_state = np.array([[0., 0., 0.]]).T):
        # steady state wind defined in the inertial frame
        self._steady_state = steady_state
        ##### TODO #####

        #   Dryden gust model parameters (pg 61 UAV book)
        # low altitude light turbulance
        # self.altitude = 50 # m
        # self.L_u = 200 # m
        # self.L_v = 200 # m
        # self.L_w = 50 # m
        # self.sigma_u = 1.06 # m
        # self.sigma_v = 1.06 # m
        # self.sigma_w = 0.7 # m
        

        L_u = 200 # m
        L_v = 200 # m
        L_w = 50 # m
        sigma_u = 1.06 # m
        sigma_v = 1.06 # m
        sigma_w = 0.7 # m
        Va = MAV.Va0 
        

        C1 = sigma_u*np.sqrt(2*Va/(np.pi*L_u))
        C2 = sigma_v*np.sqrt(3*Va/(np.pi*L_v))
        C3 = sigma_w*np.sqrt(3*Va/(np.pi*L_w))



        # Dryden transfer functions (section 4.4 UAV book) - Fill in proper num and den
        self.u_w = TransferFunction(num=np.array([[ C1 ]]), den=np.array([[1,Va/L_u]]),Ts=Ts)
        self.v_w = TransferFunction(num=np.array([[ C2, C2*Va/(np.sqrt(3)*L_v)]]), den=np.array([[1, 2*Va/L_v, (Va/L_v)**2]]),Ts=Ts)
        self.w_w = TransferFunction(num=np.array([[ C3, C3*Va/(np.sqrt(3)*L_w)]]), den=np.array([[1, 2*Va/L_w, (Va/L_w)**2]]),Ts=Ts)
        self._Ts = Ts

    def update(self):
        # returns a six vector.
        #   The first three elements are the steady state wind in the inertial frame
        #   The second three elements are the gust in the body frame
        gust = np.array([[self.u_w.update(np.random.randn())],
                         [self.v_w.update(np.random.randn())],
                         [self.w_w.update(np.random.randn())]])
        return np.concatenate(( self._steady_state, gust ))

