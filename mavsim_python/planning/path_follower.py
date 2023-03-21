import numpy as np
from math import sin, cos
import sys

sys.path.append('..')
from message_types.msg_autopilot import MsgAutopilot
from tools.wrap import wrap


class PathFollower:
    def __init__(self):
        ##### TODO #####
        self.chi_inf = 90*np.pi/180  # approach angle for large distance from straight-line path
        self.k_path = 0.02  # path gain for straight-line path following
        self.k_orbit = 5.0  # path gain for orbit following
        self.gravity = 9.82
        self.autopilot_commands = MsgAutopilot()  # message sent to autopilot

    def update(self, path, state):
        if path.type == 'line':
            self._follow_straight_line(path, state)
        elif path.type == 'orbit':
            self._follow_orbit(path, state)
        return self.autopilot_commands

    def _follow_straight_line(self, path, state):
        ##### TODO #####
        #airspeed command
        self.autopilot_commands.airspeed_command = path.airspeed
        chi_q = np.arctan2(path.line_direction.item(1),path.line_direction.item(0))
        
        # smallest angle turn logic
        while ((chi_q-state.chi) < -np.pi):
            chi_q += 2*np.pi
        while ((chi_q-state.chi) > np.pi):
            chi_q -= 2*np.pi

        # transform from inertial to path frame
        R_pi = np.array([[np.cos(chi_q) , np.sin(chi_q), 0],
                         [-np.sin(chi_q), np.cos(chi_q), 0],
                         [0             ,0             ,1]])
        p = np.array([[state.north,state.east,state.altitude]]).T
        r = path.line_origin
        ep = R_pi @ (p - r)
        # epy = -np.sin(chi_q)*(state.north-path.line_origin.item(0)) + cos(chi_q)*(state.east-path.line_origin.item(1))
        # print(epy-ep.item(1))
        # course command
        self.autopilot_commands.course_command = chi_q - self.chi_inf*2/np.pi*np.arctan(self.k_path*ep.item(1))
        
        q = path.line_direction
        qn, qe, qd = q
        rn, re, rd = r
        # qd = -qd # !!! unsure about this sign
        epi = p - r
        ki = np.array([[0,0,1]])
        n = np.linalg.norm(np.cross(ki,q.T)) 
        s = epi -(epi*n)*n
        sn,se,sd = s
        # altitude command
        self.autopilot_commands.altitude_command = (-rd - np.sqrt(sn**2 + se**2)*(qd/(np.sqrt(qn**2 +qe**2)))).item(0)

        # feedforward roll angle for straight line is zero
        self.autopilot_commands.phi_feedforward = 0

    def _follow_orbit(self, path, state):
        ##### TODO #####
        # airspeed command
        self.autopilot_commands.airspeed_command = path.airspeed
        cn, ce, cd = path.orbit_center
        p = path.orbit_radius 
        pn = state.north
        pe = state.east
        pd = state.altitude
        if path.orbit_direction == "CCW":
            gam = -1
        else:
            gam = 1
        
        d = np.sqrt((pn-cn)**2 + (pe-ce)**2)
        ro = np.arctan2(pe-ce,pn-cn) #+ 2*np.pi*m pretty sure this m is picked by the following smallest angle turn logic.
        
        # smallest angle turn logic
        while ((ro - state.chi) < -np.pi):
            ro += 2*np.pi
        while ((ro  - state.chi) > np.pi):
            ro  -= 2*np.pi

        # course command
        self.autopilot_commands.course_command = (ro + gam*(np.pi/2 + np.arctan(self.k_orbit*(d-p)/p))).item(0)

        # altitude command
        self.autopilot_commands.altitude_command = -cd.item(0)
        
        # roll feedforward command
        self.autopilot_commands.phi_feedforward = (gam*np.arctan(state.Vg**2/self.gravity/p/np.cos(state.chi-state.psi))).item(0) # coordinated turn and ramping changes and gamma




