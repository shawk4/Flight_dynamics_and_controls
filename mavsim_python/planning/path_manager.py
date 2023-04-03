"""
mavsim_python: drawing tools
    - Beard & McLain, PUP, 2012
    - Update history:
        4/15/2019 - RWB
        3/30/2022 - RWB
"""

import numpy as np
import sys
sys.path.append('..')
from planning.dubins_parameters import DubinsParameters
from message_types.msg_path import MsgPath


class PathManager:
    def __init__(self):
        # message sent to path follower
        self.path = MsgPath()
        # pointers to previous, current, and next waypoints
        self.ptr_previous = 0
        self.ptr_current = 1
        self.ptr_next = 2
        self.num_waypoints = 0
        self.halfspace_n = np.inf * np.ones((3,1))
        self.halfspace_r = np.inf * np.ones((3,1))
        # state of the manager state machine
        self.manager_state = 1
        self.manager_requests_waypoints = True
        self.dubins_path = DubinsParameters()

    def update(self, waypoints, radius, state):
        if waypoints.num_waypoints == 0:
            self.manager_requests_waypoints = True
        if self.manager_requests_waypoints is True and waypoints.flag_waypoints_changed is True:
            self.manager_requests_waypoints = False
        if waypoints.type == 'straight_line':
            self.line_manager(waypoints, state)
        elif waypoints.type == 'fillet':
            self.fillet_manager(waypoints, radius, state)
        elif waypoints.type == 'dubins':
            self.dubins_manager(waypoints, radius, state)
        else:
            print('Error in Path Manager: Undefined waypoint type.')
        return self.path

    def line_manager(self, waypoints, state):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # if the waypoints have changed, update the waypoint pointer
        ##### TODO ######
        # Use functions - self.initialize_pointers(), self.construct_line()
        # self.inHalfSpace(mav_pos), self.increment_pointers(), self.construct_line()

        # Use variables - self.ptr_current, self.manager_requests_waypoints,
        # waypoints.__, radius

        if waypoints.flag_waypoints_changed:
            self.initialize_pointers()
            self.construct_line(waypoints)
            waypoints.flag_waypoints_changed = False

        if self.inHalfSpace(mav_pos):
             self.increment_pointers()
             self.construct_line(waypoints)



    def fillet_manager(self, waypoints, radius, state):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # if the waypoints have changed, update the waypoint pointer
        ##### TODO ######
        # Use functions - self.initialize_pointers(), self.construct_fillet_line(),
        # self.inHalfSpace(), self.construct_fillet_circle(), self.increment_pointers()

        # Use variables self.num_waypoints, self.manager_state, self.ptr_current
        # self.manager_requests_waypoints, waypoints.__, radius
        self.num_waypoints = waypoints.num_waypoints
        if self.num_waypoints >= 3:
            if waypoints.flag_waypoints_changed:
                self.initialize_pointers()
                self.construct_fillet_line(waypoints, radius)
                waypoints.flag_waypoints_changed = False
                self.manager_state = 1

                
            
            if self.inHalfSpace(mav_pos):
                if self.manager_state == 1:
                    self.construct_fillet_circle(waypoints, radius)
                else:
                    self.increment_pointers()
                    self.construct_fillet_line(waypoints, radius)
        else:
            print("Need more waypoints for fillet path manager.")


    def dubins_manager(self, waypoints, radius, state):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        R = radius 
        # if the waypoints have changed, update the waypoint pointer
        ##### TODO #####
        # Use functions - self.initialize_pointers(), self.dubins_path.update(),
        # self.construct_dubins_circle_start(), self.construct_dubins_line(),
        # self.inHalfSpace(), self.construct_dubins_circle_end(), self.increment_pointers(),

        # Use variables - self.num_waypoints, self.dubins_path, self.ptr_current,
        # self.ptr_previous, self.manager_state, self.manager_requests_waypoints,
        # waypoints.__, radius

        self.num_waypoints = waypoints.num_waypoints
        if self.num_waypoints >= 3:
            if waypoints.flag_waypoints_changed:
                self.initialize_pointers()
                waypoints.flag_waypoints_changed = False
                self.manager_state = 1                



            # find L, c_s, lam_s, c_e, lam_e, z1, q1, z2, z3, q3 ## 
            self.dubins_path.update(waypoints.ned[:,self.ptr_previous:self.ptr_previous+1], waypoints.course[self.ptr_previous], 
                                    waypoints.ned[:, self.ptr_current:self.ptr_current+1 ], waypoints.course[self.ptr_current ], R) 

            if self.manager_state == 1:
                # flag = 2; c = cs; rho = R; lam = lam_s
                self.construct_dubins_circle_start(start_sign=-1)
                if self.inHalfSpace(mav_pos):
                    self.manager_state = 2

            elif self.manager_state == 2:
                self.construct_dubins_circle_start(start_sign=1)
                if self.inHalfSpace(mav_pos):
                    self.manager_state = 3

            elif self.manager_state == 3: # make and follow a line in state 3
                # flag = 1; r = z1; q = q1
                self.construct_dubins_line()
                if self.inHalfSpace(mav_pos):
                    self.manager_state = 4

            elif self.manager_state == 4: # make and follow an end orbit in state 4
                # flag = 2; c = ce; rho = R; lam = lam_e
                self.construct_dubins_circle_end(end_sign=-1)
                if self.inHalfSpace(mav_pos):
                    self.manager_state = 5

            elif self.manager_state == 5: 
                self.construct_dubins_circle_end(end_sign=1)
                if self.inHalfSpace(mav_pos):
                    self.manager_state = 1  
                    self.increment_pointers()


            # find L, c_s, lam_s, c_e, lam_e, z1, q1, z2, z3, q3 ## 
            self.dubins_path.update(waypoints.ned[:,self.ptr_previous:self.ptr_previous+1], waypoints.course[self.ptr_previous], 
                                    waypoints.ned[:, self.ptr_current:self.ptr_current+1 ], waypoints.course[self.ptr_current ], R) 
        else:
            print("Need more waypoints for fillet path manager.")



    def initialize_pointers(self):
        if self.num_waypoints >= 3:
            ##### TODO #####
            self.ptr_previous = 0
            self.ptr_current = 1
            self.ptr_next = 2
        else:
            print('Error Path Manager: need at least three waypoints')

    def increment_pointers(self):
        ##### TODO #####
        # print("before: " + str(self.ptr_previous) +";"+ str(self.ptr_current) +";"+ str(self.ptr_next))
        self.ptr_previous = self.ptr_current
        self.ptr_current = self.ptr_next # !!! this makes sense but I am not sure... 
        self.ptr_next += 1

        if self.ptr_next == self.num_waypoints:
            self.ptr_next = 0
        # print("after: " + str(self.ptr_previous) +";"+ str(self.ptr_current) +";"+ str(self.ptr_next))


    def construct_line(self, waypoints):
        previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous+1]
        current = waypoints.ned[:, self.ptr_current:self.ptr_current+1]
        next = waypoints.ned[:, self.ptr_next:self.ptr_next+1]
        ##### TODO #####

        qi_min_1 = (current-previous)/ np.linalg.norm(current-previous)
        qi = (next - current) / np.linalg.norm(next - current)
        # update halfspace variables
        self.halfspace_n = (qi_min_1 + qi) / np.linalg.norm(qi_min_1 - qi)
        self.halfspace_r = current
        
        # Update path variables
        self.path.line_direction = qi_min_1
        self.path.line_origin = previous
        # self.path.orbit_center =
        # self.path.orbit_radius = 
        # self.path.orbit_direction = 'CW'
        self.path.plot_updated = False
        self.path.type = "line"

    def construct_fillet_line(self, waypoints, radius):
        previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous+1]
        current = waypoints.ned[:, self.ptr_current:self.ptr_current+1]
        next = waypoints.ned[:, self.ptr_next:self.ptr_next+1]
        ##### TODO #####
        
        self.manager_state = 1 # Flag
        R = radius
        qi_min_1 = (current-previous)/ np.linalg.norm(current-previous)
        qi = (next - current) / np.linalg.norm(next - current)
        var_rho = np.arccos(-qi_min_1.T@qi)#np.dot(qi_min_1[:,0], qi[:,0])
        # update halfspace variables
        self.halfspace_n = qi_min_1
        self.halfspace_r = current - (R/np.tan(var_rho/2))*qi_min_1
        
        # Update path variables
        self.path.line_direction = qi_min_1
        self.path.line_origin = previous
        # self.path.orbit_center =
        # self.path.orbit_radius = 
        # self.path.orbit_direction = 'CW'
        self.path.plot_updated = False
        self.path.type = "line"

    def construct_fillet_circle(self, waypoints, radius):
        previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous+1]
        current = waypoints.ned[:, self.ptr_current:self.ptr_current+1]
        next = waypoints.ned[:, self.ptr_next:self.ptr_next+1]
        ##### TODO #####
        
        self.manager_state = 2 # Flag
        R = radius
        qi_min_1 = (current-previous)/ np.linalg.norm(current-previous)
        qi = (next - current) / np.linalg.norm(next - current)
        var_rho = np.arccos(-qi_min_1.T@qi) # np.dot(qi_min_1[:,0], qi[:,0]) previous method is not sensitve to signs
        c = current - (R/np.sin(var_rho/2))*(qi_min_1-qi)/np.linalg.norm(qi_min_1-qi)
        rho = R
        lam = np.sign(qi_min_1.item(0)*qi.item(1)-qi_min_1.item(1)*qi.item(0))
        z = current + (R/np.tan(var_rho/2))*qi

        # update halfspace variables
        self.halfspace_n = qi
        self.halfspace_r = z
        
        if lam > 0:
            orbit_dir = "CW"
        else:
            orbit_dir = "CCW"

        # Update path variables
        # self.path.line_direction = 
        # self.path.line_origin = 
        self.path.orbit_center = c
        self.path.orbit_radius = R
        self.path.orbit_direction = orbit_dir
        self.path.plot_updated = False
        self.path.type = "orbit"

    def construct_dubins_circle_start(self, start_sign):
        ##### TODO #####
        # update halfspace variables

        self.halfspace_n = start_sign*self.dubins_path.get_q1() # a normal vector to the point  n <-> q ???
        self.halfspace_r =  self.dubins_path.get_z1() # a point  r <-> z ???
        
        if self.dubins_path.get_dir_s() > 0:
            orbit_dir = "CW"
        else:
            orbit_dir = "CCW"

        # Update path variables
        # return flag, r, q, c, rho, lam
        # self.path.line_direction = 
        # self.path.line_origin = 
        self.path.orbit_center = self.dubins_path.get_center_s()
        self.path.orbit_radius = self.dubins_path.get_radius()
        self.path.orbit_direction = orbit_dir
        self.path.plot_updated = False
        self.path.type = "orbit"

    def construct_dubins_line(self):
        # update halfspace variables

        self.halfspace_n = self.dubins_path.get_q1()  # a normal vector to the point  n <-> q ???
        self.halfspace_r =  self.dubins_path.get_z2() # a point  r <-> z ???
        
        # Update path variables
        # return flag, r, q, c, rho, lam
        self.path.line_direction = self.dubins_path.get_q1()
        self.path.line_origin = self.dubins_path.get_z2()
        # self.path.orbit_center = self.dubins_path.get_center_s()
        # self.path.orbit_radius = self.dubins_path.get_radius()
        # self.path.orbit_direction = self.dubins_path.get_dir_s()
        self.path.plot_updated = False
        self.path.type = "line"

    def construct_dubins_circle_end(self, end_sign):
        ##### TODO #####
        self.halfspace_n = end_sign*self.dubins_path.get_q3() # a normal vector to the point r... n <-> q ???
        self.halfspace_r = self.dubins_path.get_z3() # a point  r <-> z ???

        if self.dubins_path.get_dir_e() > 0:
            orbit_dir = "CW"
        else:
            orbit_dir = "CCW"
        
        # Update path variables
        # return flag, r, q, c, rho, lam
        # self.path.line_direction = 
        # self.path.line_origin = 
        self.path.orbit_center = self.dubins_path.get_center_e()
        self.path.orbit_radius = self.dubins_path.get_radius()
        self.path.orbit_direction = orbit_dir
        self.path.plot_updated = False
        self.path.type = "orbit"
        

    def inHalfSpace(self, pos):
        if (pos-self.halfspace_r).T @ self.halfspace_n >= 0:
            print("enter half space")
            return True
        else:
            return False

