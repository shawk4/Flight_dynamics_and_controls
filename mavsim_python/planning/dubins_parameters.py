# dubins_parameters
#   - Dubins parameters that define path between two configurations
#
# mavsim_matlab 
#     - Beard & McLain, PUP, 2012
#     - Update history:  
#         3/26/2019 - RWB
#         4/2/2020 - RWB
#         3/30/2022 - RWB

import numpy as np
import sys
sys.path.append('..')


class DubinsParameters:

    def update(self, ps, chis, pe, chie, R):
         self.p_s = ps
         self.chi_s = chis
         self.p_e = pe
         self.chi_e = chie
         self.radius = R
         self.compute_parameters()

    def compute_parameters(self):
        ps = self.p_s
        pe = self.p_e
        chis = self.chi_s
        chie = self.chi_e
        R = self.radius
        ell = np.linalg.norm(ps[0:2] - pe[0:2])

        ##### TODO #####
        if ell < 2 * R:
            print('Error in Dubins Parameters: The distance between nodes must be larger than 2R.')
        else:
            # compute start and end circles
            crs = ps + R * rotz( np.pi/2) @ np.array([[np.cos(chis), np.sin(chis), 0]]).T
            cls = ps + R * rotz(-np.pi/2) @ np.array([[np.cos(chis), np.sin(chis), 0]]).T
            cre = pe + R * rotz( np.pi/2) @ np.array([[np.cos(chie), np.sin(chie), 0]]).T
            cle = pe + R * rotz(-np.pi/2) @ np.array([[np.cos(chie), np.sin(chie), 0]]).T

            # compute L1
            var_t_L1 = np.arctan2(cre.item(1)-crs.item(1), cre.item(0)-crs.item(0))
            L1 = np.linalg.norm(crs - cre) + R*mod(2*np.pi + mod(var_t_L1-np.pi/2) - mod(chis-np.pi/2)) + \
                                             R*mod(2*np.pi + mod(chie-np.pi/2)  - mod(var_t_L1-np.pi/2))

            # compute L2
            var_t_L2 = np.arctan2(cle.item(1)-crs.item(1), cle.item(0)-crs.item(0))
            l2 = np.linalg.norm(cle-crs)
            var_t2L2 = var_t_L2 - np.pi/2 + np.arcsin(2*R/l2)
            L2 = np.sqrt(l2**2 - 4*R**2) + R*mod(2*np.pi + mod(var_t2L2)       - mod(chis-np.pi/2)) + \
                                           R*mod(2*np.pi + mod(var_t2L2+np.pi) - mod(chie+np.pi/2))

            # compute L3
            var_t_L3 = np.arctan2(cre.item(1)-cls.item(1), cre.item(0)-cls.item(0))
            l3 = np.linalg.norm(cre-cls)
            var_t2L3 = np.arccos(2*R/l3) # !!! this seemed wrong in the book it was var_t2 again may cause issues seems there are four of them var_t_l1 var_t_l2...
            L3 = np.sqrt(l3**2 - 4*R**2) + R*mod(2*np.pi + mod(chis+np.pi/2) - mod(var_t_L3 + var_t2L3)) + \
                                           R*mod(2*np.pi + mod(chie-np.pi/2)   - mod(var_t_L3 + var_t2L3 - np.pi))

            # compute L4
            var_t_L4 = np.arctan2(cle.item(1)-cls.item(1), cle.item(0)-cls.item(0))
            L4 = np.linalg.norm(cls - cle) + R*mod(2*np.pi + mod(chis + np.pi/2) - mod(var_t_L4 + np.pi/2)) + \
                                             R*mod(2*np.pi + mod(var_t_L4 + np.pi/2) - mod(chie+np.pi/2))

            # L is the minimum distance
            L = np.min([L1, L2, L3, L4])
            min_idx = np.argmin([L1, L2, L3, L4])

            e1 = np.array([[1,0,0]]).T
            # L1
            if min_idx == 0:
                cs=crs; lam_s=1; ce=cre; lam_e=1
                q1 = (ce-cs)/np.linalg.norm(ce-cs)
                z1 = cs + R * rotz(-np.pi/2)@q1
                z2 = ce + R * rotz(-np.pi/2)@q1
            # L2
            elif min_idx == 1:
                cs=crs; lam_s=1; ce=cle; lam_e=-1
                l = np.linalg.norm(ce-cs)
                # var_t = angle(ce-cs) # !!! going on a hunch here that var_t 1 and 2 from above are the same as these but need to be different for each case so var_tL2
                cse = (ce - cs) / np.linalg.norm(ce - cs)
                var_t_rl = np.arccos(e1.T@cse).item(0)# angle(ce-cs)

                var_t2_rl = var_t_rl - np.pi/2 + np.arcsin(2*R/l) 
                q1 = rotz(var_t2_rl + np.pi/2)@e1 
                z1 = cs + R * rotz(var_t2_rl)@e1
                z2 = ce + R * rotz(var_t2_rl + np.pi)@e1
            # L3
            elif min_idx == 2:
                cs=cls; lam_s=-1; ce=cre; lam_e=1
                l = np.linalg.norm(ce-cs)
                # var_t = angle(ce-cs) # !!! still going on the hunch above so var_tL3 here probably wrong definityly makes a difference trying var_t
                cse = (ce - cs) / np.linalg.norm(ce - cs)
                var_t_lr = np.arccos(e1.T@cse).item(0) # angle(ce-cs)
                
                var_t2_lr = np.arccos(2*R/l)
                q1 = rotz(var_t_lr + var_t2_lr - np.pi/2)@e1
                z1 = cs + R * rotz(var_t_lr + var_t2_lr)@e1
                z2 = ce + R * rotz(var_t_lr + var_t2_lr - np.pi)@e1
            # L4
            elif min_idx == 3:
                cs=cls; lam_s=-1; ce=cle; lam_e=-1
                q1 = (ce-cs)/np.linalg.norm(ce-cs)
                z1 = cs + R * rotz(np.pi/2)@q1
                z2 = ce + R * rotz(np.pi/2)@q1
            z3 = pe
            q3 = rotz(chie)@e1
            
            # L, cs, lam_s, ce, lam_e, z1, q1, z2, z3, q3
            self.length = L
            self.center_s = cs
            self.dir_s = lam_s
            self.center_e = ce
            self.dir_e = lam_e
            self.r1 = z1
            self.r2 = z2
            self.r3 = z3
            self.n1 = q1
            self.n3 = q3

    def compute_points(self):
        ##### TODO ##### - uncomment lines and remove last line
        Del = 0.1  # distance between point

        # points along start circle
        th1 = np.arctan2(self.p_s.item(1) - self.center_s.item(1),
                         self.p_s.item(0) - self.center_s.item(0))
        th1 = mod(th1)
        th2 = np.arctan2(self.r1.item(1) - self.center_s.item(1),
                         self.r1.item(0) - self.center_s.item(0))
        th2 = mod(th2)
        th = th1
        theta_list = [th]
        if self.dir_s > 0:
            if th1 >= th2:
                while th < th2 + 2*np.pi - Del:
                    th += Del
                    theta_list.append(th)
            else:
                while th < th2 - Del:
                    th += Del
                    theta_list.append(th)
        else:
            if th1 <= th2:
                while th > th2 - 2*np.pi + Del:
                    th -= Del
                    theta_list.append(th)
            else:
                while th > th2 + Del:
                    th -= Del
                    theta_list.append(th)

        points = np.array([[self.center_s.item(0) + self.radius * np.cos(theta_list[0]),
                            self.center_s.item(1) + self.radius * np.sin(theta_list[0]),
                            self.center_s.item(2)]])
        for angle in theta_list:
            new_point = np.array([[self.center_s.item(0) + self.radius * np.cos(angle),
                                   self.center_s.item(1) + self.radius * np.sin(angle),
                                   self.center_s.item(2)]])
            points = np.concatenate((points, new_point), axis=0)

        # points along straight line
        sig = 0
        while sig <= 1:
            new_point = np.array([[(1 - sig) * self.r1.item(0) + sig * self.r2.item(0),
                                   (1 - sig) * self.r1.item(1) + sig * self.r2.item(1),
                                   (1 - sig) * self.r1.item(2) + sig * self.r2.item(2)]])
            points = np.concatenate((points, new_point), axis=0)
            sig += Del

        # points along end circle
        th2 = np.arctan2(self.p_e.item(1) - self.center_e.item(1),
                         self.p_e.item(0) - self.center_e.item(0))
        th2 = mod(th2)
        th1 = np.arctan2(self.r2.item(1) - self.center_e.item(1),
                         self.r2.item(0) - self.center_e.item(0))
        th1 = mod(th1)
        th = th1
        theta_list = [th]
        if self.dir_e > 0:
            if th1 >= th2:
                while th < th2 + 2 * np.pi - Del:
                    th += Del
                    theta_list.append(th)
            else:
                while th < th2 - Del:
                    th += Del
                    theta_list.append(th)
        else:
            if th1 <= th2:
                while th > th2 - 2 * np.pi + Del:
                    th -= Del
                    theta_list.append(th)
            else:
                while th > th2 + Del:
                    th -= Del
                    theta_list.append(th)
        for angle in theta_list:
            new_point = np.array([[self.center_e.item(0) + self.radius * np.cos(angle),
                                   self.center_e.item(1) + self.radius * np.sin(angle),
                                   self.center_e.item(2)]])
            points = np.concatenate((points, new_point), axis=0)
        # points = np.zeros((5,3))
        return points

    def get_L(self):
        return self.length
    def get_center_s(self):
        return self.center_s
    def get_center_e(self):
        return self.center_e
    def get_dir_s(self):
        return self.dir_s
    def get_dir_e(self):
        return self.dir_e    
    def get_q1(self): # q1= q2
        return self.n1
    def get_q3(self):
        return self.n3
    def get_z1(self):
        return self.r1
    def get_z2(self):
        return self.r2
    def get_z3(self):
        return self.r3
    def get_radius(self):
        return self.radius
    
def rotz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])


def mod(x):
    while x < 0:
        x += 2*np.pi
    while x > 2*np.pi:
        x -= 2*np.pi
    return x


