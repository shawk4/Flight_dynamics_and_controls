"""
point_gimbal
    - point gimbal at target
part of mavsim
    - Beard & McLain, PUP, 2012
    - Update history:  
        3/31/2022 - RWB
"""
import sys
sys.path.append('..')
import numpy as np
from tools.rotations import Euler2Rotation
import parameters.camera_parameters as CAM


class Gimbal:
    def pointAtGround(self, mav):
        ###### TODO #######
        # desired inertial frame vector points down
        
        # rotate line-of-sight vector into body frame and normalize
        
        ell = np.array([[0],[0],[0]])
        return( self.pointAlongVector(ell, mav.camera_az, mav.camera_el) )

    def pointAtPosition(self, mav, target_position):
        ###### TODO #######
        # line-of-sight vector in the inertial frame
        
        # rotate line-of-sight vector into body frame and normalize
        ell = np.array([[0],[0],[0]])
        return( self.pointAlongVector(ell, mav.camera_az, mav.camera_el) )

    def pointAlongVector(self, ell, azimuth, elevation):
        # point gimbal so that optical axis aligns with unit vector ell
        # ell is assumed to be aligned in the body frame
        # given current azimuth and elevation angles of the gimbal

        ##### TODO #####
        # compute control inputs to align gimbal
        azimuth_c = np.arctan2(ell.item(1), ell.item(0))
        elevation_c = np.arcsin(ell.item(2))
        # proportional control for gimbal

        u_az = CAM.k_az*(azimuth_c - azimuth)
        u_el = CAM.k_el*(elevation_c - elevation)
        return( np.array([[u_az], [u_el]]) )




