"""
mavsim_python: mav viewer (for chapter 2)
    - Beard & McLain, PUP, 2012
    - Update history:
        1/15/2019 - RWB
        4/15/2019 - BGM
        3/31/2020 - RWB
"""
import sys
sys.path.append("..")
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import pyqtgraph.Vector as Vector
from viewers.draw_mav import DrawMav


class MavViewer():
    def __init__(self, app):
        self.scale = 100
        # initialize Qt gui application and window
        self.app = app  # initialize QT, external so that only one QT process is running
        self.window = gl.GLViewWidget()  # initialize the view object
        self.window.setWindowTitle('MAV Viewer')
        grid = gl.GLGridItem() # make a grid to represent the ground
        grid.scale(20, 20, 20) # set the size of the grid (distance between each line)
        self.window.addItem(grid) # add grid to viewer
        self.window.setCameraPosition(distance=200) # distance from center of plot to camera
        self.window.setBackgroundColor('k')  # set background color to black
        self.window.setGeometry(0, 0, 750, 750)  # args: upper_left_x, upper_right_y, width, height
        # center = self.window.cameraPosition()
        # center.setX(250)
        # center.setY(250)
        # center.setZ(0)
        # self.window.setCameraPosition(pos=center, distance=self.scale, elevation=50, azimuth=-90)
        self.window.show()  # display configured window
        # self.window.raise_() # bring window to the front

        self.plot_initialized = False # has the mav been plotted yet?
        self.mav_plot = []

    def update(self, state):
        # initialize the drawing the first time update() is called
        if not self.plot_initialized:
            self.mav_plot = DrawMav(state, self.window)
            self.plot_initialized = True
        # else update drawing on all other calls to update()
        else:
            self.mav_plot.update(state)
        # update the center of the camera view to the mav location
        view_location = Vector(state.east, state.north, state.altitude)  # defined in ENU coordinates
        self.window.opts['center'] = view_location
        # redraw
        self.process_app()
    
    def process_app(self):
        self.app.processEvents()

    def clear_viewer(self):
        self.window.clear()

