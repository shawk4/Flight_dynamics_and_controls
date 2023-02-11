"""
scsim_python: spacecraft viewer (for appendix C)
    - Beard & McLain, PUP, 2012
    - Update history:
        1/13/2021 - TWM
"""
import sys
sys.path.append("..")
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import pyqtgraph.Vector as Vector
from viewers.draw_spacecraft import DrawSpacecraft

class SpaceCraftViewer():
    def __init__(self):
        # initialize Qt gui application and window
        self.scale = 100
        self.app = pg.QtWidgets.QApplication([])  # initialize QT
        self.window = gl.GLViewWidget()  # initialize the view object
        self.window.setWindowTitle('Spacecraft Viewer')
        self.window.setGeometry(0, 0, 1000, 1000)  # args: upper_left_x, upper_right_y, width, height
        grid = gl.GLGridItem() # make a grid to represent the ground
        grid.scale(20, 20, 20) # set the size of the grid (distance between each line)
        self.window.addItem(grid) # add grid to viewer
        self.window.setCameraPosition(distance=200) # distance from center of plot to camera
        self.window.setBackgroundColor('k')  # set background color to black
        self.window.show()  # display configured window
        self.window.raise_() # bring window to the front
        self.plot_initialized = False # has the spacecraft been plotted yet?
        self.sc_plot = []

    def update(self, state):
        # initialize the drawing the first time update() is called
        if not self.plot_initialized:
            self.sc_plot = DrawSpacecraft(state, self.window)
            self.plot_initialized = True
        # else update drawing on all other calls to update()
        else:
            self.sc_plot.update(state)
        # update the center of the camera view to the spacecraft location
        view_location = Vector(state.east, state.north, state.altitude)  # defined in ENU coordinates
        self.window.opts['center'] = view_location
        # redraw

    def process_app(self):
        self.app.processEvents()
