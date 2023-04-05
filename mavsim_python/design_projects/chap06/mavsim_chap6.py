"""
mavsim_python
    - Chapter 6 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        2/5/2019 - RWB
        2/24/2020 - RWB
"""
import sys
sys.path.append('../..')
import numpy as np
import pyqtgraph as pg
import parameters.simulation_parameters as SIM
from tools.signals import Signals
from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from MAV_control.autopilot import Autopilot
# from MAV_control.autopilot_lqr import Autopilot
# from MAV_control.autopilot_tecs import Autopilot
from viewers.mav_viewer import MavViewer
from viewers.data_viewer import DataViewer
from tools.quit_listener import QuitListener

from models.trim import compute_trim
from models.compute_models import compute_model
COMPUTE_MODEL = True

quitter = QuitListener()

VIDEO = False
PLOTS = True
ANIMATION = True
SAVE_PLOT_IMAGE = False

# video initialization
if VIDEO is True:
    from viewers.video_writer import VideoWriter
    video = VideoWriter(video_name="chap6_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

#initialize the visualization
if ANIMATION or PLOTS:
    app = pg.QtWidgets.QApplication([]) # use the same main process for Qt applications
if ANIMATION:
    mav_view = MavViewer(app=app)  # initialize the mav viewer
if PLOTS:
    # initialize view of data plots
    data_view = DataViewer(app=app,dt=SIM.ts_simulation, plot_period=SIM.ts_plot_refresh, 
                           data_recording_period=SIM.ts_plot_record_data, time_window_length=30)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)

# autopilot commands
from message_types.msg_autopilot import MsgAutopilot
commands = MsgAutopilot()
Va_command = Signals(dc_offset=25.0,
                     amplitude=3.0,
                     start_time=2.0,
                     frequency=0.01)
altitude_command = Signals(dc_offset=100.0,
                           amplitude=20.0, # 20.0
                           start_time=0.0,
                           frequency=0.02) # 0.02
course_command = Signals(dc_offset=np.radians(180), # 180
                         amplitude=np.radians(45), # 45
                         start_time=5.0,
                         frequency=0.015) # 0.015


# use compute_trim function to compute trim state and trim input
Va = 15.
gamma = 0.*np.pi/180.
trim_state, trim_input = compute_trim(mav, Va, gamma)
mav._state = trim_state  # set the initial state of the mav to the trim state
delta_trim = trim_input  # set input to constant constant trim input

# # compute the state space model linearized about trim
if COMPUTE_MODEL:
    compute_model(mav, trim_state, trim_input)


# initialize the simulation time
sim_time = SIM.start_time
end_time = 100

# main simulation loop
print("Press 'Esc' to exit...")
while sim_time < end_time:

    # -------autopilot commands-------------
    commands.airspeed_command = Va_command.square(sim_time)
    commands.course_command = course_command.square(sim_time)
    commands.altitude_command = altitude_command.square(sim_time)

    # -------autopilot-------------
    estimated_state = mav.true_state  # uses true states in the control
    delta, commanded_state = autopilot.update(commands, estimated_state)
    # delta.aileron = delta_trim.aileron
    # delta.elevator = -.1 #delta_trim.elevator
    delta.rudder = delta_trim.rudder
    # delta.throttle += delta_trim.throttle

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    # ------- animation -------
    if ANIMATION:
        mav_view.update(mav.true_state)  # plot body of MAV
    if PLOTS:
        plot_time = sim_time
        data_view.update(mav.true_state,  # true states
                            None,  # estimated states
                            commanded_state,  # commanded states commanded_state, None
                            delta)  # inputs to aircraft
    if ANIMATION or PLOTS:
        app.processEvents()
    if VIDEO is True:
        video.update(sim_time)
        
    # -------Check to Quit the Loop-------
    if quitter.check_quit():
        break

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if SAVE_PLOT_IMAGE:
    data_view.save_plot_image("ch6_plot")

if VIDEO is True:
    video.close()




