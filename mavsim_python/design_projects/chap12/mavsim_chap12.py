"""
mavsim_python
    - Chapter 12 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        4/3/2019 - BGM
        2/27/2020 - RWB
        3/30/2022 - RWB
"""
import sys
sys.path.append('../..')
import pyqtgraph as pg
import parameters.simulation_parameters as SIM
import parameters.planner_parameters as PLAN
from models.mav_dynamics_sensors import MavDynamics
from models.wind_simulation import WindSimulation
from MAV_control.autopilot import Autopilot
from estimation.observer import Observer
from planning.path_follower import PathFollower
from planning.path_manager import PathManager
from planning.path_planner import PathPlanner
from viewers.mav_world_viewer import MAVWorldViewer
from viewers.data_viewer import DataViewer
from tools.quit_listener import QuitListener
from message_types.msg_world_map import MsgWorldMap
quitter = QuitListener()

VIDEO = False
DATA_PLOTS = False
ANIMATION = True
PLANNING_VIEWER = True

# video initialization
if VIDEO is True:
    from viewers.video_writer import VideoWriter
    video = VideoWriter(video_name="chap12_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)
    
#initialize the visualization
if ANIMATION or DATA_PLOTS:
    app = pg.QtWidgets.QApplication([]) # use the same main process for Qt applications
if ANIMATION:
    world_view = MAVWorldViewer(app=app) # initialize the viewer
if DATA_PLOTS:
    data_view = DataViewer(app=app,dt=SIM.ts_simulation, plot_period=SIM.ts_plot_refresh, 
                           data_recording_period=SIM.ts_plot_record_data, time_window_length=30)

# initialize elements of the architecture
mav = MavDynamics(SIM.ts_simulation)
wind = WindSimulation(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)
observer = Observer(SIM.ts_simulation)
path_follower = PathFollower()
path_manager = PathManager()
planner_flag = 'simple_straight'  # return simple waypoint path
# planner_flag = 'simple_dubins'  # return simple dubins waypoint path
# planner_flag = 'rrt_straight'  # plan path through city using straight-line RRT
# planner_flag = 'rrt_dubins'  # plan path through city using dubins RRT
path_planner = PathPlanner(app=app, planner_flag=planner_flag, show_planner=PLANNING_VIEWER)
world_map = MsgWorldMap()

# initialize the simulation time
sim_time = SIM.start_time
end_time = 200

# main simulation loop
print("Press 'Esc' to exit...")
while sim_time < SIM.end_time:
    # -------observer-------------
    measurements = mav.sensors()  # get sensor measurements
    estimated_state = observer.update(measurements)  # estimate states from measurements
    # Observer occasionally gives bad results, true states always work.
    #estimated_state = mav.true_state
    # -------path planner - ----
    if path_manager.manager_requests_waypoints is True:
        waypoints = path_planner.update(world_map, estimated_state, PLAN.R_min)

    # -------path manager-------------
    path = path_manager.update(waypoints, PLAN.R_min, estimated_state)

    # -------path follower-------------
    autopilot_commands = path_follower.update(path, estimated_state)

    # -------autopilot-------------
    delta, commanded_state = autopilot.update(autopilot_commands, estimated_state)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics

        # -------update viewer-------------
    if ANIMATION:
        world_view.update(mav.true_state, path, waypoints, world_map)  # plot path and MAV
    if DATA_PLOTS:
        plot_time = sim_time
        data_view.update(mav.true_state,  # true states
                         estimated_state,  # estimated states
                         commanded_state,  # commanded states
                         delta)  # inputs to aircraft
    if ANIMATION or DATA_PLOTS or PLANNING_VIEWER:
        app.processEvents()

    # -------increment time-------------
    sim_time += SIM.ts_simulation

    # -------Check to Quit the Loop-------
    if quitter.check_quit():
        break

if VIDEO is True:
    video.close()



