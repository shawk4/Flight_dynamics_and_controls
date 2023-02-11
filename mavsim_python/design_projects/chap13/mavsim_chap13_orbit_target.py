"""
mavsim_python
    - Chapter 13 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        4/1/2022 - RWB
"""
import sys
sys.path.append('../..')
import pyqtgraph as pg
import parameters.simulation_parameters as SIM
import parameters.planner_parameters as PLAN
from models.wind_simulation import WindSimulation
from models.camera import Camera
from models.target_dynamics import TargetDynamics
from models.mav_dynamics_camera import MavDynamics
from models.gimbal import Gimbal
from control.autopilot import Autopilot
from estimation.observer import Observer
from estimation.geolocation import Geolocation
from viewers.geolocation_viewer import GeolocationViewer
from planning.path_follower import PathFollower
from planning.path_manager_follow_target import PathManager
from viewers.data_viewer import DataViewer
from viewers.mav_world_camera_viewer import MAVWorldCameraViewer
from viewers.camera_viewer import CameraViewer
from message_types.msg_world_map import MsgWorldMap
from message_types.msg_waypoints import MsgWaypoints
from tools.quit_listener import QuitListener


quitter = QuitListener()

VIDEO = False
DATA_PLOTS = False
ANIMATION = True
GEO_PLOTS = True
CAMERA_VIEW = True

# video initialization
if VIDEO is True:
    from viewers.video_writer import VideoWriter
    video = VideoWriter(video_name="chap13_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize the visualization
if ANIMATION or DATA_PLOTS or GEO_PLOTS:
    app = pg.QtWidgets.QApplication([]) # use the same main process for Qt applications
if ANIMATION:
    world_view = MAVWorldCameraViewer(app=app)  # initialize the viewer
if DATA_PLOTS:
    data_view = DataViewer(app=app,dt=SIM.ts_simulation, plot_period=SIM.ts_plot_refresh, 
                           data_recording_period=SIM.ts_plot_record_data, time_window_length=30)
if GEO_PLOTS:
    geo_viewer = GeolocationViewer(app=app,dt=SIM.ts_simulation, plot_period=SIM.ts_plot_refresh, 
                           data_recording_period=SIM.ts_plot_record_data, time_window_length=30)
if CAMERA_VIEW:
    camera_view = CameraViewer()

# initialize elements of the architecture
world_map = MsgWorldMap()
mav = MavDynamics(SIM.ts_simulation)
gimbal = Gimbal()
camera = Camera()
target = TargetDynamics(SIM.ts_simulation, world_map)
wind = WindSimulation(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)
observer = Observer(SIM.ts_simulation)
path_follower = PathFollower()
path_manager = PathManager()
waypoints = MsgWaypoints()
geolocation = Geolocation(SIM.ts_simulation)

# initialize the simulation time
sim_time = SIM.start_time
end_time = 200

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < SIM.end_time:
    # -------observer-------------
    measurements = mav.sensors()  # get sensor measurements
    camera.updateProjectedPoints(mav.true_state, target.position())
    pixels = camera.getPixels()
    #estimated_state = observer.update(measurements)  # estimate states from measurements
        # I occasionally get bad results with the observer.  true states seem to always work.
    estimated_state = mav.true_state

    # -------camera control-------------
    #gimbal_cmd = gimbal.pointAtGround(estimated_state)  # point gimbal at ground
    estimated_target_position = geolocation.update(estimated_state, pixels)
    #gimbal_cmd = gimbal.pointAtPosition(estimated_state, estimated_target_position) # point gimbal at target position
    gimbal_cmd = gimbal.pointAtPosition(estimated_state, target.position())  # point gimbal at target position

    # -------path manager-------------
    path = path_manager.update(target.position())
    #path = path_manager.update(estimated_target_position)

    # -------path follower-------------
    autopilot_commands = path_follower.update(path, estimated_state)

    # -------autopilot-------------
    delta, commanded_state = autopilot.update(autopilot_commands, estimated_state)
    delta.gimbal_az = gimbal_cmd.item(0)
    delta.gimbal_el = gimbal_cmd.item(1)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics
    target.update()  # propagate the target dynamics

    # -------update viewer-------------
    if ANIMATION:
        world_view.update(mav.true_state,
                          target.position(),
                          path, waypoints, world_map)  # plot world
    if DATA_PLOTS:
        plot_time = sim_time
        data_view.update(mav.true_state,  # true states
                         estimated_state,  # estimated states
                         commanded_state,  # commanded states
                         delta)  # inputs to aircraft
    if GEO_PLOTS:
        geo_viewer.update(estimated_target_position - target.position())
    if CAMERA_VIEW:
        camera_view.updateDisplay(camera.getProjectedPoints())
    if ANIMATION or DATA_PLOTS or GEO_PLOTS:
        app.processEvents()
    
    # -------Check to Quit the Loop-------
    if quitter.check_quit():
        break

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if VIDEO is True:
    video.close()




