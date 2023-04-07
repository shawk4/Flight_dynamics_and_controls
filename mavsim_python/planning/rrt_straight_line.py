# rrt straight line path planner for mavsim_python
#
# mavsim_python
#     - Beard & McLain, PUP, 2012
#     - Last updated:
#         4/3/2019 - Brady Moon
#         4/11/2019 - RWB
#         3/31/2020 - RWB
import numpy as np
from message_types.msg_waypoints import MsgWaypoints
from viewers.planner_viewer import PlannerViewer

class RRTStraightLine:
    def __init__(self, app, show_planner=True):
        self.segment_length = 300 # standard length of path segments
        self.show_planner = show_planner
        if show_planner:
            self.planner_viewer = PlannerViewer(app)

    def update(self, start_pose, end_pose, Va, world_map, radius):
        tree = MsgWaypoints()
        tree.type = 'straight_line'
        # tree.type = 'fillet'

        ###### TODO ######
        # add the start pose to the tree
        tree.add(start_pose, Va, np.inf, np.inf, 0, 0)
        
        # check to see if start_pose connects directly to end_pose
        pose_i = start_pose
        # pose_i = np.array([tree.ned[:,0]]).T
        # pose_i = column(tree.ned, i)
        while (distance(pose_i, end_pose) > radius): # !!! not sure if this should be radius or not 
            self.extend_tree(tree, end_pose, Va, world_map)
        
            pose_i = column(tree.ned, tree.num_waypoints-1)
        waypoints = smooth_path(tree, world_map) # !!! this won't work when I have branches
        # find path with minimum cost to end_node
        # waypoints_not_smooth = find_minimum_path()
        return waypoints

    def extend_tree(self, tree, end_pose, Va, world_map):
        # extend tree by randomly selecting pose and extending tree toward that pose
        
        ###### TODO ######
        # random_pose(world_map, Va) .... hmmmmm

        flag = True
        while flag: # loop until a valid point is found ... or path is extended???? !!!
            # Generate a random step
            rand_x = np.random.uniform(0, 1)
            if np.random.uniform(0, 1) < 0.5:
                sign = 1
            else:
                sign = -1
            rand_y = 1-rand_x
            rand_x *= self.segment_length*sign
            rand_y *= self.segment_length*sign
            prev_point = tree.ned[:,tree.num_waypoints-1]
            rand_step =  prev_point + np.array([rand_x, rand_y, 0])
            # check to see if random point is close to end pose
            if (distance(rand_step, end_pose) <= self.segment_length):
                rand_step = end_pose
            # check for collision
            flag = collision(prev_point, rand_step, world_map)
            # if there is no collision add point to the tree
            if (flag == False):
                tree.add(np.array([rand_step]).T, Va, np.inf, np.inf, 0, 0) # !!! serious weakness this tree does not have any branches
        return flag
        
    def process_app(self):
        self.planner_viewer.process_app()

def smooth_path(waypoints, world_map):

    # construct smooth waypoint path
    smooth_waypoints = MsgWaypoints()
    ##### TODO #####
    # smooth the waypoint path
    smooth_waypoints.add(column(waypoints.ned,0), 0, np.inf, np.inf, 0, 0)  # add the first waypoint
    for i in range(1, len(waypoints)-1):
        # check for collision
        flag = collision(waypoints[i-1], waypoints[i], world_map)
        if (flag == True):
            smooth_waypoints.add(column(waypoints,i-1), 0, np.inf, np.inf, 0, 0)

    return smooth_waypoints
    


    return smooth_waypoints


def find_minimum_path(tree, end_pose):
    # find the lowest cost path to the end node

    ##### TODO #####
    # find nodes that connect to end_node
    connecting_nodes = []
    
    # find minimum cost last node
    idx = 0

    # construct lowest cost path order
    path = []  # last node that connects to end node
    
    # construct waypoint path
    waypoints = MsgWaypoints()
    return waypoints


def random_pose(world_map, pd):
    # generate a random pose

    ##### TODO #####
    pn = 0
    pe = 0
    pose = np.array([[pn], [pe], [pd]])
    return pose


def distance(start_pose, end_pose):
    # compute distance between start and end pose

    ##### TODO #####
    d = np.linalg.norm(end_pose-start_pose)
    
    return d


def collision(start_pose, end_pose, world_map):
    # check to see of path from start_pose to end_pose colliding with map
    
    ###### TODO ######
    collision_flag = True
    path_points = points_along_path(start_pose, end_pose, 30)
    # for point in path_points:
    #     if world_map[int(point[0]), int(point[1])] == 0:
    #         collision_flag = False
    #     else:
    #         collision_flag = True !!! Still need to do this
    collision_flag = False

    return collision_flag


def height_above_ground(world_map, point):
    # find the altitude of point above ground level
    
    ##### TODO #####
    h_agl = 0
    return h_agl

def points_along_path(start_pose, end_pose, N):
    # returns points along path separated by Del
    unit_vector = (end_pose - start_pose)/(np.linalg.norm(end_pose-start_pose))
    points = []
    for i in range(N):
        points.append(start_pose + i*unit_vector)
    return points


def column(A, i):
    # extracts the ith column of A and return column vector
    tmp = A[:, i]
    col = tmp.reshape(A.shape[0], 1)
    return col