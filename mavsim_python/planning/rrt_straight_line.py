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
        waypoints = tree
        # find path with minimum cost to end_node
        waypoints = find_minimum_path(tree, end_pose)
        waypoints = smooth_path(waypoints, world_map) 
        return waypoints

    def extend_tree(self, tree, end_pose, Va, world_map):
        # extend tree by randomly selecting pose and extending tree toward that pose
        
        ###### TODO ######
        # random_pose(world_map, Va) .... hmmmmm

        flag = True
        while flag: # loop until a valid point is found
            # Generate a random step
            rand_p = random_pose(world_map, end_pose.item(2))
            # loop over all points in the tree to find closest point to rand_p_step
            dist_list = []
            for i in range(tree.num_waypoints):
                dist_list.append(distance(column(tree.ned, i), rand_p))
                idx = np.argmin(dist_list)
            
            prev_point = column(tree.ned, idx)
            rand_p_step = (rand_p-prev_point)/np.linalg.norm(rand_p-prev_point) * self.segment_length
            rand_step =  prev_point + rand_p_step            
            dist = distance(rand_step, prev_point)
            # check to see if random point is close to end pose
            if (distance(rand_step, end_pose) <= self.segment_length/2):
                dist = distance(rand_step, end_pose)
                rand_step = end_pose
                connect_to_goal = 1
            else:
                connect_to_goal = 0

            # check for collision
            flag = collision(prev_point, rand_step, world_map)
            # if there is no collision add point to the tree
            if (flag == False):
                tree.add(rand_step, airspeed=Va, course=np.inf, cost=dist, parent=idx, connect_to_goal=connect_to_goal)
        return flag
        
    def process_app(self):
        self.planner_viewer.process_app()

def smooth_path(waypoints, world_map):
    Va = waypoints.airspeed.item(0)
    # construct smooth waypoint path
    smooth_waypoints = MsgWaypoints()
    ##### TODO #####
    # smooth the waypoint path
    smooth_waypoints.add(column(waypoints.ned,0), airspeed=Va, course=np.inf, cost=0, parent=0, connect_to_goal=0)  # add the first waypoint
    j = 0
    for i in range(1, waypoints.num_waypoints-1):
        # check for collision
        flag = collision(column(waypoints.ned,j), column(waypoints.ned,i), world_map)
        if (flag == True):
            j = i-1 
            i = i-1
            smooth_waypoints.add(column(waypoints.ned,j), airspeed=Va, course=np.inf, cost=0, parent=j, connect_to_goal=0)
    smooth_waypoints.add(column(waypoints.ned, waypoints.num_waypoints-1), airspeed=Va, course=np.inf, cost=0, parent=j, connect_to_goal=1)
    

    return smooth_waypoints


def find_minimum_path(tree, end_pose):
    # find the lowest cost path to the end node
    ##### TODO #####
    # find nodes that connect to end_node
    connecting_nodes = np.argmax(tree.connect_to_goal) # !!! this will break when finding more than one path
    
    # find minimum cost last node
    idx = int(connecting_nodes) # !!! this will break when finding more than one path

    # construct lowest cost path order
    path = [column(tree.ned, idx)] # last node that connects to end node
    # loop untill we find the root 
    while True:
        parent = tree.parent[idx]
        path.append(column(tree.ned, int(parent)))
        idx = int(parent)
        if idx == 0:
            break
    
    # construct waypoint path
    waypoints = MsgWaypoints()
    for i in range(len(path)-1):
        idx = len(path)-1-i
        
        waypoints.add(np.array(path[idx]), airspeed=tree.airspeed.item(idx), course=np.inf, cost=0, parent=idx, connect_to_goal=0)

    waypoints.add(np.array(path[0]), airspeed=tree.airspeed.item(0), course=np.inf, cost=0, parent=idx, connect_to_goal=1)
    return waypoints


def random_pose(world_map, pd):
    # generate a point within the city limits world_map.city_width

    ##### TODO #####
    
    pn = np.random.uniform(0, world_map.city_width)
    pe = np.random.uniform(0, world_map.city_width)
    pose = np.array([[pn], [pe], [pd]])
    return pose


def distance(start_pose, end_pose):
    # compute distance between start and end pose

    ##### TODO #####
    d = np.linalg.norm(end_pose-start_pose)
    
    return d


def collision(start_pose, end_pose, world_map):
    # return False
    # check to see of path from start_pose to end_pose colliding with map
    ###### TODO ######
    collision_flag = True
    path_points = points_along_path(start_pose, end_pose, 30)
    #loop over all points between start and end pose and check for collision
    for point in path_points:
        for ni, n_pose in enumerate(world_map.building_north[0]):
            for ei, e_pose in enumerate(world_map.building_east[0]):
                bld_pose = np.array([[n_pose, e_pose, point.item(2)]]).T
                if (distance(point, bld_pose) < world_map.building_width*1.5):
                    if (abs(point.item(2)) < abs(world_map.building_height[ni,ei])):
                        return True
    return False



def height_above_ground(world_map, point):
    # find the altitude of point above ground level
    
    ##### TODO #####
    h_agl = 0
    return h_agl

def points_along_path(start_pose, end_pose, N):
    # returns points along path separated by Del
    unit_vector = (end_pose - start_pose)/(np.linalg.norm(end_pose-start_pose))
    dist = distance(start_pose, end_pose)
    res = dist/N
    points = []
    for i in range(N):
        points.append(start_pose + i*res*unit_vector)
    points.append(end_pose)

    return points


def column(A, i):
    # extracts the ith column of A and return column vector
    tmp = A[:, i]
    col = tmp.reshape(A.shape[0], 1)
    return col