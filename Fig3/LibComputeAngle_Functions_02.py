import numpy as np
import os, sys
import matplotlib.pyplot as plt

def tangent_points_at_distance(p, d):
    """
    Given a point p = (x, y) on a circle centered at the origin,
    return two points that lie on the tangent at distance d from p.

    Function arguments:
    p = np.array([x,y]) : A point  on the circle
    d = float : Distance from p along the tangent line

    Return:
    p1 = np.array([x1,y1]), p2 = np.array([x2,y2]): Two points (x1, y1), (x2, y2) on the tangent at a distance d to p
    sketch : ---*-----*-----*---- <- Tangent
                |     |     |
                p2    p     p1
    """
    x, y = p[0], p[1]

    # Tangent direction is perpendicular to (x, y), i.e. (-y, x)
    tangent_vec = np.array([-y, x])
    tangent_unit = tangent_vec / np.linalg.norm(tangent_vec)

    # Displacement vector of length d
    displacement = d * tangent_unit

    # Two points along the tangent
    p1 = (x + displacement[0], y + displacement[1])
    p2 = (x - displacement[0], y - displacement[1])

    return p1, p2


def single_entry_exit_intersection(p1, p2, R):
    """
    Returns the unique intersection point of a line segment p1->p2
    with a circle centered at the origin of radius R,
    assuming that p1 lies outside and p2 inside the circle
    point (p2 or p1) outside the circle.

    Function arguments:
    p1 = np.array([x1,y1]) : First point of line segment (x1, y1)
    p2 = np.array([x2,y2]) : Second point of line segment (x2, y2)
    R = float : Radius of the circle centered at origin

    Return:
    intersec = np.array([x,y]) : The intersection point (x, y) or None if no intersection
    """
    assert ((np.linalg.norm(p1) > R) and (np.linalg.norm(p2) < R))
    
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]
    dx = x2 - x1
    dy = y2 - y1

    # Solving (x1 + t*dx)^2 + (y1 + t*dy)^2 = R^2
    # i.e. a*t^2 + b*t + c = 0
    a = dx**2 + dy**2
    b = 2 * (x1 * dx + y1 * dy)
    c = x1**2 + y1**2 - R**2

    discriminant = b**2 - 4 * a * c

    sqrt_disc = np.sqrt(discriminant)
    t = (-b - sqrt_disc) / (2 * a)  # first intersection from p1 -> p2

    x = x1 + t * dx
    y = y1 + t * dy
    intersec = np.array([x,y])
    return intersec



def closest_point(point1, point2, reference):
    """
    For two given points (point1, point2), function
    returns the point which is closer to a reference
    point.
    
    Function arguments:
    point1 = np.array([x1,y1]) - a point in 2 dimensions (x1, y1)
    point2 = np.array([x2,y2]) - a point in 2 dimensions (x2, y2)
    reference = np.array([x,y]) - the reference point in 2 dimensions (x, y)
    
    Return:
    point = np.array([x,y]) - the point (point1 or point2) that is closer to reference
    """
    distance1 = np.linalg.norm(point1 - reference)
    distance2 = np.linalg.norm(point2 - reference)
    
    # Return the point with the smaller distance
    return point1 if distance1 < distance2 else point2
    
def furthest_point(point1, point2, reference):
    """
    For two given points (point1, point2), function
    returns the point which is furthest to a reference
    point.
    
    Function arguments:
    point1 = np.array([x1,y1]) - a point in 2 dimensions (x1, y1)
    point2 = np.array([x2,y2]) - a point in 2 dimensions (x2, y2)
    reference = np.array([x,y]) - the reference point in 2 dimensions (x, y)
    
    Return:
    point = np.array([x,y]) - the point (point1 or point2) that is furthest to reference
    """
    distance1 = np.linalg.norm(point1 - reference)
    distance2 = np.linalg.norm(point2 - reference)
    
    # Return the point with the larger distance
    return point2 if distance1 < distance2 else point1
    
def angle_between_vectors(v1, v2):
    """
    Computes angle between two two-dimensional
    vectors. For the geometry of the
    considered problem, the angle lies in the
    interval [0,pi].
    
    Function arguments:
    v1 = np.array([x1,y1]) - vector in 2d
    v2 = np.array([x2,y2]) - vector in 2d
    
    Return:
    angle = (float) - angle between v1 and v2
    """
    
    # compute dot product and magnitudes of vectors
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    
    # Compute cosine of angle using dot product formula
    cos_theta = dot_product / (norm_v1 * norm_v2)
    
    # Ensure cosine value is within valid range [-1, 1] to avoid precision errors
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    
    # Calculate angle in radians
    angle = np.arccos(cos_theta)
    
    return angle

def shortest_vector_1d(pos1, pos2, systemSize):
    """
    Calculate the shortest vector connecting point A = pos1
    to point B = pos2 in 1D with periodic boundaries.
    Periodic boundaries are at the origin and the system
    size L = systemSize.
    
    Function arguments:
    pos1 = np.array([x1,y1]) - point in 2d
    pos2 = np.array([x2,y2]) - point in 2d
    systemSize = (float) - system size
    
    Return:
    shortest_vector = np.array([x,y]) - shortes vector A->B
    """
    
    # Calculate direct vector from pos1 to pos2
    direct_vector = pos2 - pos1
    
    # Get shortest vector under periodic boundary conditions
    shortest_vector = (direct_vector + systemSize / 2) % systemSize - systemSize / 2
    
    return shortest_vector


def GetOutAngleTrajStartsAtBoundary(Trajectory, pos, theta, HeadTailLength, BulkRegionRadius, ds_min):
    """
    Function computes outgoing angle for bulk trajectories
    with a head that lies outside the bulk region.
    
    Function arguments:
    Trajectory = np.array([i0, i1, i2, ...]) - list of position indices belonging to one trajectory
    pos = np.array([[x0,y0], [x1,y1], [x2,y2], ...]) - 2d array of positions
    theta = np.array([theta0, theta1, theta2, ...]) - the polar angle for positions in pos
    HeadTailLength = (int) - number of head (tail) positions before (after) the bulk trajectory
    BulkRegionRadius = (float) - radius defining bulk region boundary (BRB)
    ds_min = (float) - minimum trajectory length traveled after the trejectory entered bulk region
    
    Return:
    Alpha_out - out angle
    ----- remaining variables relevant only for plotting
    OutCrossPos - intersection point of trajectory with BRB
    OutTangent1 - point on tangent to the BRB at OutCrossPos, at a distance d from OutCrossPos
    OutTangent2 - same as OutTangent1, but opposite direction of OutCrossPos
    sketch : -----*----------*----------*------ <- Tangent to the BRB at OutCrossPos
                  |          |          |
         OutTangent1 <- OutCrossPos -> OutTangent2
    OutTangent - either OutTangent1 or OutTangent2
    OutTraj_Bulk - trajectory index of bulk position, ball covered a trajectory of length ~ds_min since OutCrossPos
    """

    bulk_0 = Trajectory[HeadTailLength] # index of 1st bulk trajectory position
    boundary_final = Trajectory[HeadTailLength-1] # index of last boundary position before bulk trajectory starts
    head_0 = Trajectory[0] # index of 1st head trajectory position
    
    # computer intersection point at which ball enters bulk area
    CrossPos = single_entry_exit_intersection(pos[boundary_final], pos[bulk_0], BulkRegionRadius)
    # get two points at distancd d from CrossPos that lie on tangent to BRB
    Tangent1, Tangent2 = tangent_points_at_distance(CrossPos,50)
    
    # find trajectory index for which the trajectory length since CrossPos is about ds_min
    tot_distInBulk = np.linalg.norm(CrossPos-pos[bulk_0])
    qq = 0
    while tot_distInBulk < ds_min:
        tot_distInBulk += np.linalg.norm(pos[bulk_0+qq]-pos[bulk_0+qq+1])
        qq += 1
    bulk_target = Trajectory[HeadTailLength + qq] # index of position for which the angle will be computed
    
    # make sure bulk_target position is still in bulk
    if np.linalg.norm(pos[bulk_target]) < BulkRegionRadius:
        
        # get polar angle of CrossPos, first head position and target position
        theta_cross = np.angle( CrossPos[0] + 1j * CrossPos[1]) + np.pi
        theta_0 = theta[head_0]
        theta_bulk = theta[bulk_target]
        
        # determine order of positions (see sketch below)
        if (shortest_vector_1d(theta_0, theta_cross, 2*np.pi) * shortest_vector_1d(theta_cross, theta_bulk, 2*np.pi)) > 0 :
            # case : head and target on oposite sides of BRB intersection point
            # -----*-------*-------*----
            #      |       |       |
            # theta_0 theta_cross  theta_bulk
            # consequence, choose tangent that is closest to pos[bulk_target]
            Tangent = closest_point(Tangent1, Tangent2, pos[bulk_target])
        else:
            # case : head and target on same side of BRB intersection point
            # -----*-------*-------*----
            #      |       |       |
            # theta_0 theta_bulk   theta_cross
            # consequence, choose tangent that is furthest to pos[bulk_target]
            Tangent = furthest_point(Tangent1, Tangent2, pos[bulk_target])
        TangentVec = Tangent - CrossPos # tangent vector pointing away from head
        bulkTraj = pos[bulk_target] - CrossPos
        Alpha_out = angle_between_vectors(TangentVec, bulkTraj)
        return Alpha_out, CrossPos, Tangent1, Tangent2, Tangent, bulk_target
    else:
        return -1, np.array([0,0]), np.array([0,0]), np.array([0,0]), np.array([0,0]), 0
        
def GetOutAngleTrajStartsInBulk(Trajectory, pos, theta, HeadTailLength, BulkRegionRadius, ds_min):
    """
    Function computes outgoing angle for bulk trajectories
    with a head that lies inside the bulk region.
    
    Function arguments:
    Trajectory = np.array([i0, i1, i2, ...]) - list of position indices belonging to one trajectory
    pos = np.array([[x0,y0], [x1,y1], [x2,y2], ...]) - 2d array of positions
    theta = np.array([theta0, theta1, theta2, ...]) - the polar angle for positions in pos
    HeadTailLength = (int) - number of head (tail) positions before (after) the bulk trajectory
    BulkRegionRadius = (float) - radius defining bulk region boundary (BRB)
    ds_min = (float) - minimum trajectory length traveled after the trejectory entered bulk region
    
    Return:
    Alpha_out - out angle
    ----- remaining variables relevant only for plotting
    OutCrossPos - intersection point of trajectory with BRB
    OutTangent1 - point on tangent to the BRB at OutCrossPos, at a distance d from OutCrossPos
    OutTangent2 - same as OutTangent1, but opposite direction of OutCrossPos
    sketch : -----*----------*----------*------ <- Tangent to the BRB at OutCrossPos
                  |          |          |
         OutTangent1 <- OutCrossPos -> OutTangent2
    OutTangent - either OutTangent1 or OutTangent2
    OutTraj_Bulk - trajectory index of bulk position, ball covered a trajectory of length ~ds_min since OutCrossPos
    """

    # bouncing from the bulk -> boundary -> bulk
    # start of trajectory is at boundary
    # here define the traj_0 as the min to boundary
    bulk_0 = Trajectory[HeadTailLength]  # index of 1st bulk trajectory position
    boundary_final = Trajectory[HeadTailLength-1] # index of last boundary position before bulk trajectory starts
    
    # computer intersection point at which ball enters bulk area
    CrossPos = single_entry_exit_intersection(pos[boundary_final], pos[bulk_0], BulkRegionRadius)
    # get two points at distancd d from CrossPos that lie on tangent to BRB
    Tangent1, Tangent2 = tangent_points_at_distance( CrossPos, 50)
    
    # compute corse grained angle when bulk trajectory is of length >~ ds_min
    tot_distInBulk = np.linalg.norm(CrossPos-pos[bulk_0])
    qq = 0
    while tot_distInBulk < ds_min:
        tot_distInBulk += np.linalg.norm(pos[bulk_0+qq]-pos[bulk_0+qq+1])
        qq += 1
    bulk_target = Trajectory[HeadTailLength + qq] # index of position for which the angle will be computed
    
    # find last bulk position of head
    qq = HeadTailLength-1
    head_finalBulk = Trajectory[qq] # last bulk position of head
    while (np.linalg.norm(pos[head_finalBulk]) > BulkRegionRadius) and (qq > 0):
        qq -= 1
        head_finalBulk = Trajectory[qq]
    head_1stBound = Trajectory[qq+1] # 1st position of head outside BRB, i.e. at boundary
    # computer intersection point at which head of trajectory leaves BRB
    CrossPos_head = single_entry_exit_intersection(pos[head_1stBound], pos[head_finalBulk], BulkRegionRadius)
    
    # check that pos[bulk_target] is still in bulk
    if np.linalg.norm(pos[bulk_target]) < BulkRegionRadius:
        
        # get polar angle of CrossPos, CrossPos_head and target position
        theta_cross = np.mod( np.angle( CrossPos[0] + 1j * CrossPos[1]), 2*np.pi)
        theta_crossHead = np.mod( np.angle( CrossPos_head[0] + 1j * CrossPos_head[1]), 2*np.pi)
        theta_bulk = theta[bulk_target]
        
        # determine order of positions (see sketch neöpw)
        if (shortest_vector_1d(theta_crossHead, theta_cross, 2*np.pi) * shortest_vector_1d(theta_cross, theta_bulk, 2*np.pi)) > 0 :
            # case : CrossPos_head and target on oposite sides of BRB intersection point
            # -----*----------*-------*----
            #      |          |       |
            # theta_bulk theta_cross  theta_crossHead
            # consequence, choose tangent that is closest to pos[bulk_target]
            Tangent = closest_point(Tangent1, Tangent2, pos[bulk_target])
        else:
            # case : CrossPos_head and target on same side of BRB intersection point
            # -----*-------*-----------*----
            #      |       |           |
            # theta_cross  theta_bulk  theta_crossHead
            # consequence, choose tangent that is furthest to pos[bulk_target]
            Tangent = furthest_point(Tangent1, Tangent2, pos[bulk_target])
        TangentVec = Tangent - CrossPos # tangent vector pointing away from head
        bulkTraj = pos[bulk_target] - CrossPos
        Alpha_out = angle_between_vectors(TangentVec, bulkTraj)
        return Alpha_out, CrossPos, Tangent1, Tangent2, Tangent, bulk_target
    else:
        return -1, np.array([0,0]), np.array([0,0]), np.array([0,0]), np.array([0,0]), 0

def PlotSampleTrajectory( ArenaWallRadius, BulkRegionRadius, Trajectory, pos, HeadTailLength, OutTangent, OutCrossPos, OutTraj_Bulk, N_BulkTraj, Reversed, t_i, Alpha_out, Prefix):
    fig = plt.figure(dpi = 150)
    ax = fig.add_subplot(111)
    
    ''' PLOT ARENA AND BULK REGION BOUNDARY
        Note that the arena boundary is not the
        actual wall position, becasue the distance
        between the ball and the wall is at
        least the ball's radius.
    '''
    Arena = plt.Circle((0, 0), ArenaWallRadius, facecolor = 'w', edgecolor = 'k')
    BulkRegion = plt.Circle((0, 0), BulkRegionRadius, facecolor = 'w', edgecolor = 'k', linestyle = ':')
    ax.add_patch( Arena )
    ax.add_patch( BulkRegion )
    
    ''' INDICATE THE START AND END OF
        THE TRAJECTORY
    '''
    if Reversed:
        startIndex = Trajectory[0]
        endIndex = Trajectory[-1]
    else:
        startIndex = Trajectory[-1]
        endIndex = Trajectory[0]
    ax.plot( pos[ startIndex, 0], pos[ startIndex, 1], 'o', color = 'b', markersize = 8)
    ax.plot( pos[ endIndex, 0], pos[ endIndex, 1], 's', color = 'g', markersize = 8)

    ''' DRAW TRAJECTORY:
        - gray line : total trajectory
        - gray dots: positions before and after the bulk trajectory
        - orange dots: positions of the bulk trajectory
    '''
    ax.plot(pos[Trajectory,0], pos[Trajectory,1], '-', color = 'gray')
    
    TrajectoryHead = Trajectory[:HeadTailLength] # position indices before bulk trajectory
    ax.plot(pos[ TrajectoryHead, 0], pos[ TrajectoryHead, 1], '.', color = 'gray')
    
    BulkTrajectory = Trajectory[HeadTailLength:-HeadTailLength] # position indices of bulk trajectory
    ax.plot(pos[ BulkTrajectory, 0], pos[ BulkTrajectory, 1], '.', color = 'orange')
    
    TrajectoryTail = Trajectory[-HeadTailLength:] # position indices after bulk trajectory
    ax.plot(pos[ TrajectoryTail, 0], pos[ TrajectoryTail, 1], '.', color = 'gray')
    
    ''' Plot arrows that define incoming/outgoing angle
        1) tangential to boundary and pointing in
        the direction opposite head positions (outgoing)
        or opposite to tail positions (incoming)
        
        2) following the bulk trajectory leaving the boundary
        (outgoing) or following the history of the bulk
        trajectory before the collision (incoming)
    '''
    ax.annotate( "", OutTangent, OutCrossPos, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='red', lw=1.5))
    ax.annotate( "", pos[OutTraj_Bulk] , OutCrossPos, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='red', lw=1.5))

    if Reversed:
        id = N_BulkTraj - t_i - 1
        ax.text( 0, 0, '$\\Phi_{out} = %.2f \\pi$\nID = %d \ngreen square: Start \nblue circle : End' % (Alpha_out/np.pi, t_i) )
    else:
        id = t_i
        ax.text( 0, 0, '$\\Phi_{in} = %.2f \\pi$\nID = %d \ngreen square: Start \nblue circle : End' % (Alpha_out/np.pi, t_i) )
    ax.set_aspect('equal')
    
    FigName = Prefix + '%d.png' % id
    plt.savefig(FigName)
    plt.close()
    print('open ', FigName)

