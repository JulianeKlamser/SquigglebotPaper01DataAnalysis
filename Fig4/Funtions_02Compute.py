import numpy as np
import matplotlib.pyplot as plt
import pickle
from dataclasses import dataclass

@dataclass
class CollisionData:
    inAngle: float = 0.0
    current: float = 0.0
    type: str = '' # 'wall', or 'bulk', or 'mix'


def angle_between_vectors(v1, v2):
    """
    Returns angle (0 <= angle <= pi) between two
    vectors v1 and v2.
    
    Function arguments:
    v1 = np.array([x1,y1])
    v2 = np.array([x2,y2])

    Return:
    angle
    """
    
    # Calculate the dot product and the magnitudes of the vectors
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)

    # Compute the cosine of the angle using the dot product formula
    cos_theta = dot_product / (norm_v1 * norm_v2)
    
    # Ensure the cosine value is within the valid range [-1, 1] to avoid numerical issues
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    
    # Calculate the angle in radians
    angle = np.arccos(cos_theta)
    
    return angle

def getIncomingDirection(pos, frame_n, ds_min):
    """
    Takes a particle trajectory (pos), 
    a frame number (frame_n),
    and a minimum trajectory length (ds_min) as arguments. 
    Returns vector v = pos[frame_m] - pos[frame_n] and 
    frame_m (frame_m < frame_n) such that
    ds = sum_{i = m+1}^n |dr_i| >= ds_min, with
    dr_i = pos[frame_i - 1] - pos[frame_i].
    Exception: If pos[frame_m] and pos[frame_n] are 
    identical (finite resolution in space),
    then v = pos[frame_m+1] - pos[frame_n].
    
    Function arguments:
    pos = np.array([[x0,y0], [x1,y1], [x2,y2], ...])
    frame_n (int)
    ds_min (double) > 0

    Return:
    vec = np.array([x,y]]
    frame_m (int)
    """
    assert (ds_min > 0)
    
    frame_t2 = frame_n
    frame_t1 = frame_t2 - 1
    ds = 0.0 # length of trajectory
    
    # while trajectory length shorter than ds_min
    # add previous displacement
    # (displacement is between two sucessive frames)
    while (ds < ds_min):
        dr_frame = pos[ frame_t2 ] - pos[ frame_t1 ]
        ds += np.linalg.norm( dr_frame )
        if (ds < ds_min):
            frame_t2 -= 1
            frame_t1 -= 1
    
    frame_m = frame_t1
    vec =  pos[ frame_n ] - pos[ frame_m ]
    if (np.linalg.norm(vec) == 0):
        frame_m += 1
        vec =  pos[ frame_n ] - pos[ frame_m ]
    
    return vec, frame_m
    
    
def GetProbabilityOfAngles( Collision, Bins = 30):
    binWidth = np.pi / Bins
    BinCenters = np.array([(binWidth/2 + i * binWidth) for i in range(Bins)])

    n_wallEvents, n_bulkEvents, n_mixEvents = 0, 0, 0
    HistInAngleBulk = np.zeros(Bins)
    HistInAngleWall = np.zeros(Bins)
    HistInAngleMix = np.zeros(Bins)
    HistInAngle = np.zeros(Bins)

    for event in Collision:
        inA = event.inAngle
        
        if (inA != -1):
            InBin = int(inA/binWidth)
            
            if InBin == Bins:
                InBin -= 1
                        
            HistInAngle[InBin] += 1
            
            if event.type == 'wall':
                n_wallEvents += 1
                HistInAngleWall[ InBin ] += 1
                
            elif event.type == 'bulk':
                n_bulkEvents += 1
                HistInAngleBulk[ InBin ] += 1
                
            elif event.type == 'mix':
                n_mixEvents += 1
                HistInAngleMix[ InBin ] += 1
                
        else:
            print('invalid angle')
            
    HistInAngleWall /= n_wallEvents * binWidth
    HistInAngleBulk /= n_bulkEvents * binWidth
    HistInAngleMix /= n_mixEvents * binWidth
    HistInAngle /= (n_wallEvents + n_bulkEvents + n_mixEvents) * binWidth
    return BinCenters, HistInAngleBulk, HistInAngleWall, HistInAngleMix, HistInAngle

def GetProbabilityOfCurrent( Collision, minC = 0.08, maxC = 0.19, Bins = 27):

    binWidth = (maxC - minC)/Bins
    C_BinCenters = minC + (0.5 + np.arange(Bins)) * binWidth
    
    n_wallEvents, n_bulkEvents, n_mixEvents = 0, 0, 0
    HistBulkCurrent = np.zeros(Bins)
    HistWallCurrent = np.zeros(Bins)
    HistMixCurrent = np.zeros(Bins)

    for event in Collision:
        C = event.current
        inA = event.inAngle
        
        if (inA != -1):
            cBin = int( (C - minC) / binWidth )
                
            if event.type == 'wall' and cBin < Bins:
                HistWallCurrent[ cBin ] += 1
                n_wallEvents += 1
                
            elif event.type == 'bulk' and cBin < Bins:
                HistBulkCurrent[ cBin ] += 1
                n_bulkEvents += 1
                
            elif event.type == 'mix' and cBin < Bins:
                HistMixCurrent[ cBin ] += 1
                n_mixEvents += 1
            
    HistWallCurrent /= (n_wallEvents * binWidth)
    HistBulkCurrent /= (n_bulkEvents * binWidth)
    HistMixCurrent /= (n_mixEvents * binWidth)
    
    return C_BinCenters, HistWallCurrent, HistBulkCurrent, HistMixCurrent

def ShowCollision(FigFile, coll_start_frame, posB, posR, nDist_t, PreCollFrameB, PreCollFrameR, d_coll, ArenaWallRadius, BulkRegionRadius, nD_max_contact, collision_type, deltaFrame = 15):
    """ visualize the collision """

    ''' define figure'''
    fig, axs = plt.subplots(1, 2, dpi = 150, figsize = (10,5))
    ax, ax01 = axs
    ax.title.set_text('---> collision type : %s <---' % collision_type)
    
    ''' define time frame to visualize '''
    MinFrame = coll_start_frame - deltaFrame
    MaxFrame = coll_start_frame + deltaFrame
    
    ''' plot arena and bulk region boundary '''
    Arena = plt.Circle((0, 0), ArenaWallRadius, facecolor = 'none', edgecolor = 'gray')
    BulkRegion = plt.Circle((0, 0), BulkRegionRadius, facecolor = 'none', edgecolor = 'gray', linestyle = ':')
    ax.add_patch( Arena )
    ax.add_patch( BulkRegion )
            
    ''' plot balls at first contact '''
    BallB = plt.Circle(posB[ coll_start_frame ], d_coll/2, facecolor = 'none', edgecolor = 'k')
    BallR = plt.Circle(posR[ coll_start_frame ], d_coll/2, facecolor = 'none', edgecolor = 'r')
    ax.add_patch( BallB )
    ax.add_patch( BallR )
    
    ''' plot ball trajectories '''
    ax.plot(posB[ MinFrame: MaxFrame, 0], posB[ MinFrame: MaxFrame, 1], '-', color = 'k', zorder = 0, linewidth = 1)
    ax.plot(posR[ MinFrame: MaxFrame, 0], posR[ MinFrame: MaxFrame, 1], '-', color = 'r', zorder = 0, linewidth = 1)
    
    colors = np.linspace(0, 1, MaxFrame-MinFrame)
    ax.scatter(posB[ MinFrame: MaxFrame, 0], posB[ MinFrame: MaxFrame, 1], c = colors, cmap = 'jet', s = 2)
    ax.scatter(posR[ MinFrame: MaxFrame, 0], posR[ MinFrame: MaxFrame, 1], c = colors, cmap = 'jet', s = 2)
    
    ''' highlight last position before contact '''
    ax.plot(posB[coll_start_frame,0], posB[coll_start_frame,1], 'k.', label = 'last frame before contact')
    ax.plot(posR[coll_start_frame,0], posR[coll_start_frame,1], 'k.')
    
    ''' first frame '''
    ax.plot(posB[MinFrame,0], posB[MinFrame,1], '.', color = 'grey', label = 'start')
    ax.plot(posR[MinFrame,0], posR[MinFrame,1], '.', color = 'grey')
    
    ''' visualize approximation of velocity direction '''
    ax.plot([posB[coll_start_frame,0], posB[PreCollFrameB,0]], [posB[coll_start_frame,1], posB[PreCollFrameB,1]], '-', color = 'fuchsia', label = 'incoming velocity direction')
    ax.plot([posR[coll_start_frame,0], posR[PreCollFrameR,0]], [posR[coll_start_frame,1], posR[PreCollFrameR,1]], '-', color = 'fuchsia')
    
    drB = posB[PreCollFrameB] - posB[coll_start_frame]
    drBScaled = d_coll/2 * drB / np.linalg.norm(drB)
    EndPoint = posB[coll_start_frame] + drBScaled
    ax.plot([posB[coll_start_frame,0], EndPoint[0]], [posB[coll_start_frame,1], EndPoint[1]], '-', color = 'fuchsia', alpha = 0.3)
    
    drR = posR[PreCollFrameR] - posR[coll_start_frame]
    drRScaled = d_coll/2 * drR / np.linalg.norm(drR)
    EndPoint = posR[coll_start_frame] + drRScaled
    ax.plot([posR[coll_start_frame,0], EndPoint[0]], [posR[coll_start_frame,1], EndPoint[1]], '-', color = 'fuchsia', alpha = 0.3)
    
    ''' set plot range '''
    minYB, minYR = min(posB[ MinFrame: MaxFrame, 1]), min(posR[ MinFrame: MaxFrame, 1])
    minY = min(minYB, minYR)
    maxYB, maxYR = max(posB[ MinFrame: MaxFrame, 1]), max(posR[ MinFrame: MaxFrame, 1])
    maxY = max(maxYB, maxYR)
    ax.set_ylim( minY-d_coll/4, maxY+d_coll/4)
    
    
    minXB, minXR = min(posB[ MinFrame: MaxFrame, 0]), min(posR[ MinFrame: MaxFrame, 0])
    minX = min(minXB, minXR)
    maxXB, maxXR = max(posB[ MinFrame: MaxFrame, 0]), max(posR[ MinFrame: MaxFrame, 0])
    maxX = max(maxXB, maxXR)
    ax.set_xlim( minX-d_coll/4, maxX+d_coll/4)
    
    ax.set_aspect('equal')
    ax.legend()
    
    ''' plot separation distance 
        between both ballas as a
        function of time
    '''
    
    timePlot = np.arange(MinFrame, MaxFrame)
    ax01.plot(timePlot, nDist_t[timePlot], 'o-')
    ax01.plot(coll_start_frame, nDist_t[coll_start_frame], 'o', color = 'k', label = 'last frame before contact')
    ax01.plot([MinFrame, MaxFrame], [nD_max_contact, nD_max_contact])
    ax01.set_xlabel('time')
    ax01.set_ylabel('$\\sqrt{dx^2+dy^2}/(2r_{ball})$')
    ax01.legend()
    
    plt.tight_layout()
    plt.savefig(FigFile)
    print('open ', FigFile)
    plt.close()

