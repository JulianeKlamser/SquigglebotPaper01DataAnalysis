import numpy as np
import os, sys
import pickle
script_dir = os.path.dirname(os.path.abspath(__file__))

from Funtions_02Compute import angle_between_vectors
from Funtions_02Compute import getIncomingDirection
from Funtions_02Compute import ShowCollision
from Funtions_02Compute import CollisionData
from Funtions_02Compute import GetProbabilityOfAngles
from Funtions_02Compute import GetProbabilityOfCurrent


HowToUse = "\nto compute histogram incoming collision angle launch script as\npython3 02_ComputeAngle.py F\n\nto compute outgoing histogram incoming collision angle of time reversed trajectory launch script as\npython3 02_ComputeAngle.py R\n\n"
if len(sys.argv) != 2:
    print(HowToUse)
    exit()
if sys.argv[1] == 'F':
    Reversed = False
elif sys.argv[1] == 'R':
    Reversed = True
else:
    print(HowToUse)
    exit()
    
''' definition of constants
    d_coll - inter-particle distance at collision
    ds_min - minimum trajectory length to estimate velocity direction
    bulkDef - percentage of the total arena radius defining the bulk region
'''
d_coll = 87.5118
ds_min = 0.1 * (d_coll/2)
bulkDef = 0.90

''' read trajectory data and current data
    for both balls
'''
Base = script_dir + '/Data/'
FileName = Base + 'Ball1_TrajectoryAndCurrent.txt'
Data = np.loadtxt(FileName)
posB = Data[:,:2] # x y of black ball
CurrentB = Data[:,2] # current of black ball
n_framesB = len(posB) # number of frames of total black trajectory

FileName = Base + 'Ball2_TrajectoryAndCurrent.txt'
Data = np.loadtxt(FileName)
posR = Data[:,:2] # x y of red ball
CurrentR = Data[:,2] # current of red ball
n_framesR = len(posR) # number of frames of total red trajectory

assert( n_framesB == n_framesR )
n_frames = n_framesB # number of frames of the experiment

''' reverse trajecotires 
'''
if Reversed:
    posB = posB[::-1]
    posR = posR[::-1]
    CurrentB = CurrentB[::-1]
    CurrentR = CurrentR[::-1]
    FilePath = script_dir + '/Data/Out_Backward'
    Pref = script_dir + '/CollBackward'
else:
    FilePath = script_dir + '/Data/Out_Forward'
    Pref = script_dir + '/CollForward'

''' move center of arena to origin
'''
PosTot = np.concatenate((posB, posR), axis=0)

centerOfMass = np.mean(PosTot, axis = 0)
PosTot -= centerOfMass
posR -= centerOfMass
posB -= centerOfMass

''' segment arena in azimuthal bins.
    determin the average maximum distance of ball positions
    per azimuthal bins.
'''
thetaTot, rTot = np.mod(np.angle(PosTot[:,0] + 1j * PosTot[:,1]), 2*np.pi), np.linalg.norm(PosTot, axis = 1)
bins = 200
max_r = np.zeros(bins)
Theta_binBoundaries = np.linspace(0, 2*np.pi, bins+1)
Theta_BinWidth = (Theta_binBoundaries[1]-Theta_binBoundaries[0])
for i, dummy in enumerate(max_r):
    # mask = particle indices with positions inside the azimuthal bin i
    mask = ( (thetaTot > Theta_binBoundaries[i]) & (thetaTot <= Theta_binBoundaries[i+1]))
    sorted_r = np.sort( rTot[mask] ) # sort distance of particle positions to origin
    max_r[i] = np.mean( sorted_r[-50:] ) # average over largest 50 distances to origin
ArenaWallRadius = np.mean(max_r)
BulkRegionRadius = bulkDef * ArenaWallRadius


''' Identify collision events.
    d_i = |r_black(t_i) - r_red(t_i)| - distance between both particles at time t_i corresponding to frame i.
    d_coll = inter-particle distance at collision
    
    The frame i where particles enter into contact satisfies
    d_i < d_coll, d_{i-1} > d_coll, d_{i-2} > d_coll, d_{i-3} > d_coll
    
    Therefore frame i-1 is the last frame before particles enter into contact.
    
    CollStartFrames : array of all frames (i - 1)  
'''

distance_t = np.linalg.norm( posB - posR, axis = 1) # distance between particles per frame
nDist_t = distance_t / d_coll # distance normalized by collision distance
nD_max_contact = 1.0 # normalized contact distance

frame_indices = np.arange(n_frames)

# condition : (nDist_t < nD_max_contact) and the previous 3 (nDist_t > nD_max_contact)
Mask = (nDist_t < nD_max_contact) & (np.roll(nDist_t,1) > nD_max_contact) & (np.roll(nDist_t,2) > nD_max_contact) & (np.roll(nDist_t,3) > nD_max_contact)
CollStartFrames = frame_indices[Mask] - 1
n_collisions = len(CollStartFrames)


''' Choose some collision IDs to
    visualise the bulk trajectories
    with the correponing vectors
    between which angle is computed.
'''
ShowTrajectory = np.array([10, 20, 40, 80]) # IDs of collsions to visualize
if Reversed:
    ShowTrajectory = n_collisions - 1 - ShowTrajectory

''' Main loop over all collisions
'''
Collision = [CollisionData() for _ in range(n_collisions)]
for i, coll_start_frame in enumerate(CollStartFrames):
    
    ''' identify collision type
    '''
    rR, rB = np.linalg.norm(posR[coll_start_frame]), np.linalg.norm(posB[coll_start_frame])
    if (rR < BulkRegionRadius) and (rB < BulkRegionRadius):
        # both particles are in the bulk region at the beginning of collision
        Collision[i].type = 'bulk'
    elif (rR < BulkRegionRadius) or (rB < BulkRegionRadius):
        # only one particle is in the bulk region at the beginning of collision
        Collision[i].type = 'mix'
    else:
        # both particles are at the wall at the beginning of collision
        Collision[i].type = 'wall'
    
    ''' get approximate direction of trajectories
        before balls enter into contact.
        v_b - incoming direction of black ball
        v_r - incoming direction of red ball
        v_b = posB[coll_start_frame] - posB[PreCollFrameB] 
        v_r = posR[coll_start_frame] - posR[PreCollFrameR]
        Trajectory length between PreCollFrameB (PreCollFrameR) and coll_start_frame is at least ds_min.
    '''
    v_b, PreCollFrameB = getIncomingDirection(posB, coll_start_frame, ds_min)
    v_r, PreCollFrameR = getIncomingDirection(posR, coll_start_frame, ds_min)
    
    ''' Finite resolution of position space can cause 
        |v_r| = 0 or |v_b| = 0
        In that case, the angle is not evaluated.
    '''
    if (np.linalg.norm(v_r) > 0) and (np.linalg.norm(v_b) > 0):
        Collision[i].inAngle = angle_between_vectors(v_b, v_r)
    else:
        # collision could not be evaluated
        print('invalid angle')
        Collision[i].inAngle = -1 # unvalid angle
    
    ''' get current.
        C = 1/9 sum_{i = first_contact_frame - 4}^{first_contact_frame + 4} 1/2 (c_red(t_i) + c_black(t_i)) 
    '''
    timeLapse = 5
    first_contact_frame = coll_start_frame + 1
    Current = 0.5 * ( np.mean(CurrentB[ first_contact_frame - (timeLapse-1) : first_contact_frame + timeLapse ]) + np.mean(CurrentR[ first_contact_frame - (timeLapse-1) : first_contact_frame + timeLapse ]) )
    Collision[i].current = Current

    
    if i in ShowTrajectory:
        if Reversed:
            id = n_collisions - i - 1
        else:
            id = i
        FigFile = Pref + '%d.png' % id
        ShowCollision(FigFile, coll_start_frame, posB, posR, nDist_t, PreCollFrameB, PreCollFrameR, d_coll, ArenaWallRadius, BulkRegionRadius, nD_max_contact, Collision[i].type, deltaFrame = 15)

Angle, P_of_Angle_bulk, P_of_Angle_wall, P_of_Angle_mix, P_of_Angle = GetProbabilityOfAngles( Collision, Bins = 20)

current, P_of_current_wall, P_of_current_bulk, P_of_current_mix = GetProbabilityOfCurrent( Collision, minC = 0.08, maxC = 0.15, Bins = 30)


Dictionary = {'Angle' : Angle, 'P_of_Angle_bulk' : P_of_Angle_bulk, 'P_of_Angle_wall' : P_of_Angle_wall, 'P_of_Angle_mix' : P_of_Angle_mix, 'P_of_Angle' : P_of_Angle, 'current' : current, 'P_of_current_wall' : P_of_current_wall, 'P_of_current_bulk' : P_of_current_bulk, 'P_of_current_mix' : P_of_current_mix}
FilePath += '%.2f.pkl' % bulkDef
with open(FilePath, 'wb') as f:
    pickle.dump(Dictionary, f)
    
print('dump ', FilePath)

