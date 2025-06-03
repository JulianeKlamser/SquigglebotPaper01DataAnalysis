import numpy as np
import os, sys
import matplotlib.pyplot as plt
script_dir = os.path.dirname(os.path.abspath(__file__))

from LibComputeAngle_Functions_02 import tangent_points_at_distance
from LibComputeAngle_Functions_02 import single_entry_exit_intersection
from LibComputeAngle_Functions_02 import closest_point
from LibComputeAngle_Functions_02 import furthest_point
from LibComputeAngle_Functions_02 import angle_between_vectors
from LibComputeAngle_Functions_02 import shortest_vector_1d
from LibComputeAngle_Functions_02 import GetOutAngleTrajStartsAtBoundary
from LibComputeAngle_Functions_02 import GetOutAngleTrajStartsInBulk
from LibComputeAngle_Functions_02 import PlotSampleTrajectory


''' script computes outgoing angle only !
    To obtain in the incoming angle, the
    same script is used, but the particle
    trajectory is reversed.
'''
HowToUse = "\nto compute incoming angle launch script as\npython3 02_ComputeAngle.py in\n\nto compute outgoing angle launch script as\npython3 02_ComputeAngle.py out\n\n"
if len(sys.argv) != 2:
    print(HowToUse)
    exit()
if sys.argv[1] == 'in':
    Reversed = True
elif sys.argv[1] == 'out':
    Reversed = False
else:
    print(HowToUse)
    exit()
Prefix = 'TimeForward'

''' Visualise the bulk trajectories
    with the correponing vectors
    between which angle is computed.
'''
ShowTrajectory = np.array([3, 4, 5, 40, 111]) # list the trajectories you want to see
# some trajectories are unvalid because they are too short
ShowTrajOut = np.array([201,800])

''' read data '''
DataPath = script_dir + '/Data/ParticleTrajectoryDiscArena.txt'
pos = np.loadtxt(DataPath)
N_pos = len(pos)

''' shift center of mass to zero '''
centerOfMass = np.mean(pos, axis = 0)
pos -= centerOfMass
FileName = script_dir + '/Data/F3a_TrajectoryInDisk.txt'
np.savetxt(FileName, pos, header=' x and y position of ball in disc arena centered at origin' )
print('save ', FileName)
if Reversed:
    pos = pos[::-1]
    Prefix = 'TimeReversed'

''' define arena wall radius from positions.
    theta, r are trajectory positions in polar coordinates
'''
theta, r = np.mod(np.angle(pos[:,0] + 1j * pos[:,1]), 2*np.pi), np.linalg.norm(pos, axis = 1)
NAzimuthalBins = 200
max_r = np.zeros(NAzimuthalBins)
Azimuth_binBoundaries = np.linspace(0, 2*np.pi, NAzimuthalBins+1)
for i, dummy in enumerate(max_r):
    mask = (theta > Azimuth_binBoundaries[i]) & (theta <= Azimuth_binBoundaries[i+1])
    sorted_r = np.sort(r[mask])
    max_r[i] = np.mean(sorted_r[-10:])
ArenaWallRadius = np.mean(max_r)

'''
 define bulk and boundary area
 bulk : positions with radii < 0.97 * ArenaWallRadius
 boundary : positions with radii >= 0.97 * ArenaWallRadius
'''
BulkRegionRadius = 0.97*ArenaWallRadius


'''
 collision angle is computed from a vector
 pointing to the position that the ball
 has reached after traveling a minimum
 trajectory length of ds_min from the
 intersection point of the bulk region boundary
'''
RadiusOfBall = 148/2
ds_min = 0.5898 * RadiusOfBall

'''
 find bulk trajectories
 BulkTrajectories = lists of indeces of positions that correspond to one bulk trajectory
 BulkTrajectories = [[ i, i+1, i+2, i+3, ...], [ j, j+1, j+2, j+3, ...], ....]
 Each bulk trajectory is composed of :
 (10 positions preceding the bulk trajectory) +
 (the bulk trajectory) +
 (10 positions following the bulk trajectory)
 We call (10 positions following the bulk trajectory) tail positions
 We call (10 positions preceding the bulk trajectory) head positions
 We call (the bulk trajectory) bulk positions
'''
BulkTrajectories = []
i = 0
HeadTailLength = 10 # also minimum length of trajectory
# init
while r[i] < BulkRegionRadius:
    # particle is in bulk
    i += 1
# get bulk trajectories
while i < N_pos:
    Trajectory = []
    while (i < N_pos) and (r[i] < BulkRegionRadius):
        # particle is in bulk
        Trajectory += [i]
        i += 1
    if (len(Trajectory) > HeadTailLength) and (i < (N_pos - HeadTailLength - 1 )):
        Head = [ Trajectory[0] - HeadTailLength + j for j in range(HeadTailLength)]
        Tail = [ Trajectory[-1] + 1 + j for j in range(HeadTailLength)]
        Trajectory = Head + Trajectory + Tail
        BulkTrajectories += [Trajectory]
    i += 1
N_BulkTraj = len(BulkTrajectories)
if Reversed:
    ShowTrajectory = N_BulkTraj - 1 - ShowTrajectory
    ShowTrajOut = N_BulkTraj - 1 - ShowTrajOut
    
''' Get outgoing angles of bulk trajectories
    by looping over each bulk trajectory.
    Initially outgoing angles are initialiued
    to -1 which is used later to identify valid
    angles.
'''
OutAngleList = np.ones(len(BulkTrajectories))*(-1)
for t_i, Trajectory in enumerate(BulkTrajectories):
    ''' distinguish two cases:
        1) head is outside the bulk area
        2) part of the head is in bulk
    '''
    traj_index_Start = Trajectory[:HeadTailLength]
    if (np.all(r[traj_index_Start] > BulkRegionRadius)):
        ''' case 1
            Alpha_out - out angle
            ----- remaining variables relevant only for plotting
            OutCrossPos - intersection point of trajectory with bulk region boundary (brb)
            OutTangent1 - point on tangent to the brb at OutCrossPos, at a distance d from OutCrossPos
            OutTangent2 - same as OutTangent1, but opposite direction of OutCrossPos
            sketch : -----*----------*----------*------ <- Tangent to the brb at OutCrossPos
                          |          |          |
                 OutTangent1 <- OutCrossPos -> OutTangent2
            OutTangent - either OutTangent1 or OutTangent2
            OutTraj_Bulk - trajectory index of bulk position, ball covered a trajectory of length ~ds_min since OutCrossPos (angle will be computed for this position)
        '''
        Alpha_out, OutCrossPos, OutTangent1, OutTangent2, OutTangent, OutTraj_Bulk = GetOutAngleTrajStartsAtBoundary(Trajectory, pos, theta, HeadTailLength, BulkRegionRadius, ds_min)
        if Alpha_out >= 0:
            # Alpha_out = -1 if bulk trajectory does not exceed minimum length
            OutAngleList[t_i] = Alpha_out
            
    else:
        ''' case 2
            see case 1 for variable definition
        '''
        Alpha_out, OutCrossPos, OutTangent1, OutTangent2, OutTangent, OutTraj_Bulk = GetOutAngleTrajStartsInBulk(Trajectory, pos, theta, HeadTailLength, BulkRegionRadius, ds_min)
        if Alpha_out >= 0:
            # Alpha_out = -1 if bulk trajectory does not exceed minimum length
            OutAngleList[t_i] = Alpha_out
            
    ''' plotting sample trajectories '''
    if t_i in ShowTrajectory:
        PlotSampleTrajectory( ArenaWallRadius, BulkRegionRadius, Trajectory, pos, HeadTailLength, OutTangent, OutCrossPos, OutTraj_Bulk, N_BulkTraj, Reversed, t_i, Alpha_out, Prefix)
    
    
    if t_i in ShowTrajOut:
        if Reversed:
            ID = N_BulkTraj - 1 - t_i
            FileName = script_dir + '/Data/F3a_InVectors_ID%d.txt' % ID
            Data = np.vstack((OutCrossPos, OutTangent, pos[OutTraj_Bulk]))
            np.savetxt(FileName, Data, header='InCrossPos, InTangent, InBulk' )
            print('save ', FileName)
            
        else:
            FileName = script_dir + '/Data/F3a_Trajectory_ID%d.txt' % t_i
            np.savetxt(FileName, Trajectory, fmt = '%d', header='sample trajectory indices' )
            print('save ', FileName)
            
            FileName = script_dir + '/Data/F3a_OutVectors_ID%d.txt' % t_i
            Data = np.vstack((OutCrossPos, OutTangent, pos[OutTraj_Bulk]))
            np.savetxt(FileName, Data, header='OutCrossPos, OutTangent, OutTraj_Bulk' )
            print('save ', FileName)


FileName = script_dir + '/Data/F3a_Parameters.txt'
with open(FileName, "w") as f:
    head = "# HeadTailLength, BulkRegionRadius, ArenaWallRadius\n"
    f.write(head)
    line = "%d %.5f %.5f" % (HeadTailLength, BulkRegionRadius, ArenaWallRadius)
    f.write(line)
    
''' compute probability of angles
    and plot
'''
# remove invalid trajectories
CleanOutAngle = OutAngleList[OutAngleList>0]

# generate histogram of angles
Bins = 40
BinWidth = np.pi/Bins
BinCenters = BinWidth/2 + np.arange(Bins) * BinWidth
HistOut = np.zeros(Bins)
for inO in CleanOutAngle:
    index = int(inO/BinWidth)
    HistOut[index] += 1
HistOut /= (sum(HistOut) * BinWidth)
# assert histogram is normalized
assert( abs((sum(HistOut) * BinWidth) - 1.) < 10**(-5) )


plt.ylim(0,1.8)
plt.plot(BinCenters, HistOut, 'r-')
if Reversed:
    FigName = script_dir + '/ProbOfIn.png'
    DataFile = script_dir + '/Data/F3b_ProbOfInAngle.txt'
else:
    FigName = script_dir + '/ProbOfOut.png'
    DataFile = script_dir + '/Data/F3b_ProbOfOutAngle.txt'
plt.savefig(FigName)
print('open ', FigName)
print('save ', DataFile)

Data = np.vstack((BinCenters,HistOut)).T
np.savetxt( DataFile, Data, header='angle P(angle)' )

