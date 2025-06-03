import numpy as np
from matplotlib.collections import EllipseCollection
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

print('\n This may take a while. \n')
    
def AnimateFct( frame ):
    """
    The MD loop including update of frame for animation.
    """
    global WallRadius, gamma, rho
    ID, t_physical = frame
    
    ''' read configuration of frame '''
    Path_ToConf = SimuPath + '/Particles%d.txt' % ID
    Data = np.loadtxt(Path_ToConf)
    pos = Data[:,:2]
    NP = len(pos)
    
    ''' shift positions such that center of confining
        wall is at the origin '''
    pos -= WallRadius

    ''' update positions '''
    collection.set_offsets( pos )
    collection2.set_offsets( pos[0] )
    
    ''' update title '''
    plt.title( '$\\phi = %.3f$, $\\gamma = %.2f$, $N = %d$\n$t = %.2f$, ID = %06d' % (rho, gamma, NP, t_physical, ID) )
    return collection, collection2


BasePath = script_dir + '/Data/ConfinedCircular_N000256L28.000000/N000256L28.000000/v_o1.000000D_r1.000000'
PathList = [BasePath + '/gamma1.000000alignCutOffBySigma2.000000/00MovieData', BasePath + '/gamma0.100000alignCutOffBySigma2.000000/00MovieData', BasePath + '/gamma0.010000alignCutOffBySigma2.000000/00MovieData']
FPS_gamma = np.array([[10, 0.01], [20, 0.1], [20, 1.0]]) # fps, gamma


for SimuPath in PathList:

    ''' read gamma from the simulation parameter file '''
    Path_ToParams = SimuPath + '/Parameters.txt'
    Data = np.loadtxt(Path_ToParams)
    gamma = Data[-1]

    ''' read size of simulaiton box '''
    Path_ToPeriods = SimuPath + '/Periods.txt'
    Periods = np.loadtxt(Path_ToPeriods)
    assert( np.all(Periods[0, :] == Periods[1, :]) )
    assert( np.all(Periods[:, 0] == 0) )
    WallRadius = (Periods[0,1] - 5.) / 2.

    ''' read time file of the simulation '''
    Path_ToTime = SimuPath + '/Time.txt'
    Time = np.loadtxt(Path_ToTime, usecols=(0, 1), skiprows = 1)
    # simulation was running with large time differences between
    # individual frames. find first ID in movie directory from
    # which frames were generated with small time differences
    # between individual frames.
    i, NumberOfFrames = int(Time[0,0]), len(Time)
    while i < NumberOfFrames and not os.path.isfile(SimuPath + '/Particles%d.txt' % Time[i,0]):
        i += 1
    Time = Time[i:]
    Time[:,1] -= Time[0,1]
    start_ID, t_physical = Time[0]
    NumberOfFrames = len(Time)

    ''' read initial configuration '''
    Path_ToFinalConf = SimuPath + '/Particles%d.txt' % Time[0,0]
    Data = np.loadtxt(Path_ToFinalConf)
    NP = len(Data) # number of particles
    pos_init = Data[:,:2] - WallRadius # shift center of wall to origin
    radii = Data[:,5] # sigma_i / 2 = particle 'radii'
    assert( np.all(radii == radii[0])) # all particles have the same sigma
    sigmaP_Half = radii[0]

    ''' compute density'''
    rho = NP * sigmaP_Half**2 / WallRadius**2

    ''' define figure '''
    fig, ax = plt.subplots(dpi = 25, figsize=(4, 4))
    ax.set_xlim(-WallRadius, +WallRadius)
    ax.set_ylim(-WallRadius, +WallRadius)
    ax.set_aspect('equal')
    plt.title( '$\\phi = %.3f$, $\\gamma = %.2f$, $N = %d$\n$t = %.2f$, ID = %06d' % (rho, gamma, NP, t_physical, start_ID) )
    plt.tight_layout()

    ''' plot arena wall '''
    circle = plt.Circle((0, 0), WallRadius, color='k', linestyle=':', fill=False, linewidth=1)
    ax.add_patch(circle)

    ''' plot particle positions '''
    Width, Hight, Angle = 2*sigmaP_Half, 2*sigmaP_Half, 0
    collection = EllipseCollection( Width, Hight, Angle, units='x', offsets=pos_init, transOffset=ax.transData, edgecolor = 'k', facecolor = 'none')
    ax.add_collection(collection)

    collection2 = EllipseCollection( Width, Hight, Angle, units='x', offsets=pos_init[0], transOffset=ax.transData, edgecolor = 'k', facecolor = 'k')
    ax.add_collection(collection2)

    ''' create animation '''
    if not gamma in FPS_gamma[:,1]:
        Frames_Per_s = 20
    else:
        k = np.where(FPS_gamma[:,1] == gamma)[0][0]
        Frames_Per_s = FPS_gamma[k,0]
    MoviePath = SimuPath + '/00Movie'
    os.makedirs(MoviePath, exist_ok=True)
    MovieFile = MoviePath + '/Gamma%.2f.mp4' % gamma
    ani = FuncAnimation(fig, AnimateFct, frames=Time, blit=False, repeat=False)
    ani.save( MovieFile, fps=Frames_Per_s, extra_args=['-vcodec', 'libx264'])
    print('open', MovieFile, '\n')
