import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.pyplot import cm
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

def circular_text_polar(text, ax, radius=1.0, start_angle=0, end_angle=np.pi, clockwise=False, color = 'red'):
    """
    Writes a string of characters along an arc in polar coordinates.
    """
    
    n = len(text)
    # Determine direction and angle span
    angle_range = end_angle - start_angle
    if clockwise:
        angle_range = -angle_range
    
    if n > 1:
        angle_step = angle_range / (n - 1)
    else:
        angle_step = 0  # single character

    for i, char in enumerate(text):
        theta = start_angle + i * angle_step  # theta in radians
        # Convert polar to cartesian for text rotation
        rot = np.degrees(theta) + (90 if clockwise else -90)
        
        ax.text(theta, radius, char,
                rotation=rot,
                rotation_mode='anchor',
                ha='center', va='bottom',
                transform=ax.transData, color = color)

def computePolarHistogram(radius, azimut, N_radialBins):
    '''
    radius & azimut define particle positions (circular arena centered at origin)
    N_radialBins : number of radial bins
    '''
    N = len(radius)
    maxR = np.max(radius)
    
    # number of radial and azimutal bins
    N_rbins = N_radialBins # radial number of bins
    N_abins_per_rbin = np.array([4] + [4*n_r for n_r in range(1,N_rbins)]) # per radial bin, there is a growing number of azimutal bins

    # radial and azimutal bin with
    rBinWidth = maxR / N_rbins # constant radial bin bidth
    aBinWidth_perRbin = [ 2*np.pi/abins for abins in N_abins_per_rbin] # azimutal bin width depends on radial bin

    # initialize histogram to zeros
    Histogram = [[0] * N_abins_per_rbin[n_r] for n_r in range(0,N_rbins)] # dimension = N_rbins x N_abins_per_rbin

    # fill the histogram
    for i, r in enumerate(radius):
        n_r = int( radius[i] / rBinWidth )
        if n_r == N_radialBins:
            n_r -= 1

        aBinWidth = aBinWidth_perRbin[n_r]
        n_a = int( azimut[i] / aBinWidth )
        
        Histogram[n_r][n_a] += 1
        
    # area of bins depends on radial bin
    # compute area ob bins per radial bin
    maxR_of_rBin = np.linspace(0, maxR, N_rbins+1)
    RingAreas = np.pi * np.array([maxR_of_rBin[n_r]**2 - maxR_of_rBin[n_r-1]**2 for n_r in range(1,N_rbins+1)])
    binArea_perRadialBin = RingAreas / N_abins_per_rbin

    # normalization
    totalProbability = 0
    for n_r in range(N_rbins):
        N_abins = N_abins_per_rbin[n_r]
        binArea = binArea_perRadialBin[n_r]
        for n_a in range(N_abins):
            Histogram[n_r][n_a] /= ( N * binArea )
            totalProbability += ( Histogram[n_r][n_a] * binArea )
    assert ( abs(1.-totalProbability) < 10**(-4) )

    # For graphical representation use pcolormesh
    # pcolormesh works with a constant number of abins
    PseudoN_abins = 10 * max(N_abins_per_rbin)
    Pseudo_aBinWidth = 2*np.pi / PseudoN_abins
    HistogramForPlot = np.zeros((N_rbins, PseudoN_abins))

    # fill the HistogramForPlot with the Histogram data
    for n_r in range(N_rbins):
        aBinWidth = aBinWidth_perRbin[n_r]
        for n_a_pseudo in range(PseudoN_abins):
            n_a = int( ( n_a_pseudo+0.5 ) * Pseudo_aBinWidth / aBinWidth )
            HistogramForPlot[n_r,n_a_pseudo] = Histogram[n_r][n_a]

    # Define objects for plotting
    Pseudo_abins = Pseudo_aBinWidth * np.arange(PseudoN_abins+1)
    rbins = rBinWidth * np.arange(N_rbins+1)
    NonZeroHist = HistogramForPlot[ HistogramForPlot > 0 ]
    return Pseudo_abins, rbins, NonZeroHist, HistogramForPlot

'''
-----------------------------
-----------------------------
-----------------------------
----------- middle pannel (b)
- probability distribution of ball positions in circular arena
- sample trajectory of ball in circular arena
'''
fig, ax = plt.subplots(figsize=(3, 3), subplot_kw={'projection': 'polar'})
''' define color or objects '''
PathColor = 'turquoise' # ball color
DiscColor = cm.autumn(np.linspace(1,0, 2))[1] # Outer wall
AnulusColor = cm.winter(np.linspace(1,0, 2))[1] # Inner wall
colorMotor = 'lime'
colorWeight = 'gray'

''' define size of ball in the same units as the position '''
RadiusOfBall = 148/2
RadiusOuterWall = 956/2
RadiusInnerWall = 646/2

''' read data'''
DataPath = script_dir + '/Data/'
fileID = 'ParticleTrajectoryDiscArena.txt'
File = DataPath + fileID
pos = np.loadtxt(File)

''' shift center of cicular arena to origin '''
centerOfMass = np.mean(pos, axis=0)
pos = pos-centerOfMass
radius = np.linalg.norm(pos, axis = 1)
azimut = np.angle( pos[:,0] + 1j * pos[:,1]) + np.pi
Min = [min(pos[:,0]), min(pos[:,1])]
Max = [max(pos[:,0]), max(pos[:,1])]

''' compute probability distribution '''
N_radialBins = 20
Pseudo_abins, rbins, NonZeroHist, HistogramForPlot = computePolarHistogram(radius, azimut, N_radialBins)
A, R = np.meshgrid(Pseudo_abins, rbins)
pc = ax.pcolormesh(A, R, HistogramForPlot, cmap=cm.turbo_r, vmin=0, vmax=NonZeroHist.max())
circular_text_polar('arena', ax, radius=RadiusOuterWall, start_angle=0.84*np.pi, end_angle=0.7*np.pi, clockwise=False, color = DiscColor)

''' plot arena wall '''
for spine in ax.spines.values():
        spine.set_edgecolor(DiscColor)
        spine.set_linestyle('--')

ax.set_ylim(0, RadiusOuterWall)
ax.set_yticks([])
ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
ax.set_xticklabels(['0', '$\\pi/2$', '$\\pi$', '$3\\pi/2$'])
ax.tick_params(axis='x', which='major', pad=-6)
ax.grid(linestyle = '--', color='k')

''' sample trajectory and ball '''
start = 330
end = start+195
ax.plot(azimut[start:end],radius[start:end], '-', color = PathColor)
Ball = plt.Circle([azimut[start],radius[start]], RadiusOuterWall-radius[end] ,transform=ax.transData._b, facecolor = PathColor, edgecolor = None)
ax.add_patch(Ball)

''' Create the colorbar in that new axis '''
fig.colorbar(pc, label = '$P(r,\\phi)$')

''' safe image '''
plt.tight_layout()
ImName = script_dir + '/F1_PropOfPos.png'
plt.savefig(ImName, dpi = 300)
print('open ', ImName)

''' write out data '''
DataFilePath = script_dir + '/Data/F1b_HistogramForPlot.txt'
np.savetxt( DataFilePath, HistogramForPlot, header = ' probability distribution of position ' )
print('save', DataFilePath)

DataFilePath = script_dir + '/Data/F1b_AzimutBins.txt'
np.savetxt( DataFilePath, Pseudo_abins, header = ' azimutal bins of probability distribution of position ' )
print('save', DataFilePath)

DataFilePath = script_dir + '/Data/F1b_RadiusBins.txt'
np.savetxt( DataFilePath, rbins, header = ' radial bins of probability distribution of position ' )
print('save', DataFilePath)

DataFilePath = script_dir + '/Data/F1b_SampleTrajectory.txt'
np.savetxt( DataFilePath, np.column_stack((azimut[start:end], radius[start:end])) , header = ' ball trajectory in polar coordinates (azimut, radius)')
print('save', DataFilePath)
