import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import os
script_dir = os.path.dirname(os.path.abspath(__file__))


fig, ax = plt.subplots(figsize=(4, 4))
ax.set_xlabel("$|\\vec{v}|$ in $cm/s$")
ax.set_ylabel("$P\\left(\\vert\\vec{v}\\vert\\right)$")
DiscColor = cm.autumn(np.linspace(1,0, 2))[1] # Outer wall
AnulusColor = cm.winter(np.linspace(1,0, 2))[1] # Inner wall
'''
-----------------------------
-----------------------------
-----------------------------
----------- lower right pannel (e)
- probability of speed for Disc and annulus arena
'''

''' Start Disc arena, read data '''
DataPath = script_dir + '/Data/'
File = DataPath + 'ParticleTrajectoryDiscArena.txt'
pos = np.loadtxt(File)
NFrames = len(pos)

''' define time and length units & compute speed'''
FramesPerSecond = 25
timeUnitPerFrame = 1. / FramesPerSecond #0.04 s
lengthUnit = 8.2 / 148 #cm/pixel
pos = pos*lengthUnit
speed = np.linalg.norm( np.diff(pos, axis = 0), axis = 1) /timeUnitPerFrame

''' compute histogram of speed '''
maxSpeed = max(speed)
NBins = 100
binWidth = maxSpeed/NBins
SpeedHistogram_Disc = np.zeros(NBins)
BinCenters_Disc = ( 0.5 + np.arange(NBins) ) * binWidth
for s in speed:
    if s == maxSpeed:
        SpeedHistogram_Disc[-1] += 1
    else:
        i = int(s/binWidth)
        SpeedHistogram_Disc[i] += 1
SpeedHistogram_Disc /= ( NFrames * binWidth )
assert ( ( ( sum(SpeedHistogram_Disc) * binWidth ) - 1 ) < 10**(-4))

ax.plot(BinCenters_Disc, SpeedHistogram_Disc, 'o-', color = DiscColor, label = '$\\mathcal{D}$')

''' End Disc arena '''

''' Start annulus arena, read data '''
File = DataPath + 'ParticleTrajectoryAnnulusArena.txt'
pos = np.loadtxt(File)

''' positions in cm and compute speed'''
pos = pos*lengthUnit
speed = np.linalg.norm( np.diff(pos, axis = 0), axis = 1) /timeUnitPerFrame

''' compute histogram '''
maxSpeed = max(speed)
NBins = 100
binWidth = maxSpeed/NBins
SpeedHistogram_Annulus = np.zeros(NBins)
BinCenters_Annulus = ( 0.5 + np.arange(NBins) ) * binWidth
for s in speed:
    if s == maxSpeed:
        SpeedHistogram_Annulus[-1] += 1
    else:
        i = int(s/binWidth)
        SpeedHistogram_Annulus[i] += 1
SpeedHistogram_Annulus /= ( NFrames * binWidth )
assert ( ( ( sum(SpeedHistogram_Annulus) * binWidth ) - 1 ) < 10**(-4))
ax.plot(BinCenters_Annulus, SpeedHistogram_Annulus, 'o-', color = AnulusColor, markeredgecolor = AnulusColor, markerfacecolor = 'w', label = '$\\mathcal{A}$')


''' save image '''
ax.set_xlim(0)
ax.legend(loc = 'upper right', frameon=False)
plt.tight_layout()
ImName = script_dir + '/F1_Speed.png'
plt.savefig(ImName, dpi = 300)
print('open ', ImName)

''' write out data '''
DataFilePath = DataPath + 'F1e_SpeedDistributionDisc.txt'
np.savetxt( DataFilePath, np.column_stack((BinCenters_Disc, SpeedHistogram_Disc)), header = 'norm(v) P(norm(v))' )
print('save', DataFilePath)

DataFilePath = DataPath + 'F1e_SpeedDistributionAnnulus.txt'
np.savetxt( DataFilePath, np.column_stack((BinCenters_Annulus, SpeedHistogram_Annulus)), header = 'norm(v) P(norm(v))' )
print('save', DataFilePath)
