import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

print('\nThis may take some time.\n')


'''
-----------------------------
-----------------------------
-----------------------------
----------- upper right pannel (c)
- MSD
'''
fig, ax = plt.subplots(figsize=(4, 4))
ax.set_xlabel("$t$ in $s$")
ax.set_xlim(0.04, 1.4*10**3)
ax.set_ylim(0.001, 10**5)


''' define color or objects '''
DiscColor = cm.autumn(np.linspace(1,0, 2))[1] # Outer wall
AnulusColor = cm.winter(np.linspace(1,0, 2))[1] # Inner wall

''' read trajectory data '''
DataPath = script_dir + '/Data/'
File = DataPath + 'ParticleTrajectoryDiscArena.txt'
pos = np.loadtxt(File)


''' shift center of cicular arena to origin '''
centerOfMass = np.mean(pos, axis=0)
pos = pos-centerOfMass
Min = [min(pos[:,0]), min(pos[:,1])]
Max = [max(pos[:,0]), max(pos[:,1])]

''' define time and length units '''
FramesPerSecond = 25
timeUnitPerFrame = 1. / FramesPerSecond #0.04 s
lengthUnit = 8.2 / 148 #cm/pixel

''' get trajectory length and define number of data points for MSD '''
FinalTrajTimeInFrames = len(pos) # units in number of frames
N_MSD_DataPoints = 120 # final time for MSD = N_MSD_DataPoints * timeUnitPerFrame

''' init MSD and physical time '''
MSD = np.zeros(N_MSD_DataPoints)
TimeIndex_List = np.arange(N_MSD_DataPoints) + 1
time = TimeIndex_List * timeUnitPerFrame

''' rescale positions to physical units '''
pos = pos*lengthUnit

''' compute MSD with sliding average '''
NValuesPerDataPoint = FinalTrajTimeInFrames - N_MSD_DataPoints
for time_i in range(NValuesPerDataPoint):
    slidingTimeIndex_List = TimeIndex_List + time_i
    DisplacementOverTime = pos[slidingTimeIndex_List] - pos[time_i] # array - float
    MSD += np.sum( DisplacementOverTime**2, axis = 1)
MSD /= NValuesPerDataPoint

''' plot MSD and t^2 scaling'''
lns1 = ax.loglog( time, MSD, '-', label = "$\\mathcal{D}$: MSD in $\\rm{cm}^2$", color = DiscColor)
theroy = time**2
pref = MSD[1]/theroy[1]
lns2 = ax.loglog(time, pref*theroy, 'k--', label = '$\\propto t^2$')

''' create legend '''
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='lower right', ncol = 2, frameon = False, borderaxespad = 0)

''' save image '''
plt.tight_layout()
ImName = script_dir + '/F1_MSD.png'
plt.savefig(ImName, dpi = 300)
print('open ', ImName)

''' write out data '''
DataFilePath = DataPath + 'F1c_MSD.txt'
np.savetxt( DataFilePath, np.column_stack((time, MSD)), header = ' time MSD(time)' )
print('save', DataFilePath)

DataFilePath = DataPath + 'F1c_MSD_BallisticTheory.txt'
np.savetxt( DataFilePath, np.column_stack((time, pref*theroy)), header = ' time a*time**2)' )
print('save', DataFilePath)
