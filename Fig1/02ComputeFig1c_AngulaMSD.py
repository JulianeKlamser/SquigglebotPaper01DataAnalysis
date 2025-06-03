import numpy as np
import matplotlib.pyplot as plt
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

print('\nThis may take some time.\n')

''' read data '''
DataPath = script_dir + '/Data/'
fileID = 'ParticleTrajectoryAnnulusArena.txt'
Pos = np.loadtxt(DataPath + fileID)
N_pos = len(Pos)


''' shift center of mass to origin '''
center = np.mean(Pos, axis = 0)
Pos -= center

''' compute Azimute '''
Azimute = np.angle( Pos[:,0] + 1j * Pos[:,1], deg = False) + np.pi # range in [0,2pi]
UnfoldedAzimute = np.copy(Azimute)
borderCorssing = 0
for i, a_1 in enumerate(Azimute[:-1]):
    a_2 = Azimute[i+1]
    diff = a_2 - a_1
    if np.abs(diff) > np.pi:
        #crossed periodic boundary
        if diff > 0:
            borderCorssing -= 1
        else:
            borderCorssing += 1
    UnfoldedAzimute[i+1] = a_2 + borderCorssing * (2.*np.pi)
        
''' define time and length units '''
FramesPerSecond = 25
timeUnitPerFrame = 1. / FramesPerSecond #0.04 s
lengthUnit = 8.2 / 148 # cm/pixel

''' get trajectory length and define number of data points for MSD '''
FinalTrajTimeInFrames = N_pos # units in number of frames
N_MSD_DataPoints = int(3.5*10**4) # final time for MSD = N_MSD_DataPoints * timeUnitPerFrame

''' init MSD and physical time '''
MSD_angular = np.zeros(N_MSD_DataPoints)
TimeIndex_List = np.arange(N_MSD_DataPoints) + 1
time = TimeIndex_List * timeUnitPerFrame

''' compute MSD with sliding average '''
NValuesPerDataPoint = FinalTrajTimeInFrames - N_MSD_DataPoints
for time_i in range(NValuesPerDataPoint):
    slidingTimeIndex_List = TimeIndex_List + time_i
    dTheta = UnfoldedAzimute[slidingTimeIndex_List] - UnfoldedAzimute[time_i] # array - float
    MSD_angular += dTheta**2
MSD_angular /= NValuesPerDataPoint

''' plot MSD and scalings '''
plt.loglog(time, MSD_angular, ':')
theroy = time**2
pref = MSD_angular[1]/theroy[1]
plt.loglog(time, pref*time**2, '-')
theroy = time
pref = 1.1*MSD_angular[-1]/theroy[-1]
plt.loglog(time, pref*time, '-')

plt.xlabel('$t$')
plt.ylabel('$\\Delta\\theta^2$ in $rad^2$')

''' save image '''
plt.tight_layout()
ImName = script_dir + '/F1_AngularMSD.png'
plt.savefig(ImName, dpi = 300)
print('open ', ImName)

''' safe angular MSD '''
OutFile = script_dir + '/Data/F1c_AngularMSD.txt'
np.savetxt(script_dir + '/Data/F1c_AngularMSD.txt', [time, MSD_angular], header = ' time MSD_angular(time)')
print('save', OutFile)
