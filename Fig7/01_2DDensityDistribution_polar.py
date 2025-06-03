import numpy as np
import os
import pickle
script_dir = os.path.dirname(os.path.abspath(__file__))

from Functions import computePolarHistogram
from Functions import PlotDensityDistribution

''' list of paths :
    one path per simulaiton parameter
'''
BasePath = script_dir + '/Data/ConfinedCircular_N000256L28.000000/N000256L28.000000/v_o1.000000D_r1.000000/'
PathList = [BasePath + 'gamma1.000000alignCutOffBySigma2.000000', BasePath + 'gamma0.100000alignCutOffBySigma2.000000', BasePath + 'gamma0.010000alignCutOffBySigma2.000000']

''' loop over paths to evaluate each
    parameter set '''
    
NSample = 1000 # number of configuration to average over
for i, SimuPath in enumerate(PathList):
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
    
    ''' read final configuration '''
    Path_ToFinalConf = SimuPath + '/Particles.txt'
    Data = np.loadtxt(Path_ToFinalConf)
    NP = len(Data) # number of particles
    pos_final = Data[:,:2]
    radii = Data[:,5] # sigma_i / 2 = particle 'radii'
    assert( np.all(radii == radii[0])) # all particles have the same sigma
    sigmaP_Half = radii[0]
    RadialBinWidth = 2*sigmaP_Half
    
    ''' compute and print density '''
    rho = NP * sigmaP_Half**2 / WallRadius**2
    print('rho = N (sigma/2)^2 / R^2 = %.4f <----------' % rho)
    
    ''' read time file of the simulation '''
    Path_ToTime = SimuPath + '/Time.txt'
    Time = np.loadtxt(Path_ToTime, usecols=(0, 1), skiprows = 1)
    Time = Time[-NSample:] # choose the last NSample configurations
    
    ''' loop over each configuration and compute histogram '''
    for time in Time:
        ID, t_physical = time
        
        ''' read configuration at t_physical '''
        Path_ToConf = SimuPath + '/Particles%d.txt' % ID
        Data = np.loadtxt(Path_ToConf)
        pos = Data[:,:2]
        
        ''' shift positions such that center of confining
            wall is at the origin '''
        pos -= WallRadius
        
        ''' get center of mass position (X, Y) and apply
            a coordinate transformation (rotation)
            such that (X, Y) -> (X', 0) '''
        center_of_mass_position = np.mean(pos, axis = 0)
        COM_azimut = np.angle(center_of_mass_position[0] + 1j * center_of_mass_position[1])
        r_coord, azimut = np.linalg.norm( pos, axis = 1), np.angle( pos[:,0] + 1j * pos[:,1]) # polar coordinates of particles
        azimut = np.mod(azimut - COM_azimut, 2*np.pi) # rotate such that (X, Y) -> (X', 0)
        
        ''' compute histogram of particle positions 
            in polar coordinates '''
        N_radialBins = int( WallRadius // RadialBinWidth )
        Pseudo_abins, rbins, DensityDist = computePolarHistogram(r_coord, azimut, N_radialBins, WallRadius)
        
        ''' average over time '''
        if not 'averagedDensityDist' in locals():
            averagedDensityDist = np.copy(DensityDist)
        else:
            averagedDensityDist += DensityDist
    averagedDensityDist /= NSample
    
    ''' visualize density distribution together with 
        particle pistions of several configurations'''
    for time in Time[::100]:
        # Plot the density distribution
        PlotDensityDistribution(Pseudo_abins, rbins, averagedDensityDist, SimuPath, time, WallRadius, sigmaP_Half, gamma)
    
    ''' save density distribution as dictionary in
        a pickle '''
    OutDir = SimuPath + '/00EvaluatedData'
    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    OutDictionary = {"Pseudo_abins" : Pseudo_abins, "rbins" : rbins, "averagedDensityDist" : averagedDensityDist}
    OutFileName = OutDir + '/DensityDistribution.pkl'
    file = open(OutFileName, 'wb')
    pickle.dump(OutDictionary, file)
    file.close()
    print('dump', OutFileName)
    
    del averagedDensityDist
