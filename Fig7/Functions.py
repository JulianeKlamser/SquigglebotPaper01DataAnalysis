import numpy as np
import os
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
from matplotlib.ticker import ScalarFormatter


def computePolarHistogram(radius, azimut, N_radialBins, maxR):
    '''
    radius & azimut define particle positions (circular arena centered at origin)
    N_radialBins : number of radial bins
    '''
    N = len(radius) # number of positions
    azimut = np.mod(azimut, 2*np.pi)
    
    # number of radial and azimutal bins
    N_rbins = N_radialBins # radial number of bins
    N_abins_per_rbin = np.array([4] + [4*(n_r) for n_r in range(1,N_rbins)]) # per radial bin, there is a growing number of azimutal bins

    # radial and azimutal bin with
    rBinWidth = maxR / N_rbins # constant radial bin bidth
    aBinWidth_perRbin = [ 2*np.pi/abins for abins in N_abins_per_rbin] # azimutal bin width depends on radial bin

    # initialize histogram to zeros
    Histogram = [[0] * N_abins_per_rbin[n_r] for n_r in range(0,N_rbins)] # dimension = N_rbins x N_abins_per_rbin

    # fill the histogram
    for i, r in enumerate(radius):
        n_r = int( radius[i] / rBinWidth )
        if n_r >= N_radialBins:
            n_r -= 1
        if n_r > N_radialBins:
            print('particles get closer to the wall than you think')
            exit()

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

    return Pseudo_abins, rbins, HistogramForPlot


def PlotDensityDistribution(Pseudo_abins, rbins, averagedDensityDist, SimuPath, time, WallRadius, sigmaP_Half, gamma):

    fig = plt.figure()
    width, hight = 5, 5
    fig.set_size_inches(width, hight)
    ax_polar = fig.add_axes([0.07, 0.05, 0.75, 0.75], polar=True)
    ax_carte = fig.add_axes([0.07, 0.05, 0.75, 0.75])
    ax_carte.axis('off')
    ax_carte.set_facecolor('none')
    ax_cBar = fig.add_axes([0.87, 0.05, 0.1, 0.8])
    
    ''' plot average density distribution '''
    A, R = np.meshgrid(Pseudo_abins, rbins)
    c = ax_polar.pcolormesh(A, R, averagedDensityDist, cmap=cm.turbo, vmin=0, vmax=averagedDensityDist.max())
    cBar = plt.colorbar(c, ax=ax_cBar, fraction = 1, shrink = 0.8, label = '$P(r, \\theta)$')
    cBar.ax.yaxis.labelpad = 0
    ax_cBar.axis('off')
    
    ''' format the color bar '''
    # Format ticks in scientific notation
    formatter = ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))  # Always use scientific notation
    cBar.ax.yaxis.offsetText.set_position((2, 1))
    cBar.ax.yaxis.set_major_formatter(formatter)
    
    ''' read configuration at t_physical '''
    ID, t_physical = time
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
    azimut -= COM_azimut # rotate such that (X, Y) -> (X', 0)
    pos[:,0] = r_coord * np.cos(azimut)
    pos[:,1] = r_coord * np.sin(azimut)
    
    ''' plot particle positions '''
    Width, Hight, Angle = 2*sigmaP_Half, 2*sigmaP_Half, 0
    collection = EllipseCollection( Width, Hight, Angle, units='x', offsets=pos, transOffset=ax_carte.transData, edgecolor = 'w', facecolor = 'none')
    ax_carte.add_collection(collection)
    
    collection2 = EllipseCollection( Width, Hight, Angle, units='x', offsets=pos[0], transOffset=ax_carte.transData, edgecolor = 'w', facecolor = 'w')
    ax_carte.add_collection(collection2)
    ax_carte.set_xlim(-WallRadius,WallRadius)
    ax_carte.set_ylim(-WallRadius,WallRadius)
    
    fig.suptitle('Two-Dimensional density distribution in polar coords\n Particle positions of configuration %d for $\\gamma = %.2f$\n' % (ID, gamma))
    FileDir = SimuPath + '/00Snapshots'
    os.makedirs(FileDir, exist_ok=True)
    FileName = FileDir + '/PolarConf_%d.png' % ID
    plt.savefig(FileName)
    print('open', FileName)
    plt.close()
