import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

DataPath = script_dir + '/Data/'
PlotPath = script_dir + '/Plots/'

''' devided time series of current in chuncks 
    of N values '''
N = 5000

''' defines the time frame within the current 
    is stable '''
StableLimits = [ 6*N, 46*N]

''' list of cases to monitor '''
FileNames = [ 'reference', 'boundary', 'bulk']
color = ['gray', 'palevioletred', 'skyblue']

fig, axs = plt.subplots(1, 2, dpi = 100, figsize = (9,4))
ax1, ax2 = axs
''' loop over cases '''
for i, FileID in enumerate(FileNames):

    ''' read current '''
    filePath = DataPath + 'Current_%s.txt' % FileID
    current = np.loadtxt(filePath)
    
    ''' Ensure the total number of frames is a multiple of N '''
    N_frames = len(current)
    M = N_frames // N # number of chunks
    N_frames = int( M * N )
    current = current[:N_frames]
    
    ''' get minimum and maximum current '''
    MinCurrent, MaxCurrent = current.min(), current.max()
    
    ''' reshape current into blocks of size N and 
        compute average and variance per block '''
    current_inBlocksOfN = current.reshape(-1, N) # [[c_0,c_1,c_2,...,c_T-1],[c_2T, c_2T+1, ..., c_3T-1], ...]
    AveragedCurrent = np.mean( current_inBlocksOfN, axis=1 ) #averave per block
    VarianceCurrent = np.var( current_inBlocksOfN, axis=1 ) #variance per block
    
    ''' plot average and variance '''
    n_list = np.arange(M)
    ax1.plot( n_list, AveragedCurrent, 'o-', label = FileID, color = color[i])
    ax2.plot( n_list, VarianceCurrent, 'o-', label = FileID, color = color[i])
    
    ''' define a stable region of the current '''
    StableCurrent = current[ StableLimits[0] : StableLimits[1] ]
    
    ''' reshape stable current into blocks of size N and 
    compute average and variance per block '''
    SCurrent_inBlocksOfN = StableCurrent.reshape(-1, N) # [[c_0,c_1,c_2,...,c_T-1],[c_2T, c_2T+1, ..., c_3T-1], ...]
    Averaged_SCurrent = np.mean( SCurrent_inBlocksOfN, axis=1 ) #averave per block
    Variance_SCurrent = np.var( SCurrent_inBlocksOfN, axis=1 ) #variance per block
    
    N_frames = len(StableCurrent)
    M_stable = N_frames // N # number of chunks
    n_list_stable =  StableLimits[0] // N + np.arange(M_stable)
    if i == 2:
        ax1.plot(n_list_stable, Averaged_SCurrent, '.', color = 'k', label = 'stable range')
    ax1.plot(n_list_stable, Averaged_SCurrent, '.', color = 'k')
    ax2.plot(n_list_stable, Variance_SCurrent, '.', color = 'k')
    
''' define lables, title and limits '''
ax1.legend()
ax1.set_xlim(0, np.max(n_list))
ax1.set_ylim(0.04, 0.08)

ax2.set_xlim(0, np.max(n_list))
VarLim = [0.0003, 0.0010]
ax2.set_ylim(VarLim)

ax1.set_ylabel('$\\langle\\, j \\,\\rangle(n)$')
ax1.set_xlabel('$n$')

ax2.set_ylabel('$\\langle\\, j^2 \\,\\rangle(n) - \\langle\\, j \\,\\rangle^2(n)$')
ax2.set_xlabel('$n$')

plt.suptitle('$j(i)$ is the current at frame i, $\\langle\\, j \\,\\rangle(n) = \\frac{1}{N} \\sum_{i=0}^{N-1} j(i + n \\cdot N)$,\n $\\langle\\, j^2 \\,\\rangle(n)$ defined analogously, $N = %d$, $n\\in{0,1,2,...,%d}$' % ( N, M-1))
plt.tight_layout()

ImageName = 'StabilityMonitor.png'
FileName = PlotPath + ImageName
print('open', FileName)
plt.savefig(FileName, dpi = 300)

FineName = DataPath + 'StableLimits.txt'
np.savetxt(FineName, StableLimits, header = 'minimum and maximum frame defining the region of stable current', fmt = '%d')
print('write ', FineName)
