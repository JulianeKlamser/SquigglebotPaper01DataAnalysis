import numpy as np
import matplotlib.pyplot as plt
from LibPlotFig3_Class_03 import AngleAnnotation
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

fig_width_pt = 246.0        # Get this from LaTeX using \the\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = 1.2     # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          #'figure.facecolor': None,
          'axes.facecolor': 'None',
          'axes.labelsize': 8,
          'axes.titlesize': 8,
          'font.size': 8,
          'axes.titlepad': 0,
          'axes.labelpad': 0,
          'lines.markersize': 5,
          'lines.linewidth' : 1,
          'legend.fontsize': 8,
          'legend.handlelength' : 1.,
          'legend.labelspacing' : 0.1,
          'legend.handletextpad' : 0.1,
          'legend.columnspacing' : 0.3,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'xtick.major.size' : 3,
          'xtick.minor.size' : 2,
          'ytick.major.size' : 3,
          'ytick.minor.size' : 2,
          'xtick.major.pad' : 1,
          'ytick.major.pad' : 0.5,
          'mathtext.fontset' : 'cm',
          'font.family': 'serif', #'serif' (e.g., Times),
          'figure.figsize': fig_size}
plt.rcParams.update(params)

''' read Data '''

DataFile = script_dir + '/Data/F3b_ProbOfInAngle.txt'
Data = np.loadtxt(DataFile)
InAngle, probOfInAngle = Data[:,0], Data[:,1]

DataFile = script_dir + '/Data/F3b_ProbOfOutAngle.txt'
Data = np.loadtxt(DataFile)
OutAngle, probOfOutAngle = Data[:,0], Data[:,1]

DataPath = script_dir + '/Data/F3a_TrajectoryInDisk.txt'
pos = np.loadtxt(DataPath)

DataFile = script_dir + '/Data/F3a_Parameters.txt'
Data = np.loadtxt(DataFile)
HeadTailLength, BulkRegionRadius, ArenaWallRadius = Data
HeadTailLength = int(HeadTailLength)

DataFile = script_dir + '/Data/F3a_InVectors_ID201.txt'
Data = np.loadtxt(DataFile)
InCrossPos_201, InTangent_201, InBulk_201 = Data

DataFile = script_dir + '/Data/F3a_OutVectors_ID201.txt'
Data = np.loadtxt(DataFile)
OutCrossPos_201, OutTangent_201, OutBulk_201 = Data

DataFile = script_dir + '/Data/F3a_Trajectory_ID201.txt'
Data = np.loadtxt(DataFile, dtype=int)
x_201, y_201 = pos[Data,0], pos[Data,1]

DataFile = script_dir + '/Data/F3a_InVectors_ID800.txt'
Data = np.loadtxt(DataFile)
InCrossPos_800, InTangent_800, InBulk_800 = Data

DataFile = script_dir + '/Data/F3a_OutVectors_ID800.txt'
Data = np.loadtxt(DataFile)
OutCrossPos_800, OutTangent_800, OutBulk_800 = Data

DataFile = script_dir + '/Data/F3a_Trajectory_ID800.txt'
Data = np.loadtxt(DataFile, dtype=int)
x_800, y_800 = pos[Data,0], pos[Data,1]

''' def figure '''
fig = plt.figure()

'''
---------------upper pannel a)------------------
------------------------------------------------
------------------------------------------------
------------------------------------------------
---------------upper pannel a)------------------
'''
L = 0.25
W = 1.-L
H = W / golden_mean
B = 1. - H
ax = fig.add_axes([L/1.9, B, W, H])
#ax.set_ylim(-250,250)
#ax.set_xlim(-250,250)
ax.set_ylim(-1.03*ArenaWallRadius,1.02*ArenaWallRadius)
ax.set_xlim(-1.02*ArenaWallRadius,1.02*ArenaWallRadius)
ax.axis('off')

''' PLOT ARENA AND BULK REGION BOUNDARY
    Note that the arena boundary is not the
    actual wall position, becasue the distance
    between the ball and the wall is at
    least the ball's radius.
'''
Arena = plt.Circle((0, 0), ArenaWallRadius, facecolor = 'w', edgecolor = 'k')
BulkRegion = plt.Circle((0, 0), BulkRegionRadius, facecolor = 'w', edgecolor = 'k', linestyle = ':')
ax.add_patch( Arena )
ax.add_patch( BulkRegion )

''' DRAW TRAJECTORY 204:
    - gray line : total trajectory
    - gray dots: positions before and after the bulk trajectory
    - orange dots: positions of the bulk trajectory
    - black square: start
    - black circle: end
'''
ax.plot( x_201, y_201, '-', color = 'gray')
ax.plot( x_201[0], y_201[0], 's', color = 'k', markersize = 6)
ax.plot(x_201[-1], y_201[-1], 'o', color = 'k', markersize = 7)
ax.text(x_201[0], 0.9*y_201[1], "START", ha = 'left', va = 'top')

x_head, y_head = x_201[:HeadTailLength], y_201[:HeadTailLength]
ax.plot( x_head, y_head, '.', color = 'gray')

x_bulk, y_bulk = x_201[HeadTailLength:-HeadTailLength], y_201[HeadTailLength:-HeadTailLength]
ax.plot( x_bulk, y_bulk , '.', color = 'orange')

x_tail, y_tail = x_201[-HeadTailLength:], y_201[-HeadTailLength:]
ax.plot(x_tail, y_tail, '.', color = 'gray')


''' Plot arrows that define incoming/outgoing angle
    1) tangential to boundary and pointing in
    the direction opposite head positions (outgoing)
    or opposite to tail positions (incoming)
    
    2) following the bulk trajectory leaving the boundary
    (outgoing) or following the history of the bulk
    trajectory before the collision (incoming)
'''
ax.annotate( "", InTangent_201, InCrossPos_201, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='lightseagreen', lw=1.5))
ax.annotate( "", InBulk_201 , InCrossPos_201, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='lightseagreen', lw=1.5))

ax.annotate( "", OutTangent_201, OutCrossPos_201, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='slateblue', lw=1.5))
ax.annotate( "", OutBulk_201 , OutCrossPos_201, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='slateblue', lw=1.5))

 
AngleAnnotation(InCrossPos_201, InTangent_201, InBulk_201, ax=ax, size=100, color = 'lightseagreen')
ax.text(0.98*InCrossPos_201[0], 0.8*InCrossPos_201[1], "$\\varphi_i$", ha = 'left', color = 'lightseagreen')
AngleAnnotation(OutCrossPos_201, OutTangent_201, OutBulk_201, ax=ax, size=120, color = 'slateblue')
ax.text(-13, 0.95*OutCrossPos_201[1], "$\\varphi_o$", ha = 'left', color = 'slateblue')

''' DRAW TRAJECTORY 800:
    - gray line : total trajectory
    - gray dots: positions before and after the bulk trajectory
    - orange dots: positions of the bulk trajectory
    - black square: start
    - black circle: end
'''
ax.plot( x_800, y_800, '-', color = 'gray')
ax.plot( x_800[0], y_800[0], 's', color = 'k', markersize = 6)
ax.plot(x_800[-1], y_800[-1], 'o', color = 'k', markersize = 7)
ax.text(1.3*x_800[0], 1.02*y_800[1], "START", ha = 'right', va = 'top')

x_head, y_head = x_800[:HeadTailLength], y_800[:HeadTailLength]
ax.plot( x_head, y_head, '.', color = 'gray')

x_bulk, y_bulk = x_800[HeadTailLength:-HeadTailLength], y_800[HeadTailLength:-HeadTailLength]
ax.plot( x_bulk, y_bulk , '.', color = 'orangered')

x_tail, y_tail = x_800[-HeadTailLength:], y_800[-HeadTailLength:]
ax.plot(x_tail, y_tail, '.', color = 'gray')


''' Plot arrows that define incoming/outgoing angle
    1) tangential to boundary and pointing in
    the direction opposite head positions (outgoing)
    or opposite to tail positions (incoming)
    
    2) following the bulk trajectory leaving the boundary
    (outgoing) or following the history of the bulk
    trajectory before the collision (incoming)
'''
ax.annotate( "", InTangent_800, InCrossPos_800, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='lightseagreen', lw=1.5))
ax.annotate( "", InBulk_800 , InCrossPos_800, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='lightseagreen', lw=1.5))

ax.annotate( "", OutTangent_800, OutCrossPos_800, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='slateblue', lw=1.5))
ax.annotate( "", OutBulk_800 , OutCrossPos_800, arrowprops=dict( arrowstyle='->', shrinkA=0, shrinkB=0, connectionstyle="arc3,rad=0", color='slateblue', lw=1.5))

 
AngleAnnotation(InCrossPos_800, InBulk_800, InTangent_800, ax=ax, size=100, color = 'lightseagreen')
ax.text(InCrossPos_800[0], 0.75*InCrossPos_800[1], "$\\varphi_i$", ha = 'left', color = 'lightseagreen')
AngleAnnotation(OutCrossPos_800, OutTangent_800, OutBulk_800, ax=ax, size=120, color = 'slateblue')
ax.text(50, 0.93*OutCrossPos_800[1], "$\\varphi_o$", ha = 'left', color = 'slateblue')


'''
---------------lower pannel b)------------------
------------------------------------------------
------------------------------------------------
------------------------------------------------
---------------lower pannel b)------------------
'''

L = 0.07
B = 0.07
W = 0.99 * (1. - L)
H = 0.95 * (1. - H - B)
ax = fig.add_axes([L, B, W, H])


''' plot probability distribution
    of incoming and outgoing angles
'''
ax.plot( InAngle, probOfInAngle, 'o-', label = '$P(\\varphi_i)$', color = 'lightseagreen')
ax.plot( OutAngle, probOfOutAngle, 's-', label = '$P(\\varphi_o)$', color = 'slateblue')

''' define plot makeup
'''
ax.set_xlim(0,np.pi)
ax.set_ylim(0,2)
ax.set_xticks([i*np.pi/4 for i in range(5)])
ax.set_xticklabels(['$0$', '$\\pi/4$', '$\\pi/2$', '$3\\pi/4$', '$\\pi$'])
ax.set_xlabel('$\\varphi$')
ax.legend(frameon = False)

''' create frame for pannel labels
'''
L = 0.0
B = 0.0
W = 1
H = 1
ax = fig.add_axes([L, B, W, H])
ax.axis('off')
ax.text(0.9, 0.45, '(a)')
ax.text(0.9, 0.09, '(b)')

''' save figure '''
FName = script_dir + '/F3.png'
plt.savefig(FName, dpi = 300)
print('open ', FName)

''' compute and print kullback leibler divergence
'''
mask = (probOfOutAngle > 0) & (probOfInAngle > 0)
print('KLD = ', np.sum(probOfOutAngle[mask] * np.log( probOfOutAngle[mask]/probOfInAngle[mask] ) ))
