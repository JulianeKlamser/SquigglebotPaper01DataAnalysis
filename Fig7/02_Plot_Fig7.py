import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.pyplot import cm
import pickle
from matplotlib.collections import EllipseCollection
from matplotlib.ticker import ScalarFormatter
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

fig_width_pt = 246.0*2        # Get this from LaTeX using \the\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = 0.2     # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.facecolor': 'None',
          'axes.labelsize': 8,
          'axes.titlesize': 8,
          'font.size': 8,
          'axes.titlepad': 0,
          'axes.labelpad': 0,
          'lines.markersize': 2,
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
fig = plt.figure()

''' pannel a) of Fig 7 - experiment '''
L=0.0
B=0.0
H = 1.
W=H*golden_mean
ax = fig.add_axes([L, B, W, H])
ax.set_axis_off()

# import photo of experiment
img = Image.open(script_dir + '/Config2_20balls.png')
ax.imshow(img)
# highlight the transparent circular wall
width, height = img.size
Center = [0.97*width/2, 0.98*height/2]
Radius = 0.89*(width/2)
circle = plt.Circle(Center, Radius, color='k', linestyle = ':', fill=False, linewidth=1)
ax.add_patch(circle)


''' pannel b), c), d) of Fig 7 - simulations '''
#define data path
BasePath = script_dir + '/Data/ConfinedCircular_N000256L28.000000/N000256L28.000000/v_o1.000000D_r1.000000/'
PathList = [BasePath + 'gamma1.000000alignCutOffBySigma2.000000', BasePath + 'gamma0.100000alignCutOffBySigma2.000000', BasePath + 'gamma0.010000alignCutOffBySigma2.000000']

#------- START : define axis for each pannel ------
WCBar = 0.1 # width of color bar axis
L = 0.01 # horizontal spacing between individual axis
HCBar = 0.9*H # hight of color bar
HLabel = 1- B+0.9*HCBar # hight of axis containing the scintific scale
WLabel = (1 - 4*W - 3*L)/3 # width of axis containing the scintific scale

LNew = W+L
ax1_polar = fig.add_axes([LNew, B, W, H], polar=True)
ax1_carte = fig.add_axes([LNew, B, W, H])
LNew += W + L
ax1_cbar = fig.add_axes([LNew, B, WCBar, HCBar])
ax1_label = fig.add_axes([LNew, B+0.91*HCBar, WLabel, HLabel])

LNew += WLabel
ax2_polar = fig.add_axes([LNew, B, W, H], polar=True)
ax2_carte = fig.add_axes([LNew, B, W, H])
LNew += W + L
ax2_cbar = fig.add_axes([LNew, B, WCBar, HCBar])
ax2_label = fig.add_axes([LNew, B+0.91*HCBar, WLabel, HLabel])

LNew += WLabel
ax3_polar = fig.add_axes([LNew, B, W, H], polar=True)
ax3_carte = fig.add_axes([LNew, B, W, H])
LNew += W + L
ax3_cbar = fig.add_axes([LNew, B, WCBar, HCBar])
ax3_label = fig.add_axes([LNew, B+0.91*HCBar, WLabel, HLabel])
#------- End : define axis for each pannel ------

# create list of axis to loop over them
axs_polar = [ax1_polar, ax2_polar, ax3_polar]
axs_carte = [ax1_carte, ax2_carte, ax3_carte]
axs_cbar = [ax1_cbar, ax2_cbar, ax3_cbar]
axs_label = [ax1_label, ax2_label, ax3_label]

# loop over pannel b), c), d)
for i, DataPath in enumerate(PathList):
    ax_p = axs_polar[i]
    ax_c = axs_carte[i]
    ax_cbar = axs_cbar[i]
    ax_label = axs_label[i]
    
    ''' start - plot density distribution '''
    # read density distribution
    path = DataPath + '/00EvaluatedData'
    DensityDistFile = path + '/DensityDistribution.pkl'
    DensityDistDict = pickle.load( open( DensityDistFile, "rb" ) )

    # plot density distribution
    A, R = np.meshgrid(DensityDistDict['Pseudo_abins'], DensityDistDict['rbins'])
    averagedDensityDist = DensityDistDict['averagedDensityDist']
    MaxP = averagedDensityDist.max()
    c = ax_p.pcolormesh(A, R, averagedDensityDist, cmap=cm.turbo, vmin=0, vmax=MaxP)
    cBar = plt.colorbar(c, ax=ax_cbar, fraction = 1, shrink = 0.8, label = '$P(r, \\theta)$')
    # Format ticks in scientific notation
    formatter = ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))  # Always use scientific notation
    cBar.ax.yaxis.offsetText.set_position((7, 1))
    cBar.ax.yaxis.set_major_formatter(formatter)
    # get exponent of scientific notion
    exponent = int(np.floor(np.log10(abs(MaxP))))
    
    # Hide grid, ticks, and labels
    ax_p.grid(False)
    ax_p.set_xticks([])  # Remove angular ticks
    ax_p.set_yticks([])  # Remove radial ticks
    ax_p.set_xticklabels([])
    ax_p.set_yticklabels([])
    ax_p.spines['polar'].set_color('k')
    ax_p.spines['polar'].set_linestyle(':')
    ax_p.spines['polar'].set_linewidth(1.)
    ax_cbar.set_axis_off()
    
    # overwrite automatic scientrific formate scale
    # on top of the color bar as latex style
    ax_label.text(0,0, '$\\times 10^{%d}$' % exponent, ha = 'left', va = 'bottom')
    ax_label.set_facecolor('w')
    ax_label.set_xlim(0,1)
    ax_label.set_ylim(0,1)
    for spine in ax_label.spines.values():
        spine.set_visible(False)
    ax_label.set_xticks([])
    ax_label.set_yticks([])
    ax_cbar.set_zorder(1)
    ax_label.set_zorder(2)
    
    
    ''' start - plot paricle positions '''
    # read the size of the simulaiton box
    Periods = np.loadtxt(DataPath + '/Periods.txt')
    WallRadius = (Periods[0,1] - 5.) / 2
    
    # read the ID of the final configuration
    Time = np.loadtxt(DataPath + '/Time.txt', usecols=(0, 1), skiprows = 1)
    FinalTime = Time[-1]
    
    # read final configuration
    Conf = np.loadtxt(DataPath + '/Particles%d.txt' % FinalTime[0])
    N = len(Conf) # number of particles
    Radii = Conf[:,5] # sigma_i / 2 in WCA potential (radii of particles)
    Pos = Conf[:,:2] - WallRadius # shift positions such that the center of the wall is at the origin
    r, azimut = np.linalg.norm(Pos, axis = 1), np.angle(Pos[:, 0] + 1j * Pos[:, 1]) # polar coordinates of positions
    
    # get center of mass position on its polar angle
    CenterOfMass = np.mean(Pos, axis = 0)
    azimutCM = np.angle(CenterOfMass[0] + 1j * CenterOfMass[1])
    
    # coordinate transformation (rotation) such that
    # the center of mass lies on the positive x-axis
    # i.e. azimutCM -> 0
    azimut = azimut - azimutCM
    Pos[:,0] = r * np.cos(azimut)
    Pos[:,1] = r * np.sin(azimut)
        
    # plot all particles as empty circles with balck boundary
    Width, Hight, Angle = 2*Radii, 2*Radii, np.zeros(N)
    collection = EllipseCollection( Width, Hight, Angle, units='x', offsets=Pos, transOffset=ax_c.transData, edgecolor = 'k', facecolor = 'none', linewidth = 0.4)
    ax_c.add_collection(collection)
    
    # plot the active particle as whit circle
    Width, Hight, Angle = 2*Radii[0], 2*Radii[0], 0.0
    collection2 = EllipseCollection( Width, Hight, Angle, units='x', offsets=Pos[0], transOffset=ax_c.transData, edgecolor = 'w', facecolor = 'w', linewidth = 0.4)
    ax_c.add_collection(collection2)
    
    # set plot limits
    ax_c.set_ylim(-WallRadius,WallRadius)
    ax_c.set_xlim(-WallRadius,WallRadius)
    ax_c.set_axis_off()
    
    
''' add pannel labels '''
ax = fig.add_axes([0, 0, 1., 1.])
ax.axis('off')
ax.text(0.002, 0.05, '(a)')
ax.text(0.21, 0.05, '(b)')
ax.text(0.47, 0.05, '(c)')
ax.text(0.74, 0.05, '(d)')

''' save figure '''
FName = script_dir + '/Fig7_new.png'
plt.savefig(FName, dpi = 300)
print('open ', FName)
