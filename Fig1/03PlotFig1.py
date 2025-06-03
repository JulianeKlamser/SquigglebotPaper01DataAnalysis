import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from PIL import Image
from matplotlib.pyplot import cm
# Only import FigureCanvasAgg if necessary
try:
    from matplotlib.backends.backend_agg import FigureCanvasAgg
except ImportError:
    FigureCanvasAgg = None
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

def circular_text_polar(text, ax, radius=1.0, start_angle=0, end_angle=np.pi, clockwise=False, color = 'red'):
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

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)

fig_width_pt = 246.0*2        # Get this from LaTeX using \the\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = 0.3     # Aesthetic ratio
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

''' define Figure '''
fig = plt.figure()
ax = fig.add_axes([0, 0, 1., 1.])
ax.axis('off')
ax.text(0.23, 0.9, '(a)')
ax.text(0.45, 0.9, '(b)')
ax.text(0.65, 0.9, '(c)')
ax.text(0.65, 0.4, '(e)')
ax.text(0.53, 0.05, '(d)')

''' define color or objects '''
PathColor = 'turquoise' # ball color
DiscColor = cm.autumn(np.linspace(1,0, 2))[1] # Outer wall
AnulusColor = cm.winter(np.linspace(1,0, 2))[1] # Inner wall
colorMotor = 'lime'
colorWeight = 'gray'

'''
-----------------------------
-----------------------------
-----------------------------
----------- middle pannel (b)
- probability distribution of ball positions in circular arena
- sample trajectory of ball in circular arena
'''
L=0.31
B=0.34
W=.18
H = W / golden_mean
ax = fig.add_axes([L , B, W, H], polar=True)

''' define size of ball in the same units as the position '''
RadiusOfBall = 148/2
RadiusOuterWall = 956/2
RadiusInnerWall = 646/2

''' read data'''
DataPath = script_dir + '/Data/'
DataFilePath = DataPath + 'F1b_HistogramForPlot.txt'
HistogramForPlot = np.loadtxt( DataFilePath)
DataFilePath = DataPath + 'F1b_AzimutBins.txt'
Pseudo_abins = np.loadtxt( DataFilePath )
DataFilePath = DataPath + 'F1b_RadiusBins.txt'
rbins = np.loadtxt( DataFilePath )
DataFilePath = DataPath + 'F1b_SampleTrajectory.txt'
Data = np.loadtxt( DataFilePath )
Azimut, Radius = Data[:,0], Data[:,1]

''' plot histogram '''
A, R = np.meshgrid(Pseudo_abins, rbins)
pc = ax.pcolormesh(A, R, HistogramForPlot, cmap=cm.turbo_r, vmin=0, vmax=HistogramForPlot.max())

''' plot arena wall '''
for spine in ax.spines.values():
        spine.set_edgecolor(DiscColor)
        spine.set_linestyle('--')
ax.set_ylim(0, RadiusOfBall+Radius[-1])
ax.set_yticks([])
ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2])
ax.set_xticklabels(['0', '$\\pi/2$', '$\\pi$', '$3\\pi/2$'])
ax.tick_params(axis='x', which='major', pad=-6)
ax.grid(linestyle = '--', color='k')
circular_text_polar('arena', ax, radius=RadiusOfBall+Radius[-1], start_angle=0.84*np.pi, end_angle=0.7*np.pi, clockwise=False, color = DiscColor)

''' sample trajectory and ball '''
ax.plot(Azimut, Radius, '-', color = PathColor)
Ball = plt.Circle([Azimut[0], Radius[0]], RadiusOfBall,transform=ax.transData._b, facecolor = PathColor, edgecolor = None)
ax.add_patch(Ball)

''' add the color bar in a new axis '''
cbaxes = fig.add_axes([L+1.1*W, 1.15*B, 0.01, 0.9*H])
cb = plt.colorbar(pc, cax = cbaxes)
cb.set_label('      $P(r,\\phi)$')
cb.ax.yaxis.offsetText.set_position((4, 1))

totH = 1.15*B+0.91*H
scaleAxis = fig.add_axes([L+1*W, totH, 0.06, 1-totH])
scaleAxis.set_facecolor('w')
scaleAxis.set_xlim(0,1)
scaleAxis.set_ylim(0,1)
scaleAxis.text(0.25,0,'$\\times 10^{-5}$')
for spine in scaleAxis.spines.values():
    spine.set_visible(False)
scaleAxis.set_xticks([])
scaleAxis.set_yticks([])

'''
-----------------------------
-----------------------------
-----------------------------
----------- lower middle pannel (d)
- pictogram of experimental setup
'''
L=0.42
B=0.00
W=.18/1.5
H = W / golden_mean
ax = fig.add_axes([L , B, W, H])

''' plot disc arena wall = outer wall of annulus (circle1) & plot inner wall of annulus (circle2) '''
circle1 = plt.Circle((0, 0), RadiusOuterWall, facecolor = 'w', edgecolor = DiscColor, linestyle = '--')
circle2 = plt.Circle((0, 0), RadiusInnerWall, facecolor = AnulusColor, edgecolor = None)
ax.add_patch(circle1)
ax.add_patch(circle2)

''' plot ball in inside the annulus setup '''
pos2 = (RadiusInnerWall+RadiusOfBall,0)
circle3 = plt.Circle(pos2, RadiusOfBall, facecolor = PathColor, edgecolor = None)
ax.add_patch(circle3)

''' indicate diameters of outer wall with arrows '''
ax.arrow(0, 0, 0, RadiusOuterWall, head_width=91, head_length=91, color=DiscColor, length_includes_head=True)
ax.arrow(0, 0, 0, -RadiusOuterWall, head_width=91, head_length=91, color=DiscColor, length_includes_head=True)
ax.text(0.1*RadiusOuterWall,0.95*RadiusOuterWall,'$d_\\mathcal{D}$', color = DiscColor, va = 'top')

''' indicate diameters of inner wall with arrows '''
phi=2*np.pi/3
ax.arrow(0, 0, RadiusInnerWall*np.cos(phi), RadiusInnerWall*np.sin(phi), head_width=91, head_length=91, color='w', length_includes_head=True)
ax.arrow(0, 0, -RadiusInnerWall*np.cos(phi), -RadiusInnerWall*np.sin(phi), head_width=91, head_length=91, color='w', length_includes_head=True)
ax.text(0.5*RadiusInnerWall*np.cos(phi), 0.5*RadiusInnerWall*np.sin(phi),'$d_\\mathcal{A}$', color = 'w', ha = 'right', va = 'top')

''' no axis '''
ax.set_xlim(-1.03*RadiusOuterWall,1.03*RadiusOuterWall)
ax.set_ylim(-1.03*RadiusOuterWall,1.03*RadiusOuterWall)
ax.axis('off')


'''
-----------------------------
-----------------------------
-----------------------------
----------- upper right pannel (c)
- MSD
'''
L2=0.105
B2=0.15
W2= .353
H2=.345-B
ax = fig.add_axes([L+W+L2 , 2*B2+H2, W2, H2])
ax.set_xlabel("$t$ in $s$")
ax.set_xlim(0.04, 1.4*10**3)
ax.set_ylim(0.001, 10**5)

''' read MSD data'''
DataFilePath = DataPath + 'F1c_MSD.txt'
Data = np.loadtxt( DataFilePath )
time, MSD = Data[:,0], Data[:,1]
DataFilePath = DataPath + 'F1c_MSD_BallisticTheory.txt'
Data = np.loadtxt( DataFilePath )
time, BallisticTheory = Data[:,0], Data[:,1]

''' plot MSD and t^2 scaling'''
lns1 = ax.loglog( time, MSD, '-', label = "$\\mathcal{D}$: MSD in $\\rm{cm}^2$", color = DiscColor)
ax.loglog(time, BallisticTheory, 'k--', label = '$\\propto t^2$')

''' read angular MSD data for anulus geometry '''
DataFilePath = DataPath + 'F1c_AngularMSD.txt'
Data = np.loadtxt( DataFilePath )
time = Data[0]
MSD_angular = Data[1]

''' plot angular MSD for anulus geometry & ballistic and diffusive scaling '''
lns2 = ax.loglog(time, MSD_angular, '-', label = '$\\mathcal{A}$: $\\Delta \\theta^2$ in $rad^2$', color = 'blue')

theroy = time**2
pref = MSD_angular[1]/theroy[1]
lns3 = ax.loglog(time, pref*theroy, 'k--', label = '$\\propto t^2$')
theroy = time
pref = MSD_angular[-1]/theroy[-1]
mask = time > 5
lns4 = ax.loglog(time[mask], pref*theroy[mask], 'k:', label = '$\\propto t$')

''' create legend '''
lns = lns1+lns2+lns3+lns4
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='lower right', ncol = 2, frameon = False, borderaxespad = 0)


'''
-----------------------------
-----------------------------
-----------------------------
----------- lower right pannel (e)
- probability of speed
'''
ax = fig.add_axes([L+W+L2 , B2, W2, H2])
ax.set_xlabel("$|\\vec{v}|$ in $cm/s$")
ax.set_ylabel("$P\\left(\\vert\\vec{v}\\vert\\right)$")

DataFilePath = DataPath + 'F1e_SpeedDistributionDisc.txt'
Data = np.loadtxt( DataFilePath )
BinCenters_Disc, SpeedHistogram_Disc = Data[:,0], Data[:,1]
DataFilePath = DataPath + 'F1e_SpeedDistributionAnnulus.txt'
Data = np.loadtxt( DataFilePath )
BinCenters_Annulus, SpeedHistogram_Annulus = Data[:,0], Data[:,1]

''' plot speed distributions '''
ax.plot(BinCenters_Disc, SpeedHistogram_Disc, 'o-', color = DiscColor, label = '$\\mathcal{D}$')
ax.plot(BinCenters_Annulus, SpeedHistogram_Annulus, 'o-', color = AnulusColor, markeredgecolor = AnulusColor, markerfacecolor = 'w', label = '$\\mathcal{A}$')

ax.set_xlim(0, 62)
ax.legend(loc = 'upper right', frameon=False)


'''
-----------------------------
-----------------------------
-----------------------------
----------- left most pannel (a)
- schematic of the ball
'''
ax = fig.add_axes([-0.66, -0.29, 1.6, 1.6], projection='3d')

ax.set_box_aspect((2., 1.5, 1.3))
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

''' define vertical and horizontal rotation of the ball '''
elev = 8.0 # horizontal rotation
rot = 60.0 / 180 * np.pi #vertical rotation

''' plot turquoise ball '''
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = 1 * np.outer(np.cos(u), np.sin(v))
y = 1 * np.outer(np.sin(u), np.sin(v))
z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color=PathColor, linewidth=0, alpha=0.3)

''' plot the stator elements inside the ball '''
ax.plot( 0, 0, 0, '*', markersize = 5, color = 'r', zorder = 3) # center of stator : red star
ax.text( 0., 0.04, 0, 'center', color = 'r', va = 'center') # the red word 'center'
ax.plot(np.sin(np.pi-rot), np.cos(np.pi-rot), 0, 'x', markersize = 3, color = 'k', zorder = 11) # fixation of stator to shell : balck x in the front
ax.plot(np.sin(-rot),np.cos(-rot), 0, 'x', markersize = 3, color = 'k' ) # fixation of stator to shell : balck x in the back
ax.plot([np.sin(np.pi-rot),np.sin(-rot)],[np.cos(np.pi-rot),np.cos(-rot)],[0,0], color = 'k') # the stator axis : black line

''' plot rotor '''
offCenterX, offCenterY = np.sin(np.pi-rot)*0.5, np.cos(np.pi-rot)*0.5
ax.plot_surface(0.06*x+offCenterX, 0.06*y+offCenterY, 0.06*z,  rstride=4, cstride=4, color=colorMotor, linewidth=0, alpha=1) # rotor : green ball
ax.text( offCenterX, offCenterY+0.05, -0.04, 'rotor', color = colorMotor, va = 'center') # the green word 'rotor'

''' plot load '''
ax.plot_surface(-0.2+0.1*x+offCenterX, -0.2+0.1*y+offCenterY, 0.2+0.1*z,  rstride=4, cstride=4, color=colorWeight, linewidth=0, alpha=1) # load : gray ball
ax.text( -0.2+offCenterX, -0.2+offCenterY+0.05, 0.35, 'load', color = colorWeight, va = 'center', ha = 'right') # the gray word 'load'
ax.plot([offCenterX,-0.2+offCenterX],[offCenterY,-0.2+offCenterY],[0,0.2], color = colorWeight) # connection load -> rotor : gray line

''' indicate trajectory of load (sens of rotation : two gray arrows) '''
a = Arrow3D([offCenterX, 0.2+offCenterX], [offCenterY, 0.2+offCenterY], [0.35, 0.35], mutation_scale=10, arrowstyle="->", color=colorWeight)
ax.add_artist(a)
a = Arrow3D([offCenterX, -0.2+offCenterX], [offCenterY, -0.2+offCenterY], [-0.36, -0.35], mutation_scale=10, arrowstyle="->", color=colorWeight)
ax.add_artist(a)

''' plot as additional guidance for the horizontal belt for ball '''
ax.plot(np.sin(u),np.cos(u),0,color=PathColor, linestyle = 'dashed') # dashed full horizontal belt for shell
horiz_front = np.linspace(0, np.pi, 100)
ax.plot(np.sin(horiz_front),np.cos(horiz_front),0,color=PathColor, zorder = 10)# solid front horizontal belt for shell

''' plot as additional guidance for the vertical belt for ball '''
a = np.array([-np.sin(elev / 180 * np.pi), 0, np.cos(elev / 180 * np.pi)])
b = np.array([0, 1, 0])
b = b * np.cos(rot) + np.cross(a, b) * np.sin(rot) + a * np.dot(a, b) * (1 - np.cos(rot))
vert_front = np.linspace(np.pi / 2, 3 * np.pi / 2, 100)
ax.plot(a[0] * np.sin(u) + b[0] * np.cos(u), b[1] * np.cos(u), a[2] * np.sin(u) + b[2] * np.cos(u),color=PathColor, linestyle = 'dashed')# vertical belt back
ax.plot(a[0] * np.sin(vert_front) + b[0] * np.cos(vert_front), b[1] * np.cos(vert_front), a[2] * np.sin(vert_front) + b[2] * np.cos(vert_front),color=PathColor, zorder = 10)# vertical belt front

''' indicate trajectory of load (trajectory : gray dashed cicle) '''
b2 = np.array([0, 1, 0])
b2 = b2 * np.cos(rot+np.pi/2) + np.cross(a, b2) * np.sin(rot+np.pi/2) + a * np.dot(a, b2) * (1 - np.cos(rot+np.pi/2))
ax.plot(0.35*a[0] * np.sin(u) + b2[0] * np.cos(u) + offCenterX, 0.35*b2[1] * np.cos(u) + offCenterY, 0.35*a[2] * np.sin(u) + b2[2] * np.cos(u),color=colorWeight, linestyle = '--', linewidth = 0.7) # trajectory of load : dashed gray line

''' adjust view and remove lables '''
ax.view_init(elev = elev, azim = 0)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

''' include png of more detailed schematic showing battery, sensor etc. '''
L=0.0
B=0.01
W=.17-L
H = W / golden_mean
ax = fig.add_axes([L+W, B, W, H])
ImPath = script_dir + '/OpenBall.png'
img = Image.open(ImPath)
ax.imshow(img, )
ax.axis('off')


''' safe Figure '''
FigPath = script_dir + '/F1.png'
plt.savefig( FigPath, dpi = 300)
print('open', FigPath)

