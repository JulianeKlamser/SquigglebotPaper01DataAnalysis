import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
script_dir = os.path.dirname(os.path.abspath(__file__))

fig_width_pt = 246.0        # Get this from LaTeX using \the\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = 0.8     # Aesthetic ratio
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
          'legend.borderaxespad' : 0.05,
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

''' read data '''
BulkDef = 0.90
DataPath = script_dir + '/Data/Out_Forward%.2f.pkl' % BulkDef
with open(DataPath, 'rb') as f:
    DatF = pickle.load(f)
    
DataPath = script_dir + '/Data/Out_Backward%.2f.pkl' % BulkDef
with open(DataPath, 'rb') as f:
    DatB = pickle.load(f)


fig = plt.figure()
L = 0.09
B = 0.106
W = 0.99*(0.492 - L/2)
H = 0.99*(0.5 - B/2)


''' lower left pannel of Fig 4 '''
ax = fig.add_axes([L, B, W, H])
ax.plot( DatF['Angle'], DatF['P_of_Angle_wall'], 'o-', label = '$P(\\varphi_i)$', color = 'lightseagreen')
ax.plot( DatB['Angle'], DatB['P_of_Angle_wall'], 's-', label = '$P(\\varphi_o)$', color = 'slateblue')
ax.set_xlabel('$\\varphi$')
ax.set_ylim(0,1.8)
ax.set_xlim(0,np.pi)
ax.set_xticks([i*np.pi/4 for i in range(5)])
ax.set_xticklabels(['$0$', '$\\pi/4$', '$\\pi/2$', '$3\\pi/4$', '$\\pi$'])
ax.legend(frameon = False)
mask = (DatF['P_of_Angle_wall'] > 0) & (DatB['P_of_Angle_wall'] > 0)
KLdWall = np.sum( DatB['P_of_Angle_wall'][mask] * np.log(DatB['P_of_Angle_wall'][mask]/DatF['P_of_Angle_wall'][mask]))
ax.text(np.pi/2, 1., '$KLd(\\mathcal{W}) = %.3f$' % (KLdWall), ha = 'center')

''' upper left pannel of Fig 4 '''
ax = fig.add_axes([L, B + H, W, H])
ax.plot( DatF['Angle'], DatF['P_of_Angle_bulk'], 'o-', label = '$P(\\varphi_i)$', color = 'lightseagreen')
ax.plot( DatB['Angle'], DatB['P_of_Angle_bulk'], 's-', label = '$P(\\varphi_o)$', color = 'slateblue')
ax.set_ylim(0,1.1)
ax.set_xlim(0,np.pi)
ax.set_xticks([])
ax.set_xticklabels([])
ax.legend(frameon = False)
mask = (DatF['P_of_Angle_bulk'] > 0) & (DatB['P_of_Angle_bulk'] > 0)
KLdBulk = np.sum( DatB['P_of_Angle_bulk'][mask] * np.log(DatB['P_of_Angle_bulk'][mask]/DatF['P_of_Angle_bulk'][mask]))
ax.text(np.pi/2, 0.7, '$KLd(\\mathcal{B}) = %.3f$' % (KLdBulk), ha = 'center')

''' upper right pannel of Fig 4 '''
ax = fig.add_axes([1.15*L + W, B + H, W, H])
ax.plot( DatF['Angle'], DatF['P_of_Angle_mix'], 'o-', label = '$P(\\varphi_i)$', color = 'lightseagreen')
ax.plot( DatB['Angle'], DatB['P_of_Angle_mix'], 's-', label = '$P(\\varphi_o)$', color = 'slateblue')
ax.set_xlabel('$\\varphi$')
ax.set_ylim(0,1.1)
ax.set_xlim(0,np.pi)
ax.set_xticks([i*np.pi/4 for i in range(5)])
ax.set_xticklabels(['$0$', '$\\pi/4$', '$\\pi/2$', '$3\\pi/4$', '$\\pi$'])
ax.set_yticks([0.2*i for i in range(6)])
ax.set_yticklabels([])
ax.legend(frameon = False)
mask = (DatF['P_of_Angle_mix'] > 0) & (DatB['P_of_Angle_mix'] > 0)
KLdMix = np.sum( DatB['P_of_Angle_mix'][mask] * np.log(DatB['P_of_Angle_mix'][mask]/DatF['P_of_Angle_mix'][mask]))
ax.text(np.pi/2, 0.7, '$KLd(\\mathcal{B}\\text{\\&}\\mathcal{W}) = %.3f$' % (KLdMix), ha = 'center')


''' lower right pannel of Fig 4 '''
B = 0.092
W = 0.36
H = 0.35
L = (1-W)
ax = fig.add_axes([L, B, 0.99*W, H])
ax.semilogy( DatF['current'], DatF['P_of_current_wall'], 's-', label = '$\\mathcal{W}$', color = 'k')
ax.semilogy( DatF['current'], DatF['P_of_current_bulk'], 'o-', label = '$\\mathcal{B}$', color = 'gray')
ax.semilogy( DatF['current'], DatF['P_of_current_mix'], 'P-', label = '$\\mathcal{B}\\text{\\&}\\mathcal{W}$', color = 'cadetblue')
ax.legend(frameon = False)
ax.set_xlabel('$J$')
ax.set_ylabel('$P(J)$', labelpad = -1)
ax.set_ylim(1,100)
ax.set_xlim(0.09,0.14)
ax.set_xticks([0.09, 0.1, 0.11, 0.12, 0.13, 0.14])
ax.set_xticklabels(['0.09', '', '0.11', '', '0.13', ''])

L = 0.0
B = 0.0
W = 1
H = 1
ax = fig.add_axes([L, B, W, H])
ax.axis('off')
ax.text(0.1, 0.95, '(a)')
ax.text(0.6, 0.95, '(b)')
ax.text(0.1, 0.48, '(c)')
ax.text(0.65, 0.4, '(d)')

FName = script_dir + '/InOutAngle_new_logy_bulkDef%.2f.png' % BulkDef
plt.savefig(FName, dpi = 300)
print('open ', FName)
