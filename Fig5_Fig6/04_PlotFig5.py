#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from matplotlib.pyplot import cm
from scipy.optimize import curve_fit
script_dir = os.path.dirname(os.path.abspath(__file__))

def Gaussian( x, sigma, mu ):
    return 1/np.sqrt(sigma**2 * 2 * np.pi) * np.exp(-0.5* ((x-mu) / sigma )**2 )


fig_width_pt = 246.0        # Get this from LaTeX using \the\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = 0.9      # Aesthetic ratio
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
          'legend.columnspacing' : 0.6,
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

refDataFile = script_dir + '/Data/Reference_Analysed.p'
bulkDataFile = script_dir + '/Data/bulk_Analysed.p'
boundaryDataFile = script_dir + '/Data/boundary_Analysed.p'

refDict = pickle.load( open( refDataFile, "rb" ) )
bulkDict = pickle.load( open( bulkDataFile, "rb" ) )
boundaryDict = pickle.load( open( boundaryDataFile, "rb" ) )


colorAnu = cm.winter(np.linspace(1,0,8))
fig = plt.figure()
L=0.115
B=0.093
W=.495-L
H=.497-B

#lower left pannel of Fig 5
ax = fig.add_axes([0.75*L, B, W, H])
ax.set_xlabel('$|\\Delta j_T|$')
ax.set_ylabel('$\\log\\left[ P(|\\Delta j_T|) / P(-|\\Delta j_T|)\\right]$', labelpad = 2)
for i, T in enumerate(boundaryDict['T']):
    if T in [20, 40, 60, 80]:
        x, y = boundaryDict['Ratio_ProbPosC_PropNegC_perT'][i]
        y = np.log(y)
        ax.plot(x, y, 'x:', color = colorAnu[i] , linewidth = 0.5)
        ex = boundaryDict['exludeFromFit'][i]
        if ex > 0:
            x, y = x[:-ex], y[:-ex]
        ax.plot( x, y, 'o-', color = colorAnu[i], label = "$T = %d$" % T)
ax.legend(loc = 'upper left', frameon=False, bbox_to_anchor=(-0.05,1.08))
ax.set_xlim(0,0.0151)
ax.set_ylim(0,6)



#lower right pannel of Fig 5
ax = fig.add_axes([2.05*L+W, B, W, H])
ax.set_xlabel('$|\\Delta j_T|$')
ax.set_ylabel('$T^{-1} \\,\\log\\left[ P(|\\Delta j_T|) / P(-|\\Delta j_T|)\\right]$', labelpad = 1.7)
for i, T in enumerate(boundaryDict['T']):
    if T in [20, 40, 60, 80]:
        x, y = boundaryDict['Ratio_ProbPosC_PropNegC_perT'][i]
        y = 1/T * np.log(y)
        ex = boundaryDict['exludeFromFit'][i]
        if ex > 0:
            x, y = x[:-ex], y[:-ex]
        ax.plot( x, y, 'o-', color = colorAnu[i], label = "$T = %d$" % T)
ax.legend(loc = 'upper left', frameon=False, bbox_to_anchor=(-0.05,1.08))
ax.set_xlim(0,0.011)
ax.set_ylim(0,0.1)
ax.set_xticks([0, 0.005, 0.01])
ax.set_xticklabels(['0.0', '0.005', '0.01'])
ax.set_yticks([0, 0.05, 0.1])
ax.set_yticklabels(['0.0', '0.05', '0.1'])


# upper left pannel of Fig 5
ax = fig.add_axes([0.91*L, 2*B+H, W, H])
ax.set_xlabel('$j_T$')
ax.set_ylabel('$P\\left(j_T\\right)$')
NumOfT = len(refDict['T'])
colors=cm.gray(np.linspace(0.65,0,NumOfT))
for i, T in enumerate(refDict['T']):
    if T in [75, 175]:
        avCurrentBins, P_of_avC = refDict['P_of_avCurrents'][i]
        ax.semilogy(avCurrentBins, P_of_avC, '.-', color = colors[i], label = '$T = %d$' % T)

ax.legend(loc = 'upper right', frameon=False, bbox_to_anchor=(1.07,1.07))
popt, pcov = curve_fit( Gaussian, avCurrentBins, P_of_avC)
gaussianFit = Gaussian(avCurrentBins, *popt)
ax.plot(avCurrentBins, gaussianFit, '-', color = 'r', linewidth = 0.5 )
Mean = refDict['MeanCurrent_perT'][-1]
ax.plot( [Mean,Mean] , [0,400] , ':', color = colors[-1])
ax.set_ylim(0.6,400)
ax.set_xlim(0.04,0.062)


# upper right pannel of Fig 5
ax = fig.add_axes([2.05*L+W, 2*B+H, W, H])
ax.set_xlabel('$\\Delta j_T$')
ax.set_ylabel('$P\\left(\\Delta j_T\\right)$', labelpad = -2)
for i, T in enumerate(boundaryDict['T']):
    if T in [20, 80]:
        refCorr_avCurrentBins, P_of_corrAvC = boundaryDict['P_of_CorrAvCurrents'][i]
        ax.semilogy(refCorr_avCurrentBins, P_of_corrAvC, '.-', color = colorAnu[i], label = "$T = %d$" % T)
popt, pcov = curve_fit( Gaussian, refCorr_avCurrentBins, P_of_corrAvC)
gaussianFit = Gaussian(refCorr_avCurrentBins, *popt)
ax.plot(refCorr_avCurrentBins, gaussianFit, '-', color = 'r', linewidth = 0.5 )
ax.legend(loc = 'upper right', frameon=False, bbox_to_anchor=(1.06,1.07))
ax.set_xlim(-0.015,0.03)
ax.set_ylim(0.1,400)


# add the frame lables
ax = fig.add_axes([0,0,1,1])
ax.axis('off')
ax.text(0.11,0.95, '(a)')
ax.text(0.62,0.95, '(b)')
ax.text(0.42,0.12, '(c)')
ax.text(0.95,0.12, '(d)')

FigName = script_dir + '/Fig5.png'
plt.savefig(FigName, facecolor=fig.get_facecolor(), dpi=300)
print('open', FigName)
