#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib.pyplot import cm
from scipy.optimize import curve_fit
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

def Gaussian( x, sigma, mu ):
    return 1/np.sqrt(sigma**2 * 2 * np.pi) * np.exp(-0.5* ((x-mu) / sigma )**2 )

fig_width_pt = 246.0        # Get this from LaTeX using \the\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = 0.4      # Aesthetic ratio
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
          'legend.columnspacing' : -0.5,
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

fig = plt.figure()
L=0.11
B=0.2
W=.495-L
H=.95-B

''' left pannel of Fig 6 '''
ax = fig.add_axes([L, B, 0.95*W, H])
ax.set_xlabel('$\\Delta j_T$')
ax.set_ylabel('$P\\left(\\Delta j_T\\right)$')

# reference current distribution
assert (refDict['T'][-1] == 175)
x, y = refDict['P_of_mean0_avCurrents'][-1]
ax.semilogy(x, y, '.', color = 'k', label = "$\\mathcal{R}$")
popt, pcov = curve_fit( Gaussian, x, y)
gaussianFit = Gaussian(x, *popt)
ax.plot(x, gaussianFit, '-', color = 'k', linewidth = 0.5 )

# annulus arena current distribution
winterColor = cm.winter(np.linspace(1,0, 2))[1]
i = np.where(boundaryDict['T'] == 80)[0][0]
x, y = boundaryDict['P_of_CorrAvCurrents'][i]
ax.semilogy(x, y, '.', color = winterColor, label = "$\\mathcal{A}$")
popt, pcov = curve_fit( Gaussian, x, y)
gaussianFit = Gaussian(x, *popt)
ax.plot(x, gaussianFit, '-', color = winterColor, linewidth = 0.5 )

# disk arena current distribution
autumnColor = cm.autumn(np.linspace(1,0, 2))[1]
i = np.where(bulkDict['T'] == 80)[0][0]
x, y = bulkDict['P_of_CorrAvCurrents'][i]
ax.semilogy(x, y, '.', color = autumnColor, label = "$\\mathcal{D}$")
popt, pcov = curve_fit( Gaussian, x, y)
gaussianFit = Gaussian(x, *popt)
ax.plot(x, gaussianFit, '-', color = autumnColor, linewidth = 0.5 )

# limits and legend
ax.legend(loc = 'upper right', frameon=False, bbox_to_anchor=(1.05,1.08))
ax.set_xlim(-0.009,0.022)
ax.set_ylim(0.6, 400)



''' right pannel of Fig 6 '''
xLim = np.array([ 0, 0.012])
ax = fig.add_axes([2.05*L+W, B, W, H])
ax.set_xlabel('$\\Delta j_T$')
ax.set_ylabel('$\\frac{1}{T}\\log\\left[ P(\\Delta j_T) / P(-\\Delta j_T)\\right]$')

PlotMarkers = ["o", "d", "v", "P"]
markerSize = [ 2., 2., 2., 2.]

# annulus arena
winterColor = cm.winter(np.linspace(1,0, len(boundaryDict['T'])))
k = 0
for i, T in enumerate(boundaryDict['T']):
    if T in [20, 40, 60, 80]:
        x, y = boundaryDict['Ratio_ProbPosC_PropNegC_perT'][i]
        y = 1/T * np.log(y)
        ex = boundaryDict['exludeFromFit'][i]
        if ex > 0:
            x, y = x[:-ex], y[:-ex]
        ax.plot( x, y, PlotMarkers[k], color = winterColor[i], label = "$\\!$")
        k += 1
i = np.where(boundaryDict['T'] == 80)[0][0]
slope = boundaryDict['Slope'][i]
ax.plot( xLim, slope*xLim, '-', color = winterColor[-1], zorder = 0)

# disk arena
autumnColor = cm.autumn(np.linspace(1,0, len(bulkDict['T'])))
k = 0
for i, T in enumerate(bulkDict['T']):
    if T in [20, 40, 60, 80]:
        x, y = bulkDict['Ratio_ProbPosC_PropNegC_perT'][i]
        y = 1/T * np.log(y)
        ex = bulkDict['exludeFromFit'][i]
        if ex > 0:
            x, y = x[:-ex], y[:-ex]
        ax.plot( x, y, PlotMarkers[k], color = autumnColor[i], label = "$\\!$" % T)
        k += 1
i = np.where(bulkDict['T'] == 80)[0][0]
slope = bulkDict['Slope'][i]
ax.plot( xLim, slope*xLim, '-', color = autumnColor[-1], zorder = 0)

# reference
grayColor = cm.binary(np.linspace(0.1,1, len(refDict['T'])))
k = 0
for i, T in enumerate(refDict['T']):
    if T in [100, 125, 150, 175]:
        x, y = refDict['Ratio_ProbPosC_PropNegC_perT'][i]
        y = 1/T * np.log(y)
        ax.plot( x, y, PlotMarkers[k], color = grayColor[i], label = "$T_%d$" % (k+1))
        k += 1
ax.plot([0,1], [0,0], 'k-', zorder = 0)

ax.legend(loc = 'upper left', frameon=False, bbox_to_anchor=(-0.05,1.08), ncol = 3)
ax.set_xlim(0,0.012)
ax.set_ylim(-0.05,0.3)

''' pannel labels '''
ax = fig.add_axes([0,0,1,1])
ax.axis('off')
ax.text(0.43,0.22, '(a)')
ax.text(0.95,0.22, '(b)')

FigFile = script_dir + '/Fig6.png'
plt.savefig(FigFile, facecolor=fig.get_facecolor(), dpi=300)
print('open', FigFile)
