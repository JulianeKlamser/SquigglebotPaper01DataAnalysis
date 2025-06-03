import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from scipy.optimize import curve_fit
import pickle
import numpy as np

def Gaussian( x, sigma, mu ):
    return 1/np.sqrt( sigma**2 * 2 * np.pi ) * np.exp( -(x-mu)**2 / ( 2 * sigma**2) )

def GenerateHistogram(MIN , MAX, nBins, data):
    binWidth = ( MAX - MIN  ) / nBins
    Histogram = np.zeros(nBins)
    
    for dat in data:
        bin = int( (dat - MIN ) / binWidth )
        if bin == nBins:
            bin -= 1
        Histogram[bin] += 1
    Histogram /= sum(Histogram) * binWidth
    assert( (sum(Histogram) * binWidth - 1) < 10**(-5) )
    
    BinCenters = MIN  + (np.arange(nBins) + 0.5) * binWidth
    return BinCenters, Histogram
    
def GenerateCenteredHistogram(nBins, data):
    if nBins % 2 != 0:
        nBins += 1
    AbsMAX = max( max(data), -min(data)  )
    MIN   = -AbsMAX
    binWidth = 2*AbsMAX / nBins
    Histogram = np.zeros(nBins)
    
    for dat in data:
        bin = int( (dat - MIN  ) / binWidth )
        if bin == nBins:
            bin -= 1
        Histogram[bin] += 1
    Histogram /= sum(Histogram) * binWidth
    assert( (sum(Histogram) * binWidth - 1) < 10**(-5) )
    
    BinCenters = MIN   + (np.arange(nBins) + 0.5) * binWidth
    return BinCenters, Histogram

def VisualizeDataAnalysis(SaveFigAs, T, current, rawCurrentBins, P_of_rawCurrent, P_of_avCurrents, P_of_mean0_avCurrents, Ratio_ProbPosC_PropNegC_perT):
    
    NumOfT = len(T)
    color=cm.nipy_spectral(np.linspace(0,0.95,NumOfT))

    fig, axs = plt.subplots(3, 2, dpi = 100, figsize = (10,8.5))
    ax11, ax12 = axs[0]
    ax21, ax22 = axs[1]
    ax31, ax32 = axs[2]
    
    ''' show time evolution of the raw current
        first 150 and last 150 values '''
    ax11.plot(current[:150], '.-', linewidth = 0.5, label = 'first 150')
    ax11.plot(current[-150:], '.-', linewidth = 0.5, label = 'last 150')
    ax11.set_xlabel('time')
    ax11.set_ylabel('current')
    ax11.set_title('raw referece current')
    ax11.set_xlim(0,150)
    ax11.legend()

    ''' show probability distribution of the 
        raw reference current '''
    ax12.step( rawCurrentBins, P_of_rawCurrent, where='mid')
    ax12.set_xlabel('current $\\Delta j_i$')
    ax12.set_ylabel('probability')
    
    ''' show probability distribution of the
        sliding averaged currents per T '''
    for i, dat in enumerate(P_of_avCurrents):
        avCurrentBins, P_of_avCurrent = dat
        
        ax21.semilogy( avCurrentBins, P_of_avCurrent, '.', color = color[i], label = '$T = %d$' % T[i])
        
        popt, pcov = curve_fit( Gaussian, avCurrentBins, P_of_avCurrent)
        Gaussian_fit = Gaussian(avCurrentBins, *popt)
        ax21.semilogy(avCurrentBins, Gaussian_fit, '-' , color = color[i], alpha = 0.5)
    ax21.set_title('solid lines - gaussian fit')
    ax21.set_ylim(10**(-1),400)
    ax21.set_xlabel('sliding averaged current $j_T$')
    ax21.set_ylabel('probability $P(j_T)$')
    ax21.legend()

    ''' show probability distribution of the mean-substracted
        sliding averad currents per T '''
    for i, dat in enumerate(P_of_mean0_avCurrents):
        mean0_avCurrentBins, P_of_mean0avC = dat
        
        ax22.semilogy(mean0_avCurrentBins, P_of_mean0avC, '.', color = color[i], label = '$T = %d$' % T[i] )
        
        popt, pcov = curve_fit( Gaussian, mean0_avCurrentBins, P_of_mean0avC)
        Gaussian_fit = Gaussian(mean0_avCurrentBins, *popt)
        ax22.semilogy(mean0_avCurrentBins, Gaussian_fit, '-' , color = color[i], alpha = 0.5)
        
        #    https://en.wikipedia.org/wiki/Coefficient_of_determination
        SS_res = sum( (P_of_mean0avC - Gaussian_fit)**2 )
        mean = np.mean(P_of_mean0avC)
        SS_tot = sum( (P_of_mean0avC - mean)**2 )
        rSq = 1 - SS_res/SS_tot
        
    yLim = [0.01,400]
    MarkerC = 0.007303
    ax22.semilogy( [ -MarkerC, -MarkerC], yLim, 'k-')
    ax22.semilogy( [ MarkerC, MarkerC], yLim, 'k-', label = '$x = %.3e$' % MarkerC)
    ax22.set_title('solid lines - gaussian fit')
    ax22.set_ylim(yLim)
    ax22.set_xlim(-0.015,0.015)
    ax22.set_xlabel('$j_T - \\langle j_T \\rangle$')
    ax22.set_ylabel('probability $P(j_T - \\langle j_T \\rangle)$')
    ax22.legend(frameon = False, handlelength = 0.5)

    ''' show log ( P(j)/P(-j) ) vs j '''
    for i, dat in enumerate(Ratio_ProbPosC_PropNegC_perT):
        current, Ratio_ProbPosC_PropNegC = dat
        ax31.plot(current, np.log(Ratio_ProbPosC_PropNegC), 'o-', color = color[i], label = '$T = %d$' % T[i])
    ax31.grid(True)
    ax31.set_xlim(0)
    ax31.set_xlabel('$|j_T - \\langle j_T \\rangle|$')
    ax31.set_ylabel('$\\log\\left[ P(|j_T - \\langle j_T \\rangle|) / P(-|j_T - \\langle j_T \\rangle| )\\right]$')
    ax31.legend()
    
    ''' show 1/T log ( P(j)/P(-j) ) vs j '''
    for i, dat in enumerate(Ratio_ProbPosC_PropNegC_perT):
        current, Ratio_ProbPosC_PropNegC = dat
        ax32.plot(current, 1/T[i] * np.log(Ratio_ProbPosC_PropNegC), 'o-', color = color[i], label = '$T = %d$' % T[i])
    yLim = [-0.1,0.1]
    ax32.plot( [ MarkerC, MarkerC], yLim, 'k-', label = '$x = %.3e$' % MarkerC)
    ax32.grid(True)
    ax32.set_xlim(0)
    ax32.set_ylim(-0.1, 0.1)
    ax32.set_xlabel('$|j_T - \\langle j_T \\rangle|$')
    ax32.set_ylabel('$\\frac{1}{T}\\log\\left[ P(|j_T - \\langle j_T \\rangle|) / P(-|j_T - \\langle j_T \\rangle| )\\right]$')
    ax32.legend(ncol = 2, frameon = False, handlelength = 0.5)

    plt.tight_layout()
    plt.savefig(SaveFigAs, dpi = 300)
    print('open', SaveFigAs)


def VisualizeDataAnalysis_nonRef(SaveFigAs, case, T, current, rawCurrentBins, P_of_rawCurrent, P_of_CorrAvCurrents, Ratio_ProbPosC_PropNegC_perT, exludeFromFit, Slope, rSqList):
    NumOfT = len(T)
    color=cm.nipy_spectral(np.linspace(0,0.95,NumOfT))

    fig, axs = plt.subplots(3, 2, dpi = 100, figsize = (10,8.5))
    ax11, ax12 = axs[0]
    ax21, ax22 = axs[1]
    ax31, ax32 = axs[2]
    
    ''' show time evolution of the raw current
        first 150 and last 150 values '''
    ax11.plot(current[:150], '.-', label = 'first 150')
    ax11.plot(current[-150:], '.-', label = 'last 150')
    ax11.set_xlabel('time')
    ax11.set_ylabel('current')
    ax11.set_title(case)
    ax11.set_xlim(0)
    ax11.legend()
    
    
    ''' show probability distribution of the 
        raw reference current '''
    ax12.step( rawCurrentBins, P_of_rawCurrent, where='mid')
    ax12.set_xlabel('current $\\Delta j_i$')
    ax12.set_ylabel('probability')

    ''' '''
    for i, dat in enumerate(P_of_CorrAvCurrents):
        refCorr_avCurrentBins, P_of_corrAvC = dat
        ax21.semilogy( refCorr_avCurrentBins, P_of_corrAvC, '.', color = color[i], label = '$T = %d$' % T[i])
        
        popt, pcov = curve_fit( Gaussian, refCorr_avCurrentBins, P_of_corrAvC)
        Gaussian_fit = Gaussian(refCorr_avCurrentBins, *popt)
        ax21.semilogy(refCorr_avCurrentBins, Gaussian_fit, '-' , color = color[i], alpha = 0.5)
    ax21.set_title('solid lines - gaussian fit')
    ax21.set_ylim(10**(-1),400)
    ax21.set_xlabel('$\\Delta j_T$')
    ax21.set_ylabel('probability $P(\\Delta j_T)$')
    ax21.legend(frameon = False, handlelength = 0.5)

    ''' '''
    for i, dat in enumerate(Ratio_ProbPosC_PropNegC_perT):
        PositiveCurrent, Ratio_ProbPosC_PropNegC = dat
        
        x, y = PositiveCurrent, np.log(Ratio_ProbPosC_PropNegC)
        ax22.plot(x, y, '.:', color = color[i])
        ex = exludeFromFit[i]
        if ex > 0:
            x, y = x[:-ex], y[:-ex]
        ax22.plot(x, y, '-', color = color[i], label = '$T = %d$' % T[i])
    ax22.set_xlim(0)
    ax22.set_xlabel('$|\\Delta j_T|$')
    ax22.set_ylabel('$\\log\\left[ P(|\\Delta j_T|) / P(-|\\Delta j_T|)\\right]$')
    ax22.legend()
    
    
    ''' '''
    MAX_i = 0
    for i, dat in enumerate(Ratio_ProbPosC_PropNegC_perT):
        PositiveCurrent, Ratio_ProbPosC_PropNegC = dat
        
        x, y = PositiveCurrent, 1/T[i] * np.log(Ratio_ProbPosC_PropNegC)
        MAX_i = max( MAX_i, y.max())
        ax31.plot(PositiveCurrent, y, '.:', color = color[i])
        ex = exludeFromFit[i]
        if ex > 0:
            x, y = x[:-ex], y[:-ex]
        ax31.plot( x, y, '-', color = color[i], label = '$T = %d$' % T[i])
        
    ax31.set_xlim(0)
    ax31.set_ylim(0, MAX_i)
    ax31.set_xlabel('$|\\Delta j_T|$')
    ax31.set_ylabel('$\\frac{1}{T}\\log\\left[ \\frac{P(|\\Delta j_T|) }{P(-|\\Delta j_T|)}\\right]$')
    ax31.legend()
    
    ''' '''
    ax32Right = ax32.twinx()
    
    ax32.plot( T, Slope, 'o-', label = '$dy / dx$' )
    ax32.set_ylim(0, 30)
    ax32.set_xlabel( '$T$' )
    ax32.set_ylabel( '$dy/dx$' )
    ax32.legend(loc = 'upper left')
    
    ax32Right.plot( T, rSqList, 'rs-', label = '$R^2$')
    ax32Right.set_ylabel( '$R^2$' )
    ax32Right.set_ylim( 0.97, 1 )
    ax32Right.legend( loc = 'upper right')
    
    
    plt.tight_layout()
    plt.savefig(SaveFigAs, dpi = 300)
    print('open', SaveFigAs)
