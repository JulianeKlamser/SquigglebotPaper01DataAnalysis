import os, sys
import pickle
import numpy as np
from scipy.optimize import curve_fit

from Functions import VisualizeDataAnalysis_nonRef
from Functions import GenerateHistogram
from Functions import GenerateCenteredHistogram

script_dir = os.path.dirname(os.path.abspath(__file__))

DataPath = script_dir + '/Data/'
PlotPath = script_dir + '/Plots/'


def LinearFunction(x, slope):
    return slope * x

HowToUse = "\nto analyze boundary current, launch script as\npython3 03_EvaBCurrent.py boundary\n\nto analyze bulk current, launch  script as\npython3 03_EvaBCurrent.py bulk\n\n"
if len(sys.argv) != 2:
    print(HowToUse)
    exit()
if sys.argv[1] == 'boundary':
    case = sys.argv[1]
elif sys.argv[1] == 'bulk':
    case = sys.argv[1]
else:
    print(HowToUse)
    exit()
Prefix = 'TimeForward'


''' read reference current data '''
refFilePath = DataPath + 'Reference_Analysed.p'
print('load ', refFilePath)
RefDict = pickle.load( open( refFilePath, "rb" ) )
ReferenceAverageCurrent = RefDict['MeanCurrent_perT'][-1] # mean of the largest T

''' read time range of stable current '''
Limits = np.loadtxt(DataPath + 'StableLimits.txt', dtype=np.int32)

''' read bulk current '''
SaveFigAs = PlotPath + 'Analysis_%s.png' % case
FileName = 'Current_%s.txt' % case
filePath = DataPath + FileName
current = np.loadtxt(filePath)
current = current[Limits[0]:Limits[1]]
MinCurrent, MaxCurrent = current.min(), current.max()

''' compute probability distribution of the raw
    current '''
secondSmallestMaxCurrent = min( current[current > MinCurrent] )
binWidth = (secondSmallestMaxCurrent - MinCurrent)
nBins = np.round( (MaxCurrent - MinCurrent) / binWidth )
nBins = int(nBins) + 1
rawCurrentBins, P_of_rawCurrent = GenerateHistogram(MinCurrent - 0.5*binWidth, MaxCurrent + 0.5*binWidth, nBins, current)

''' define bins '''

if case == 'bulk':
    BinList = [48, 54, 54, 54, 60, 64, 66, 66]#54
    exludeFromFit = [2, 6, 5, 0, 2, 0, 4, 5]
elif case == 'boundary':
    BinList = [54, 72, 70, 80, 90, 92, 90, 102, 114]
    exludeFromFit = [9, 9, 6, 7, 5, 3, 4, 5, 7]
else:
    exit()
 
NumOfT = len(BinList)
T = np.array([10 + 10*(i) for i in range(NumOfT)])

''' Compute sliding averages of the current over different time windows T
    For each value of T:
      - Compute the sliding (moving) average of the current with window size T_i
      - Store the result
      - Compute the mean of that sliding average and store it
      - Also store a mean-subtracted version of the sliding average (centered around 0) 
      
    finally safe data'''
avCurrent_perT = [] # sliding averages of current for each T
refCorr_avCurrents_perT = [] # store zero-mean (mean-subtracted) sliding averages
P_of_avCurrents = []
P_of_CorrAvCurrents = []
for i, T_i in enumerate(T):
    # Compute sliding average with window size T_i
    AveragedCurrent_i = np.convolve(current, np.ones(T_i)/T_i, mode='valid')
    avCurrent_perT += [ AveragedCurrent_i ]
        
    # compute probability distrubution
    MinCurrent, MaxCurrent = min(AveragedCurrent_i), max(AveragedCurrent_i)
    avCurrentBins, P_of_avC = GenerateHistogram(MinCurrent, MaxCurrent, BinList[i], AveragedCurrent_i)
    P_of_avCurrents += [[avCurrentBins, P_of_avC]]
    
    # Subtract the mean to center the data and store it
    referenceCorrected_current = AveragedCurrent_i - ReferenceAverageCurrent
    refCorr_avCurrents_perT += [ AveragedCurrent_i - ReferenceAverageCurrent]
    
    # compute probability distribution
    refCorr_avCurrentBins, P_of_corrAvC = GenerateCenteredHistogram(BinList[i], referenceCorrected_current)
    P_of_CorrAvCurrents += [[refCorr_avCurrentBins, P_of_corrAvC]]


''' Compute the ratio of probabilities for positive vs. negative mean-zeroed currents.
    The bins for P(c) are symmetric around 0 and evenly spaced. '''
Ratio_ProbPosC_PropNegC_perT = []
Slope = np.zeros(NumOfT)
rSqList = np.zeros(NumOfT)
for i, data in enumerate(P_of_CorrAvCurrents):
    refCorr_avCurrentBins, P_of_corrAvC = data
    
    # Create masks to separate positive and negative current values
    MaskPositiveCurrent = refCorr_avCurrentBins > 0
    MaskNegativeCurrent = refCorr_avCurrentBins < 0
    
    # Extract and align probabilities for positive and negative current
    Prob_PositiveCurrent = P_of_corrAvC[ MaskPositiveCurrent ]
    Prob_NegativeCurrent = P_of_corrAvC[ MaskNegativeCurrent ]
    Prob_NegativeCurrent = np.flip( Prob_NegativeCurrent )
    
    # Extract and align current values
    PositiveCurrent = refCorr_avCurrentBins[MaskPositiveCurrent]
    NegativeCurrent = refCorr_avCurrentBins[MaskNegativeCurrent]
    NegativeCurrent = np.flip(NegativeCurrent)
    # Ensure symmetry of bin centers (check: c / -(-c) ≈ 1)
    assert ( np.mean(PositiveCurrent / - NegativeCurrent) - 1 < 10**(-6) )
    
    # Keep only bins where the negative-side probability is nonzero to avoid division by zero
    NonZeroProb = Prob_NegativeCurrent > 0
    Prob_PositiveCurrent = Prob_PositiveCurrent[ NonZeroProb ]
    Prob_NegativeCurrent = Prob_NegativeCurrent[ NonZeroProb ]
    PositiveCurrent = PositiveCurrent[ NonZeroProb ]
    
    # Compute the ratio of probabilities for +c vs. -c
    Ratio_ProbPosC_PropNegC = Prob_PositiveCurrent / Prob_NegativeCurrent
    
    # Store results for this T
    Ratio_ProbPosC_PropNegC_perT += [[ PositiveCurrent, Ratio_ProbPosC_PropNegC]]
    
    if exludeFromFit[i] > 0:
        ex = exludeFromFit[i]
        xList, yList = PositiveCurrent[:-ex], 1/T[i] * np.log(Ratio_ProbPosC_PropNegC[:-ex])
    else:
        xList, yList = PositiveCurrent, 1/T[i] * np.log(Ratio_ProbPosC_PropNegC)

    popt, pcov = curve_fit(LinearFunction, xList, yList)
    LinearFit = LinearFunction(xList, popt[0] )
    #ax31.plot(xList, LinearFit, '-', lw = 2, color = color[i], label = '$T=%d$' % T[i])
    Slope[i] = popt[0]
    
    #    https://en.wikipedia.org/wiki/Coefficient_of_determination
    SS_res = sum( (yList - LinearFit)**2 )
    mean = sum(yList)/len(yList)
    SS_tot = sum( (yList - mean)**2 )
    rSq = 1 - SS_res/SS_tot
    rSqList[i] = rSq
        
VisualizeDataAnalysis_nonRef(SaveFigAs, case, T, current, rawCurrentBins, P_of_rawCurrent, P_of_CorrAvCurrents, Ratio_ProbPosC_PropNegC_perT, exludeFromFit, Slope, rSqList)

# store data
OutDir = {'T' : T,
        'P_of_avCurrents' : P_of_avCurrents,
        'P_of_CorrAvCurrents' : P_of_CorrAvCurrents,
        'Ratio_ProbPosC_PropNegC_perT' : Ratio_ProbPosC_PropNegC_perT,
        'Slope' : Slope,
        'exludeFromFit' : exludeFromFit}
FileName = DataPath + "%s_Analysed.p" % case
pickle.dump( OutDir, open( FileName, 'wb' ) )
print('dump data in ', FileName)


