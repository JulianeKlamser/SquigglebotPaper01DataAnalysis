import pickle
import numpy as np
import os
from Functions import GenerateHistogram
from Functions import GenerateCenteredHistogram
from Functions import VisualizeDataAnalysis
script_dir = os.path.dirname(os.path.abspath(__file__))

DataPath = script_dir + '/Data/'
PlotPath = script_dir + '/Plots/'

''' read time range of stable current '''
Limits = np.loadtxt(DataPath + 'StableLimits.txt', dtype=np.int32)

''' read reference current 
    and apply stable limits '''
FileID = 'Current_reference.txt'
filePath = DataPath + FileID
current = np.loadtxt(filePath)
current = current[Limits[0]:Limits[1]]
MinCurrent, MaxCurrent = current.min(), current.max()

''' compute probability distribution of the raw
    current '''
secondLargestMaxCurrent = max( current[current < MaxCurrent] )
binWidth = (MaxCurrent - secondLargestMaxCurrent)
nBins = np.round( (MaxCurrent - MinCurrent) / binWidth )
nBins = int(nBins)+1
rawCurrentBins, P_of_rawCurrent = GenerateHistogram(MinCurrent - 0.5*binWidth, MaxCurrent + 0.5*binWidth, nBins, current)

''' define time window T for the sliding average
    and the number of histogram bins per T '''
BinList = [60, 74, 80, 80, 90]
NumOfT = len(BinList)
T = np.array([75+25*i for i in range(NumOfT)])

''' Compute sliding averages of the current over different time windows T
    For each value of T:
      - Compute the sliding (moving) average of the current with window size T_i
      - Store the result
      - Compute the mean of that sliding average and store it
      - Also store a mean-subtracted version of the sliding average (centered around 0) 
      
    finally safe data'''
avCurrent_perT = [] # sliding averages of current for each T
MeanCurrent_perT = np.zeros(NumOfT)  # mean value of each sliding average
mean0_avCurrents_perT = [] # store zero-mean (mean-subtracted) sliding averages
for i, T_i in enumerate(T):
    # Compute sliding average with window size T_i
    AveragedCurrent_i = np.convolve(current, np.ones(T_i)/T_i, mode='valid')
    avCurrent_perT += [ AveragedCurrent_i ]
    
     # Compute and store the mean of the sliding average
    meanC = np.mean(AveragedCurrent_i)
    MeanCurrent_perT[i] = meanC
    
    # Subtract the mean to center the data and store it
    mean0_avCurrents_perT += [ AveragedCurrent_i - meanC]

''' compute probability distributions of sliding averaged
    and mean-subtracted version of the sliding averaged
    current per T '''
P_of_avCurrents = []
P_of_mean0_avCurrents = []
for i, AveragedCurrent_i in enumerate(avCurrent_perT):

    MinCurrent, MaxCurrent = min(AveragedCurrent_i), max(AveragedCurrent_i)
    avCurrentBins, P_of_avC = GenerateHistogram(MinCurrent, MaxCurrent, BinList[i], AveragedCurrent_i)
    P_of_avCurrents += [[avCurrentBins, P_of_avC]]
    
    mean0_avCurrent_i = mean0_avCurrents_perT[i]
    mean0_avCurrentBins, P_of_mean0avC = GenerateCenteredHistogram(BinList[i], mean0_avCurrent_i)
    P_of_mean0_avCurrents += [[mean0_avCurrentBins,P_of_mean0avC]]
        
''' Compute the ratio of probabilities for positive vs. negative mean-zeroed currents.
    The bins for P(c) are symmetric around 0 and evenly spaced. '''
Ratio_ProbPosC_PropNegC_perT = []
for i, data in enumerate(P_of_mean0_avCurrents):
    mean0_avCurrentBins, P_of_mean0avC = data
    
    # Create masks to separate positive and negative current values
    MaskPositiveCurrent = mean0_avCurrentBins > 0
    MaskNegativeCurrent = mean0_avCurrentBins < 0
    
    # Extract and align probabilities for positive and negative current
    Prob_PositiveCurrent = P_of_mean0avC[ MaskPositiveCurrent ]
    Prob_NegativeCurrent = P_of_mean0avC[ MaskNegativeCurrent ]
    Prob_NegativeCurrent = np.flip( Prob_NegativeCurrent )
    
    # Extract and align current values
    PositiveCurrent = mean0_avCurrentBins[MaskPositiveCurrent]
    NegativeCurrent = mean0_avCurrentBins[MaskNegativeCurrent]
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
    
# store data
OutDir = {'T' : T,
        'MeanCurrent_perT' : MeanCurrent_perT,
        'P_of_avCurrents' : P_of_avCurrents,
        'P_of_mean0_avCurrents' : P_of_mean0_avCurrents,
        'Ratio_ProbPosC_PropNegC_perT' : Ratio_ProbPosC_PropNegC_perT }
FileName = DataPath + "Reference_Analysed.p"
pickle.dump( OutDir, open( FileName, 'wb' ) )
print('dump data in ', FileName)

SaveFigAs = PlotPath + 'DataAnalysisReferenceCurrent.png'
VisualizeDataAnalysis(SaveFigAs, T, current, rawCurrentBins, P_of_rawCurrent, P_of_avCurrents, P_of_mean0_avCurrents, Ratio_ProbPosC_PropNegC_perT)
