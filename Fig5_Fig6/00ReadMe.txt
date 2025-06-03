This directory contains everything to create Fig 5 and 6.

######################################
Underlying data:
./Data/Current_boundary.txt./Data/Current_bulk.txt./Data/Current_reference.txt

The data files contain the time series of the current of one sqigglebot in annulus arena (boundary), disk arena (bulk), and in the reference state.

Content of data file: current recorded with a rate of 25 frames per second


######################################
The python script 01_StableC.py :

- identifies stable region of the current based on the average current and its variance

Input :
./Data/Current_boundary.txt./Data/Current_bulk.txt./Data/Current_reference.txt

Output :

'./Data/StableLimits.txt'
- minimum and maximum frame ID defining the region of stable current

'./Plots/StabilityMonitor.png'
- illustration

######################################
The python script 02_AnalyseReferenceCurrent.py : 
- analysis of the reference current within the stable limits

Input :
./Data/Current_reference.txt
./Data/StableLimits.txt

Output :
./Plots/DataAnalysisReferenceCurrent.png
- illustration of the analysis

./Data/Reference_Analysed.p
   Type : python pickle
   Pickles Object : python dictionary 
   Dictionary content:
 	DataPath = './Data/Reference_Analysed.p'
	with open(DataPath, 'rb') as f:
	    Dir = pickle.load(f)
	Dir['T'] = np.array([T_0, T_1, T_2, ...]) : time window of the sliding average in numbers of frames
	Dir['MeanCurrent_perT'] = [M0, M1, M2, ...] : mean current per sliding average with T
       	Dir['P_of_avCurrents'] = [ [[C0_0, C0_1, C0_2, ...], [PofC0_0, PofC0_1, ...]], [[C1_0, C1_1, ...], [PofC1_0, PofC1_1, ...]] ] Probability distribution of sliding averaged current per T
        Dir['P_of_mean0_avCurrents'] : same as P_of_avCurrents for mean-subtracted sliding averaged current P( c - <c> )
        Dir['Ratio_ProbPosC_PropNegC_perT'] : same data structure as P_of_avCurrents but for P( c - <c> ) / P( -{c - <c>} )

######################################
The python script 03_EvaBCurrent.py : 
- analysis the current of the annulus and disk arena with stable limits

----> run as 
python3 03_EvaBCurrent.py bulk
to analyse the current of the disk arena

----> run as 
python3 03_EvaBCurrent.py boundary
to analyse the current of the annulus arena

Input :
./Data/Current_reference.txt
./Data/StableLimits.txt
./Data/Reference_Analysed.p

Output :
./Plots/Analysis_boundary.png
./Plots/Analysis_bulk.png
- illustration of data analysis

./Data/boundary_Analysed.p
./Data/bulk_Analysed.p
   Type : python pickle
   Pickles Object : python dictionary 
   Dictionary content:
 	DataPath = './Data/boundary_Analysed.p'
	with open(DataPath, 'rb') as f:
	    Dir = pickle.load(f)

	Dir['T'] = np.array([T_0, T_1, T_2, ...]) : time window of the sliding average in numbers of frames
        Dir['P_of_avCurrents'] = [ [[C0_0, C0_1, C0_2, ...], [PofC0_0, PofC0_1, ...]], [[C1_0, C1_1, ...], [PofC1_0, PofC1_1, ...]] ] Probability distribution of sliding averaged current per T
        Dir['P_of_CorrAvCurrents'] : same as P_of_avCurrents, but for P( c - <c>_ref )
        Dir['Ratio_ProbPosC_PropNegC_perT'] : same data structure as P_of_avCurrents but for P( c - <c>_ref ) / P( -{c - <c>_ref} )
        Dir['Slope'] = [s_0, s_1, s_2, ...] : slope of 1/T log( P( c - <c>_ref ) / P( -{c - <c>_ref} ) ) per T
        Dir['exludeFromFit'] = [e_0, e_1, e_2, ...] : number of data points in Ratio_ProbPosC_PropNegC_perT that were excluded from linear fit


######################################
The python script 04_PlotFig5.py : 
- generates figure 5

Input :
./Data/Reference_Analysed.p
./Data/boundary_Analysed.p
./Data/bulk_Analysed.p

######################################
The python script 05_PlotFig6.py : 
- generates figure 6

Input :
./Data/Reference_Analysed.p
./Data/boundary_Analysed.p
./Data/bulk_Analysed.p

######################################
The python script Functions.py
- contains function definitions called in 02_AnalyseReferenceCurrent.py & 03_EvaBCurrent.py