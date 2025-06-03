This directory contains everything to create Fig 1.

###################################### 
Underlying data:

./Data/ParticleTrajectoryAnnulusArena.txt
./Data/ParticleTrajectoryDiscArena.txt

Data files contain trajectories of the sqigglebot in annulus arena and Disc arena.
Content of data files: 1st column x position, 2nd column y position.
Units of positions are in pixel and can be converted to cm through the relation 8.2 cm = 148 pix.
Frames (individual positions) are recorded at a rate of 25 frames per second.

######################################
The python scripts 02ComputeFig1* evaluate the data and generate curves for Fig1.

--------------------------------------
02ComputeFig1b_SpartialPropDensity.py : compute quantities of Fig 1b
- creates the histogram (probability distribution) of positions in Disc arena in polar coordinates (radius, azimuth)
- extracts sample section of total particle trajectories 

Input : ./Data/ParticleTrajectoryDiscArena.txt

Output :

./Data/F1b_RadiusBins.txt
list of radial bin boundaries [r_0, r_1, r_2, ..., r_{N_r+1}]
1st bin : r_0 to r_1
N_r'th bin: r_{N_r} to r_{N_r+1}

./Data/F1b_AzimutBins.txt
list of azimuthal bin boundaries [theta_0, theta_1, theta_2, ..., theta_{N_a+1}]
1st bin: theta_0 to theta_1
N_a'th bin: theta_{N_a} to theta_{N_a+1}

./Data/F1b_HistogramForPlot.txt
(N_r x N_a) matrix with histogram entries
Histogram is normalised probability, i.e. sum_i H_i A_i = 1, where the sum runs over all histogram bins H_i and A_i is the corresponding area of each bin

./Data/F1b_SampleTrajectory.txt
a sample section of the total trajectory in ./Data/ParticleTrajectoryDiscArena.txt

--------------------------------------
02ComputeFig1c_AngulaMSD.py : compute quantities of Fig 1c
- computes the angular MSD of sqigglebot in the annulus arena

Input : ./Data/ParticleTrajectoryAnnulusArena.txt

Output : 
./Data/F1c_AngularMSD.txt
1st column physical time in seconds
2nd column angular MSD in rad^2

--------------------------------------
02ComputeFig1c_MSD.py : compute quantities of Fig 1c
- computes the MSD of sqigglebot in the disc arena

Input : ./Data/ParticleTrajectoryDiscArena.txt

Output : 
./Data/F1c_MSD.txt
1st column physical time in seconds
2nd column MSD in cm^2

./Data/F1c_MSD_BallisticTheory.txt
1st column physical time in seconds
2nd column ballistic scaling, i.e. time^2

--------------------------------------
02ComputeFig1e_SpeedDistribution.py : compute quantities of Fig 1e
- computes the probability distribution of speed in Annulus and Disc arenas

Input : 
./Data/ParticleTrajectoryAnnulusArena.txt
./Data/ParticleTrajectoryDiscArena.txt

Output:
./Data/F1e_SpeedDistributionDisc.txt
1st column speed in Disc arena
2nd column probability of speed

./Data/F1e_SpeedDistributionAnnulus.txt
1st column speed in Annulus arena
2nd column probability of speed

######################################
03PlotFig1.py generates Figure 1 with the date generated with scripts 02ComputeFig1*
