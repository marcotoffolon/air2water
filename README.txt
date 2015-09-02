       .__       ________                  __                
_____  |__|______\_____  \__  _  _______ _/  |_  ___________ 
\__  \ |  \_  __ \/  ____/\ \/ \/ /\__  \\   __\/ __ \_  __ \
 / __ \|  ||  | \/       \ \     /  / __ \|  | \  ___/|  | \/
(____  /__||__|  \_______ \ \/\_/  (____  /__|  \___  >__|   
     \/                  \/             \/          \/       
A model to predict Lake Surface Temperature (LST) using air temperature.
Version 1.0.0 - August 2015

Provided by Sebastiano Piccolroaz and Marco Toffolon

Department of Civil, Environmental, and Mechanical Engineering, University of Trento (Italy)
email contacts: s.piccolroaz@unitn.it, marco.toffolon@unitn.it

How to cite:

- Piccolroaz S., M. Toffolon, and B. Majone (2013), A simple lumped model to convert air temperature into surface water temperature in lakes, Hydrol. Earth Syst. Sci., 17, 3323-3338, doi:10.5194/hess-17-3323-2013

- Toffolon M., S. Piccolroaz, B. Majone, A.M. Soja, F. Peeters, M. Schmid and A. Wüest 2014), Prediction of surface water temperature from air temperature in lakes with different morphology, Limnology and Oceanography, 59(6), 2185-2202, doi: 10.4319/lo.2014.59.6.2185
-----------------------------------------------------------------------------------------------------------------------------------------
 
How to use air2water 

0. Preamble
The model is presented and discussed in Piccolroaz et al., 2013.
The version of air2water released here is that described in Toffolon et al., 2014, which coincides with the original version (Piccolroaz et al., 2013), except for some changes in the order/definition of parameters. The mathematical meaning is the same.

Precompiled executable files are available for Linux (air2water_1.0.0.out) and Windows (air2water_1.0.0.exe) systems.

The code is provided with an example test case. The study site is Lake Superior and the data have been downloaded from (downloaded data have been pre-processed):
- National Oceanic and Atmospheric Administration’s (NOAA) National Data Buoy Center (NDBC, webpage: http://www.ndbc.noaa.gov/) --> daily air temperature (STDM4 – Stannard Rock station);
- National Oceanic and Atmospheric Administration’s (NOAA) Great Lakes Environmental Research Laboratory (GLERL, webpage: http://www.glerl.noaa.gov/) --> daily lake surface temperature (satellite-derived LST averaged over the whole lake).


1. Input files
-----> 'input.txt' to be located in the same folder of the executable file of air2water. This file contains the input information:
! Main input
Superior		! name of the lake
stndrck   		! name/ID of the air temperature station
sat				! name/ID of the water temperature station
c				! type of series: c=continuous series, m=mean year
1d        		! time resolution: 1d=daily, nw=n weeks (n=1,2,...), 1m=monthly
8          		! version: 4, 6, and 8 parameters
0				! Threshold temperature for ice formation
RMS				! objective function: RMS (Root Mean Square Error), NSE (Nash-Sutcliffe Efficiency Index), KGE (Kling-Gupta Efficiency Index)
RK4				! mod_num : numerical method used to solve the model equation EUL (Euler Explicit), RK2 (Runge-Kutta 2), RK4 (Runge-Kutta 4)
PSO        		! run mode: PSO (Particle Swarm Optimization), LATHYP (Latin Hypercube), RANSAM (Random Sampling), or FORWARD
0.60			! minimum percentage of data: 0...1. E.g., when using 1m time resolution, the monthly average is set to NaN if available data during one month are less than this percentage (60% in this case) 
500				! n_run: number of iterations
-999			! minimum efficiency index (i.e., RMS, NSE, KGE). All model realization with efficiency index greater than this threshold are saved in file 0_...

-----> 'PSO.txt'   to be located in the same folder of the executable file of air2water. This file contains the parameters of PSO:
see e.g.:	
- Kennedy, J., and R. Eberhart (1995), Particle swarm optimization, in Proceedings of IEEE International Conference on Neural Networks, pp. 1942-1948, doi:10.1109/ICNN.1995.488968.
- Robinson, J., and Y. Rahmat-Samii (2004), Particle swarm optimization in electromagnetics, IEEE T. Antenn. Propag., 52, 397-407, doi:10.1109/TAP.2004.823969.

! PSO parameters
500			! number of particles
2 2	   		! c1, c2: constants known as cognitive and social learning factors (used in the equation of motion of particles)
0.9  0.4    ! inertia max and min (used in the equation of motion of particles)

-----> A folder named according to the first line in file input.txt (e.g., Superior), located in the same folder of the executable file of air2water and containing:

-----> input files for calibration 'IDair_IDwater_typeofseriesc.txt' (e.g., stndrck_sat_cc.txt) and validation 'IDair_IDwat_typeofseriesv.txt' (e.g., stndrck_sat_cv.txt)
	The two files have 5 columns (year, month, day, air temperature, water temperature):

	2001		1		1		-0.125		1.000		
	2001		1		2		 0.158		1.155		
	2001		1		3		 1.155		-999		
	2001		1		4		 3.458		-999		
	2001		1		5		10.125		3.254		
	...
		
	NOTE: 	The series of observed air temperature must be complete. It cannot have gaps or no data. 
		The series of observed water temperature can contain no-data (-999). 
		Both series are always at daily time scale, as the equation of the model is solved with daily time step. The model automatically evaluates weekly, multi-weekly, or monthly averages (of water temperature) when using different time scales for model calibration. 
	
-----> file 'parameters.txt'
	The file contains the range of variation of each parameter of the model. The range of variation of the parameters can be defined on the basis of physical considerations and is strictly lake dependent, thus its reasonable a priori estimation is crucial to obtain reliable results. The range of parameters can be estimated using the equations presented in Appendix A of Piccolroaz et al., 2013 (notice that the version of air2water released here is that described in Toffolon et al., 2014, thus there are some differences compared to the parameters defined in Piccolroaz et al., 2013). A pre-processing script in matlab is available to evaluate the a priori physical range of model's parameter (see point 3. Pre-processing).
	The structure of the file 'parameters.txt' is as follows:

	-0.017   0.00047   0.00069   2   0.005   0   0     0
	 0.100   0.00850   0.02000   6   0.150   1   150   0.5

	The first line contains the minimum value of each parameters, while the second contains the maximum values. 
		
-----> file 'parameters_forward.txt'
	The file contains a set of parameters to be used to run the model in forward mode (required only if the model is run	in forward), see 10th line of file input.txt . The structure of the file is as follows:

	 0.000284   0.006668   0.006719   2.888324   0.017322   0.212055 148.709752   0.499891  -1.303755

	
2. Output files
The model writes the output in a folder called 'output_numberofparameters' (e.g., output_8) located in the same folder of the input files (e.g., Superior).
The folder contains:

-----> file '0_PSO_..........   .out'
	Binary file that contains a matrix with 9 columns. Each row contains the set of parameters (columns 1-8) and the 	associated efficiency index (column 9) of each iteration performed by the optimization algorithm. Values are saved in double precision.
	This file allows: a) drawing the dotty plots of parameters to evaluate whether the a priori range of variations of 	parameters have been reasonably defined (i.e., not too narrow, not too large), b) evaluating whether the optimization (searching) algorithm converged towards the best solution, c) perform uncertainty analyses (when using RANSAM or LATHYP).
	
-----> file '1_PSO_..........   .out'
	ASCII file that contains the best set of parameters (1st line), the value of the efficiency index during the calibration period (2nd line), and the value of the efficiency index during the validation period (if any, 3nd line)

-----> file '2_PSO_..........   .out'
	ASCII file containing the results of the calibration period in which the columns are:
	year
	month
	day
	observed air temperature
	observed water temperature
	simulated water temperature
	observed water temperature aggregated at the time scale chosen in file 'input.txt' (e.g., 1m) 
	simulated water temperature aggregated at the time scale chosen in file 'input.txt' (e.g., 1m) 
	
	Note that the first year is replicated and is used as warm up year to mitigate the influence of initial conditions. 
	During the warm up year, year, month and day are equal to -999.
	
-----> file '3_PSO_..........   .out'
	ASCII file with the same structure of file '2_PSO... .out' but referred to the validation period.

-----> file '4_PSO_..........   .out'
     	Delta (dimensionless thickness of the well-mixed surface layer) during the calibration period.

-----> file '5_PSO_..........   .out'
     	Delta (dimensionless thickness of the well-mixed surface layer) during the validation period.
	 

3. Pre-processing
A pre-processing script in matlab is available to evaluates the a priori physical range of model's parameter. The script ('pre_processing.m') is located in the folder 'pre_processing'. The required data to run the script are:
- estimates of the annual minimum and maximum solar radiation at the lake site (W/m2);
- estimated range of variation of the maximum reactive volume of the lake (in terms of depth, m). 

Definition: the reactive volume is the volume participating to the heat exchange with the atmosphere. The reactive volume is maximum when the lake is not stratified. 

The output of the pre-processing script is the file 'parameters_physical_range.txt' and is saved in the folder 'pre_processing'. The estimated range of variation of parameters may require a narrowing in order to get a more accurate calibration of model's parameters. This is especially the case when the chosen range of variation of the maximum reactive volume depth is too large. In any case, the range of parameters should be carefully checked during the post-processing step through the analysis of the dotty plots (the cloud of points obtained by performing random sampling or latin hypercube should be centred within the searching domain, and in any case the best set of parameters should be sufficiently far from the upper and lower bounds. Exceptions are represented by parameters p6 and p7, see Piccolroaz et al., 2013 and Toffolon et al., 2014). 
  
	 
4. Post-processing
A matlab script is available for post-processing. The script ('post_processing.m') is located in the folder 'post_processing'. The output of the script is saved in the same folder where the output of air2water is located. 
The figures produced by the post-processing script are:
- dotty plots of parameters;
- comparison between observed air temperature, observed water temperature, and simulated water temperature for the calibration period;
- comparison between observed air temperature, observed water temperature, and simulated water temperature for the validation period.
