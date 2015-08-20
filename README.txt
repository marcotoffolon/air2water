# air2water
Model to predict Lake Surface Temperature (LST) using air temperature.

       .__       ________                  __                
_____  |__|______\_____  \__  _  _______ _/  |_  ___________ 
\__  \ |  \_  __ \/  ____/\ \/ \/ /\__  \\   __\/ __ \_  __ \
 / __ \|  ||  | \/       \ \     /  / __ \|  | \  ___/|  | \/
(____  /__||__|  \_______ \ \/\_/  (____  /__|  \___  >__|   
     \/                  \/             \/          \/       
Version 1.0 - August 2015

Provided by sebastiano Piccolroaz and Marco Toffolon

Department of Civil, Environmental, and Mechanical Engineering, University of Trento (Italy)
email contacts: s.piccolroaz@unitn.it, marco.toffolon@unitn.it

How to cite:

- Piccolroaz S., M. Toffolon, and B. Majone (2013), A simple lumped model to convert air temperature into surface water temperature in lakes, Hydrol. Earth Syst. Sci., 17, 3323-3338, doi:10.5194/hess-17-3323-2013

- Toffolon M., S. Piccolroaz, B. Majone, A.M. Soja, F. Peeters, M. Schmid and A. Wüest 2014), Prediction of surface water temperature from air temperature in lakes with different morphology, Limnology and Oceanography, 59(6), 2185-2202, doi: 10.4319/lo.2014.59.6.2185
-----------------------------------------------------------------------------------------------------
 
How to use air2water 

1. Input files

-----> 'input.txt' to be located in the same folder of the executable file of air2water. This file contains the input information:

! Main input
Superior		! name of the lake
stndrck   		! name/ID of the air temperature station
mrqtt			! name/ID of the water temperature station
c				! type of series: c=continuous series, m=mean year
1d        		! time resolution: 1d=daily, nw=n weeks (n=1,2,...), 1m=monthly
8           	! version: 4, 6, and 8 parameters
0				! Threshold temperature for ice formation
RMS				! objective function: RMS (Root Mean Square Error), NSE (Nash-Sutcliffe Efficiency Index), KGE (Kling-Gupta Efficiency Index)
RK4				! mod_num : numerical method used to solve the model equation EUL (Euler Explicit), RK2 (Runge-Kutta 2), RK4 (Runge-Kutta 4)
PSO            	! run mode: PSO (Particle Swarm Optimization), LATHYP (Latin Hypercube), RANSAM (Random Sampling), or FORWARD
0.60			! minimum percentage of data: 0...1. E.g., when using 1m time resolution, the monthly average is set to NaN if available data during one month are less than this percentage (60% in this case) 
500				! n_run: number of iterations
-999			! minimum efficiency index (i.e., RMS, NSE, KGE). All model realization with efficiency index greater than this threshold are saved in file 0_...


-----> 'PSO.txt'   to be located in the same folder of the executable file of air2water. This file contains the parameters of PSO:
		see e.g.:	
		- Kennedy, J., and R. Eberhart (1995), Particle swarm optimization, in Proceedings of IEEE International Conference on 	Neural Networks, pp. 1942-1948, doi:10.1109/ICNN.1995.488968.
		- Robinson, J., and Y. Rahmat-Samii (2004), Particle swarm optimization in electromagnetics, IEEE T. Antenn. Propag., 52, 397-407, doi:10.1109/TAP.2004.823969.

! PSO parameters
500			! number of particles
2 2	   		! c1, c2: constants known as cognitive and social learning factors (used in the equation of motion of particles)
0.9  0.4    ! inertia min and max (used in the equation of motion of particles)


-----> A folder named according to the first line in file input.txt (e.g., Superior), located in the same folder of the executable file of air2water and containing:

	-----> input files for calibration 'IDair_IDwater_typeofseriesc.txt' (e.g., stndrck_mrqtt_cc.txt) and validation 'IDair_IDwat_typeofseriesv.txt' (e.g., stndrck_mrqtt_cv.txt)
		The two files have 5 columns (year, month, day, air temperature, water temperature):

		2001		1		1		-0.125		1.000		
		2001		1		2		 0.158		1.155		
		2001		1		3		 1.155		-999		
		2001		1		4		 3.458		-999		
		2001		1		5		10.125		3.254		
		...
		
		NOTE: 	the series of observed air temperature must be complete. It cannot have gaps or no data. 
				the series of observed water temperature can contain no-data (-999). 
				Both series are always at daily time scale, as the equation of the model is solved with daily time step. The model automatically evaluates weekly, multi-weekly, or monthly averages (of water temperature) when using different time scales for model calibration. 
	
	-----> file 'parameters.txt'
		The file contains the range of variation of each parameter of the model. The range of variation of the parameters is lake dependent, and can be estimated using the equations presented in Piccolroaz et al., 2013. 
		The structure of the file is as follows:

		-0.017   0.00047   0.00069   2   0.005   0   0     0
		 0.100   0.00850   0.02000   6   0.150   1   150   0.5

		The first line contains the minimum value of each parameters, while the second contains the maximum values. 
		
	-----> file 'parameters_forward.txt'
		The file contains a set of parameters to be used to run the model in forward mode (required only if the model is run in forward), see 10th line of file input.txt . The structure of the file is as follows:

		 0.000284   0.006668   0.006719   2.888324   0.017322   0.212055 148.709752   0.499891  -1.303755

		
2. Output files

The model writes the output in a folder called 'output_numberofparameters' (e.g., output_8) located in the same folder of the input files (e.g., Superior).
The folder contains:


-----> file '0_PSO_..........   .out'
	Binary file that contains a matrix with 9 columns. Each row contains the set of parameters (columns 1-8) and the associated efficiency index (column 9) of each iteration performed by the optimization algorithm. Values are saved in double precision.
	This file allows: a) drawing the dotty plots of parameters to evaluate whether the a priori range of variations of parameters have been reasonably defined (i.e., not too narrow, not too large), b) evaluating wheter the optimization (searching) algorithm converged towards the best solution, c) perform uncertainty anayses (when using RANSAM or LATHYP).
	

-----> file '1_PSO_..........   .out'
	ASCII file that contains the best set of parameters (1st line) and the value of the efficiency index (2nd line)

	
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
	
	Note that the first year is replicated and is used as warm up year to mitigate the influence of initial conditions. During the warm up year, year, month and day are equal to -999.
	
	
-----> file '3_PSO_..........   .out'
     ASCII file with the same structure of file '2_PSO... .out' but referred to the validation period.

	 
-----> file '4_PSO_..........   .out'
     Delta (dimensionless thickness of the well-mixed surface layer) during the calibration period.

	 
-----> file '5_PSO_..........   .out'
     Delta (dimensionless thickness of the well-mixed surface layer) during the validation period.
	 
	 
3. Output files

Coming soon...

