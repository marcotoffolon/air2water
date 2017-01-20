PROGRAM air2water

USE commondata
IMPLICIT NONE
INTEGER :: i,j,k,status
REAL(KIND=8):: T1,T2

WRITE(*,*) '       .__       ________                  __                '
WRITE(*,*) '_____  |__|______\_____  \__  _  _______ _/  |_  ___________ '
WRITE(*,*) '\__  \ |  \_  __ \/  ____/\ \/ \/ /\__  \\   __\/ __ \_  __ \'
WRITE(*,*) ' / __ \|  ||  | \/       \ \     /  / __ \|  | \  ___/|  | \/'
WRITE(*,*) '(____  /__||__|  \_______ \ \/\_/  (____  /__|  \___  >__|   '
WRITE(*,*) '     \/                  \/             \/          \/       '
WRITE(*,*) 'Version 2.0.0 - January 2017'
WRITE(*,*) 
WRITE(*,*) 
WRITE(*,*) 
!-------------------------------------------------------------------------------------
!
! Provided by Sebastiano Piccolroaz and Marco Toffolon
!
! Department of Civil, Environmental, and Mechanical Engineering, University of Trento (Italy)
! email contacts: s.piccolroaz@unitn.it, marco.toffolon@unitn.it
!
! How to cite:
!
! Piccolroaz S., M. Toffolon, and B. Majone (2013), A simple lumped model to convert
! air temperature into surface water temperature in lakes, Hydrol. Earth Syst. Sci., 
! 17, 3323-3338, doi:10.5194/hess-17-3323-2013
!
! Toffolon M., S. Piccolroaz, B. Majone, A.M. Soja, F. Peeters, M. Schmid and A. Wüest 
! (2014), Prediction of surface water temperature from air temperature in lakes with 
! different morphology, Limnology and Oceanography, 59(6), 2185-2202, doi: 10.4319/lo.2014.59.6.2185
! 
! Piccolroaz S., M. Toffolon, and B. Majone (2015), The role of stratification on lakes’ thermal
! response: The case of Lake Superior, Water Resources Research, 51(10):7878–7894,
! DOI:10.1002/2014WR016555
!
! Piccolroaz S. (2016), Prediction of lake surface temperature using the air2water model: 
! guidelines, challenges, and future perspectives, Advances in Oceanography and Limnology,
! 7:36-50, DOI: http://dx.doi.org/10.4081/aiol.2016.5791
!-------------------------------------------------------------------------------------

CALL CPU_TIME(T1)

! allocation of parameter matrices
ALLOCATE(parmin(n_par),stat=status) 
ALLOCATE(parmax(n_par),stat=status) 
ALLOCATE(flag_par(n_par),stat=status) 
ALLOCATE(par(n_par),stat=status) 
ALLOCATE(par_best(n_par),stat=status) 

! lettura dati in input
CALL read_calibration

! aggregate calibration data on the basis of time_res
CALL aggregation

! evaluate mean, TSS e standard deviation of data
CALL statis
WRITE(*,*) 'mean, TSS and standard deviation (calibration)'
WRITE(*,*)  SNGL(mean_obs),SNGL(TSS_obs),SNGL(std_obs)

OPEN(UNIT=11,FILE=TRIM(folder)//'/1_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')

ii=0
IF (run .eq. 'FORWARD') THEN
    CALL forward_mode
ELSE IF (run .eq. 'PSO') THEN
    CALL PSO_mode
ELSE IF (run .eq. 'LATHYP') THEN
    CALL LH_mode
END IF

OPEN(UNIT=99,FILE=TRIM(folder)//'/00_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
WRITE(99,*) ii
CLOSE(99)

! lancio in forward con il set di parametri migliore  
CALL forward

CALL CPU_TIME(T2)
print *, 'Computation time was ', T2-T1, 'seconds.'

STOP
END PROGRAM air2water
