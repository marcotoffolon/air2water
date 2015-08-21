SUBROUTINE read_calibration

USE commondata
USE ifport              ! necessario per il comando makedirqq

IMPLICIT NONE
INTEGER:: i, j, status
CHARACTER(LEN=1) :: string
LOGICAL result          ! necessario per il comando makedirqq

! read input information
OPEN(unit=1,file='input.txt',status='old',action='read')
READ(1,*)		! header
READ(1,*) name
READ(1,*) air_station
READ(1,*) water_station
READ(1,*) series
READ(1,*) time_res
READ(1,*) version
READ(1,*) Tice_cover
READ(1,*) fun_obj
READ(1,*) num_mod
READ(1,*) run
READ(1,*) prc
READ(1,*) n_run
READ(1,*) mineff_index
CLOSE(1)

station=TRIM(air_station)//'_'//TRIM(water_station)

WRITE(string,'(i1)' ) version

folder = TRIM(name)//'/output_'//string//'/'
result=makedirqq(folder)
    
WRITE(*,*) 'Objective function ',fun_obj

IF (run .eq. 'FORWARD') THEN
    OPEN(unit=1,file=TRIM(name)//'/parameters_forward.txt',status='old',action='read')    
    READ(1,*) (par(i), i=1,n_par)  
ELSE IF (run .eq. 'PSO') THEN    
    ! read PSO parameters
    OPEN(unit=1,file='PSO.txt',status='old',action='read')
    READ(1,*)		    ! header
    READ(1,*) n_particles
    READ(1,*) c1,c2
    READ(1,*) wmax,wmin
    CLOSE(1)
END IF

IF (run .eq. 'PSO' .or. run .eq. 'LATHYP') THEN    
    ! read model parameters
    OPEN(unit=1,file=TRIM(name)//'/parameters.txt',status='old',action='read')

    READ(1,*) (parmin(i),i=1,n_par);	
    READ(1,*) (parmax(i),i=1,n_par);	

    ! parameters that are not used are zeroed
    flag_par=.true.
    IF (version<=6 .and. version .ne. 1) THEN                    ! air2water 6 parameters
	    parmin(7)=0;    parmax(7)=0;    flag_par(7)=.false.;
	    parmin(8)=0;	parmax(8)=0;    flag_par(8)=.false.;
	    IF (version==4) THEN                ! air2water 4 parameters
		    parmin(5)=0;    parmax(5)=0;        flag_par(5)=.false.;
		    parmin(6)=0;	parmax(6)=0;        flag_par(6)=.false.;
	    END IF
    END IF
    
    IF (version==1) THEN
    		parmin(5)=0;    parmax(5)=0;        flag_par(5)=.false.;
		    parmin(6)=0;	parmax(6)=0;        flag_par(6)=.false.;
    END IF
    CLOSE(1)
    	
    ! write parameters
    OPEN(unit=2,file=TRIM(folder)//'/parameters.txt',status='unknown',action='write')
    WRITE(2,'(I2,A)') n_par, '   !numero parametri'
    WRITE(2,'(<n_par>(F10.5,1x))') (parmin(i),i=1,n_par)
    WRITE(2,'(<n_par>(F10.5,1x))') (parmax(i),i=1,n_par)
    CLOSE(2)
END IF

! read T series (calibration)
CALL read_Tseries('c')

RETURN
END

!-------------------------------------------------------------------------------
!				LETTURA PERIODO VALIDAZIONE
!-------------------------------------------------------------------------------
SUBROUTINE read_validation

USE commondata

IMPLICIT NONE

DEALLOCATE(date, tt, Tair, Twat_obs, Twat_obs_agg, Twat_mod, Twat_mod_agg, delta)
DEALLOCATE(I_pos, I_inf)

! read T series (validation)
CALL read_Tseries('v')

RETURN
END

!-------------------------------------------------------------------------------
!				LETTURA FILE TEMPERATURA
!-------------------------------------------------------------------------------
SUBROUTINE read_Tseries(p)

USE commondata

IMPLICIT NONE

INTEGER :: i, j, k, status
INTEGER :: leap, year_ini
CHARACTER(LEN=1),INTENT(IN) :: p
CHARACTER(LEN=10) :: period

n_tot=0;

IF (p=='c') THEN
    period='calibration'
ELSE
    period='validation'
END IF

OPEN(unit=3,file=TRIM(name)//'/'//TRIM(station)//'_'//series//p//'.txt',status='unknown',action='read', iostat=status)
openif3: IF (status==0) THEN
	readloop3: DO
		READ(3,*,iostat=status)
		IF (status/=0) EXIT
		n_tot=n_tot+1
	END DO readloop3
	readif3: IF(status>0) THEN
	END IF readif3	
END IF openif3
REWIND(3)

! allocation + replication of the 1st year
WRITE(*,1001)  n_tot/365.25,TRIM(period)
1001 FORMAT('There are ',f4.1,' years for ', a12)

IF (p=='v' .and. n_tot .lt. 365) THEN
    WRITE(*,*) 'Validation period < 1 year --> validation is skipped'
    GO TO 100
END IF

n_year=CEILING(n_tot/365.25)
n_tot=n_tot+365             ! the 1st year is replicated. The 1st year is always considered 365 days long
ALLOCATE(date(n_tot,3),stat=status)
ALLOCATE(Tair(n_tot),stat=status)
ALLOCATE(Twat_obs(n_tot),stat=status) 
ALLOCATE(Twat_obs_agg(n_tot),stat=status) 
ALLOCATE(Twat_mod(n_tot),stat=status) 
ALLOCATE(Twat_mod_agg(n_tot),stat=status) 
ALLOCATE(tt(n_tot),stat=status)
ALLOCATE(delta(n_tot),stat=status) 

DO i=366,n_tot
	READ(3,*) (date(i,j),j=1,3),Tair(i),Twat_obs(i)
END DO
date(1:365,:)=-999
Tair(1:365)=Tair(366:730)
Twat_obs(1:365)=Twat_obs(366:730)

CLOSE(3)

year_ini=date(366,1)

! check leap years + define tt
k=0
DO j=1,365
    tt(k+j)=REAL(j)/365.0d0
END DO
k=365
DO i=1,n_year
    CALL leap_year(year_ini+i-1,leap)
    IF(leap==0) THEN
        DO j=1,365
            IF (k+j .gt. n_tot) THEN
                EXIT
            END IF
            tt(k+j)=REAL(j)/365.0d0
        END DO
        k=k+365
    ELSE
        DO j=1,366
            IF (k+j .gt. n_tot) THEN
                EXIT
            END IF
            tt(k+j)=REAL(j)/366.0d0
        END DO
        k=k+366
    END IF
END DO


100 RETURN 
END
