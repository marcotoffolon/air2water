!-------------------------------------------------------------------------------
!				sub_1
!-------------------------------------------------------------------------------
SUBROUTINE sub_1(ei)
USE commondata

! Dichiarazione variabili
implicit none
INTEGER:: i, j, k
REAL(KIND=8) :: ei
REAL(KIND=8) :: ind, TSS, mean_mod_ymo,mean_obs_ymo

! air2water
CALL model
	
! Evaluate the objective function
CALL funcobj(ei)

RETURN
END

!-------------------------------------------------------------------------------
!				AIR2WATER
!-------------------------------------------------------------------------------
SUBROUTINE model
USE commondata

! Variables declaration
IMPLICIT NONE

INTEGER :: j, k
INTEGER :: nsub
REAL(KIND=8) :: DD
REAL(KIND=8) :: K1, K2, K3, K4
REAL(KIND=8) :: pp, lambda, dt
REAL(KIND=8) :: dTair , ttt
REAL(KIND=8) :: Tairk, Tairk1, Twatk, ttk
REAL(KIND=8) :: alpha, d1, d2, A, B

Tmin=MAX(0.0d0,MINVAL(Twat_obs))
Tmin=MAX(4.d0,Tmin)	! Tmin

IF (Twat_obs(1)==-999) THEN	! Condizione iniziale
	Twat_mod(1)=Tmin
ELSE
	Twat_mod(1)=Twat_obs(1)
END IF

IF (num_mod .eq. 'CRN') THEN
    dt = 1.0d0
    ttt = DoY
END IF
        
DO j=1,n_tot-1

    IF (num_mod .eq. 'CRN') THEN
        CALL a2w(Tair(j), Twat_mod(j), tt(j), K1, DD, pp)  
        delta(j)=DD
        Twatk=( Twat_mod(j)*2.0d0*DD +                     &
        dt*( pp + par(1) + par(2)*Tair(j+1) + par(5)*cos(2.d0*pi*(tt(j)+ttt-par(6))) ) ) &
        /(2.0d0*DD+dt*par(3)) ;                                    	
        Twatk=MAX(Twatk,Tice_cover)                
    ELSE
        ! Substepping procedure
        CALL a2w(Tair(j), Twat_mod(j), tt(j), K1, DD, pp)
        DD=MAX(DD,0.01d0)
        lambda=(pp/par(4) - par(3))/DD
        pp=-LIM/lambda    
        IF (lambda .le. 0.0d0 .and. pp .lt. 1.0d0) THEN
            nsub=CEILING(1.0d0/pp)
            nsub=MIN(100,nsub)
            dt=1.0d0/nsub           
        ELSE
            dt=1.0d0
            nsub=1 
        END IF  
        dTair=(Tair(j+1)-Tair(j))/DBLE(nsub)
        ttt = DoY/DBLE(nsub)
            
	    Twatk=Twat_mod(j)
	    DO k=1,nsub ! Substepping cycle		
	        Tairk = Tair(j) + dTair*DBLE(k-1)
	        Tairk1= Tairk + dTair
	        ttk=tt(j) + ttt*DBLE(k-1)
    	    	    		
	        IF (num_mod .eq. 'RK4') THEN
	            CALL a2w(Tairk, Twatk, ttk, K1, DD, pp)
	            delta(j)=DD
	            CALL a2w(0.5d0*(Tairk + Tairk1), Twatk + 0.5d0*K1, ttk + 0.5*ttt, K2, DD, pp)
	            CALL a2w(0.5d0*(Tairk + Tairk1), Twatk + 0.5d0*K2, ttk + 0.5*ttt, K3, DD, pp)
	            CALL a2w(Tairk1, Twatk + K3, ttk + ttt, K4, DD, pp)

                Twatk=Twatk + 1.0d0/6.0d0*(K1 + 2.0d0*K2 + 2.0d0*K3 + K4 )*dt
            ELSEIF (num_mod .eq. 'RK2') THEN
	            CALL a2w(Tairk, Twatk, ttk, K1, DD, pp)
	            delta(j)=DD
	            CALL a2w(Tairk1, Twatk + K1, ttk + ttt, K2, DD, pp)
            	
                Twatk=Twatk + 0.5d0*(K1 + K2)*dt
            ELSEIF (num_mod .eq. 'EUL') THEN
	            CALL a2w(0.5d0*(Tairk + Tairk1), Twatk, ttk, K1, DD, pp)
	            delta(j)=DD
            	
	            Twatk=Twatk + K1*dt                    	    
	        ELSE
	            WRITE(*,*) 'Error in the choice of the numerical model'
	            STOP
            END IF
            
            Twatk=MAX(Twatk,Tice_cover)
            
        END DO
    END IF
    
    Twat_mod(j+1)=Twatk 
    Tmin=MAX(4.d0,MIN(Tmin,Twat_mod(j+1)))
            
END DO

delta(n_tot) = DD

RETURN
END
!-------------------------------------------------------------------------------
!				NUMERICAL INTEGRATION
!-------------------------------------------------------------------------------
SUBROUTINE a2w(Ta, Tw, time, K, DD, pp)

USE commondata

IMPLICIT NONE

REAL(KIND=8) :: lambda
REAL(KIND=8), INTENT(OUT) :: K, DD, pp
REAL(KIND=8), INTENT(IN) :: Ta, Tw, time


IF (Tw>=Tmin) THEN
	DD=DEXP( -(Tw-Tmin)/par(4) );	
ELSE
	IF (version=='3') THEN
		DD=DEXP( (Tw-Tmin)/par(7) ) + DEXP( -Tw/par(8) );	
    ELSE
		DD=1.0d0
	END IF
END IF

pp = par(1) + par(2)*Ta -  par(3)*Tw + par(5)*COS(2.d0*pi*(time-par(6))) ! Note that if version == '1' par(5)==0

!!! lower bound for numerical stability
!!lambda=pp/par(4) - par(3)
!!IF (lambda .le. 0.0d0) THEN
!!    IF (num_mod .eq. 'RK4' .and. DD .lt. -lambda/LIM) THEN
!!        DD=-lambda/LIM
!!    ELSEIF ((num_mod .eq. 'RK2' .or. num_mod .eq. 'EUL') .and. DD .lt. -lambda/LIM ) THEN
!!        DD=-lambda/LIM
!!    END IF
!!END IF
    
K = pp/DD    
       
RETURN
END SUBROUTINE

!-------------------------------------------------------------------------------
!				OBJECTIVE FUNCTION
!-------------------------------------------------------------------------------
SUBROUTINE funcobj(ind)
! Subroutine per il calcolo della funzione obiettivo della simulazione
! Dichiarazione delle variabili
USE commondata
IMPLICIT NONE

INTEGER :: i, j
REAL(KIND=8),INTENT(OUT) :: ind
REAL(KIND=8) :: TSS, mean_mod, TSS_mod, std_mod, covar_mod
REAL(KIND=8) :: tmp, max_, max_err

max_err=-9999
Twat_mod_agg=-999
DO i=1,n_dat            ! 1st year = warm-up period
    tmp=0.0d0
    DO j=I_inf(i,1),I_inf(i,2)
        tmp=tmp+Twat_mod(I_pos(j)) 
    END DO
    Twat_mod_agg(I_inf(i,3))=tmp/REAL(I_inf(i,2)-I_inf(i,1)+1)    
END DO
    
IF (fun_obj=='NSE') THEN
    TSS=0.d0
    DO i=1,n_dat		
	    TSS=TSS+(Twat_mod_agg(I_inf(i,3))-Twat_obs_agg(I_inf(i,3)))**2
    END DO
    ind=1.d0-TSS/TSS_obs
ELSEIF  (fun_obj=='KGE') THEN
    mean_mod=0.0d0
    DO i=1,n_dat		
	    mean_mod=mean_mod+Twat_mod_agg(I_inf(i,3))
    END DO
    mean_mod=mean_mod/REAL(n_dat)
    
    covar_mod=0.0d0
    TSS_mod=0.0d0
    DO i=1,n_dat
	    TSS_mod=TSS_mod+(Twat_mod_agg(I_inf(i,3))-mean_mod)**2
	    covar_mod=covar_mod+(Twat_mod_agg(I_inf(i,3))-mean_mod)*(Twat_obs_agg(I_inf(i,3))-mean_obs)
    END DO
    std_mod=DSQRT(TSS_mod/REAL(n_dat-1))
    covar_mod=covar_mod/REAL(n_dat-1)
    ind=1.0d0-DSQRT((std_mod/std_obs-1.0d0)**2 + (mean_mod/mean_obs-1.0d0)**2 + (covar_mod/(std_mod*std_obs)-1.0d0)**2)
ELSEIF  (fun_obj=='RMS') THEN
    TSS=0.d0
    DO i=1,n_dat
	    TSS=TSS+(Twat_mod_agg(I_inf(i,3))-Twat_obs_agg(I_inf(i,3)))**2
    END DO
    ind=-DSQRT(TSS/n_dat)        ! sign - --> the calibration procedure maximizes ind.    
ELSEIF  (fun_obj=='XXX') THEN
	WRITE(*,*) 'XXX non definito'
ELSE
	WRITE(*,*) 'Errore nella scelta della f. obiettivo'
END IF

RETURN
END
!-------------------------------------------------------------------------------
!				AGGREGATION
!-------------------------------------------------------------------------------
SUBROUTINE aggregation

! Dichiarazione delle variabili
USE commondata
IMPLICIT NONE

INTEGER :: i, j, k, status, pp
INTEGER :: n_inf, n_pos, n_units, n_days, count, pos_tmp
INTEGER :: month, month_curr
INTEGER, DIMENSION(:,:), ALLOCATABLE :: A
INTEGER, DIMENSION(:), ALLOCATABLE :: B
REAL(KIND=8) :: tmp

pp=LEN_TRIM(time_res)
IF (pp==2) THEN
    unit=time_res(2:)
    READ(time_res,'(i1)') qty
ELSEIF (pp==3) THEN
    unit=time_res(3:)
    READ(time_res,'(i2)') qty
END IF

ALLOCATE(I_pos(n_tot),stat=status)
I_pos=-999
Twat_obs_agg=-999
n_inf=1   
n_pos=1
IF (time_res=='1d') THEN    ! daily resolution
    n_units=n_tot-365
    ALLOCATE(I_inf(n_units,3))
    I_inf=-999
    DO i=366,n_tot    ! 1st year = warm-up period
        IF (Twat_obs(i) .ne. -999) THEN
            I_inf(n_inf,2)=n_pos
            I_inf(n_inf,3)=i
            I_pos(n_pos)=i
            Twat_obs_agg(I_inf(n_inf,3))=Twat_obs(i)  
            n_inf=n_inf+1
            n_pos=n_pos+1
        END IF  
    END DO
ELSEIF (unit=='w') THEN      ! weekly resolution
    n_days=qty*7      
    n_units=CEILING(REAL(n_tot-365)/REAL(n_days))
    ALLOCATE(I_inf(n_units,3))
    I_inf=-999
    DO i=366,n_tot,n_days
        tmp=0.0d0
        count=0
        pos_tmp=i+CEILING(0.5*n_days)-1
        DO j=0,n_days-1
            k=i+j
            IF (k .gt. n_tot) THEN
                EXIT
            END IF
            IF (Twat_obs(k) .ne. -999) THEN
                tmp=tmp+Twat_obs(k)
                I_pos(n_pos)=k
                n_pos=n_pos+1
                count=count+1
            END IF
        END DO
        IF (count .ge. n_days*prc) THEN           
            I_inf(n_inf,2)=n_pos-1
            I_inf(n_inf,3)=pos_tmp  
            Twat_obs_agg(I_inf(n_inf,3))=tmp/count
            n_inf=n_inf+1
        ELSE
            I_pos(n_pos-count:n_pos)=-999
            n_pos=n_pos-count
        END IF
    END DO
ELSEIF (unit=='m') THEN      ! monthly resolution
    n_units=CEILING(n_tot/(30.5))
    ALLOCATE(I_inf(n_units,3))
    I_inf=-999
    n_days=0
    month_curr=-999    
    count=0
    DO i=366,n_tot
        month=date(i,2)
        IF (month .ne. month_curr) THEN
            IF (count .ge. n_days*prc .and. i .ne. 366) THEN  
                I_inf(n_inf,2)=n_pos-1
                I_inf(n_inf,3)=i-FLOOR(0.5*n_days)-1
                Twat_obs_agg(I_inf(n_inf,3))=tmp/count
                n_inf=n_inf+1
            ELSE
                I_pos(n_pos-count:n_pos)=-999
                n_pos=n_pos-count
            END IF        
            month_curr=month
            count=0           
            n_days=1
            tmp=0.0d0            
        ELSE
            n_days=n_days+1
        END IF  
        IF (Twat_obs(i) .ne. -999) THEN
            tmp=tmp+Twat_obs(i)
            I_pos(n_pos)=i
            n_pos=n_pos+1
            count=count+1
        END IF
    END DO
    ! Last month
    IF (count .ge. n_days*prc) THEN  
        I_inf(n_inf,2)=n_pos-1
        I_inf(n_inf,3)=i-FLOOR(0.5*n_days)-1
        Twat_obs_agg(I_inf(n_inf,3))=tmp/count
        n_inf=n_inf+1
    ELSE
        I_pos(n_pos-count:n_pos)=-999
        n_pos=n_pos-count
    END IF                    
ELSE
    WRITE(*,*) 'Error: variable time_res'
END IF

n_dat=n_inf-1
n_pos=n_pos-1

I_inf(1,1)=1
I_inf(2:n_dat,1)=I_inf(1:n_dat-1,2)+1
ALLOCATE(A(n_units,3),B(n_tot))
A=I_inf;    
B=I_pos
DEALLOCATE(I_inf,I_pos)
ALLOCATE(I_inf(n_dat,3),I_pos(n_pos))
I_inf=A(1:n_dat,:)
I_pos=B(1:n_pos)
DEALLOCATE(A,B)

RETURN
END
!-------------------------------------------------------------------------------
!				STATIS
!-------------------------------------------------------------------------------
SUBROUTINE statis
! Subroutine per il calcolo di media, somma degli scarti quadratici e std dei dati

! Dichiarazione delle variabili
USE commondata
IMPLICIT NONE

INTEGER :: i, k, d, status

mean_obs=0.d0
TSS_obs=0.d0
DO i=1,n_dat		    
	mean_obs=mean_obs+Twat_obs_agg(I_inf(i,3))
END DO
mean_obs=mean_obs/REAL(n_dat)

DO i=1,n_dat		
	TSS_obs=TSS_obs+(Twat_obs_agg(I_inf(i,3))-mean_obs)**2
END DO

std_obs=DSQRT(TSS_obs/REAL(n_dat-1))

RETURN
END

!-------------------------------------------------------------------------------
!				FORWARD
!-------------------------------------------------------------------------------
SUBROUTINE forward
USE commondata
IMPLICIT NONE

INTEGER :: i, j
REAL(KIND=8) :: ei_check,ei

par=par_best		! uso miglior set di paramteri

CALL model

! Controllo: calcolo dell'indice di efficienza
CALL funcobj(ei_check)
WRITE(*,*) 'Calibration: objective function', ABS(ei_check)

IF (ABS(ei_check - finalfit) .gt. 0.0001) THEN
	    WRITE(*,*) 'Errore efficienza in forward'
	    WRITE(*,*) ei_check, finalfit
	    PAUSE
ELSE
    WRITE(*,*) 'Check completed'
END IF
    
WRITE(11,'(<n_par>(f10.6,1x))') (par_best(i),i=1,n_par)
WRITE(11,'(f10.6)') ABS(ei_check)
    	
OPEN(UNIT=12,FILE=TRIM(folder)//'/2_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'c_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
DO i=1,n_tot
	WRITE(12,1004) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Twat_mod(i),Twat_obs_agg(i),Twat_mod_agg(i)
END DO
CLOSE(12)

OPEN(UNIT=14,FILE=TRIM(folder)//'/4_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'c_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
DO i=1,n_tot
    WRITE(14,*) delta(i)
END DO
CLOSE(14)

CALL read_validation

IF (n_tot .lt. 365) THEN
    ei=-999
    GO TO 200
END IF

! aggregate calibration data on the basis of time_res
CALL aggregation

CALL statis
WRITE(*,*) 'mean, TSS and standard deviation (validation)'
WRITE(*,*)  SNGL(mean_obs),SNGL(TSS_obs),SNGL(std_obs)

CALL model
CALL funcobj(ei)
WRITE(11,'(f10.6)') ABS(ei)

CLOSE(11)

WRITE(*,*) 'Validation: objective function', ABS(ei)

OPEN(UNIT=13,FILE=TRIM(folder)//'/3_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'v_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
DO i=1,n_tot
	WRITE(13,1004) (date(i,j),j=1,3),Tair(i),Twat_obs(i),Twat_mod(i),Twat_obs_agg(i),Twat_mod_agg(i)
END DO
CLOSE(13)

OPEN(UNIT=15,FILE=TRIM(folder)//'/5_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'v_'//TRIM(time_res)//'.out',STATUS='unknown',ACTION='write')
DO i=1,n_tot
    WRITE(15,*) delta(i)
END DO
CLOSE(15)  
    
1004 FORMAT(i4,1x,i4,1x,i4,1x,5(1x,f10.5))

200 RETURN
END

!-------------------------------------------------------------------------------
!				BEST
!-------------------------------------------------------------------------------
SUBROUTINE best(fit,part,foptim)
USE commondata
IMPLICIT NONE

INTEGER,INTENT(OUT)::part
INTEGER:: k
REAL(KIND=8),INTENT(IN),DIMENSION(n_particles):: fit
REAL(KIND=8),INTENT(OUT):: foptim

foptim=-1e30				! valore molto piccolo
DO k=1,n_particles
    IF(fit(k).gt.foptim) then
	    foptim=fit(k)
	    part=k
    END IF
END DO

RETURN
END

!-------------------------------------------------------------------------------
!				LEAP YEAR
!-------------------------------------------------------------------------------
SUBROUTINE leap_year(Y,I)
USE commondata
IMPLICIT NONE

INTEGER,INTENT(IN) :: Y
INTEGER,INTENT(OUT) :: I

IF(MOD(Y,100).NE.0.AND.MOD(Y,4).EQ.0) THEN 
    I=1
ELSEIF(MOD(Y,400).EQ.0) THEN
    I=1
ELSE 
    I=0
END IF

RETURN
END


!-------------------------------------------------------------------------------
!				LINEAR REGRESSION
!-------------------------------------------------------------------------------
SUBROUTINE linreg(n,X,Y,m,b,r2)
IMPLICIT NONE

! adapted from http://www.pgccphy.net/Linreg/linreg_f90.txt
! Dr. David G. Simpson

INTEGER :: i
INTEGER, INTENT(IN):: n
REAL(KIND=8), INTENT(IN):: X(n),Y(n)
REAL(KIND=8), INTENT(OUT):: m, b, r2
REAL(KIND=8):: sumx, sumx2, sumxy, sumy, sumy2

sumx=0.0d0; sumx2=0.0d0; sumxy=0.0d0; sumy=0.0d0; sumy2=0.0d0;
DO i=1,n
    sumx  = sumx + x(i)                     ! compute sum of x
    sumx2 = sumx2 + x(i) * x(i)             ! compute sum of x**2
    sumxy = sumxy + x(i) * y(i)             ! compute sum of x * y
    sumy  = sumy + y(i)                     ! compute sum of y
    sumy2 = sumy2 + y(i) * y(i)             ! compute sum of y**2
END DO

m = (n * sumxy  -  sumx * sumy) / (n * sumx2 - sumx**2)                          ! compute slope
b = (sumy * sumx2  -  sumx * sumxy) / (n * sumx2  -  sumx**2)                    ! compute y-intercept
r2 = (sumxy - sumx * sumy / n) /                                     &            ! compute correlation coefficient
                 sqrt((sumx2 - sumx**2/n) * (sumy2 - sumy**2/n))
r2=r2**2
RETURN
END
 