SUBROUTINE forward_mode

USE commondata

IMPLICIT NONE 

REAL(KIND=8):: eff_index

CALL sub_1(eff_index)

par_best=par
finalfit=eff_index

RETURN
END SUBROUTINE


!-------------------------------------------------------------------------------
!				PSO
!-------------------------------------------------------------------------------


SUBROUTINE PSO_mode

USE commondata

IMPLICIT NONE

INTEGER :: i, j, k, status, count
REAL(KIND=8):: eff_index

!!! PSO parameters
REAL(KIND=8) ::  norm
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)::x,v,pbest
REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::fit,r,gbest,fitbest
REAL(KIND=8)::w,dw,dxmax,dvmax,foptim

WRITE(*,*) 'N. particles = ',n_particles,', N. run = ', n_run
! allocazione + inizializzazione per PSO
ALLOCATE (x(n_par,n_particles),v(n_par,n_particles),pbest(n_par,n_particles))
ALLOCATE (r(2*n_par),gbest(n_par),fit(n_particles),fitbest(n_particles))

! open file for the writing of all parameter set + efficiency index
OPEN(unit=10,file=TRIM(folder)//'/0_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'_'//TRIM(time_res)//'.out',status='unknown',action='write',form='binary')
 
x=0; v=0;   ! e velocità
r=0
pbest=0; gbest=0; 
fit=0; fitbest=0
dw=(wmax-wmin)/n_run
w=wmax

! generazione random set parametri 1° giro
CALL random_seed()
CALL random_number(x)
CALL random_number(v)
DO j=1,n_par
    DO k=1,n_particles
        dxmax=(parmax(j)-parmin(j))     ! range di ciascun parametro
        dvmax=1.0*dxmax                 ! massima velocità particella
        x(j,k)=x(j,k)*dxmax+parmin(j)   ! valore random del parametro per il parametro j per la particella k
        v(j,k)=v(j,k)*dvmax             ! valore random di velocità per il parametro j per la particella k
        pbest(j,k)=x(j,k)               ! inizializzo partial best
    END DO
END DO
DO k=1,n_particles
    par=x(:,k)                          ! assegno set di parametri della particella k-esima
    CALL sub_1(eff_index)        
    fitbest(k)=eff_index                ! assegno valore best fit (miglior indice efficienza)
END DO
CALL best(fit,k,foptim)                 ! k è in output=particella più performante (indice performance=foptim)
DO j=1,n_par
    gbest(j)=x(j,k)                     ! global best
END DO

! ciclo n_run
DO i=1,n_run
    CALL random_seed()      
    DO k=1,n_particles
	    CALL random_number(r)
	    status=0
        !update the velocity and the position of the particles
        DO j=1,n_par
	        v(j,k)=w*v(j,k)+c1*r(j)*(pbest(j,k)-x(j,k))+c2*r(j+n_par)*(gbest(j)-x(j,k))
            x(j,k)=x(j,k)+v(j,k)
		
	        ! muro assorbente
	        IF (x(j,k).gt.parmax(j)) THEN
	            x(j,k)= parmax(j)
	            v(j,k)=0.0          
                status=1            
	        END IF
            IF (x(j,k).lt.parmin(j)) THEN
	            x(j,k)= parmin(j)
	            v(j,k)=0.0          
                status=1            
            END IF
        END DO

        ! Calcolo nuova performance
        IF (status.eq.0) THEN	
            par=x(:,k)                          ! assegno set di parametri della particella k-esima	
            CALL sub_1(eff_index) 
		    fit(k)=eff_index
		    ! Scrittura su file se eff_index maggiore di soglia indbound
            IF (eff_index .ge. mineff_index) THEN
                WRITE(10)(x(j,k),j=1,n_par),eff_index
            END IF
        ELSE
	        fit(k)=-1e30
        ENDIF
        
        ! Valuto se la particella k ha migliorato la sua performance
        IF (fit(k).gt.fitbest(k)) THEN
		    fitbest(k)=fit(k)
		    DO j=1,n_par
			    pbest(j,k)=x(j,k)
		    END DO
	    END IF
    END DO      !ciclo k (particelle)
    
    ! valuto qual'è la particella più performante
    CALL best(fitbest,k,foptim)
    DO j=1,n_par
	    gbest(j)=pbest(j,k)
    END DO

    w=w-dw

    IF (i>=10) THEN
	    IF (MOD(i,INT(REAL(n_run)/10))==0 ) THEN
		    WRITE(*,1003) REAL(i)/REAL(n_run)*100.
	    END IF
    END IF
    
    ! If the norm between pbest and gbest is less then tol for the perc percentage of the particles --> exit the cycle
    count=0
    DO k=1,n_particles
        norm=0.
        DO j=1,n_par
            IF (flag_par(j)) THEN
                norm=norm+( (pbest(j,k)-gbest(j))/(parmax(j)-parmin(j)) )**2
            END IF    
        END DO
        norm=SQRT(norm)
        IF (norm .lt. 0.01) THEN           ! 0.001 --> 1 °/oo
            count=count+1
        END IF
    END DO
    IF (count .ge. (0.9*n_particles)) THEN           ! 1 --> 100 °/o of npc
        WRITE(*,*)  '- Warning:  PSO has been exit'
        EXIT
    END IF
	
END DO

par_best=gbest
finalfit=foptim
   
WRITE(*,*) 'Indice efficienza calibrazione', finalfit

!    CLOSE(10)

1003 FORMAT('Calcolo al', 1x, f5.1 ,1x,'%')

RETURN
END SUBROUTINE


!-------------------------------------------------------------------------------
!				Latin Hypercube
!-------------------------------------------------------------------------------


SUBROUTINE LH_mode

USE commondata

IMPLICIT NONE

INTEGER :: i, j
REAL(KIND=8):: fit, foptim, eff_index
INTEGER, DIMENSION(n_run,n_par):: permut
REAL(KIND=8):: r

!!! LH parameters
REAL(KIND=8),ALLOCATABLE,DIMENSION(:)::gbest

WRITE(*,*) 'N. run = ', n_run

ALLOCATE(gbest(n_par))
  
foptim=-999

! open file for the writing of all parameter set + efficiency index
OPEN(unit=10,file=TRIM(folder)//'/0_'//TRIM(run)//'_'//fun_obj//'_'//TRIM(station)//'_'//series//'_'//TRIM(time_res)//'.out',status='unknown',action='write',form='binary')

CALL random_seed()

! Initialization of matrix permut + permutation (shuffle)
DO j=1,n_par
    permut(:,j)= (/ (i, i=1,n_run) /) 
    CALL Shuffle(permut(:,j),n_run) 
END DO     

DO i=1,n_run
    DO j=1,n_par
        CALL random_number(r)
        r=r + (REAL(permut(i,j))-1.0)
        r=r/REAL(n_run) 

        par(j)=parmin(j) + (parmax(j)-parmin(j))*r
    END DO
        
    CALL sub_1(eff_index)
    fit=eff_index
    
    ! Scrittura su file se eff_index maggiore di soglia indbound
    IF (eff_index .ge. mineff_index) THEN
        WRITE(10)(par(j),j=1,n_par),eff_index
    END IF
                
    IF (fit .gt. foptim) THEN
        foptim=fit
        gbest=par
    END IF        
        
    IF (i>=10) THEN
	    IF (MOD(i,INT(REAL(n_run)/10))==0 ) THEN
		    WRITE(*,1003) REAL(i)/REAL(n_run)*100.
	    END IF
    END IF
           
END DO            

par_best=gbest
finalfit=foptim

WRITE(*,*) 'Indice efficienza calibrazione', finalfit

!    CLOSE(10)

1003 FORMAT('Calcolo al', 1x, f5.1 ,1x,'%')

RETURN
END SUBROUTINE


!----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------

 
SUBROUTINE Shuffle(a,n)

IMPLICIT NONE
INTEGER, INTENT(IN) :: n
INTEGER, DIMENSION(n), INTENT(INOUT) ::  a
INTEGER :: i, randpos, temp
REAL :: r

DO i = n, 2, -1
    CALL random_number(r)
    randpos = INT(r * i) + 1
    temp = a(randpos)
    a(randpos) = a(i)
    a(i) = temp
END DO

END SUBROUTINE Shuffle


!----------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------

