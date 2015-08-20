MODULE commondata
IMPLICIT NONE
SAVE 

INTEGER, PARAMETER :: n_par = 8
REAL(KIND=8), PARAMETER :: pi = ACOS(0.d0)*2.d0
REAL(KIND=8), PARAMETER :: ttt = 1.0d0/365.0d0

INTEGER :: n_tot, n_dat
INTEGER :: version, qty
INTEGER :: n_year
INTEGER, ALLOCATABLE, DIMENSION(:) :: I_pos
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I_inf
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: date
REAL(KIND=8) :: Tice_cover, prc
REAL(KIND=8) :: mean_obs, TSS_obs, std_obs
REAL(KIND=8) :: Tmin
REAL(KIND=8) :: mineff_index,finalfit
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: tt, Tair, Twat_obs_agg, Twat_obs
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: Twat_mod, Twat_mod_agg, delta

REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: parmin
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: parmax
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: par, par_best
REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: m, q, r2
LOGICAL,ALLOCATABLE,DIMENSION(:) :: flag_par

CHARACTER(len=100) :: folder
CHARACTER(LEN=30) :: name, air_station, water_station, station, run
CHARACTER(LEN=1) :: series, unit
CHARACTER(LEN=3) :: time_res
CHARACTER(LEN=3) :: fun_obj, num_mod


INTEGER :: n_run,n_particles      ! PSO
REAL(KIND=8) ::c1,c2,wmin,wmax    ! PSO

END MODULE commondata

