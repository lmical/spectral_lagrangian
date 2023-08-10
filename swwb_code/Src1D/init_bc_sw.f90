MODULE init_bc

  USE param2d
  USE overloading
  USE precision
  USE utils
  USE preprocessing


  IMPLICIT NONE

CONTAINS

  !---------------------------------------
  ! Setup domain
  !---------------------------------------
  SUBROUTINE init_geom(DATA)
    TYPE(donnees), INTENT(inout):: DATA

    SELECT CASE(DATA%test)
    CASE(-2)
       !*Solitary wave
       DATA%domain_left = -100_dp
       DATA%Length      =  200_dp
    CASE(0,40,41,42)
       ! Isentropic
       DATA%domain_left = -1.0_dp
       DATA%Length      =  2.0_dp

    CASE(1,11,12,13,32,500)
       ! Sod
       DATA%domain_left = 0.0_dp
       DATA%Length      = 1.0_dp
       
    CASE(2)
       ! Shu-Osher
       DATA%domain_left = -5.0_dp
       DATA%Length      = 10.0_dp

    CASE(3)
       ! Woodward-Colella
       DATA%domain_left = 0.0_dp
       DATA%Length      = 1.0_dp

    CASE(4)
       ! channel
       DATA%domain_left = 0.0_dp
       DATA%Length      = 1.0_dp

    CASE(5)
       ! Isentropic vortex
       DATA%domain_left = -1.0_dp
       DATA%Length      =  2.0_dp

!       DATA%domain_left = 0.0_dp
!       DATA%Length      = 1.0_dp

    CASE(6)
       ! Woodward-Colella left
       DATA%domain_left = 0.0_dp
       DATA%Length      = 1.0_dp
       ! LeBlanc
    CASE(21, 20,24,25,19,14,26,27) !convection param = 2. for 20/21, conv = 5.5 for 24,conv = 3 for 25
       !simetrie
      DATA%domain_left=0.0_dp
      DATA%Length  = 25._dp
       ! 1-2-3
    CASE(22,23) !convection bigger 8.
      DATA%domain_left=0.0_dp
      DATA%Length  = 25._dp

    CASE(29,30) !convection bigger than 4.
      DATA%domain_left=0.0_dp
      DATA%Length  = 4._dp

    CASE(31) !convection bigger than 4.
      DATA%domain_left=0.0_dp
      DATA%Length  = 15._dp
    CASE(501) !*bump
       DATA%domain_left = -2.0_dp
       DATA%Length      =  4.0_dp
       !*A New Approach for Designing Moving-Water Equilibria Preserving Schemes for the Shallow Water Equations - Yuanzhen Cheng, Alina Chertock, Michael Herty, Alexander Kurganov, Tong Wu, section 3 example 2
       !*100 supercritical transitory
       !*101 subcritical transitory
       !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 3.1.3 subcritical
       !*111 subcritical but with IC taken from scratch
       !*120 supercritical with IC taken from scratch but smooth
       !*121 subcritical with IC taken from scratch but smooth
       !*122 transcritical with IC taken from scratch but smooth

    CASE(100,101,111,120,121,122) 
       DATA%domain_left = 0.0_dp
       DATA%Length      =  25.0_dp


    !*FOR CONVERGENCE ANALYSIS 
    !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
    !*3.2.2 Short channel: 100 m
    !*150 Supercritical case 
    CASE(150) 
       DATA%domain_left = 0.0_dp
       DATA%Length      =  100.0_dp





    !*200 lake at rest - periodic BCs (perturbed with 10**(-2) and 10**(-4))
    !*201 discrete equilibrium of 200 from scratch 
    !*210 lake at rest - outflow BCs
    !*211 discrete equilibrium of 210 from scratch 
    CASE(200,201,210,211) 
       DATA%domain_left = 0.0_dp
       DATA%Length      =  25.0_dp

    !*230 lake at rest with C^0 bathymetry - outflow BCs
    !*231 discrete equilibrium of 230 from scratch 
    CASE(230,231) 
       DATA%domain_left = 0.0_dp
       DATA%Length      =  25.0_dp

    !*Well-Balancing Via Flux Globalization: Applications to Shallow Water Equations with Wet/Dry Fronts
    !*Alina Chertock, Alexander Kurganov, Xin Liu, Yongle Liu, and Tong Wu
    CASE(300)
       DATA%domain_left = -1.0_dp
       DATA%Length      =  2.0_dp


    !*NB: THE FINAL TIME IN THE NEXT TESTS IS GIVEN HERE AND NOT IN IC
    !*NBB: IT IS IN CASE MODIFIED WHEN I INTRODUCE THE PERTURBATION
    !*SWASHES WITH BATHYMETRY NON SMOOTH
    !*400 Lake at rest with an immersed bump 3.1.1
    !*410 Supercritical
    !*420 Subcritical
    !*430 Transcritical
    CASE(400,401,410,411,415,416,420,421,425,426,430,431,435)
       DATA%domain_left = 0._dp
       DATA%Length      =  25.0_dp
       DATA%tmax=10._DP**(20)
       DATA%ktmax=50000 !*50000
       IF (DATA%test==415 .OR. DATA%test==425) THEN !*Many many many iterations to reach the steady state
           DATA%ktmax=1000000 
       END IF
       IF (DATA%test==400) THEN !*10 seconds
           DATA%tmax=10._DP 
       END IF

    !*IC taken from scratch
    !*401 Lake at rest with an immersed bump 3.1.1 
    !*411 Supercritical
    !*421 Subcritical
    !*431 Transcritical


    !*Same tests with friction
    !*415 Supercritical for the generation of the steady state with friction
    !*425 Subcritical for the generation of the steady state with friction
    !*435 Transcritical for the generation of the steady state with friction !*ALERT, DO NOT RELY


    !*IC taken from scratch with friction
    !*416 Supercritical
    !*426 Subcritical


    !*NB: THE FINAL TIME IN THE NEXT TESTS IS GIVEN HERE AND NOT IN IC
    !*NB: THIS IS LIKE 150 BUT FOR PERTURBATION ANALYSIS
    !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
    !*3.2.2 Short channel: 100 m
    !*440 Supercritical case 
    CASE(440,441) 
       DATA%domain_left = 0.0_dp
       DATA%Length      =  100.0_dp
       DATA%tmax=10._DP**(20)
       DATA%ktmax=100000 !*50000
       
		
    CASE default
       PRINT*, "Wrong test number for geom(), test = ", DATA%test
       STOP
    END SELECT

  END SUBROUTINE init_geom

  !---------------------------------------
  ! Set initial and boundary conditions
  !---------------------------------------
  FUNCTION IC(x,DATA) RESULT(Var)
    REAL(dp),       INTENT(in)   :: x
    TYPE(donnees), INTENT(inout):: DATA
    !*REAL(dp), DIMENSION(n_dim) :: eta !*NOT USED
    TYPE(PVar) :: Var
    REAL(dp)    :: alpha,r,x0
    REAL(DP) :: A, k, eta, h_ref, u


    !*eta = 0._dp ! profile_bathymetry( DATA%test, (/x/) ) !*NOT USED

    !---------------
    ! for Euler
    SELECT CASE(DATA%test)
    CASE(-2)
        !*Solitary wave
        DATA%tmax = 100._DP
        h_ref=10._DP
        A = 0.35_DP
        k=SQRT(3._DP*A/(4._DP*h_ref**2*(A+h_ref)))
        eta=A*(1._DP/COSH(   k*x   ))**2
        u=eta*SQRT(grav*(A+h_ref))/(h_ref+eta)
        Var%u(1) = eta+h_ref
        Var%u(2) = Var%u(1)*u
    CASE(0)
        ! Convergence test: isentropic flow
        DATA%tmax = 1._DP
        alpha = 0.1_DP
        Var%u(1) = 1._DP + alpha*SIN(PI*x) ! + 1._DP 
        Var%u(2) = 0._DP

    CASE(40,41)
        ! Convergence test: isentropic flow
        DATA%tmax = 1._DP
        alpha = 0.1_DP
        Var%u(1) = 1._DP + alpha*SIN(PI*x) - bathymetry(DATA%test,x) ! + 1._DP 
        Var%u(2) = 0.05_DP

    CASE(42)
        ! Convergence test: isentropic flow
        DATA%tmax = 1.0_dp
        alpha = 0.01_dp
        Var%u(1) = 1._dp + alpha*SIN(PI*x) - bathymetry(DATA%test,x) ! + 1. 
        Var%u(2) = 0._dp

    CASE(1)
        ! Sod
        DATA%tmax = 0.13_DP
        IF (x < 0.5_DP) THEN
            Var%u(1) = 10._DP-bathymetry(DATA%test,x) ! [kg/m^3]
            Var%u(2) = 0._DP ! [m/s] !*NEIN, [kg/m^3*m/s]
        ELSE
            Var%u(1) = 5._DP -bathymetry(DATA%test,x)
            Var%u(2) = 0._DP
        END IF

    CASE(32)
        ! Dam break
        DATA%tmax = 0.05_DP
        IF (x < 0.5_DP) THEN
            Var%u(1) = 1._DP-bathymetry(DATA%test,x) ! [kg/m^3]
            Var%u(2) = 0._DP ! [m/s] !*NEIN, [kg/m^3*m/s]
        ELSE
            Var%u(1) = 0.5_DP  !*IT WAS 0._DP 
            Var%u(2) = 0._DP
        END IF


    CASE(11)
        ! Sod
        DATA%tmax = 5.3_DP
        Var%u(1) = 0.01_DP+ 0.3_DP*EXP(-100._DP* (x-0.7_DP)**2) -bathymetry(DATA%test,x)
        Var%U(2) = 0._DP*Var%U(1)

    CASE(12,13)
      DATA%tmax = 1._DP
      Var%U(1) = 1._DP- bathymetry(DATA%test,x)
      Var%U(2) = 0._DP
    CASE(500,501)
      DATA%tmax = 100._DP
      Var%U(1) = 1._DP- bathymetry(DATA%test,x)
      Var%U(2) = 0._DP
    CASE(21)
      DATA%tmax = 5._DP
      Var%U(1) = 0.5_DP -bathymetry(DATA%test,x)
      Var%U(2) = 0._DP

    CASE(20)
      DATA%tmax = 20._DP
      Var%U(1) = 0.5_DP -bathymetry(DATA%test,x)
      Var%U(2) = 0._DP
      IF (x<18._DP .AND. x>16._DP) THEN
        Var%U(1) = Var%U(1) + 0.1_DP
      END IF

    CASE(14)
      DATA%tmax = 60._DP
      Var%U(1) = 0.205_DP -bathymetry(DATA%test,x)
      Var%U(2) = 1._DP*Var%U(1)

    CASE(19)
      DATA%tmax = 12.5_DP
      Var%U(1) = 0.5_DP -bathymetry(DATA%test,x)
      Var%U(2) = -1._DP
      IF (x<18._DP .AND. x>16._DP) THEN
        Var%U(1) = Var%U(1) + 0.05_DP*EXP(1._dp-1._dp/(1._DP-(x-17._DP)**2))
      END IF

    CASE(22)
      DATA%tmax = 100._DP
      Var%U(1) = 2._DP -bathymetry(DATA%test,x)
      Var%U(2) = 0._DP
    CASE(23)
      DATA%tmax = 100._DP
      Var%U(1) = 2._DP -1.5_DP*bathymetry(DATA%test,x)
      Var%U(2) = 4.42_dp

    CASE(24)
      DATA%tmax = 60._DP
      IF (x<5._DP) THEN
        Var%U(1) = 1.01445_DP -bathymetry(DATA%test,x)
        Var%U(2) = 1.53_dp
      ELSE IF (x<15._DP) THEN
        Var%U(1) = (1.01445_DP)*(15._DP-x)/10._DP + 0.405781_DP*(x-5._DP)/10._DP -bathymetry(DATA%test,x)
        Var%U(2) = 1.53_dp
      ELSE
        Var%U(1) = 0.405781_DP -bathymetry(DATA%test,x)
        Var%U(2) = 1.53_dp
      END IF

      ! Var%U(1) = 0.66_DP -bathymetry(DATA%test,x)
      ! Var%U(2) = 0._DP

    CASE(25)
      DATA%tmax = 100._DP
      IF (x<8._DP) THEN
        Var%U(1) = 0.4_DP -bathymetry(DATA%test,x)
      ELSE
        Var%U(1) = 0.33_DP -bathymetry(DATA%test,x)
      END IF
      Var%U(2) = 0.14_DP

    CASE(26)
      DATA%tmax = 60._DP
      Var%U(1) = 0.3_DP -bathymetry(DATA%test,x)
      Var%U(2) = 0.2_DP*Var%U(1)

   CASE(30) ! thucker oscillations in parabola
       DATA%tmax = 10.0303_dp
       IF (x<2.5_DP .AND. x>0.5_DP) THEN
          Var%U(1) = -0.5_DP*((x-1.5_DP)**2-1._DP)
       ELSE
          Var%U(1) = 0._dp
       END IF
       Var%U(2) = 0._dp
   CASE(29) !lake at rest in parabola
       DATA%tmax = 3._dp
       Var%U(1) = MAX(0._dp, 0.5_DP-bathymetry(DATA%test,x))
       Var%U(2) = 0._dp

    CASE(31) !waves on the shore
        DATA%tmax = 20.0_dp
        Var%U(1) = MAX(0.37_DP-bathymetry(DATA%test,x),0._dp)
        Var%U(2) = 0._dp

    CASE(2)
        ! Shu-Osher
        DATA%tmax = 1.8_DP
        alpha = 0.2_DP
        IF (x <= -4._DP) THEN
            Var%u(1) = 3.857143_DP ! [kg/m^3]
            Var%u(2) = Var%u(1)*2.629369_DP ! [m/s] !*NEIN, [kg/m^3*m/s]
        ELSE
            Var%u(1) = 1._DP + alpha*SIN(5._DP*x)
            Var%u(2) = 0._DP
        END IF

    CASE(3)
        ! Woodward-Colella
        DATA%tmax = 0.038_DP
        IF (x <= 0.1_DP) THEN
            Var%u(1) = 1._DP ! [kg/m^3]
            Var%u(2) = 0._DP ! [m/s] !*NEIN, [kg/m^3*m/s]
        ELSE
             IF (x > 0.1_DP .AND. x < 0.9_DP) THEN
                Var%u(1) = 1._DP ! [kg/m^3]
                Var%u(2) = 0._DP ! [m/s] !*NEIN, [kg/m^3*m/s]
             ELSE
                Var%u(1) = 1._DP ! [kg/m^3]
                Var%u(2) = 0._DP ! [m/s] !*NEIN, [kg/m^3*m/s]
             END IF
        END IF

    CASE(4)
        ! 1D version of 2D channel
        Var%u(1) = 0.5_dp ! [kg/m^3]
        Var%u(2) = 0.0_dp ! [m/s] !*NEIN, [kg/m^3*m/s]

    CASE(5)
        ! Convergence test: isentropic vortex
        DATA%tmax = 0.5_DP
!        DATA%tmax = 0.25_DP
        Var%u(1) = 1._DP + 0.1_DP*EXP(-(x)**2/(2._DP*0.01_DP))
!        Var%u(1) = 1._DP+EXP(-80.0_dp*(x-0.4_dp)**2)
        Var%u(2) = Var%u(1)*1._DP

    CASE(6)
        ! Woodward-Colella left
        DATA%tmax = 0.012_DP
       IF (x <= 0.5_DP) THEN
          Var%u(1) = 1._DP ! [kg/m^3]
          Var%u(2) = 0._DP ! [m/s] !*NEIN, [kg/m^3*m/s]
       ELSE
          Var%u(1) = 1._DP ! [kg/m^3]
          Var%u(2) = 0._DP ! [m/s] !*NEIN, [kg/m^3*m/s]
       END IF
           
    CASE(10)
       DAta%tmax=0.5_DP
       Var%u(1) = 0.001_DP
       IF(x< 0._DP) THEN
         Var%u(2)= -Var%u(1)
        ELSE
         Var%u(2) = Var%u(1)
        END IF

    CASE (27)
      DATA%tmax = 15.0_dp
      r=0.5_DP
      x0=3_DP
      alpha=0.1_DP
      IF (x<x0+r .AND. x>x0-r) THEN
        var%u(1) = 1._DP +alpha*EXP(1._DP-1._DP/(1._DP-((x-x0)/r)**2))!(4.0_DP/9.81_DP) **(1._DP/3._DP) !*(1._DP-x) ! + 1._DP 
      ELSE
        var%u(1)= 1._dp ! (4._dp/9.81_dp) **(1._dp/3._dp) 
      ENDIF
      Var%u(2) = 4._dp



       !*A New Approach for Designing Moving-Water Equilibria Preserving Schemes for the Shallow Water Equations - Yuanzhen Cheng, Alina Chertock, Michael Herty, Alexander Kurganov, Tong Wu, section 3 example 2
       !*100 supercritical transitory
       !*101 subcritical transitory
       !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 3.1.3 subcritical
       !*111 subcritical but with IC taken from scratch
       !*120 supercritical with IC taken from scratch but smooth
       !*121 subcritical with IC taken from scratch but smooth
       !*122 transcritical with IC taken from scratch but smooth
    CASE(100,101) 
      DATA%tmax = 100._DP
      Var%U(1) = 2._DP -bathymetry(DATA%test,x)
      Var%U(2) = 0._DP
    CASE(111,120,121,122)
      DATA%tmax = 100._DP !*100._DP for convergence analyses
      Var%U=0._DP
      DATA%tmax = 0.3_DP !*ALERT: COMPARISON LAGRANGIAN


    !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
    !*3.2.2 Short channel: 100 m
    !*150 Supercritical case 
    !*440 like 150 but for perturbation analyisis
    CASE(150) !*But I also take from scratch 
		DATA%tmax=100._DP
      Var%U(1)=(4._DP/grav)**(1._DP/3._DP) * (1._DP-0.25_DP*EXP( -4._DP*(x/100._DP-0.5_DP)**2 ) )
      Var%U(2)=2._DP
      !*PRINT*, grav, 9.81_DP

    !*NB: IN 440 AND 441 THE DATA%tmax IS GIVEN IN init_geom AND NOT HERE
    CASE(440) !*But I also take from scratch the same initial condition (because I need to enter the bathymetry)
      Var%U(1)=(4._DP/grav)**(1._DP/3._DP) * (1._DP-0.25_DP*EXP( -4._DP*(x/100._DP-0.5_DP)**2 ) )
      Var%U(2)=2._DP
    CASE(441)  
      Var%U(1)=0._DP
      Var%U(2)=0._DP



    !*200 lake at rest - periodic BCs (perturbed with 10**(-2) and 10**(-4))
    !*201 discrete equilibrium of 200 from scratch 
    !*210 lake at rest - strong BCs
    !*211 discrete equilibrium of 210 from scratch 
    CASE(200,210) 
      DATA%tmax = 10._DP !*100._DP !*Long time for discrete equilibrium
      Var%U(1) = 2._DP -bathymetry(DATA%test,x)
      Var%U(2) = 0._DP
    CASE(201,211) !*From scratch
      DATA%tmax = 8._DP
      Var%U=0._DP

    !*230 lake at rest with C^0 bathymetry - outflow BCs (perturbed with 10**(-2) and 10**(-4))
    !*231 discrete equilibrium of 230 from scratch 
    CASE(230) 
      DATA%tmax = 100._DP !*100._DP !*Long time for discrete equilibrium
      Var%U(1) = 2._DP -bathymetry(DATA%test,x)
      Var%U(2) = 0._DP
    CASE(231) !*From scratch
      DATA%tmax = 8._DP
      Var%U=0._DP

    !*Well-Balancing Via Flux Globalization: Applications to Shallow Water Equations with Wet/Dry Fronts
    !*Alina Chertock, Alexander Kurganov, Xin Liu, Yongle Liu, and Tong Wu
    CASE(300)
      DATA%tmax = 0.08_DP 
      Var%U(1) = 3._DP -bathymetry(DATA%test,x)
      Var%U(2) = 0._DP       



    !*NB: THE FINAL TIME IN THE NEXT TESTS IS GIVEN IN init_geom AND NOT HERE
    !*SWASHES WITH BATHYMETRY NON SMOOTH
    !*400 Lake at rest with an immersed bump 3.1.1
    !*410 Supercritical
    !*420 Subcritical
    !*430 Transcritical
    CASE(400)
      Var%U(1) = 0.5_DP -bathymetry(DATA%test,x)
      Var%U(2) = 0._DP       
    CASE(401,410,411,415,416,420,421,425,426,430,431,435) !*Fake initialization, later taken from scratch
      Var%U=0._DP       


    !*IC taken from scratch
    !*401 Lake at rest with an immersed bump 3.1.1 
    !*411 Supercritical
    !*421 Subcritical
    !*431 Transcritical

    !*Same tests with friction
    !*415 Supercritical for the generation of the steady state with friction
    !*425 Subcritical for the generation of the steady state with friction
    !*435 Transcritical for the generation of the steady state with friction


    !*IC taken from scratch with friction
    !*416 Supercritical
    !*426 Subcritical
     







    CASE default
       PRINT*, "Wrong test number for SW_IC, test = ", DATA%test
       STOP
    END SELECT
  END FUNCTION IC

 !*---------------------------------------------------------------------------------------
 !*Initialize from scratch when it is needed
 !*---------------------------------------------------------------------------------------
 SUBROUTINE InitializeFromScratch(DATA,Mesh,Ucoeff)
    TYPE(donnees), INTENT(IN) :: DATA
    TYPE(maillage), INTENT(IN) :: Mesh
    TYPE(PVar), DIMENSION(:), INTENT(INOUT) :: UCoeff

    INTEGER :: indi !*Loop on the DoFs
    INTEGER :: indt !*Loop on the elements

    !*Variables read from the .dat
    REAL(DP), DIMENSION(Mesh%nDoFs) :: h, hu, b
    CHARACTER(LEN = 64) :: fakevariabletoreadchar    
    INTEGER :: fakevariabletoreadint
    REAL(DP) :: fakevariabletoreadreal

    !*To deal with the numeration issues
    INTEGER, DIMENSION(Mesh%nDoFs) :: GlobalNumeration !*GobalNumeration(indi) is the global index of the indi DoF with increasing order of abscissa   
    INTEGER :: orderpolynomials 
    INTEGER :: indsupp !*Loop on the DoFs

    TYPE(PVar), DIMENSION(Mesh%nDoFs) :: UVal
    TYPE(PVar), DIMENSION(Mesh%nDoFs) :: BVal, BCoeff !*Used (improperly) for the bathymetry

    TYPE(element) :: e

    IF(DATA%test==111 .OR. DATA%test==120 .OR. DATA%test==121 .OR. &
       & DATA%test==122  .OR. DATA%test==201  .OR. DATA%test==211 .OR. DATA%test==231 .OR. &
       & DATA%test==150 .OR. DATA%test==401 .OR. DATA%test==410 .OR. DATA%test==411 .OR. &
       & DATA%test==420 .OR. DATA%test==421 .OR. DATA%test==430 .OR. DATA%test==431 .OR. &
       & DATA%test==440 .OR. DATA%test==441 .OR. DATA%test==415 .OR. DATA%test==416 .OR. &
       & DATA%test==425 .OR. DATA%test==426 .OR. DATA%test==435) THEN !*Only for some cases


       !*0) Open the file
       SELECT CASE(DATA%test)
       !*A New Approach for Designing Moving-Water Equilibria Preserving Schemes for the Shallow Water Equations - Yuanzhen Cheng, Alina Chertock, Michael Herty, Alexander Kurganov, Tong Wu, section 3 example 2
       !*100 supercritical transitory
       !*101 subcritical transitory
       !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 3.1.3 subcritical
       !*111 subcritical but with IC taken from scratch
       !*120 supercritical with IC taken from scratch but smooth
       !*121 subcritical with IC taken from scratch but smooth
       !*122 transcritical with IC taken from scratch but smooth
       CASE(111) 
          !*subcritical 3.1.3
          OPEN(100, file="subcritical.dat")
       CASE(120) 
          !*supercritical
          OPEN(100, file="supercriticalsmooth.dat")
       CASE(121) 
          !*subcritical 3.1.3
          OPEN(100, file="subcriticalsmooth.dat")
       CASE(122) 
          !*transcritical 3.1.4
          OPEN(100, file="transcriticalsmooth.dat")


		 !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
		 !*3.2.2 Short channel: 100 m
		 !*150 Supercritical case
       !*440 like 150 but for perturbation analyisis 
       !*441 discrete steady state of 440
		 CASE(150,440) 
          !*3.2.2 Short channel: 100 m supercritical case
          OPEN(100, file="friction_bathymetry_supercritical.dat")
       CASE(441)
          OPEN(100, file="SupercriticalFrictionShortChannelSwashesDiscreteEquilibrium.dat")

       CASE(201) 
          !*Periodic
          OPEN(100, file="LakeAtRestDiscreteEquilibriumPeriodicBCs.dat")
       CASE(211) 
          !*BCs
          OPEN(100, file="LakeAtRestDiscreteEquilibriumOutflowBCs.dat")
       CASE(231) 
          !*BCs
          OPEN(100, file="LakeAtRestDiscreteEquilibriumOutflowBCsC0bathymetry.dat")

		 !*SWASHES WITH BATHYMETRY NON SMOOTH
       !*DISCRETE EQUILIBRIUM PERTURBED
		 !*400+1 Lake at rest with an immersed bump 3.1.1
		 !*410+1 Supercritical
		 !*420+1 Subcritical
		 !*430+1 Transcritical
       !*Same test with friction
		 !*415+1 Supercritical
		 !*425+1 Subcritical
		 !*435+1 Transcritical
       CASE(401)
          OPEN(100, file="LakeAtRestSwashesDiscreteEquilibrium.dat")
       CASE(410)
          OPEN(100, file="supercriticalnotsmooth.dat")
       CASE(411)
          OPEN(100, file="SupercriticalNotSmoothSwashesDiscreteEquilibrium.dat")
       CASE(415) !*Open the same as 410 to generate the steady state with friction
          OPEN(100, file="supercriticalnotsmooth.dat")
       CASE(416)
          OPEN(100, file="SupercriticalNotSmoothFrictionSwashesDiscreteEquilibriumInterpolated.dat")
       CASE(420)
          OPEN(100, file="subcriticalnotsmooth.dat")
       CASE(421)
          OPEN(100, file="SubcriticalNotSmoothSwashesDiscreteEquilibrium.dat")
       CASE(425) !*Open the same as 420 to generate the steady state with friction
          OPEN(100, file="subcriticalnotsmooth.dat")
       CASE(426)
          OPEN(100, file="SubcriticalNotSmoothFrictionSwashesDiscreteEquilibriumInterpolated.dat")
       CASE(430)
          OPEN(100, file="transcriticalnotsmooth.dat")
       CASE(431)
          OPEN(100, file="TranscriticalNotSmoothSwashesDiscreteEquilibrium.dat")
       CASE(435) !*Open the same as 430 to generate the steady state with friction
          OPEN(100, file="transcriticalnotsmooth.dat")
       CASE default
          PRINT*, "Wrong test number for InitializeFromScratch, test = ", DATA%test
          STOP
       END SELECT
 
       READ(100,*) fakevariabletoreadchar, fakevariabletoreadchar, fakevariabletoreadchar, fakevariabletoreadchar

       !*1) Copying in the wrong order
       DO indi=1,Mesh%nDoFs !*Loop on the DoFs
		    SELECT CASE(DATA%test)
		    CASE(111,120,121,122,201,211,231,401,410,411,415,416,420,421,425,426,430,431,435) !*READ THE CONSERVED VARIABLES ONLY
		       READ(100,*) fakevariabletoreadint, fakevariabletoreadreal, h(indi), hu(indi) !*Read the line indi
		       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		       !*SAFETY CHECK to assess whether what I have read is fine
		       !*PRINT*, fakevariabletoreadint, fakevariabletoreadreal, h(indi), hu(indi)
		       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    CASE(150,440,441) !*READ ALSO THE BATHYMETRY
		       READ(100,*) fakevariabletoreadint, fakevariabletoreadreal, h(indi), hu(indi), b(indi) !*Read the line indi
		       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		       !*SAFETY CHECK to assess whether what I have read is fine
		       !*PRINT*, fakevariabletoreadint, fakevariabletoreadreal, h(indi), hu(indi), b(indi)
		       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          CASE default
             PRINT*, "Wrong test number for InitializeFromScratch, test = ", DATA%test
             STOP
		    END SELECT

       END DO
       CLOSE(100)

       !*2)GlobalNumeration, structure needed for the reordering
       !*GobalNumeration(indi) is the global index of the indi DoF with increasing order of abscissa    

       GlobalNumeration=0

       !*RMK:
       !*itype 
       !*1: P1, 
       !*2: B2, 
       !*3: P2, 
       !*4: P3, 
       !*5: B3
       !*6: B4
       !*7: P4
       !*11, 12, 13, 14: P1, P2, P3, P4 in Gauss-Lobatto


       SELECT CASE(DATA%itype)
       CASE(1,11) !*P1 (eventually in Gauss-Lobatto points)
          orderpolynomials=1 
       CASE(2,3,12) !*B2, P2 (eventually in Gauss-Lobatto points)
          orderpolynomials=2
       CASE(4,5,13) !*P3, B3 (P3 eventually in Gauss-Lobatto points)
          orderpolynomials=3
       CASE(6,7,14) !*B4, P4 (PGL4 eventually in Gauss-Lobatto points)
          orderpolynomials=4
       CASE default
          PRINT*, "Error in Initialize from Scratch: this element is not yet defined", DATA%itype
          STOP
       END SELECT

       !*Vertices
       DO indt=1,Mesh%Nt+1
          GlobalNumeration((indt-1)*orderpolynomials+1)=indt
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, "Vertex", indt 
          !*PRINT*, "Index in the increasing numeration", indt+(indt-1)*orderpolynomials
          !*PRINT*, "Index in the Global numeration", GlobalNumeration((indt-1)*orderpolynomials+1)
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO



       !*Internal nodes
       indsupp=0
       SELECT CASE(DATA%itype)
       CASE(1,11) !*P1 (eventually in Gauss-Lobatto points)
          !*No internal nodes
       CASE(2,3,12) !*B2, P2 (eventually in Gauss-Lobatto points)
          indsupp=0
          DO indt=1,Mesh%Nt
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+1)=Mesh%Nt+1+indsupp
          END DO
       CASE(4,5,13) !*P3, B3 (P3 eventually in Gauss-Lobatto points)
          indsupp=0
          DO indt=1,Mesh%Nt
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+1)=Mesh%Nt+1+indsupp
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+2)=Mesh%Nt+1+indsupp
          END DO
       CASE(6,7,14) !*B4, P4 (PGL4 eventually in Gauss-Lobatto points) !*Not checked but it should be ok
          indsupp=0
          DO indt=1,Mesh%Nt
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+1)=Mesh%Nt+1+indsupp
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+2)=Mesh%Nt+1+indsupp
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+3)=Mesh%Nt+1+indsupp
          END DO

       CASE default
          PRINT*, "Error in Initialize from Scratch: this element is not yet defined", DATA%itype
          STOP
       END SELECT

       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*DO indi=1,Mesh%nDoFs
       !*   PRINT*, indi, GlobalNumeration(indi)
       !*END DO
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !*3)Acquisition in Uval
       Uval=0._DP
       Bval=0._DP
       !*GobalNumeration(indi) is the global index of the indi DoF with increasing order of abscissa  
       DO indi=1,Mesh%nDoFs
		    SELECT CASE(DATA%test)
		    CASE(111,120,121,122,201,211,231,401,410,411,415,416,420,421,425,426,430,431,435) !*READ THE CONSERVED VARIABLES ONLY
		       Uval(GlobalNumeration(indi))%u(1)=h(indi) !*h
		       Uval(GlobalNumeration(indi))%u(2)=hu(indi) !*hu
		       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		       !*SAFETY CHECK
		       !*PRINT*, indi, GlobalNumeration(indi), h(indi), Uval(GlobalNumeration(indi))%u(1)
		       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		    CASE(150,440,441) !*ACQUISITION ALSO OF THE BATHYMETRY
		       Uval(GlobalNumeration(indi))%u(1)=h(indi) !*h
		       Uval(GlobalNumeration(indi))%u(2)=hu(indi) !*hu
             Bval(GlobalNumeration(indi))%u(:)=b(indi) !*b
          CASE default
             PRINT*, "Wrong test number for InitializeFromScratch, test = ", DATA%test
             STOP             
		    END SELECT




       END DO

      !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !*SAFETY CHECK
      !*DO indt=1, Mesh%nt
      !*   PRINT*, indt
      !*   e=Mesh%e(indt)
      !*   DO indi=1,e%nsommets
      !*      PRINT*, "local", indi, "coor", e%coor(indi), "h", Uval(e%nu(indi))%u(1)
      !*   END DO
      !*END DO
      !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

	    SELECT CASE(DATA%test)
	    CASE(111,120,121,122,201,211,231,401,410,411,415,416,420,421,425,426,430,431,435) !*CONSERVED VARIABLES ONLY
	       Ucoeff=global_Cons_to_Control(Uval,Mesh)
		    !*Acquisition in the utils vectors
		    ALLOCATE(CoeffICUtils(Mesh%Ns),ValICUtils(Mesh%Ns))

		    CoeffICUtils=0._DP
		    ValICUtils=0._DP
		    CoeffICUtils=Ucoeff
		    ValICUtils=Uval


	    CASE(150,440,441) !*ALSO OF THE BATHYMETRY
			 Ucoeff=global_Cons_to_Control(Uval,Mesh)
          Bcoeff=global_Cons_to_Control(Bval,Mesh)
		    !*Acquisition in the utils vectors
		    ALLOCATE(CoeffICUtils(Mesh%Ns),ValICUtils(Mesh%Ns))

		    CoeffICUtils=0._DP
		    ValICUtils=0._DP
		    CoeffICUtils=Ucoeff
		    ValICUtils=Uval


			 IF (ALLOCATED(BathymetryAtDoFsUtils) .OR. ALLOCATED(CoeffBathymetryAtDoFsUtils) ) THEN !*Then there is sth wrong, it means that we have already initialized these structures. We stop the execution
				 PRINT*, "You are trying to initialize again the bathymetry and its gradient. Not possible." 
				 PRINT*, "(I'm in init_bc_sw, InitializeFromScratch)"
				 STOP
			 END IF
			 !*Otherwise continue

          ALLOCATE(BathymetryAtDoFsUtils(Mesh%Ns),CoeffBathymetryAtDoFsUtils(Mesh%Ns))

          DO indi=1,Mesh%Ns
	          BathymetryAtDoFsUtils(indi)=Bval(indi)%u(1)
   	       CoeffBathymetryAtDoFsUtils(indi)=Bcoeff(indi)%u(1)
             !*PRINT*, indi, BathymetryAtDoFsUtils(indi), CoeffBathymetryAtDoFsUtils(indi)
          END DO

       CASE default
          PRINT*, "Wrong test number for InitializeFromScratch, test = ", DATA%test
          STOP             
	    END SELECT



END IF



 END SUBROUTINE InitializeFromScratch

 !*---------------------------------------------------------------------------------------
 !*Introduce a perturbation in the IC if needed
 !*---------------------------------------------------------------------------------------
 SUBROUTINE Perturbation(DATA,Mesh,Ucoeff)
    TYPE(donnees), INTENT(INOUT) :: DATA
    TYPE(maillage), INTENT(IN) :: Mesh
    TYPE(PVar), DIMENSION(:), INTENT(INOUT) :: UCoeff

    INTEGER :: indi !*Loop on the DoFs
    INTEGER :: indt !*Loop on the elements

    TYPE(PVar), DIMENSION(Mesh%nDoFs) :: UVal

    TYPE(element) :: e !*For the quick reference to an element
    REAL(dp) ::x !*For quick reference to the coordinates of the DoFs

    REAL(dp) ::x0,r,A !*Parameters for the perturbations
    !*x0 center
    !*r radius
    !*A amplitude

    
    IF(DATA%perturbation/=0) THEN !*Only for perturbation

       !*We go back to the values
       Uval=0._DP
       Uval=global_Control_to_Cons(UCoeff,Mesh)

       !*We introduce the perturbation
       DO indt=1,Mesh%Nt !*Loop on the elements
          e=Mesh%e(indt) !*Quick reference
          DO indi=1,e%nsommets !*Loop on the DoFs of the element

             x=e%coor(indi) !*Quick reference

             SELECT CASE(DATA%test) !*Test
             CASE(120,122) !*Supercritical and transcritical
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) 
                      x0=17._dp
                      r=2._dp
                      A=5._DP*10._DP**(-1)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(2) !*Big perturbation
                      x0=17._dp
                      r=1._dp
                      A=10._DP**(-2)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(3) !*Small perturbation
                      x0=17._dp
                      r=1._dp
                      A=10._DP**(-4)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT
             CASE(200) !*Lake at Rest
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) !*Big perturbation
                      x0=17._dp
                      r=1._dp
                      A=10._DP**(-2)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=IC(x,DATA)
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(2) !*Small perturbation
                      x0=17._dp
                      r=1._dp
                      A=10._DP**(-4)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=IC(x,DATA)
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT
             CASE(201,211,231) !*Lake at Rest discrete equilibrium from scratch
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) !*Big perturbation
                      x0=17._dp
                      r=1._dp
                      A=10._DP**(-2)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(2) !*Small perturbation
                      x0=17._dp
                      r=1._dp
                      A=10._DP**(-4)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT
             !*Well-Balancing Via Flux Globalization: Applications to Shallow Water Equations with Wet/Dry Fronts
             !*Alina Chertock, Alexander Kurganov, Xin Liu, Yongle Liu, and Tong Wu
             CASE(300) 
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) !*Perturbation
                      IF ((x.GE.-0.2_DP) .AND. (x.LE.-0.1_DP)) THEN
                         UVal(e%nu(indi))=IC(x,DATA)
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + 10._DP**(-5)
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT










				 !*SWASHES WITH BATHYMETRY NON SMOOTH
				 !*400 Lake at rest with an immersed bump 3.1.1
				 !*410 Supercritical
				 !*420 Subcritical
				 !*430 Transcritical
             CASE(400) 
                DATA%tmax=1.5_DP
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) !*BIG 5*10**(-2) 
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-2)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=IC(x,DATA)
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(2) !*INTERMEDIATE 5*10**(-4)
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-4)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=IC(x,DATA)
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(3) !*Small perturbation
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-5)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=IC(x,DATA)
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT

             CASE(401) 
                DATA%tmax=1.5_DP
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) !*BIG 5*10**(-2) 
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-2)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(2) !*INTERMEDIATE 5*10**(-4)
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-4)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(3) !*Small perturbation
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-5)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT

             CASE(410,411,416) !*NB: Perturbation not allowed for 415 because used to generate the steady state with friction
                DATA%tmax=1._DP
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) !*BIG 5*10**(-2) 
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-2)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(2) !*INTERMEDIATE 5*10**(-4)
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-4)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(3) !*Small perturbation
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-5)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT

             CASE(420,421,426) !*NB: Perturbation not allowed for 425 because used to generate the steady state with friction
                DATA%tmax=1.5_DP
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) !*BIG 5*10**(-2) 
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-2)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(2) !*INTERMEDIATE 5*10**(-4)
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-4)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(3) !*Small perturbation
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-5)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT

             CASE(430,431) !*NB: Perturbation not allowed for 435 because used to generate the steady state with friction
                DATA%tmax=1.5_DP
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) !*BIG 5*10**(-2) 
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-2)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(2) !*INTERMEDIATE 5*10**(-4)
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-4)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(3) !*Small perturbation
                      x0=6._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-5)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT

             CASE(440,441) 
                DATA%tmax=5._DP !*Longer domain
                SELECT CASE(DATA%perturbation) !*Perturbation
                   CASE(1) !*BIG 5*10**(-2) 
                      x0=50._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-2)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(2) !*INTERMEDIATE 5*10**(-4)
                      x0=50._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-4)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(3) !*Small perturbation
                      x0=50._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-5)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE(4) !*Smaller perturbation
                      x0=50._dp
                      r=0.5_dp
                      A=5._DP*10._DP**(-6)
                      IF (x>x0-r .AND. x<x0+r ) THEN
                         UVal(e%nu(indi))=ValICUtils(e%nu(indi))
                         UVal(e%nu(indi))%u(1) = UVal(e%nu(indi))%u(1) + A*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
                      ELSE
                         !*Do nothing
                      END IF
                   CASE default
                      PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                      STOP
                END SELECT


             CASE default
                PRINT*, "Error in Perturbation", DATA%test, DATA%perturbation
                STOP
             END SELECT

          END DO !*End loop on the DoFs
       END DO !*End loop on the elements

       !*We go back to the coefficients
       Ucoeff=global_Cons_to_Control(Uval,Mesh)


    END IF !*End if perturbarion

 END SUBROUTINE Perturbation


  !*--------------------------------------------------------------------------------------
  !*BC_Residuals deals with the boundary residuals and it is called before the updating in fin() of the solution
  !*
  !*IMPORTANT RMK: Also in the 1D case the boundary residuals of the different subtimesteps MUST be combined through the thetas.
  !*
  !*NB: the time DATA%time=tn+alpha(k,DATA%iordret)*dt so it is better to use tn and alpha to refer to the times of the subtimesteps if needed (or at least don't forget this)
  !*
  !*In practice we can often "avoid" to perform a combination inside BC_Residuals because we have: 
  !*-> PERIODIC BC: Var%un(1) and Var%un(Mesh%nt) (node residuals of the boundary nodes (first and last) and so ALREADY combined) are given to the first and the last node. We just need to keep into account the fact that Var%un(Mesh%nt) goes also to the first node and Var%un(1) goes also to the last node
  !*-> STRONG IMPOSITION: nothing to combine at the boundary. We just impose the values of Var%ua without caring for the updating.
  !*-> OUTFLOW: no contribution of the boundary residuals to the residuals at the nodes
  !*
  !*In 2d it is more critical: the weak inflow and wall boundary conditions generate boundary residuals that must be combined through the thetas
  !*--------------------------------------------------------------------------------------
  SUBROUTINE BC_Residuals(Var,i1,iN,DATA,k, n_theta,theta,alpha,dt,tn)
    TYPE(variables), INTENT(inout):: Var
    INTEGER,         INTENT(in)   :: i1, iN, k
    TYPE(donnees),   INTENT(in)   :: DATA
    INTEGER, DIMENSION(5), INTENT(IN):: n_theta 
    REAL(dp),DIMENSION(0:4,1:4,1:5), INTENT(IN) :: theta
    REAL(dp),DIMENSION(4,2:5), INTENT(IN) :: alpha
    REAL(DP), INTENT(IN) :: dt
    REAL(DP), INTENT(IN) :: tn !*NB: DATA%time=tn+alpha(k,DATA%iordret)*dt    
    
    TYPE(Pvar):: a0,a1
    INTEGER:: p1, pN

    SELECT CASE(DATA%test)
    CASE(0,11,12,20,19,14,40,41,42,27,501,-2,200,201,211,231) !*No combination
!!$        ! periodic
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=a0+a1
        Var%un(pN)=a0+a1
    CASE(210,230) !*Steady strong
        Var%un(i1)=0._DP
        Var%un(iN)=0._DP
    CASE(1,2,6,10,13,30,29,32,500) !*No combination, anyway this is not outflow because we should'n to anything when outflow
       !*Basically strong steady
       ! outflow
       Var%un(i1)%u(:)=0._dp
       Var%un(iN)%u(:)=0._dp

    CASE(3) !*No combination
       ! reflective
       !dirichlet for velocity
       Var%un(i1)%u(2) = 0._DP
       Var%un(iN)%u(2) = 0._DP
	  
    CASE(4) !*No combination
        ! inflow left + outflow right
        !*Inflow imposed in Strong_BC after the updating
    CASE(21) !*No combination
        !dirichlet bc
        !*Imposed after the updating
    CASE(22,23) !*No combination
        !dirichlet bc
        !*Imposed after the updating
    CASE(31) !waves on the shore
        !dirichlet bc
        !*eta = 0._dp !profile_bathymetry( DATA%test, (/ DATA%domain_left /) ) !*NOT USED
        Var%un(iN)%u(:) = 0._dp
    CASE(24)
        !dirichlet bc
        !*Imposed after the updating
    CASE(25)
        !dirichlet bc
        !*Imposed after the updating
    CASE(26)
        !dirichlet bc
        !*Imposed after the updating
    CASE(5)
        ! periodic
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=a0+a1
        Var%un(pN)=a0+a1
       !*A New Approach for Designing Moving-Water Equilibria Preserving Schemes for the Shallow Water Equations - Yuanzhen Cheng, Alina Chertock, Michael Herty, Alexander Kurganov, Tong Wu, section 3 example 2
       !*100 supercritical transitory
       !*101 subcritical transitory
       !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 3.1.3 subcritical
       !*111 subcritical but with IC taken from scratch
       !*120 supercritical with IC taken from scratch but smooth
       !*121 subcritical with IC taken from scratch but smooth
       !*122 transcritical with IC taken from scratch but smooth
    CASE(100,101,111,120,121,122) 
        !*supercritical
        !*inflow left, imposed later strongly
        !*outflow right, nothing

        !*subcritical
        !*q left and h right given in strong bc

        !*transcritical
        !*q left


    !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
	 !*3.2.2 Short channel: 100 m
	 !*150 Supercritical case
    !*440 like 150 but for perturbation analyisis 
	 CASE(150,440,441) 
        !*supercritical
        !*inflow left, imposed later strongly
        !*outflow right, nothing
             




    !*Well-Balancing Via Flux Globalization: Applications to Shallow Water Equations with Wet/Dry Fronts
    !*Alina Chertock, Alexander Kurganov, Xin Liu, Yongle Liu, and Tong Wu
    CASE(300)
        !*Free BCs !*SET PERIODIC OR EXPLOSION
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=a0+a1
        Var%un(pN)=a0+a1





    !*SWASHES WITH BATHYMETRY NON SMOOTH
    !*400 Lake at rest with an immersed bump 3.1.1
    !*410 Supercritical
    !*420 Subcritical
    !*430 Transcritical
    CASE(400,401)!*Steady strong
        Var%un(i1)=0._DP
        Var%un(iN)=0._DP    

    !*IC taken from scratch
    !*401 Lake at rest with an immersed bump 3.1.1 
    !*411 Supercritical
    !*421 Subcritical
    !*431 Transcritical

    !*Same tests with friction
    !*415 Supercritical for the generation of the steady state with friction
    !*425 Subcritical for the generation of the steady state with friction
    !*435 Transcritical for the generation of the steady state with friction


    !*IC taken from scratch with friction
    !*416 Supercritical
    !*426 Subcritical


	 CASE(410,411,415,416,420,421,425,426,430,431,435)
    !*supercritical
    !*inflow left, imposed later strongly
    !*outflow right, nothing

    !*subcritical
    !*q left and h right given in strong bc

    !*transcritical
    !*q left



    CASE default
       PRINT*, "Wrong test number for BC_Residuals, test = ", DATA%test
       STOP
    END SELECT

  END SUBROUTINE BC_Residuals

  !*--------------------------------------------------------------------------------------
  !*Strong_BC imposed after the updating
  !*--------------------------------------------------------------------------------------
  SUBROUTINE Strong_BC(Var,i1,iN,DATA,k, n_theta,theta,alpha,dt,tn)
    TYPE(variables), INTENT(inout):: Var
    INTEGER,         INTENT(in)   :: i1, iN, k
    TYPE(donnees),   INTENT(in)   :: DATA
    INTEGER, DIMENSION(5), INTENT(IN):: n_theta 
    REAL(dp),DIMENSION(0:4,1:4,1:5), INTENT(IN) :: theta
    REAL(dp),DIMENSION(4,2:5), INTENT(IN) :: alpha
    REAL(DP), INTENT(IN) :: dt
    REAL(DP), INTENT(IN) :: tn !*NB: DATA%time=tn+alpha(k,DATA%iordret)*dt    
    REAL(DP) :: x 

	 REAL(DP) :: ur, ssr !*Velocity and speed of sound at the right boundary to check the characteristics for the transcritical case

    SELECT CASE(DATA%test)
       CASE(0,11,12,20,19,14,40,41,42,27,1,2,6,10,13,30,29,32,3,5,500,501,-2,200,201,210,211,230,231,300) !*NO STRONG -> Do nothing

       
       CASE(4)
          Var%ua(i1,k)%u(1) = 1._dp
          Var%ua(i1,k)%u(2) = 3._dp
        ! Var%ua(k,i1)%u(3) = 4.5_dp + 1.7857_dp
       CASE(21)
          Var%ua(i1,k)%u(1) = 0.5_dp
          Var%ua(iN,k)%u(1) = 0.5_dp
          Var%ua(i1,k)%u(2) = 0.0_dp
          Var%ua(iN,k)%u(2) = 0.0_dp
       CASE(22,23)
          Var%ua(i1,k)%u(2) = 4.42_dp
          Var%ua(iN,k)%u(1) = 2._dp!MIN(2.0_dp,Var%ua(k,iN)%u_lim(1) ) 
       CASE(31)
          Var%ua(i1,k)%u(1) = 0.37_dp+0.0615_DP*SIN(DATA%temps*12._DP*ASIN(1._DP)/10._DP)-bathymetry(DATA%test,DATA%domain_left)
          Var%ua(iN,k)%u(1) = 0.0_dp
          Var%ua(iN,k)%u(2) = 0.0_dp
          ! Var%ua(k,iN)%u_lim(1) = 2._dp!MIN(2.0_dp,Var%ua(k,iN)%u_lim(1) ) 
       CASE(24)
          Var%ua(i1,k)%u(2) = 1.53_dp
          IF (Var%ua(iN,k)%u(1) >0.66_DP ) THEN
             Var%ua(iN,k)%u(1) = 0.60_DP
          ENDIF
       CASE(25)
          Var%ua(i1,k)%u(2) = 0.18_dp
          Var%ua(iN,k)%u(1) = 0.33_dp
       CASE(26)
          Var%ua(i1,k)%u(2) = 0.3_dp
          Var%ua(iN,k)%u(1) = 0.3_dp
       CASE(100) !*A New Approach for Designing Moving-Water Equilibria Preserving Schemes for the Shallow Water Equations - Yuanzhen Cheng, Alina Chertock, Michael Herty, Alexander Kurganov, Tong Wu, section 3 example 2 supercritical
          Var%ua(i1,k)%u(1)=2._DP
          Var%ua(i1,k)%u(2)=24._DP
       !*A New Approach for Designing Moving-Water Equilibria Preserving Schemes for the Shallow Water Equations - Yuanzhen Cheng, Alina Chertock, Michael Herty, Alexander Kurganov, Tong Wu, section 3 example 2
       !*100 supercritical transitory
       !*101 subcritical transitory
       !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 3.1.3 subcritical
       !*111 subcritical but with IC taken from scratch
       !*120 supercritical with IC taken from scratch but smooth
       !*121 subcritical with IC taken from scratch but smooth
       !*122 transcritical with IC taken from scratch but smooth
       CASE(120) 
          Var%ua(i1,k)%u(1)=2._DP !*h on the left
          Var%ua(i1,k)%u(2)=24._DP !*q on the left
       CASE(101,111,121) 
          Var%ua(i1,k)%u(2)=4.42_DP !*q on the left
          Var%ua(iN,k)%u(1)=2._DP !*h on the right
       CASE(122)
          Var%ua(i1,k)%u(2)=1.53_DP !*q on the left

		 !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
		 !*3.2.2 Short channel: 100 m
		 !*150 Supercritical case
       !*440 like 150 but for perturbation analyisis 
		 CASE(150,440,441) 
        !*supercritical
        !*inflow left
        !*outflow right, nothing
         x=0._DP
         Var%ua(i1,k)%u(1)=(4._DP/grav)**(1._DP/3._DP) * (1._DP-0.25_DP*EXP( -4._DP*(x/100._DP-0.5_DP)**2 ) ) !*h on the left
         Var%ua(i1,k)%u(2)=2._DP !*q on the left



		 !*SWASHES WITH BATHYMETRY NON SMOOTH
		 !*400 Lake at rest with an immersed bump 3.1.1
		 !*410 Supercritical
		 !*420 Subcritical
		 !*430 Transcritical
		 CASE(400,401)!*Steady strong
		 		!*Given by nullifying the residual, but also strong      
            Var%ua(i1,k)%u(1)=0.5_DP-bathymetry(DATA%test,DATA%domain_left)
            Var%ua(i1,k)%u(2)=0._DP !*q on the left
            Var%ua(iN,k)%u(1)=0.5_DP-bathymetry(DATA%test,DATA%domain_left+DATA%Length)
            Var%ua(iN,k)%u(2)=0._DP !*q on the right


		 !*IC taken from scratch
		 !*401 Lake at rest with an immersed bump 3.1.1 
		 !*411 Supercritical
		 !*421 Subcritical
		 !*431 Transcritical

       !*Same tests with friction
		 !*Same tests with friction
		 !*415 Supercritical for the generation of the steady state with friction
		 !*425 Subcritical for the generation of the steady state with friction
		 !*435 Transcritical for the generation of the steady state with friction

		 !*IC taken from scratch with friction
		 !*416 Supercritical
		 !*426 Subcritical


       CASE(410,411,415,416) !*supercritical 
          Var%ua(i1,k)%u(1)=2._DP !*h on the left
          Var%ua(i1,k)%u(2)=24._DP !*q on the left

       CASE(420,421,425,426) !*subcritical
          Var%ua(i1,k)%u(2)=4.42_DP !*q on the left
          Var%ua(iN,k)%u(1)=2._DP !*h on the right

       CASE(430,431) !*transcritical
          Var%ua(i1,k)%u(2)=1.53_DP !*q on the left

		 CASE(435) !*DO NOT RELY !*ALERT
          Var%ua(i1,k)%u(2)=1.53_DP !*q on the left
          ur=Var%ua(iN,k)%u(2)/Var%ua(iN,k)%u(1) !*Velocity on the right
			 ssr=SQRT(grav*Var%ua(iN,k)%u(1))       !*Speed of sound on the right
			 IF ((ur-ssr) .LE. 0._DP) THEN !*subcritical
			    Var%ua(iN,k)%u(1)=0.66_DP
             !*!!!!!!!!!!!!!!!!!!!!!!!!!
             !*SAFETY CHECK
             !*PRINT*, "I'm here"
             !*STOP
             !*!!!!!!!!!!!!!!!!!!!!!!!!!
			 ELSE
             !*Nothing
             !*!!!!!!!!!!!!!!!!!!!!!!!!!
             !*SAFETY CHECK
				 !*PRINT*, "Look, I'm asymptotically supercritical"
             !*!!!!!!!!!!!!!!!!!!!!!!!!!
			 END IF

       CASE DEFAULT
          PRINT*, "Wrong test number for Strong_BC, test = ", DATA%test
          STOP
       END SELECT




  END SUBROUTINE Strong_BC



  SUBROUTINE isPeriodic(DATA)
    TYPE(donnees),   INTENT(inout)   :: DATA

    SELECT CASE(DATA%test)
    CASE(0,11,12,20,19,14,40,41,42,27,501,-2,200,201,211,231,300)
      DATA%periodicBC = .TRUE.
    CASE default
      DATA%periodicBC = .FALSE.
    END SELECT
  END SUBROUTINE isPeriodic


END MODULE init_bc
