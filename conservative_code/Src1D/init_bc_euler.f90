!!!  HIGH ORDER IN SPACE AND TIME DEFERRED CORRECTION (EXPLICIT) 
!!!     RESIDUAL DISTRIBUTION METHOD 
!!!  DESIGNED FOR THE SYSTEM GIVEN BY THE EULER EQUATIONS in 1D and 2D
!!!
!!!  Authors:
!!!  Remi Abgrall (University of Zurich),
!!!  Paola Bacigaluppi (University of Zurich),
!!!  Svetlana Tokareva (University of Zurich)
!!!  Institute of Mathematics and Institute of Computational Sciences
!!!  University of Zurich
!!!  July 10, 2018
!!!  Correspondance:	remi.abgrall@math.uzh.ch
!!!  ------------------------------------------
MODULE init_bc

  USE param2d
  USE overloading
  use precision
  USE utils

  IMPLICIT NONE

CONTAINS

  !---------------------------------------
  ! Setup domain
  !---------------------------------------
  !*--------------------------------------------------------------------------------------
  !*init_geom, depending on the choice of the test, sets
  !*DATA%domain_left = starting point (the left one) of the (one dimensional) domain
  !*DATA%Length = length of the domain      
  !*--------------------------------------------------------------------------------------
  SUBROUTINE init_geom(DATA)
    TYPE(donnees), INTENT(inout):: DATA

    SELECT CASE(DATA%test) !*Choice of the test
    CASE(0)
       ! Isentropic
       DATA%domain_left = -1.0_dp
       DATA%Length      =  2.0_dp

    CASE(1,11)
       ! Sod
       !*1 with gravity
       !*11 without gravity
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

    CASE(6)
       ! Woodward-Colella left
       DATA%domain_left = 0.0_dp
       DATA%Length      = 1.0_dp
    CASE(7)
       ! LeBlanc
       DATA%domain_left = 0.0_dp
       DATA%Length      = 9.0_dp
    CASE(8)
       !simetrie
       DATA%domain_left = 0.0_dp
       DATA%Length      = 10.0_dp
    CASE(9)
       ! 1-2-3
       DATA%domain_left = 0.0_dp
       DATA%Length      = 1.0_dp
    CASE(10)
       DATA%domain_left = 0.0_dp
       DATA%Length      = 1.0_dp
    CASE(12)
       DATA%domain_left = 0.0_dp
       DATA%Length      = 1.0_dp
    CASE default
       PRINT*, "Wrong test number for geom(), test = ", DATA%test
       STOP
    END SELECT

  END SUBROUTINE init_geom

  !---------------------------------------
  ! Set initial and boundary conditions
  !---------------------------------------
  !*--------------------------------------------------------------------------------------
  !*IC, depending on the choice of the test, sets 
  !*the initial conditions
  !*DATA%tmax = the final time, so the final time in  DATA is fake
  !*--------------------------------------------------------------------------------------
  !*Initial condition is given in values, not in coefficients
  !*   It is projected later in the main if we use Bernstein polynomials
  !*It is given directly in terms of the conservative variables
  !*--------------------------------------------------------------------------------------
  FUNCTION IC(x,DATA) RESULT(Var)
    REAL(dp),       INTENT(in)   :: x
    TYPE(donnees), INTENT(inout):: DATA
    TYPE(PVar) :: Var
    REAL(dp)       :: alph, p, u, ro
    REAL(dp):: ms, rhoR,uR,pR,aR,mr,S3, rho_star,u_star,p_star

    !---------------
    ! for Euler
    SELECT CASE(DATA%test) !*Choice of the test
    CASE(0)
        ! Convergence test: isentropic flow
        DATA%tmax = 0.1_dp
        gmm   = 3.0_dp
        alph = 0.9999995_dp
        Var%u(1) = 1._dp + 1._dp + alph*SIN(PI*x)
        Var%u(2) = 0._dp
        Var%u(3) = Var%u(1)*epsEOS(Var%u(1),(Var%u(1))**gmm)

    CASE(1,11)
        ! Sod
        DATA%tmax = 0.2_dp
        IF (x <= 0.5_dp) THEN
            Var%u(1) = 1._dp !* [kg/m^3]
            Var%u(2) = 0._dp !* [kg/m^3*m/s]
            Var%u(3) = Var%u(1)*epsEOS(Var%u(1),1._dp) !* [J/m^3]
            !*NB: we should add to Var%u(3) the velocity part [0.5_dp*Var%u(2)**2/Var%u(1)]. Luckily it is 0 in this case.
        ELSE
            Var%u(1) = 0.125_dp !* [kg/m^3]
            Var%u(2) = 0._dp !* [kg/m^3*m/s]
            Var%u(3) = Var%u(1)*epsEOS(Var%u(1),0.1_dp) !* [J/m^3]
            !*NB: we should add to Var%u(3) the velocity part [0.5_dp*Var%u(2)**2/Var%u(1)]. Luckily it is 0 in this case.
        END IF

    CASE(2)
        ! Shu-Osher
        DATA%tmax = 1.8_dp
        alph = 0.2_dp
        IF (x <= -4._dp) THEN
            Var%u(1) = 3.857143_dp! [kg/m^3]
            Var%u(2) = Var%u(1)*2.629369_dp !* [kg/m^3*m/s]
            Var%u(3) = Var%u(1)*epsEOS(Var%u(1),10.33333333333_dp) + 0.5_dp*Var%u(2)**2/Var%u(1) !* [J/m^3]
        ELSE
            Var%u(1) = 1._dp + alph*SIN(5._dp*x)
            Var%u(2) = 0._dp
            Var%u(3) = Var%u(1)*epsEOS(Var%u(1),1.0_dp) !*NB: we should add to Var%u(3) the velocity part [0.5_dp*Var%u(2)**2/Var%u(1)]. Luckily it is 0 in this case.
        END IF

    CASE(3)
        ! Woodward-Colella
        DATA%tmax = 0.038_dp
        IF (x <= 0.1_dp) THEN
            Var%u(1) = 1._dp !* [kg/m^3]
            Var%u(2) = 0._dp !* [kg/m^3*m/s]
            Var%u(3) = Var%u(1)*epsEOS(Var%u(1),10.0_dp**3) !* [J/m^3]
        ELSE
             IF (x > 0.1_dp .AND. x < 0.9_dp) THEN
                Var%u(1) = 1._dp !* [kg/m^3]
                Var%u(2) = 0._dp !* [kg/m^3*m/s]
                Var%u(3) = Var%u(1)*epsEOS(Var%u(1),10.0_dp**(-2)) !* [J/m^3] !*NB: we should add to Var%u(3) the velocity part [0.5_dp*Var%u(2)**2/Var%u(1)]. Luckily it is 0 in this case.
             ELSE
                Var%u(1) = 1._dp !* [kg/m^3]
                Var%u(2) = 0._dp !* [kg/m^3*m/s]
                Var%u(3) = Var%u(1)*epsEOS(Var%u(1),10.0_dp**2) !* [J/m^3] !*NB: we should add to Var%u(3) the velocity part [0.5_dp*Var%u(2)**2/Var%u(1)]. Luckily it is 0 in this case.
             END IF
        END IF

    CASE(4)
        ! 1D version of 2D channel
        Var%u(1) = 0.5_dp !* [kg/m^3]
        Var%u(2) = 0.0_dp !* [kg/m^3*m/s]
        Var%u(3) = 0.125_dp !* [J/m^3]

    CASE(5)
        ! Convergence test: isentropic vortex
        DATA%tmax = 0.5_dp
!        DATA%tmax = 0.25
        gmm   = 1.4_dp
        Var%u(1) = 1._dp + 0.1_DP*EXP(-(x)**2/(2._dp*0.01_dp))
!        Var%u(1) = 1.+EXP(-80.0d0*(x-0.4d0)**2)
        Var%u(2) = Var%u(1)*1._dp
        Var%u(3) = Var%u(1)*epsEOS(Var%u(1),1.0_dp) + 0.5_dp*Var%u(2)**2/Var%u(1)

    CASE(6)
        ! Woodward-Colella left
        DATA%tmax = 0.012_dp
       IF (x <= 0.5_dp) THEN
          Var%u(1) = 1._dp !* [kg/m^3]
          Var%u(2) = 0._dp !* [kg/m^3*m/s]
          Var%u(3) = Var%u(1)*epsEOS(Var%u(1),10.0_dp**3) !* [J/m^3]
       ELSE
          Var%u(1) = 1._dp ! [kg/m^3]
          Var%u(2) = 0._dp !* [kg/m^3*m/s]
          Var%u(3) = Var%u(1)*epsEOS(Var%u(1),10.0_dp**(-2)) !* [J/m^3]
       END IF
    CASE(7)
       ! LeBlanc
       DATA%tmax = 6._dp
       gmm   = 5._dp/3._dp
       IF (x <= 3_dp) THEN
          ro=1._dp   ; u=0._dp; p=0.1_dp*(gmm-1._dp)
       ELSE
          ro=0.001_dp; u=0._dp; p=1.e-07_DP*(gmm-1._dp) ! energie interne
       END IF
       Var%u(1)=ro
       Var%u(2)=u*ro
       Var%u(3)=ro*epsEOS(ro,p)+0.5_dp*ro*u*u
    CASE(8)
       DATA%tmax=0.003_dp
       gmm=5._dp/3._dp
       IF (x<=5._DP) THEN
          ro=1._dp; u=-1000_dp; p=0.1_dp
       ELSE
          ro=1._dp; u=1000_dp; p=0.1_dp
       ENDIF
       Var%u(1)=ro
       Var%u(2)=ro*u
       Var%u(3)=ro*epsEOS(ro,p)+0.5_dp*ro*u*u
    CASE(9) ! 1-2-3
       DATA%tmax=0.15_dp
       gmm=1.4_dp
       IF (x<=0.5_dp) THEN
          ro=1.0_dp
          u=-2.0_dp
          p=0.4_dp
       ELSE
          ro=1.0_dp
          u=2.0_dp
          p=0.4_dp
       ENDIF
       Var%u(1)=ro
       Var%u(2)=ro*u
       Var%u(3)=ro*epsEOS(ro,p)+0.5_dp*ro*u*u
    CASE(10)
       ! Pure Shock
       ms = 3._dp
       DATA%tmax = 10_dp**(-5)
       gmm=1.4_dp

       rhoR = 1.225_dp
       uR = 0._dp
       pR = 101325._dp
       
       aR = SQRT(gmm*pR/rhoR)
       mr=uR/aR
       S3 = ms*aR
       rho_star = rhoR*( (gmm+1._DP)*((mr-ms)**2) )/( (gmm-1._dp)*((mr-ms)**2)+2._dp  )
       u_star = (1._dp-(rhoR/rho_star)  )*S3+(ur*rhoR/rho_star)
       p_star = pR*( 2._dp*gmm*((mr-ms)**2)-(gmm-1._dp))/(gmm+1._dp)
    
        IF (x <= 0.5_dp) THEN
            Var%u(1) = rho_star ! [kg/m^3]
            Var%u(2) = u_star*Var%u(1)  ! [kg/m^3*m/s]
            Var%u(3) = Var%u(1)*epsEOS(Var%u(1),p_star) + 0.5_dp*Var%u(2)**2/Var%u(1) ! [J/m^3]
        ELSE
            Var%u(1) = rhoR
            Var%u(2) = uR
            Var%u(3) = Var%u(1)*epsEOS(Var%u(1),pR) + 0.5_dp*Var%u(2)**2/Var%u(1)
        END IF
    CASE(12)
       DATA%tmax = 1._dp
       Var%u(1)=EXP(-potential(x,DATA%test))
       Var%u(2)=0._DP
       Var%u(3)=Var%u(1)*epsEOS(Var%u(1),EXP(-potential(x,DATA%test)))
    CASE default
       PRINT*, "Wrong test number for Euler_IC, test = ", DATA%test
       STOP
    END SELECT

  END FUNCTION IC

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
  SUBROUTINE BC_Residuals(Var,i1,iN,DATA,k,n_theta,theta,alpha,dt,tn)
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

    CASE(0) !*PERIODIC: no combination
!!$        ! periodic
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=0._dp
        Var%un(pN)=0._dp
        Var%un(p1)=a0+a1 
        Var%un(pN)=a0+a1
        !*Clear, in both the first and the last node we only have the contribution coming from the intern of the mesh
        !*What we have in p1 is the internal contribution but also the external for pN and viceversa
        !*So we give both the contributions to the first and the last node because they are the same node in the end if we have periodic BCs
 
    CASE(1,11) !*NAIVE REFLECTION: no combination
       Var%un(i1)%u(2)=0._DP
       Var%un(iN)%u(2)=0._DP
    CASE(2,6,7,8,10,12) !*DIRICHLET CONSTANT IN TIME: no combination
       !*We nullify the residual at the boudary so we never update that value
       Var%un(i1)=0._dp
       Var%un(iN)=0.0_dp
       
    CASE(3) !*Honestly this BC doesn't make sense to me (and Davide):
        !*He's setting Var%un equal to a pointwise evaluation of the solution (up to the factor 2)
        !*I think this is not even conistent
        ! reflective
        Var%un(i1)%u(2) = 2.0_dp*Var%ua(k,i1)%u(2)
        Var%un(iN)%u(2) = 2.0_dp*Var%ua(k,iN)%u(2)
        

    CASE(4) !*STRONG DIRICHLET ON THE LEFT + NOTHING ON THE RIGHT
        ! inflow left + outflow right
        !*STRONG DIRICHLET ON THE LEFT CALLED AFTER THE UPDATING 

    CASE(5) !*PERIODIC, see above !*No combination
        ! periodic
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=0._dp!a0+a1
        Var%un(pN)=0._dp!a0+a1
        Var%un(p1)=a0+a1
        Var%un(pN)=a0+a1

    CASE default
       PRINT*, "Wrong test number for BC_Residuals, test = ", DATA%test
       STOP
    END SELECT

  END SUBROUTINE BC_Residuals

  !*--------------------------------------------------------------------------------------
  !*Strong_BC imposed after the updating
  !*--------------------------------------------------------------------------------------
  SUBROUTINE Strong_BC(Var,i1,iN,DATA,k,n_theta,theta,alpha,dt,tn)
    TYPE(variables), INTENT(inout):: Var
    INTEGER,         INTENT(in)   :: i1, iN, k
    TYPE(donnees),   INTENT(in)   :: DATA
    INTEGER, DIMENSION(5), INTENT(IN):: n_theta 
    REAL(dp),DIMENSION(0:4,1:4,1:5), INTENT(IN) :: theta
    REAL(dp),DIMENSION(4,2:5), INTENT(IN) :: alpha
    REAL(DP), INTENT(IN) :: dt
    REAL(DP), INTENT(IN) :: tn !*NB: DATA%time=tn+alpha(k,DATA%iordret)*dt

    SELECT CASE(DATA%test)
       CASE(0,1,2,5,6,7,8,10,11,12) !*NO STRONG -> Do nothing

       CASE(3) !*I do not understand
          PRINT*, "I do not understand these BCs"
          STOP
       CASE(4)
          Var%ua(k,i1)%u(1) = 1.0_dp
          Var%ua(k,i1)%u(2) = 3.0_dp
          Var%ua(k,i1)%u(3) = 4.5_dp + 1.7857_dp
       CASE DEFAULT
          PRINT*, "Wrong test number for Strong_BC, test = ", DATA%test
          STOP
       END SELECT

  END SUBROUTINE Strong_BC
   


  SUBROUTINE isPeriodic(DATA)
    TYPE(donnees),   INTENT(inout)   :: DATA

    SELECT CASE(DATA%test)
    CASE(0,5)
      DATA%periodicBC = .TRUE.
    CASE default
      DATA%periodicBC = .FALSE.
    END SELECT
  END SUBROUTINE isPeriodic


END MODULE init_bc
