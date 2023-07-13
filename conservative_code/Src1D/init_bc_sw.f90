MODULE init_bc

  USE param2d
  USE overloading
  USE precision
  USE utils

  IMPLICIT NONE

CONTAINS

  !---------------------------------------
  ! Setup domain
  !---------------------------------------
  SUBROUTINE init_geom(DATA)
    TYPE(donnees), INTENT(inout):: DATA

    SELECT CASE(DATA%test)
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

    CASE(1234) !Smooth periodic
      DATA%domain_left=0.0_dp
      DATA%Length  = 1._dp


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

    !*eta = 0._dp ! profile_bathymetry( DATA%test, (/x/) ) !*NOT USED

    !---------------
    ! for Euler
    SELECT CASE(DATA%test)
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
            Var%u(1) = 1._DP-bathymetry(DATA%test,x) ! [kg/m^3]
            Var%u(2) = 0._DP ! [m/s] !*NEIN, [kg/m^3*m/s]
        ELSE
            Var%u(1) = 0.05_DP -bathymetry(DATA%test,x)
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

    CASE(12,13,500)
      DATA%tmax = 1._DP
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

    CASE(1234) !Smooth periodic
      Var%u(1)=2._DP+COS(2._DP*Pi*x)
      Var%u(2)=1._DP*Var%u(1)

    CASE default
       PRINT*, "Wrong test number for SW_IC, test = ", DATA%test
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
    CASE(0,11,12,20,19,14,40,41,42,27,500,1234) !*No combination
!!$        ! periodic
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=a0+a1
        Var%un(pN)=a0+a1

    CASE(1,2,6,10,13,30,29,32) !*No combination, anyway this is not outflow because we should'n to anything when outflow
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


    SELECT CASE(DATA%test)
       CASE(0,11,12,20,19,14,40,41,42,27,1,2,6,10,13,30,29,32,3,5,500,1234) !*NO STRONG -> Do nothing

       
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
       CASE DEFAULT
          PRINT*, "Wrong test number for Strong_BC, test = ", DATA%test
          STOP
       END SELECT




  END SUBROUTINE Strong_BC



  SUBROUTINE isPeriodic(DATA)
    TYPE(donnees),   INTENT(inout)   :: DATA

    SELECT CASE(DATA%test)
    CASE(0,11,12,20,19,14,40,41,42,27,1234)
      DATA%periodicBC = .TRUE.
    CASE default
      DATA%periodicBC = .FALSE.
    END SELECT
  END SUBROUTINE isPeriodic


END MODULE init_bc
