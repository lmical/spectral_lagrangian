MODULE init_bc

  USE param2d
  USE overloading
  USE precision

  IMPLICIT NONE

CONTAINS

!---------------------------------------
! Setup domain
!---------------------------------------

  SUBROUTINE init_geom(DATA)
    TYPE(donnees), INTENT(inout):: DATA

    SELECT CASE(DATA%test)
    CASE(0,10,11)
       ! wave packet
       DATA%domain_left = 0.0_dp
       DATA%Length      = 1.0_dp
    CASE(1)
       ! shock moving
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

  FUNCTION IC(x,DATA) RESULT(Var)
    REAL(dp),       INTENT(in)   :: x
    TYPE(donnees), INTENT(inout):: DATA
    TYPE(PVar) :: Var

    SELECT CASE(DATA%test)
    CASE(0)
      DATA%tmax = 0.1_dp
      Var%u(1) = 0.5*sin(2._dp*ACOS(-1._dp)*x)+0.25_dp
    CASE(10,11)
      DATA%tmax = 0.1_dp
      Var%u(1) = 0.5*sin(2._dp*ACOS(-1._dp)*x)
    CASE(1)
      IF (x<0.5) THEN
        Var%u(1) = 1._dp 
      ELSE
        Var%u(1) = 0._dp
      ENDIF
    CASE default
       PRINT*, "Wrong test number for geom(), test = ", DATA%test
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
    CASE(0,11)
        ! periodic
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=0.0_dp!a0+a1
        Var%un(pN)=0.0_dp!a0+a1
        Var%un(p1)=a0+a1
        Var%un(pN)=a0+a1

    CASE(1,2)
        ! outflow
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=0._dp
        Var%un(pN)=0._dp
        Var%ua(p1,k)%u(1) = 1._dp

    CASE(3)
        ! reflective

    CASE(10)
        Var%un(i1)%u(1) = 0._dp
        Var%un(iN)%u(1) = 0._dp
        Var%ua(i1,k)%u(1) = 0._dp
        Var%ua(iN,k)%u(1) = 0._dp

    CASE default
       PRINT*, "Wrong test number for wave_1D, test = ", DATA%test
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

    !*For the moment I call BC_Residuals again
    !*In the future let's adjust it by isolating the strong imposition and removing it from BC_residuals
    CALL BC_Residuals(Var,i1,iN,DATA,k,n_theta,theta,alpha,dt,tn)

  END SUBROUTINE Strong_BC


!---------------------------------------
! Exact solution
!---------------------------------------

  FUNCTION exact_burgers(DATA,x,t) RESULT (U_ex)
    TYPE(donnees), INTENT(inout)          :: DATA
    REAL(dp),                   INTENT(in) :: x
    REAL(dp), INTENT(in)                   :: t
    TYPE(Pvar)                        :: U_ex
    U_ex = IC(x-a*t,DATA)
    print*, "exact solution to be implemented"
!    print*, 'lambda = ', lambda
!    print*, 'R = ', R
!    print*, 'L = ', L
!    print*, 'U01  = ', U01%U
!    print*, 'U02  = ', U02%U
!    print*, 'W    = ', W%U
!    print*, 'U_ex = ', U_ex%U

  END FUNCTION exact_burgers


  SUBROUTINE isPeriodic(DATA)
    TYPE(donnees),   INTENT(inout)   :: DATA

    SELECT CASE(DATA%test)
    CASE(0,11)
      DATA%periodicBC = .TRUE.
    CASE default
      DATA%periodicBC = .FALSE.
    END SELECT
  END SUBROUTINE isPeriodic


END MODULE init_bc
