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
       !*Test linear advection
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

    SELECT CASE(DATA%test) !*Choice of the test
    CASE(0)
        DATA%tmax = 1._dp
        Var%u(1)=COS(2._DP*Pi*x)
    CASE default
       PRINT*, "Wrong test number for scalar_ic, test = ", DATA%test
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
        ! periodic
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
       CASE(0) !*NO STRONG -> Do nothing
       CASE DEFAULT
          PRINT*, "Wrong test number for Strong_BC, test = ", DATA%test
          STOP
       END SELECT

  END SUBROUTINE Strong_BC
   


  SUBROUTINE isPeriodic(DATA)
    TYPE(donnees),   INTENT(inout)   :: DATA

    SELECT CASE(DATA%test)
    CASE(0)
      DATA%periodicBC = .TRUE.
    CASE default
      DATA%periodicBC = .FALSE.
    END SELECT
  END SUBROUTINE isPeriodic


END MODULE init_bc
