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
  USE preprocessing
  IMPLICIT NONE

CONTAINS

!---------------------------------------
! Setup domain
!---------------------------------------

  SUBROUTINE init_geom(DATA)
    TYPE(donnees), INTENT(inout):: DATA

    SELECT CASE(DATA%test)
    CASE(1)
       ! wave packet
       DATA%domain_left = -1.0_dp
       DATA%Length      = 3.0_dp

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

    Var%u(1) = 0.0_dp
    Var%u(2) = EXP( -bb*(x-0.5_dp)**2_dp ) * ( aa*COS(aa*x) - 2._dp*bb*(x-0.5_dp)*SIN(aa*x) )

  END FUNCTION IC

 !*---------------------------------------------------------------------------------------
 !*Initialize from scratch when it is needed
 !*---------------------------------------------------------------------------------------
 SUBROUTINE InitializeFromScratch(DATA,Mesh,Ucoeff)
    TYPE(donnees), INTENT(IN) :: DATA
    TYPE(maillage), INTENT(IN) :: Mesh
    TYPE(PVar), DIMENSION(:), INTENT(INOUT) :: UCoeff



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

             CASE default
                PRINT*, "Error in Perturbation", DATA%test
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
    CASE(1)
        ! periodic
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=0.0_dp!a0+a1
        Var%un(pN)=0.0_dp!a0+a1
        Var%un(p1)=a0+a1
        Var%un(pN)=a0+a1

    CASE(-1,-2)
        ! outflow
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=0.0_dp
        Var%un(pN)=0.0_dp

    CASE(-3)
        ! reflective
        Var%un(i1)%u(2) = 2.0_dp*Var%ua(k,i1)%u(2)
        Var%un(iN)%u(2) = 2.0_dp*Var%ua(k,iN)%u(2)

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



  SUBROUTINE isPeriodic(DATA)
    TYPE(donnees),   INTENT(inout)   :: DATA

    SELECT CASE(DATA%test)
    CASE(1)
      DATA%periodicBC = .TRUE.
    CASE default
      DATA%periodicBC = .FALSE.
    END SELECT
  END SUBROUTINE isPeriodic

!---------------------------------------
! Exact solution
!---------------------------------------

  FUNCTION exact_wave(DATA,x,t) RESULT (U_ex)
    TYPE(donnees), INTENT(inout)          :: DATA
    REAL(dp),                   INTENT(in) :: x
    REAL(dp)                               :: t
    TYPE(PVar)                            :: U_ex, U0, U01, U02, W
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: R
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: L
    REAL(dp), DIMENSION(n_dim)             :: n = (/1.0/)

    ! U0 only needed to evaluate eigenvectors (constant)
    U0 = IC(x,DATA)
    lambda = U0%evalues((/x/),n)
    R      = U0%rvectors(n)
    L      = U0%lvectors(n)

    ! calculate ICs shifted in time
    U01 = IC(x-lambda(1,1)*t,DATA)  ! evaluated at lambda(1,1)
    U02 = IC(x-lambda(2,2)*t,DATA)  ! evaluated at lambda(2,2)

    ! calculate solution in characteristic variables
    W%U(1) = L(1,1)*U01%U(1) + L(1,2)*U01%U(2)
    W%U(2) = L(2,1)*U02%U(1) + L(2,2)*U02%U(2)

    ! transter back to original variables
    U_ex%U = MATMUL(R,W%U)

!    print*, 'lambda = ', lambda
!    print*, 'R = ', R
!    print*, 'L = ', L
!    print*, 'U01  = ', U01%U
!    print*, 'U02  = ', U02%U
!    print*, 'W    = ', W%U
!    print*, 'U_ex = ', U_ex%U

  END FUNCTION exact_wave

END MODULE init_bc
