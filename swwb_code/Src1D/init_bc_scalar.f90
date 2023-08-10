MODULE init_bc

  USE param2d
  USE overloading
  USE precision
  USE preprocessing

  IMPLICIT NONE

CONTAINS

!---------------------------------------
! Setup domain
!---------------------------------------

  SUBROUTINE init_geom(DATA)
    TYPE(donnees), INTENT(inout):: DATA

    SELECT CASE(DATA%test)
    CASE(0,10)
       ! wave packet
       DATA%domain_left = 0.0_dp
       DATA%Length      = 2.0_dp
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
		CASE(0,10)
		  DATA%tmax = 0.5_dp
	    Var%u(1) = sin(ACOS(-1._dp)*x)*0.1_dp
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
    CASE(0)
        ! periodic
        p1=i1
        pN=iN
        a0=Var%un(p1)
        a1=Var%un(pN)
        Var%un(p1)=0._dp!a0+a1
        Var%un(pN)=0._dp!a0+a1
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

    CASE(10)
        ! Dirichlet
        Var%ua(i1,k)%u(1) = sin(ACOS(-1.)*(-a*DATA%temps))*0.1_dp 
        Var%ua(iN,k)%u(1) = sin(ACOS(-1.)*(2.0-a*DATA%temps))*0.1_dp 

    CASE(100)
       ! weak periodic
       p1=i1
       pN=iN
       a0=Var%un(p1)
       a1=Var%un(pN)
       Var%un(p1)%u(1)=Var%un(p1)%u(1)+Var%dt*0.5*( Var%ua(k,p1)%u(1)-Var%ua(k,pN)%u(1) )
       Var%un(pN)%u(1)=Var%un(pN)%u(1)+Var%dt*0.5*( Var%ua(k,p1)%u(1)-Var%ua(k,pN)%u(1) )


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

  FUNCTION exact_scalar(DATA,x,t) RESULT (U_ex)
    TYPE(donnees), INTENT(inout)          :: DATA
    REAL(dp),                   INTENT(in) :: x
    REAL(dp), INTENT(in)                   :: t
    TYPE(Pvar)                        :: U_ex
    U_ex = IC(x-a*t,DATA)

!    print*, 'lambda = ', lambda
!    print*, 'R = ', R
!    print*, 'L = ', L
!    print*, 'U01  = ', U01%U
!    print*, 'U02  = ', U02%U
!    print*, 'W    = ', W%U
!    print*, 'U_ex = ', U_ex%U

  END FUNCTION exact_scalar

  SUBROUTINE isPeriodic(DATA)
    TYPE(donnees),   INTENT(inout)   :: DATA

    SELECT CASE(DATA%test)
    CASE(0,100)
      DATA%periodicBC = .TRUE.
    CASE default
      DATA%periodicBC = .FALSE.
    END SELECT
  END SUBROUTINE isPeriodic


END MODULE init_bc
