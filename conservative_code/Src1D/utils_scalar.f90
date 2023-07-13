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
MODULE utils
  USE param2d
  USE PRECISION
  USE algebra
  IMPLICIT NONE


  
CONTAINS


 !*--------------------------------------------------------------------------------------------------
 !*sourcetermfunction contains the different source terms depending on a parameter which is specified in the DATA
 !*-> 1 euler+gravity with phi=x
 !*NB: You may eventually set the source to 0 to get the conservation law
 !*--------------------------------------------------------------------------------------------------
 FUNCTION sourcetermfunction(x,U,GlobalIndexDoF,choicetest) RESULT(S)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x !*NB: No dimension because 1d and we have just a single real as cooridnate
    TYPE(Pvar), INTENT(IN) :: U !*Conservative
    INTEGER, INTENT(IN) :: GlobalIndexDoF
    INTEGER, INTENT(IN) :: choicetest
    INTEGER :: choicesourceterm
    TYPE(Pvar)             :: S !*Vector of the source term
    REAL(DP), DIMENSION(1) :: gradientphi 

    SELECT CASE(choicetest) 
        CASE(0)
           choicesourceterm=1
        CASE DEFAULT
           PRINT*, "Wrong test choice in source"
           STOP
    END SELECT


    SELECT CASE(choicesourceterm) 
       CASE(1)
          S%u=0._DP
       CASE DEFAULT
          PRINT*, "Wrong source term"
          STOP
    END SELECT
 END FUNCTION sourcetermfunction



 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !*Interface to initialize the source term in the main
 !*In principle I could call InitializePotAndGrad right in the main but the code would be Euler specific so I call this general function
 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE InitializeSourceTerm(choicetest,Mesh,VecProx)
    IMPLICIT NONE
    TYPE (maillage), INTENT(IN) :: Mesh  
    INTEGER, INTENT(IN) :: choicetest
    TYPE(Proximity), DIMENSION(:), INTENT(IN) :: VecProx
    
    SELECT CASE(choicetest)
    CASE(0)
       !*No potential required
    CASE DEFAULT
          PRINT*, "Wrong test number in InitializeSourceTerm in Utils"
          STOP
    END SELECT

 END SUBROUTINE InitializeSourceTerm










  ! TYPE(Pvar) FUNCTION fonc(x,DATA)
  !   REAL(dp),       INTENT(in)   :: x
  !   TYPE(donnees), INTENT(inout):: DATA

  !   REAL(dp):: y,z
  !   INTEGER:: p
  !   y=x!/lenght
  !   fonc = IC(x,DATA)
  ! END FUNCTION fonc

  REAL(dp) FUNCTION fbon(x)
    REAL(dp), INTENT(in)::x
    REAL(dp):: s, r, pi=ACOS(-1._dp)
    fbon=SIN(2.d0*pi*x)
  END FUNCTION fbon

  REAL(dp) FUNCTION fonc_d(x)
    REAL(dp), INTENT(in):: x
    !      integer, intent(in):: k
    REAL(dp):: y, pi=ACOS(-1._dp)
    INTEGER:: p
    p=1.d0
    fonc_d=2*pi*p*COS(2*pi*x*p)
  END FUNCTION fonc_d

  ! FUNCTION interpol(e,x,u) RESULT(v)
  !   REAL(dp), INTENT(in):: x
  !   TYPE(element), INTENT(in):: e
  !   TYPE(Pvar),DIMENSION(:), INTENT(in):: u
  !   REAL(dp),DIMENSION(n_vars):: v
  !   REAL(dp), DIMENSION(e%nsommets):: base
  !   INTEGER:: i
  !   REAL(dp), DIMENSION(2):: y
  !   y(1)=x; y(2)=1.d0-x
  !   v=0.0d0
  !   DO i=1,e%nsommets
  !      base(i)=e%base(i,y)
  !   ENDDO
  !   DO i=1,n_vars
  !      v(i)=SUM(base*u%u(i))
  !   ENDDO
  ! END FUNCTION interpol

END MODULE  utils
