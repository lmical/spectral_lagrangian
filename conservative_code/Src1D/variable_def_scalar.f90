!*-----------------------------------------------------------------------
!*It contains the definition of PVar, type used for the solution (the vector u) in a single DoF (in a single subtimestep). In the code it is mostly used as  vector of coefficients (but sometimes also as values) of the conserved variables (but sometimes also of the primitive variables)
!*The type contains also some procedures like functions to compute the flux NB: They can be used only if we have the values of u (not with the coefficients)
!*-----------------------------------------------------------------------

MODULE variable_def
  USE algebra
  use precision

  IMPLICIT NONE

  INTEGER, PUBLIC, PARAMETER:: n_vars = 1   ! number of primitive variables
  INTEGER, PUBLIC, PARAMETER:: n_dim  = 1   ! number of physical dimensions
  REAL(dp),         PARAMETER:: pi=ACOS(-1._DP) ! pi
  REAL(dp), PUBLIC, PARAMETER:: a = 1._DP

  ! Type for the vector of primitive variables
  TYPE, PUBLIC :: PVar
     INTEGER                              :: NVars = n_vars
     REAL(dp), DIMENSION(n_vars)           :: U

   CONTAINS 
     PROCEDURE, PUBLIC:: flux            => flux_s
     !*PROCEDURE, PUBLIC:: evalues         => evalues_s
     PROCEDURE, PUBLIC:: spectral_radius => spectral_radius_s
     PROCEDURE, PUBLIC:: rvectors        => rvectors_s
     PROCEDURE, PUBLIC:: lvectors        => lvectors_s
     PROCEDURE, PUBLIC:: Jacobian        => Jacobian_s
     PROCEDURE, PUBLIC:: AbsJacobian     => AbsJacobian_s
  END TYPE PVar

  PRIVATE
  PUBLIC:: convert_cons2prim, convert_prim2cons


CONTAINS

!---------------------------------------------
! Convert conservative variables to primitive
!---------------------------------------------
  FUNCTION convert_cons2prim(Q) RESULT(W)
    TYPE(Pvar), INTENT(in) :: Q
    TYPE(Pvar)             :: W
    W%U = Q%U
  END FUNCTION convert_cons2prim

!---------------------------------------------
! Convert primitive variables to conservative
!---------------------------------------------
  FUNCTION convert_prim2cons(W) RESULT(Q)
    TYPE(Pvar), INTENT(in) :: W
    TYPE(Pvar)             :: Q
    Q%U = W%U
  END FUNCTION convert_prim2cons


!!!!-----------------------


  FUNCTION flux_s(Var,x) RESULT(f)
    CLASS(PVar),                  INTENT(in)  :: Var
    REAL(dp),    DIMENSION(n_dim), INTENT(in)  :: x
    TYPE(PVar), DIMENSION(n_dim)              :: f
    REAL(dp),    DIMENSION(n_vars,n_vars,n_dim):: JJ


    F(1)%u(1)=a*Var%u(1)

  END FUNCTION flux_s



  FUNCTION Jacobian_s(Var,x) RESULT(J)
    CLASS(Pvar),                 INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim),   INTENT(in) :: x
    REAL(dp), DIMENSION(n_Vars,n_Vars,n_dim) :: J

    J(1,1,1) =  a

  END FUNCTION Jacobian_s


  !*----------------------------------------------------------------
  !*Absolute value of the Jacobian
  !*The absolute value of a matrix A is defined as follows
  !*If A=RDL where D is diagonal we define |A|=R|D|L where |D| is the diagonal matrix with entries the absolute values of the entries of D (which are the eigenvalues)
  !*NB: We have L,R and D(=lambda) for the Jacobian
  !*----------------------------------------------------------------
  FUNCTION AbsJacobian_s(Var,x) RESULT(J)
    !compute abs value of jacobian matrix
    REAL(dp), DIMENSION(n_dim), PARAMETER:: nn=(/1.0_dp/),xx=(/1.0_dp/)
    CLASS(Pvar),              INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim),   INTENT(in) :: x
    REAL(dp), DIMENSION(n_Vars,n_Vars,n_dim) :: J
    !*REAL(dp) :: Vi2, u, H, eps, p !*NOT NEEDED
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: L
    REAL(dp), DIMENSION(n_Vars,n_Vars) :: R
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda
    INTEGER:: i

    R=rvectors_s(Var,nn)
    L=lvectors_s(Var,nn)
    lambda=evalues_s(var,xx,nn) !*D computed through a function down here
    lambda=ABS(lambda) !*Absolute value of D

    J(:,:,1)=MATMUL(R,MATMUL(lambda,L))

  END FUNCTION AbsJacobian_s


  FUNCTION evalues_s(Var,x,n) RESULT(lambda)
    ! eigenvalues: diagonal matrix. It is written as a matrix for ease of calculations
    CLASS(PVar),            INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim)             :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda

    REAL(dp), DIMENSION(n_vars, n_vars, n_dim) :: J
    INTEGER:: i

    lambda(1,1) = a;

  END FUNCTION evalues_s
  
!*  !!! Entropy flux for the scalar equation
!*    FUNCTION flux_scal_e(e, x) RESULT(f)
!*    CLASS(PVar),                  INTENT(in) :: e
!*    REAL(dp),    DIMENSION(n_dim), INTENT(in) :: x
!*    TYPE(PVar), dimension(n_dim)             :: f
!*! transport
!*    f(1)%U(1) =  0.5_dp*e%u(1)**2
!*  END FUNCTION flux_scal_e

  REAL(dp) FUNCTION spectral_radius_s(Var,x,n)
    ! compute the maximum value of eigenvalues:
    ! max_i {lambda_ii}
    CLASS(PVar),            INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: n
    REAL(dp),DIMENSION(n_Vars,n_Vars)      :: lambda


    lambda = evalues_s(Var,x,n)
    spectral_radius_s = MAXVAL(ABS(lambda))
    !*ABS(a)

  END  FUNCTION spectral_radius_s

  FUNCTION rvectors_s(e,n) RESULT(R)
    ! right e-vectors
    CLASS(PVar),        INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim)         :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars) :: R

    R(1,1) = 1._DP;
  END FUNCTION rvectors_s

  FUNCTION lvectors_s(e,n) RESULT(L)
    ! left e-vectors
    CLASS(PVar),         INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim)          :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: L

    L(1,1) =  1._DP; 

  END FUNCTION lvectors_s


!*  FUNCTION entropyProduct(v1,v2) RESULT(prod)
!*    CLASS(PVar), INTENT(in) :: v1,v2
!*    REAL(dp)                 :: prod
!*    prod = v1%U(1)*v2%U(1)
!*  END FUNCTION entropyProduct
  !--------------------------


  !--------------------------


END MODULE variable_def
