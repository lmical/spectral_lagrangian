!==================================================
! Definition of variables for the 1D burgers equation
! u_tt - a u^2_x = 0
!==================================================
MODULE variable_def
  USE algebra
  USE precision

  IMPLICIT NONE
  INTEGER, PUBLIC, PARAMETER:: n_vars = 1   ! number of primitive variables
  INTEGER, PUBLIC, PARAMETER:: n_dim  = 1   ! number of physical dimensions
  REAL(dp),         PARAMETER:: pi=ACOS(-1._dp) ! pi
  REAL(dp), PUBLIC, PARAMETER:: a = 1.0_dp

  ! Type for the vector of primitive variables
  TYPE, PUBLIC :: PVar
     INTEGER                              :: NVars = n_vars
     REAL(dp), DIMENSION(n_vars)           :: U

   CONTAINS 
     PROCEDURE, PUBLIC:: flux            => flux_w
     PROCEDURE, PUBLIC:: source          => source_b
     PROCEDURE, PUBLIC::  entropy_flux   => flux_scal_e
     PROCEDURE, PUBLIC:: evalues         => evalues_w
     PROCEDURE, PUBLIC:: spectral_radius => spectral_radius_w
     PROCEDURE, PUBLIC:: rvectors        => rvectors_w
     PROCEDURE, PUBLIC:: lvectors        => lvectors_w
     PROCEDURE, PUBLIC:: Jacobian        => Jacobian_w
  END TYPE PVar

  PRIVATE
  PUBLIC:: convert_cons2prim, entropyProduct


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


  FUNCTION flux_w(e,x) RESULT(f)
    CLASS(PVar),                  INTENT(in)  :: e
    REAL(dp),    DIMENSION(n_dim), INTENT(in)  :: x
    TYPE(PVar), DIMENSION(n_dim)              :: f
    REAL(dp),    DIMENSION(n_vars,n_vars,n_dim):: JJ
    JJ=Jacobian_w(e,x)
    F(1)%u(:)=a*0.5_dp *e%u*e%u

  END FUNCTION flux_w

  FUNCTION source_b(e,x,test) RESULT(s)
    CLASS(PVar),                INTENT(IN) :: e
    REAL(dp), DIMENSION(n_dim), INTENT(IN) :: x
    INTEGER ,                   INTENT(IN) :: test
    TYPE(PVar)                             :: s
    s%u = 0._dp
  END FUNCTION source_b


  !!! Entropy flux for the Burgers Equation
  
    FUNCTION flux_scal_e(e, x) RESULT(f)
    CLASS(PVar),                  INTENT(in) :: e
    REAL(dp),    DIMENSION(n_dim), INTENT(in) :: x
    TYPE(PVar), dimension(n_dim)             :: f

   f(1)%U(1) = a*1_dp/3_dp*e%u(1)**3_dp

  END FUNCTION flux_scal_e
  
  
  
  
  FUNCTION Jacobian_w(e,x) RESULT(J)
    CLASS(Pvar),                 INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim),   INTENT(in) :: x
    REAL(dp), DIMENSION(n_Vars,n_Vars,n_dim) :: J

    J(1,1,1) =  a*e%u(1)

  END FUNCTION Jacobian_w

  FUNCTION evalues_w(e,x,n) RESULT(lambda)
    ! eigenvalues: diagonal matrix. It is written as a matrix for ease of calculations
    CLASS(PVar),            INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim)             :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda

    REAL(dp), DIMENSION(n_vars, n_vars, n_dim) :: J
    INTEGER:: i

    lambda(1,1) = a*e%u(1);

  END FUNCTION evalues_w

  REAL(dp) FUNCTION spectral_radius_w(e,x,n)
    ! compute the maximum value of eigenvalues:
    ! max_i {lambda_ii}
    CLASS(PVar),            INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: n

    spectral_radius_w = a*e%u(1)

  END  FUNCTION spectral_radius_w

  FUNCTION rvectors_w(e,n) RESULT(R)
    ! right e-vectors
    CLASS(PVar),        INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim)         :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars) :: R

    R(1,1) = 1._dp;
  END FUNCTION rvectors_w

  FUNCTION lvectors_w(e,n) RESULT(L)
    ! left e-vectors
    CLASS(PVar),         INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim)          :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: L

    L(1,1) =  1.0_dp; 

  END FUNCTION lvectors_w


  FUNCTION entropyProduct(v1,v2) RESULT(prod)
    CLASS(PVar), INTENT(in) :: v1,v2
    REAL(dp)                 :: prod
    prod = v1%U(1)*v2%U(1)
  END FUNCTION entropyProduct
  !--------------------------


END MODULE variable_def
