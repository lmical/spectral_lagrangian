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
!==================================================
! Definition of variables for the 1D wave equation
! u_tt - a^2 u_xx = 0
! in the form of a 2 x 2 linear system
! v_t - a^2 w_x = 0
! w_t -     v_x = 0
! with v = u_t, w = u_x
!==================================================
MODULE variable_def
  USE algebra

  IMPLICIT NONE
  INTEGER, PUBLIC, PARAMETER:: n_vars = 2   ! number of primitive variables
  INTEGER, PUBLIC, PARAMETER:: n_dim  = 1  ! number of physical dimensions
  REAL(dp),         PARAMETER:: pi=ACOS(-1._dp) ! pi
  REAL(dp), PUBLIC, PARAMETER:: a = 1.0_dp
  REAL(dp), PUBLIC, PARAMETER:: aa = 100._dp, bb = 100._dp

  ! Type for the vector of primitive variables
  TYPE, PUBLIC :: PVar
     INTEGER                              :: NVars = n_vars
     REAL(dp), DIMENSION(n_vars)           :: U

   CONTAINS 
     PROCEDURE, PUBLIC:: flux            => flux_w
     PROCEDURE, PUBLIC:: source          => source_w
     PROCEDURE, PUBLIC:: evalues         => evalues_w
     PROCEDURE, PUBLIC:: spectral_radius => spectral_radius_w
     PROCEDURE, PUBLIC:: rvectors        => rvectors_w
     PROCEDURE, PUBLIC:: lvectors        => lvectors_w
     PROCEDURE, PUBLIC:: Jacobian        => Jacobian_w
     PROCEDURE, PUBLIC:: Nmat            => Nmat_w
  END TYPE PVar

  PRIVATE
  PUBLIC:: convert_cons2prim


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


  FUNCTION roe_w(e,u,n) RESULT (J)
    ! evaluate a roe average to estimate a rough speed for 
    ! Burman jump operator
    CLASS(Pvar),                INTENT(in):: e
    REAL(dp), DIMENSION(N_Vars), INTENT(in):: u
    REAL(dp), DIMENSION(n_dim),  INTENT(in):: n ! here it will be a normal of norme 1
    REAL(dp), DIMENSION(n_vars,n_vars):: J
    REAL(dp), DIMENSION(n_vars,n_vars,n_dim):: JJ
    REAL(dp), DIMENSION(n_dim):: v=0.
    JJ=Jacobian_w(e,v)
    J(:,:) = JJ(:,:,1)*n(1)

  END FUNCTION roe_w

  FUNCTION flux_w(e,x) RESULT(f)
    CLASS(PVar),                  INTENT(in)  :: e
    REAL(dp),    DIMENSION(n_dim), INTENT(in)  :: x
    TYPE(PVar), DIMENSION(n_dim)              :: f
    REAL(dp),    DIMENSION(n_vars,n_vars,n_dim):: JJ
    JJ=Jacobian_w(e,x)
    F(1)%u(:)=MATMUL(JJ(:,:,1),e%u)

  END FUNCTION flux_w


  FUNCTION source_w(e,x,test) RESULT(s)
    CLASS(PVar),                INTENT(IN) :: e
    REAL(dp), DIMENSION(n_dim), INTENT(IN) :: x
    INTEGER ,                   INTENT(IN) :: test
    TYPE(PVar)                             :: s
    s%u = 0._dp
  END FUNCTION source_w

  FUNCTION Jacobian_w(e,x) RESULT(J)
    CLASS(Pvar),                 INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim),   INTENT(in) :: x
    REAL(dp), DIMENSION(n_Vars,n_Vars,n_dim) :: J

    J=0.
    J(1,1,1) =  0.0_dp; J(1,2,1) = -a**2;
    J(2,1,1) = -1.0_dp; J(2,2,1) = 0.0_dp

  END FUNCTION Jacobian_w

  FUNCTION evalues_w(e,x,n) RESULT(lambda)
    ! eigenvalues: diagonal matrix. It is written as a matrix for ease of calculations
    CLASS(PVar),            INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim)             :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda

    REAL(dp), DIMENSION(n_vars, n_vars, n_dim) :: J
    INTEGER:: i

    lambda=0._dp
    lambda(1,1) = -a; lambda(2,2) = a;

  END FUNCTION evalues_w

  REAL(dp) FUNCTION spectral_radius_w(e,x,n)
    ! compute the maximum value of eigenvalues:
    ! max_i {lambda_ii}
    CLASS(PVar),            INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda

    lambda = evalues_w(e,x,n)
    spectral_radius_w = MAXVAL(ABS(lambda))

  END  FUNCTION spectral_radius_w

  FUNCTION rvectors_w(e,n) RESULT(R)
    ! right e-vectors
    CLASS(PVar),        INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim)         :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars) :: R

    R(:,:) = 0.0_dp
    R(1,1) = a;   R(1,2) = -a;
    R(2,1) = 1.0_dp; R(2,2) = 1.0_dp;

  END FUNCTION rvectors_w

  FUNCTION lvectors_w(e,n) RESULT(L)
    ! left e-vectors
    CLASS(PVar),         INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim)          :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: L

    L(:,:) = 0.0_dp
    L(1,1) =  1.0_dp; L(1,2) = a;
    L(2,1) = -1.0_dp; L(2,2) = a;
    L = L/(2.0_dp*a)

  END FUNCTION lvectors_w

  FUNCTION Nmat_w(e, n_ord, grad, x) RESULT (Nmat)
    CLASS(Pvar),                  INTENT(in) :: e
    INTEGER,                      INTENT(in) :: n_ord
    REAL(dp), DIMENSION(n_dim),       INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim,n_ord), INTENT(in) :: grad
    REAL(dp), DIMENSION(n_vars,n_vars)           :: Nmat
    REAL(dp), DIMENSION(n_vars, n_vars, n_dim)   :: J
    INTEGER:: l

    J= Jacobian_w(e,x)
    Nmat=0_dp
    DO l=1, n_ord
       Nmat = Nmat + ABS( grad(1,l)*J(1,1,1) )
    ENDDO
    Nmat =0._dp!Inverse(Nmat)

  END FUNCTION Nmat_w

  FUNCTION min_mat_w(e, x, n, alpha) RESULT (Ap)
    ! in this A must be the result of tensor*vec where tensor are
    ! the jacobian evaluated for class e. It computes the negative part of A
    CLASS(PVar),            INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp),                   INTENT(in) :: alpha
    REAL(dp), DIMENSION(n_dim), INTENT(IN) :: n
    REAL(dp), DIMENSION(N_vars, N_vars)    :: Ap
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: R
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: L

    lambda = evalues_w(e,x,n)
    R      = rvectors_w(e,n)
    L      = lvectors_w(e,n)
    Ap     = MATMUL( R, MATMUL( MIN(lambda(:,:),alpha), L) )

  END FUNCTION min_mat_w

  FUNCTION max_mat_w(e,x,n,alpha) RESULT (Ap)
    ! in this A must be the result of tensor*vec where tensor are
    ! the jacobian evaluated for class e. It computes the negative part of A
    CLASS(PVar),            INTENT(in) :: e
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp),                   INTENT(in) :: alpha
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: n
    REAL(dp), DIMENSION(N_vars, N_vars)    :: Ap
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: R
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: L

    lambda = evalues_w(e,x,n)
    R      = rvectors_w(e,n)
    L      = lvectors_w(e,n)
    Ap     = MATMUL( R, MATMUL( MAX(lambda(:,:),alpha), L) )

  END FUNCTION max_mat_w
  !--------------------------


END MODULE variable_def
