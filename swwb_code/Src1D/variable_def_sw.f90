MODULE variable_def
  !------------------------------------------------------------------------------------------------
  ! MODULE SPECIFICALLY DESIGNED FOR 1D SYSTEM OF SW EQUATIONS
  !------------------------------------------------------------------------------------------------
  ! This module collects all the information related to the following items:
  ! - equations of state: pressure, internal specific energy, speed of sound
  ! - conversion from conservative to primitive and viceversa
  ! - definition of the Jacobian of the considered system
  ! - definition of the Fluxes of the considered system
  ! - definition of the Jacobian in absolute values as: |J|=R* |lambda| *L (cf. AbsJacobian_sw)
  ! - definition of the Eigenvalues (cf. evalues_sw)
  ! - defintion of the Spectral Radius = max ( |lambda|)
  ! - definition of the Right eigenvectors
  ! - definition of the Left eigenvectors
  ! - definition of the positive and negative Jacobian


  !*Conservative variables
  !*H, Hv
  !*Primitive variables
  !*H, v


  
    USE algebra
    USE precision
    !*USE utils
  IMPLICIT NONE

  INTEGER, PUBLIC, PARAMETER:: n_vars = 2   ! number of primitive variables
  INTEGER, PUBLIC, PARAMETER:: n_dim  = 1   ! number of physical dimensions
  REAL(dp),         PARAMETER:: pi=ACOS(-1._DP) ! pi
  REAL(dp), PUBLIC           :: grav=9.81_DP      ! EOS parameters
  REAL(dp), PUBLIC           :: n_manning


  ! Type for the vector of primitive variables
  TYPE, PUBLIC :: PVar
     INTEGER                              :: NVars = n_vars
     REAL(dp), DIMENSION(n_vars)           :: U

   CONTAINS 
     PROCEDURE, PUBLIC:: flux            => flux_sw
     !     PROCEDURE, PUBLIC:: evalues         => evalues_eul
     PROCEDURE, PUBLIC:: spectral_radius => spectral_radius_sw
     PROCEDURE, PUBLIC:: rvectors        => rvectors_sw
     PROCEDURE, PUBLIC:: lvectors        => lvectors_sw
     PROCEDURE, PUBLIC:: Jacobian        => Jacobian_sw
     PROCEDURE, PUBLIC:: AbsJacobian     => AbsJacobian_sw
     PROCEDURE, PUBLIC:: AbsJacobianTricked     => AbsJacobian_Tricked_sw
     PROCEDURE, PUBLIC:: spectral_radius_Tricked => spectral_radius_Tricked_sw
     !*PROCEDURE, PUBLIC:: Nmat            => Nmat_sw !*NOT USED
  END TYPE PVar

  PRIVATE
  PUBLIC:: convert_cons2prim, convert_prim2cons, InitializeParameters

CONTAINS


  !--------------------------
  ! Speed of sound
  ! a=sqrt(gH)
  !--------------------------
  FUNCTION SoundSpeed(rho) RESULT(a)
    REAL(dp), INTENT(in) :: rho
    REAL(dp)             :: a
    a = SQRT(grav*rho)
  END FUNCTION SoundSpeed


  !---------------------------------------------
  ! Convert conservative variables to primitive
  ! H,Hv->H,v
  !---------------------------------------------
  FUNCTION convert_cons2prim(Q) RESULT(W)
    TYPE(Pvar), INTENT(in) :: Q !*Conservative H, Hv
    TYPE(Pvar)             :: W !*Primitive H, v
    !*REAL(dp)                :: eps !*Useless
    ! rho
    W%U(1) = Q%U(1)
    ! u
    W%U(2) = Q%U(2)/Q%U(1)
  END FUNCTION convert_cons2prim

  !---------------------------------------------
  ! Convert primitive variables to conservative
  ! H,v->H,Hv
  !---------------------------------------------
  FUNCTION convert_prim2cons(W) RESULT(Q)
    TYPE(Pvar), INTENT(in) :: W !*Primitive
    TYPE(Pvar)             :: Q !*Conservative
    !*REAL(dp)                :: eps
    ! q1 = rho
    Q%U(1) = W%U(1)
    ! q2 = rho*u
    Q%U(2) = W%U(1)*W%U(2)
  END FUNCTION convert_prim2cons


!!!!-----------------------

  !*!*---------------------------------------------------
  !*!*normal Jacobian !*NOT NEEDED
  !*!*NB: n=(1) because we are 1d
  !*!*The function is inheritated by 2d and can be canceled
  !*!*---------------------------------------------------
  !*FUNCTION roe_sw(e,u,n) RESULT (J)
  !*  ! evaluate a roe average to estimate a rough speed for 
  !*  ! Burman jump operator
  !*  CLASS(Pvar), INTENT(in):: e
  !*  REAL(dp),DIMENSION(N_Vars), INTENT(in):: u
  !*  REAL(dp), DIMENSION(n_dim), INTENT(in):: n ! here it will be a normal of norme 1
  !*  REAL(dp), DIMENSION(n_vars,n_vars):: J
  !*  REAL(dp), DIMENSION(n_vars,n_vars,n_dim):: JJ
  !*  REAL(dp),DIMENSION(n_dim):: v=0._DP
  !*  JJ=Jacobian_sw(e,v)
  !*  J(:,:) = JJ(:,:,1)*n(1)
  !*END FUNCTION roe_sw

  !*-------------------------------------------------------
  !* Flux (evaluated from conservative variables)
  !*-------------------------------------------------------
  FUNCTION flux_sw(Var,x) RESULT(f)
    CLASS(PVar),                  INTENT(in) :: Var    ! vector of conservative variables
    REAL(dp),       DIMENSION(n_dim), INTENT(in) :: x
    TYPE(PVar), DIMENSION(n_dim)             :: f

    F(1)%u(1) = Var%u(2) !*m
    F(1)%u(2) = Var%u(2)**2/Var%u(1) + grav*0.5_dp*Var%u(1)**2 !*m^2/rho+g*rho^2/2

  END FUNCTION flux_sw

  !*--------------------------------------------------------------------------------------
  !*Jacobian matrix df/du
  !*--------------------------------------------------------------------------------------
  FUNCTION Jacobian_sw(Var,x) RESULT(J)
    CLASS(Pvar),              INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim),   INTENT(in) :: x
    REAL(dp), DIMENSION(n_Vars,n_Vars,n_dim) :: J
    REAL(dp) :: Vi2, u

    u = Var%u(2)/Var%u(1) !*Velocity

    !*Vi2 = u**2 !*Square of the velocity

    J=0._DP !*Initialization
    J(1,:,1) = (/ 0.0_dp, 1.0_dp /) !*First row: 0, 1 
    J(2,:,1) = (/ -u**2 + grav*Var%u(1), 2.0_dp*u /) !*Second row: -m^2/rho^2+g*rho, 2*m/rho
   
  END FUNCTION Jacobian_sw

  !*----------------------------------------------------------------
  !*Absolute value of the Jacobian
  !*The absolute value of a matrix A is defined as follows
  !*If A=RDL where D is diagonal we define |A|=R|D|L where |D| is the diagonal matrix with entries the absolute values of the entries of D (which are the eigenvalues)
  !*NB: We have L,R and D(=lambda) for the Jacobian
  !*----------------------------------------------------------------
  FUNCTION AbsJacobian_sw(Var,x) RESULT(J)
    !compute abs value of jacobian matrix
    REAL(dp), DIMENSION(n_dim), PARAMETER:: nn=(/1.0_dp/),xx=(/1.0_dp/)
    CLASS(Pvar),              INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim),   INTENT(in) :: x
    REAL(dp), DIMENSION(n_Vars,n_Vars,n_dim) :: J
    !*REAL(dp) :: Vi2, u, H, eps, p !*USELESS
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: L
    REAL(dp), DIMENSION(n_Vars,n_Vars) :: R
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda
    INTEGER:: i

    R=rvectors_sw(Var,nn)
    L=lvectors_sw(Var,nn)
    lambda=evalues_sw(var,xx,nn)
    lambda=ABS(lambda)

    J(:,:,1)=MATMUL(R,MATMUL(lambda,L))

  END FUNCTION AbsJacobian_sw



  !*--------------------------------------------------------------
  !*eigenvalues of the Jacobian
  !*--------------------------------------------------------------
  FUNCTION evalues_sw(Var,x,n) RESULT(lambda)
    ! eigenvalues: diagonal matrix. It is written as a matrix for ease of calculations
    CLASS(PVar),            INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim)             :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda
    REAL(dp)    :: un, c

    c   = SoundSpeed(Var%u(1)) !*Speed of the sound

    lambda = 0._dp !*Initialization
    un = ( Var%u(2)*n(1) )/Var%u(1) !*Normal velocity
    lambda(1,1) = un-c
    lambda(2,2) = un+c

  END FUNCTION evalues_sw

  !*------------------------------------------------------------
  !*Spectral radius = Maximum absolute value of the eigenvalues
  !*------------------------------------------------------------
  REAL(dp) FUNCTION spectral_radius_sw(Var,x,n)
    ! compute the maximum value of eigenvalues:
    ! max_i {lambda_ii}
    CLASS(PVar),            INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: n
    REAL(dp),DIMENSION(n_Vars,n_Vars)      :: lambda

    lambda = evalues_sw(Var,x,n)
    spectral_radius_sw = MAXVAL(ABS(lambda))

  END  FUNCTION spectral_radius_sw

  !*-------------------------------------------------------------
  !*Right eigenvectors
  !*-------------------------------------------------------------  
  FUNCTION rvectors_sw(Q,n) RESULT(R)
    ! right e-vectors
    ! assume ||n||=1
    CLASS(PVar),           INTENT(in) :: Q
    REAL(dp), DIMENSION(n_dim)         :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars) :: R
    REAL(dp) :: rho, u,  a

    rho = Q%U(1)
    u   = Q%U(2)/Q%U(1)
    a   = SoundSpeed(rho)


    R(:,1) = (/ 1.0_dp, u-a /) !*First column
    R(:,2) = (/ 1.0_dp, u+a /) !*Second column


  END FUNCTION rvectors_sw

  !*--------------------------------------------------------------
  !*Left eigenvectors
  !*--------------------------------------------------------------
  FUNCTION lvectors_sw(Q,n) RESULT(L)
    ! left e-vectors
    ! assumes ||n||=1
    CLASS(PVar),            INTENT(in) :: Q
    REAL(dp), DIMENSION(n_dim)          :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: L
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: R
    !*REAL(dp) :: rho, u, p, a, E, H, gmm1, eps !*USELESS
        R=rvectors_sw(Q,n)
        L=inverse(R)

  END FUNCTION lvectors_sw


  !*----------------------------------------------------------------
  !*Absolute value of the Jacobian tricked not to have it singular
  !*The absolute value of a matrix A is defined as follows
  !*If A=RDL where D is diagonal we define |A|=R|D|L where |D| is the diagonal matrix with entries the absolute values of the entries of D (which are the eigenvalues)
  !*NB: We have L,R and D(=lambda) for the Jacobian
  !*----------------------------------------------------------------
  FUNCTION AbsJacobian_Tricked_sw(Var,x) RESULT(J)
    !compute abs value of jacobian matrix
    REAL(dp), DIMENSION(n_dim), PARAMETER:: nn=(/1.0_dp/),xx=(/1.0_dp/)
    CLASS(Pvar),              INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim),   INTENT(in) :: x
    REAL(dp), DIMENSION(n_Vars,n_Vars,n_dim) :: J
    !*REAL(dp) :: Vi2, u, H, eps, p !*USELESS
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: L
    REAL(dp), DIMENSION(n_Vars,n_Vars) :: R
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda
    INTEGER:: i
    REAL(dp),parameter :: epsilontrickk=1e-1_dp
    REAL(dp) :: epsilontrick

    R=rvectors_sw(Var,nn)
    L=lvectors_sw(Var,nn)
    lambda=evalues_sw(var,xx,nn)
    lambda=ABS(lambda)
    epsilontrick=maxval(lambda)*epsilontrickk
    
    DO i=1,SIZE(lambda,DIM=1)
       IF(lambda(i,i)<epsilontrick) THEN
          lambda(i,i)=(lambda(i,i)**2+epsilontrick**2)/(2._DP*epsilontrick)
       END IF
    END DO

    J(:,:,1)=MATMUL(R,MATMUL(lambda,L))

  END FUNCTION AbsJacobian_Tricked_sw


  !*------------------------------------------------------------
  !*Spectral radius tricked not to have it 0
  !*------------------------------------------------------------
  REAL(dp) FUNCTION spectral_radius_Tricked_sw(Var,x,n)
    ! compute the maximum value of eigenvalues:
    ! max_i {lambda_ii}
    REAL(dp), DIMENSION(n_dim), PARAMETER:: nn=(/1.0_dp/),xx=(/1.0_dp/)
    CLASS(PVar),            INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: n
    REAL(dp),DIMENSION(n_Vars,n_Vars)      :: lambda
    REAL(dp),parameter :: epsilontrickk=1e-1_dp
    REAL(dp) :: epsilontrick
    INTEGER:: i

    lambda=evalues_sw(var,xx,nn)
    lambda=ABS(lambda)
    epsilontrick=maxval(lambda)*epsilontrickk
    
    DO i=1,SIZE(lambda,DIM=1)
       IF(lambda(i,i)<epsilontrick) THEN
          lambda(i,i)=(lambda(i,i)**2+epsilontrick**2)/(2._DP*epsilontrick)
       END IF
    END DO

    spectral_radius_Tricked_sw = MAXVAL(ABS(lambda))

  END  FUNCTION spectral_radius_Tricked_sw





  SUBROUTINE InitializeParameters(choicetest)
     INTEGER :: choicetest

     SELECT CASE(choicetest)


		 !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
		 !*3.2.2 Short channel: 100 m
		 !*150 Supercritical case 
       CASE(150,440,441)
           n_manning=0.03_DP
           PRINT*, "Manning set", n_manning

       !*SWASHES but with friction
       !*415 SWASHES Supercritical with friction to get the steady state
       !*416 SWASHES Supercritical with friction from steady state
       !*425 SWASHES Subcritical with friction to get the steady state
       !*426 SWASHES Subcritical with friction from steady state
       !*435 SWASHES Transcritical with friction to get the steady state

       CASE(415,416,425,426,435)
           n_manning=0.03_DP
           PRINT*, "Manning set", n_manning
       CASE DEFAULT
          PRINT*, "No parameters to initialize in variable_def, but Manning set to be 0 for security reasons"
          n_manning=0._DP
     END SELECT





  END SUBROUTINE


  
!*  !*NOT USED, in fact it is set to 0 in the end
!*  FUNCTION Nmat_sw(e, n_ord, grad, x) RESULT (Nmat)
!*    CLASS(Pvar),                  INTENT(in) :: e
!*    INTEGER,                      INTENT(in) :: n_ord
!*    REAL(dp), DIMENSION(n_dim),       INTENT(in) :: x
!*    REAL(dp), DIMENSION(n_dim,n_ord), INTENT(in) :: grad
!*    REAL(dp), DIMENSION(n_vars,n_vars)           :: Nmat
!*    REAL(dp), DIMENSION(n_vars, n_vars, n_dim)   :: J
!*    INTEGER:: l
!*
!*    J= Jacobian_sw(e,x)
!*    Nmat=0._dp
!*    DO l=1, n_ord
!*       Nmat = Nmat + ABS( grad(1,l)*J(1,1,1) )
!*    ENDDO
!*    Nmat =0._dp!Inverse(Nmat)
!*  END FUNCTION Nmat_sw

!*These two last functions are drafts of the positive and the negative parts of the Jacobian but they are not public and not used

!*  FUNCTION min_mat_sw(e, x, n, alpha) RESULT (Ap)
!*    ! in this A must be the result of tensor*vec where tensor are
!*    ! the jacobian evaluated for class e. It computes the negative part of A
!*    CLASS(PVar),            INTENT(in) :: e
!*    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
!*    REAL(dp),                   INTENT(in) :: alpha
!*    REAL(dp), DIMENSION(n_dim), INTENT(IN) :: n
!*    REAL(dp), DIMENSION(N_vars, N_vars)    :: Ap
!*    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda
!*    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: R
!*    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: L
!*
!*    lambda = evalues_sw(e,x,n)
!*    R      = rvectors_sw(e,n)
!*    L      = lvectors_sw(e,n)
!*    Ap     = MATMUL( R, MATMUL( MIN(lambda(:,:),alpha), L) )
!*
!*  END FUNCTION min_mat_sw
!*
!*
!*  FUNCTION max_mat_sw(e,x,n,alpha) RESULT (Ap)
!*    ! in this A must be the result of tensor*vec where tensor are
!*    ! the jacobian evaluated for class e. It computes the negative part of A
!*    CLASS(PVar),            INTENT(in) :: e
!*    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
!*    REAL(dp),                   INTENT(in) :: alpha
!*    REAL(dp), DIMENSION(n_dim), INTENT(in) :: n
!*    REAL(dp), DIMENSION(N_vars, N_vars)    :: Ap
!*    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda
!*    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: R
!*    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: L
!*
!*    lambda = evalues_sw(e,x,n)
!*    R      = rvectors_sw(e,n)
!*    L      = lvectors_sw(e,n)
!*    Ap     = MATMUL( R, MATMUL( MAX(lambda(:,:),alpha), L) )
!*
!*  END FUNCTION max_mat_sw
!*  !--------------------------



END MODULE variable_def
