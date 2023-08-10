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

!*-----------------------------------------------------------------------
!*It contains the definition of PVar, type used for the solution (the vector u) in a single DoF (in a single subtimestep). In the code it is mostly used as  vector of coefficients (but sometimes also as values) of the conserved variables (but sometimes also of the primitive variables)
!*The type contains also some procedures like functions to compute the flux NB: They can be used only if we have the values of u (not with the coefficients)
!*-----------------------------------------------------------------------

!*A BIT OF DEFINITIONS
!*EULER EQUATION IN 1D
!*u_t+f_x(u)=0
!*where
!*   (r  )    (r*u    ) 
!* u=(r*u)  f=(r*u^2+p)
!*   (E  )    ((E+p)*u)
!*where
!*E=1/2*r*u^2+rho*eps
!*p=(gamma-1)*rho*eps
!*eps specific internal energy

!*EULER EQUATION IN 1D CONSERVED VARIABLES
!*  (r)     (m                              )
!*u=(m)   f=(m^2/r*(3-gamma)/2+(gamma-1)*E  )
!*  (E)     ([gamma*E+(1-gamma)/2*m^2/r]*m/r)
!*
!*m=r*u

!*PRIMITIVE VARIABLES
!*r, u, p


MODULE variable_def
  !------------------------------------------------------------------------------------------------
  ! MODULE SPECIFICALLY DESIGNED FOR 1D SYSTEM OF EULER EQUATIONS
  !------------------------------------------------------------------------------------------------
  ! This module collects all the information related to the following items:
  ! - equations of state: pressure, internal specific energy, speed of sound
  ! - conversion from conservative to primitive and viceversa
  ! - definition of the Jacobian of the considered system
  ! - definition of the Fluxes of the considered system
  ! - definition of the Jacobian in absolute values as: |J|=R* |lambda| *L (cf. AbsJacobian_eul)
  ! - definition of the Eigenvalues (cf. evalues_eul)
  ! - defintion of the Spectral Radius = max ( |lambda|)
  ! - definition of the Right eigenvectors
  ! - definition of the Left eigenvectors
  ! - definition of the positive and negative Jacobian

  
  USE algebra
  use precision
  !*USE utils

  IMPLICIT NONE

  INTEGER, PUBLIC, PARAMETER:: n_vars = 3   ! number of primitive variables
  INTEGER, PUBLIC, PARAMETER:: n_dim  = 1   ! number of physical dimensions
  REAL(dp),         PARAMETER:: pi=ACOS(-1._DP) ! pi !*<-ALERT: It was not ._DP 
  REAL(dp), PUBLIC           :: gmm=1.4_DP      ! EOS parameters !*<-ALERT: It was not ._DP 


  ! Type for the vector of primitive variables
  TYPE, PUBLIC :: PVar
     INTEGER                              :: NVars = n_vars
     REAL(dp), DIMENSION(n_vars)           :: U !*REAL(DP) vector containing the coefficients/values of the conserved/primitive variables

   CONTAINS 
     PROCEDURE, PUBLIC:: flux            => flux_eul
     !     PROCEDURE, PUBLIC:: evalues         => evalues_eul
     PROCEDURE, PUBLIC:: spectral_radius => spectral_radius_eul
     PROCEDURE, PUBLIC:: rvectors        => rvectors_eul
     PROCEDURE, PUBLIC:: lvectors        => lvectors_eul
     PROCEDURE, PUBLIC:: Jacobian        => Jacobian_eul
     PROCEDURE, PUBLIC:: AbsJacobian     => AbsJacobian_eul
     !*PROCEDURE, PUBLIC:: Nmat            => Nmat_eul !*NOT NEEDED AND NOT EVEN PROPERLY DEFINED AT THE END IT IS SET 0._DP
  END TYPE PVar

  PRIVATE
  PUBLIC:: pEOS, epsEOS, convert_cons2prim, convert_prim2cons, InitializeParameters

CONTAINS

  !--------------------------
  ! p = EOS(rho,eps) 
  !*pressure p=(gamma-1)*r*eps
  !--------------------------
  FUNCTION pEOS(rho,eps) RESULT(p)
    REAL(dp), INTENT(in) :: rho, eps
    REAL(dp)             :: p
    p = (gmm-1._DP)*rho*eps
  END FUNCTION pEOS

  !--------------------------
  ! eps = EOS(rho,p)
  !*specific internal energy eps=p/[(gamma-1)*r]
  !--------------------------
  FUNCTION epsEOS(rho,p) RESULT(eps)
    REAL(dp), INTENT(in) :: rho, p
    REAL(dp)             :: eps
    eps = p/((gmm-1.0_DP)*rho) !*<-ALERT: It was not ._DP
  END FUNCTION epsEOS

  !--------------------------
  ! Speed of sound
  !*a=sqrt(gamma*p/r)
  !--------------------------
  FUNCTION SoundSpeed(rho,p) RESULT(a)
    REAL(dp), INTENT(in) :: rho, p
    REAL(dp)             :: a
    a = SQRT(gmm*p/rho)
  END FUNCTION SoundSpeed

  !-----------------------------------------------
  ! Calulate pressure from conservative variables
  !*p=(gamma-1)*(E-1/2*m^2/r)=(gamma-1)*r*eps
  !*eps=p/[(gamma-1)*r] computed in epsEOS
  !-----------------------------------------------
  FUNCTION getPressureFromCons(Q) RESULT(p)
    TYPE(Pvar), INTENT(in) :: Q
    REAL(dp)                :: rho, eps, p
    rho = Q%U(1)
    eps = ( Q%U(3) - 0.5_DP*Q%U(2)**2/Q%U(1) )/Q%U(1) !*<-ALERT: It was not ._DP
    !*eps=(E-1/2*m^2/r)/r
    p = pEOS(rho,eps) !*Function defined up
  END FUNCTION getPressureFromCons

  !---------------------------------------------
  ! Convert conservative variables to primitive
  !*RMK:
  !*conservative (r,m,E)
  !*primitive (r,u,p)
  !*m=r*u->u=m/r
  !*E=1/2*m^2/r+r*eps
  !*p=(gamma-1)*r*eps
  !---------------------------------------------
  FUNCTION convert_cons2prim(Q) RESULT(W)
    TYPE(Pvar), INTENT(in) :: Q !*Conservative
    TYPE(Pvar)             :: W !*Primitive
    REAL(dp)                :: eps
    ! rho
    W%U(1) = Q%U(1) !*r=r
    ! u
    W%U(2) = Q%U(2)/Q%U(1) !*u=m/r
    ! p
    eps    = Q%U(3)/Q%U(1) - 0.5_DP*( Q%U(2)**2 )/Q%U(1)**2 !*eps=(E-1/2*m^2/r)/r=E/r-1/2*m^2/r^2
    !*ALERT, do not put dp in the exponent if it is a natural
    W%U(3) = pEOS(W%U(1),eps) 
  END FUNCTION convert_cons2prim

  !---------------------------------------------
  ! Convert primitive variables to conservative
  !*RMK:
  !*primitive (r,u,p)
  !*conservative (r,m,E)
  !*m=r*u
  !*E=1/2*m^2/r+r*eps=1/2*r*u^2+r*eps
  !*p=(gamma-1)*r*eps
  !---------------------------------------------
  FUNCTION convert_prim2cons(W) RESULT(Q)
    TYPE(Pvar), INTENT(in) :: W !*primitive
    TYPE(Pvar)             :: Q !*conservative
    REAL(dp)                :: eps
    ! q1 = rho
    Q%U(1) = W%U(1) !*r=r
    ! q2 = rho*u
    Q%U(2) = W%U(1)*W%U(2) !*m=r*v
    ! q3 = E = rho*eps + 0.5*rho*u^2.0
    eps    = epsEOS(W%U(1),W%U(3))
    Q%U(3) = W%U(1)*eps + 0.5_dp*W%U(1)*W%U(2)**2 !*ALERT, do not put dp in the exponent if it is a natural
  END FUNCTION convert_prim2cons


!!!!-----------------------

  !*!*---------------------------------------------------
  !*!*normal Jacobian !*NOT NEEDED
  !*!*NB: n=(1) because we are 1d
  !*!*The function is inheritated by 2d and can be canceled
  !*!*---------------------------------------------------
  !*FUNCTION roe_eul(e,u,n) RESULT (J) 
  !*  ! evaluate a roe average to estimate a rough speed for 
  !*  ! Burman jump operator
  !*  CLASS(Pvar), INTENT(in):: e
  !*  REAL(dp),DIMENSION(N_Vars), INTENT(in):: u
  !*  REAL(dp), DIMENSION(n_dim), INTENT(in):: n ! here it will be a normal of norme 1
  !*  REAL(dp), DIMENSION(n_vars,n_vars):: J
  !*  REAL(dp), DIMENSION(n_vars,n_vars,n_dim):: JJ
  !*  REAL(dp),DIMENSION(n_dim):: v=0._dp
  !*  JJ=Jacobian_eul(e,v)
  !*  J(:,:) = JJ(:,:,1)*n(1)
  !*END FUNCTION roe_eul

  !*------------------------------------------------------------------------------
  !*flux (evaluated from the values of the conserved variables)
  !*------------------------------------------------------------------------------
  FUNCTION flux_eul(Var,x) RESULT(f)
    CLASS(PVar),                  INTENT(in) :: Var    ! vector of conservative variables
    REAL(dp),       DIMENSION(n_dim), INTENT(in) :: x !*NOT NEEDED OR AT LEAST OPTIONAL
    TYPE(PVar), DIMENSION(n_dim)             :: f
    REAL(dp) :: eps, p

    eps = Var%u(3)/Var%u(1) - 0.5_dp*( Var%u(2)/Var%u(1))**2 !*eps=(E-1/2*m^2/r)/r=E/r-1/2*m^2/r^2
    
    p = pEOS(Var%u(1),eps) !*pressure

    f(1)%u(1) = Var%u(2) !*m
    f(1)%u(2) = Var%u(2)**2/Var%u(1) + p !*m^2/r+p
    f(1)%u(3) = Var%u(2)*( Var%u(3) + p )/Var%u(1) !*(E+p)*u=(E+p)*m/r

  END FUNCTION flux_eul

  !*------------------------------------------------------------------
  !*Jacobian matrix df/du
  !*------------------------------------------------------------------
  FUNCTION Jacobian_eul(Var,x) RESULT(J)
    CLASS(Pvar),              INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim),   INTENT(in) :: x
    REAL(dp), DIMENSION(n_Vars,n_Vars,n_dim) :: J
    REAL(dp) :: Vi2, u, H, eps, p

    eps = Var%u(3)/Var%u(1) - 0.5_DP*( Var%u(2)**2 )/Var%u(1)**2
    p = pEOS(Var%u(1),eps)

    u = Var%u(2)/Var%u(1)

    Vi2 = u**2

    H = (Var%u(3) + p)/Var%u(1)

    J=0._DP
    J(1,:,1) = (/ 0.0_dp, 1.0_dp, 0.0_dp /) !*correct
    J(2,:,1) = (/ -Vi2 + 0.5_dp*(gmm-1._dp)*Vi2, (3._dp-gmm)*u, gmm-1._dp /) !*correct
    J(3,:,1) = (/ u*( 0.5_dp*(gmm-1._DP)*Vi2 - H ), H - (gmm-1._dp)*Vi2, gmm*u /) !*correct

  END FUNCTION Jacobian_eul

  !*----------------------------------------------------------------
  !*Absolute value of the Jacobian
  !*The absolute value of a matrix A is defined as follows
  !*If A=RDL where D is diagonal we define |A|=R|D|L where |D| is the diagonal matrix with entries the absolute values of the entries of D (which are the eigenvalues)
  !*NB: We have L,R and D(=lambda) for the Jacobian
  !*----------------------------------------------------------------
  FUNCTION AbsJacobian_eul(Var,x) RESULT(J)
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

    R=rvectors_eul(Var,nn)
    L=lvectors_eul(Var,nn)
    lambda=evalues_eul(var,xx,nn) !*D computed through a function down here
    lambda=ABS(lambda) !*Absolute value of D

    J(:,:,1)=MATMUL(R,MATMUL(lambda,L))

  END FUNCTION AbsJacobian_eul

  !*--------------------------------------------------------------------
  !*Diagonal matrix containing the eigenvalues of the Jacobian
  !*--------------------------------------------------------------------
  FUNCTION evalues_eul(Var,x,n) RESULT(lambda)
    ! eigenvalues: diagonal matrix. It is written as a matrix for ease of calculations
    CLASS(PVar),            INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x !*NOT NEEDED OR AT LEAST OPTIONAL
    REAL(dp), DIMENSION(n_dim)             :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda !*Matrix

    REAL(dp), DIMENSION(n_vars, n_vars, n_dim) :: J
    REAL(dp)    :: un, c, p, eps

    eps = Var%u(3)/Var%u(1) - 0.5_dp*( Var%u(2)**2)/Var%u(1)**2
    p   = pEOS(Var%u(1),eps)
    c   = SoundSpeed(Var%u(1),p)

    lambda = 0._DP
    un = ( Var%u(2)*n(1) )/Var%u(1)
    lambda(1,1) = un-c
    lambda(2,2) = un
    lambda(3,3) = un+c

  END FUNCTION evalues_eul

  !*--------------------------------------------------------------------------
  !*Spectral radius of the Jacobian matrix = Maximum absolute value of the eigenvalues
  !*--------------------------------------------------------------------------
  REAL(dp) FUNCTION spectral_radius_eul(Var,x,n)
    ! compute the maximum value of eigenvalues:
    ! max_i {lambda_ii}
    CLASS(PVar),            INTENT(in) :: Var
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: x
    REAL(dp), DIMENSION(n_dim), INTENT(in) :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)     :: lambda1, lambda2
    REAL(dp):: vi, c, p, eps
    REAL(dp),DIMENSION(n_Vars,n_Vars)      :: lambda

    lambda = evalues_eul(Var,x,n)
    spectral_radius_eul = MAXVAL(ABS(lambda))

  END  FUNCTION spectral_radius_eul

  !*-----------------------------------------
  !*Right eigenvectors (one per column)
  !*-----------------------------------------
  FUNCTION rvectors_eul(Q,n) RESULT(R)
    ! right e-vectors
    ! assume ||n||=1
    CLASS(PVar),           INTENT(in) :: Q
    REAL(dp), DIMENSION(n_dim)         :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars) :: R
    REAL(dp) :: rho, u, p, a, E, H, eps

    rho = Q%U(1)
    u   = Q%U(2)/Q%U(1)
    eps = Q%u(3)/Q%u(1) - 0.5_DP*( Q%u(2)**2 )/Q%u(1)**2
    p   = pEOS(Q%u(1),eps)
    a   = SoundSpeed(rho,p)
    E   = Q%U(3)
    H   = (E+p)/rho

    R(:,1) = (/ 1.0_dp, u-a, H-u*a /)
    R(:,2) = (/ 1.0_dp, u, 0.5_dp*u**2 /)
    R(:,3) = (/ 1.0_dp, u+a, H+u*a /)


  END FUNCTION rvectors_eul

  !*-----------------------------------------
  !*Left eigenvectors (one per row)
  !*-----------------------------------------
  FUNCTION lvectors_eul(Q,n) RESULT(L)
    ! left e-vectors
    ! assumes ||n||=1
    CLASS(PVar),            INTENT(in) :: Q
    REAL(dp), DIMENSION(n_dim)          :: n
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: L
    REAL(dp), DIMENSION(n_Vars,n_Vars)  :: R
    !*REAL(dp) :: rho, u, p, a, E, H, gmm1, eps
        R=rvectors_eul(Q,n)
        L=inverse(R)

  END FUNCTION lvectors_eul



  SUBROUTINE InitializeParameters(choicetest)
     INTEGER :: choicetest


  END SUBROUTINE


  

!*  FUNCTION Nmat_eul(e, n_ord, grad, x) RESULT (Nmat)
!*    CLASS(Pvar),                  INTENT(in) :: e
!*    INTEGER,                      INTENT(in) :: n_ord
!*    REAL(dp), DIMENSION(n_dim),       INTENT(in) :: x
!*    REAL(dp), DIMENSION(n_dim,n_ord), INTENT(in) :: grad
!*    REAL(dp), DIMENSION(n_vars,n_vars)           :: Nmat
!*    REAL(dp), DIMENSION(n_vars, n_vars, n_dim)   :: J
!*    INTEGER:: l
!*
!*    J= Jacobian_eul(e,x)
!*    Nmat=0._dp
!*    DO l=1, n_ord
!*       Nmat = Nmat + ABS( grad(1,l)*J(1,1,1) )
!*    ENDDO
!*    Nmat =0._dp!Inverse(Nmat)
!*
!*  END FUNCTION Nmat_eul


!*These two last functions are drafts of the positive and the negative parts of the Jacobian but they are not public and not used

!*  !*----------------------------------------------------------------
!*  !*"minimum" Jacobian
!*  !*----------------------------------------------------------------
!*  FUNCTION min_mat_eul(e, x, n, alpha) RESULT (Ap)
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
!*    lambda = evalues_eul(e,x,n)
!*    R      = rvectors_eul(e,n)
!*    L      = lvectors_eul(e,n)
!*    Ap     = MATMUL( R, MATMUL( MIN(lambda(:,:),alpha), L) )
!*
!*  END FUNCTION min_mat_eul

!*  !*----------------------------------------------------------------
!*  !*"maximum" Jacobian
!*  !*----------------------------------------------------------------
!*  FUNCTION max_mat_eul(e,x,n,alpha) RESULT (Ap)
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
!*    lambda = evalues_eul(e,x,n)
!*    R      = rvectors_eul(e,n)
!*    L      = lvectors_eul(e,n)
!*    Ap     = MATMUL( R, MATMUL( MAX(lambda(:,:),alpha), L) )
!*
!*  END FUNCTION max_mat_eul
!*  !--------------------------




END MODULE variable_def
