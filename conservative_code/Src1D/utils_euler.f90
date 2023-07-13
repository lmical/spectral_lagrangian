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

  !*--------------------------------------------------------------------------------------
  !*We may want to use the potential and its gradient also "outside" and not as analytical functions.
  !*Maybe we may want to have the values at the DoFs of the potential and of the gradient (NB: with a gradient reconstruction). This happens for example when the potential is not analytical and the gradient can just be reconstructed from the values
  !*This is why I SAVE here these PUBLIC structures
  !*They're filled/initialized down here in InitializePotAndGrad and are at our disposal
  !*--------------------------------------------------------------------------------------
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: PotentialAtDoFsUtils 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: GradPotentialAtDoFsUtils


  
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
        CASE(1,12)
           choicesourceterm=1
        CASE(0,2,3,4,5,6,7,8,9,10,11)
           choicesourceterm=2
        CASE DEFAULT
           PRINT*, "Wrong test choice"
           STOP
    END SELECT


    SELECT CASE(choicesourceterm) 
       CASE(1)            
          gradientphi=0._DP
#if(1==1)
          !*Analytical gradient
          gradientphi=grad_phi(x,choicetest)
#else
          !*Reconstructed gradient
          gradientphi=GradPotentialAtDoFsUtils(GlobalIndexDoF,:)
#endif

          S%u(1) = 0._DP ! 0
          S%u(2) = -U%u(1)*gradientphi(1) ! -rho*phi_x
          S%u(3) = -U%u(2)*gradientphi(1) ! -rho*u*phi_x
       CASE(2)
          S%u=0._DP
       CASE DEFAULT
          PRINT*, "Wrong source term"
          STOP
    END SELECT
 END FUNCTION sourcetermfunction

 FUNCTION grad_phi(x,choicetest) RESULT(gradient) 
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x !*NB: No dimension because 1d and we have just a single real as cooridnate
    INTEGER, INTENT(IN) :: choicetest
    INTEGER :: choicepotential 
    REAL(DP), DIMENSION(1) :: gradient 

    SELECT CASE(choicetest)
        CASE(1,12)
           choicepotential=1
        CASE DEFAULT
           PRINT*, "Wrong test choice"
           STOP
    END SELECT

    SELECT CASE(choicepotential) 
       CASE(1) !*x 
#if(1==1)
          gradient(1)=1._DP !*<- ALERT
#else
          gradient(1)=EXP(x)
#endif
       CASE DEFAULT
          PRINT*, "Wrong potential"
          STOP
    END SELECT
    
 END FUNCTION grad_phi
 
 FUNCTION potential(x,choicetest) RESULT(phi) 
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x !*NB: No dimension because 1d and we have just a single real as cooridnate
    INTEGER, INTENT(IN) :: choicetest
    INTEGER :: choicepotential 
    REAL(DP) :: phi 

    SELECT CASE(choicetest)
        CASE(1,12)
           choicepotential=1
        CASE DEFAULT
           PRINT*, "Wrong test choice"
           STOP
    END SELECT


    SELECT CASE(choicepotential) 
       CASE(1) !*x
#if(1==1)
          phi=x
#else
          phi=EXP(x)
#endif
       CASE DEFAULT
          PRINT*, "Wrong potential"
          STOP
    END SELECT
    
 END FUNCTION potential


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
    CASE(0,2,3,4,5,6,7,8,9,10,11)
       !*No potential required
    CASE(1,12)
       CALL InitializePotAndGrad(choicetest,Mesh,VecProx)
    CASE DEFAULT
          PRINT*, "Wrong test number in InitializeSourceTerm in Utils"
          STOP
    END SELECT

 END SUBROUTINE InitializeSourceTerm

 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !*Potential at the DoFs and reconstruction of its gradient at the DoFs 
 !*NB: High-Order
 !*- Compute the gradient in every DoF
 !*- In the DoFs shared by more elements make the average between the gradients (in that DoF) in the different elements  
 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE InitializePotAndGrad(choicetest,Mesh,VecProx)
    IMPLICIT NONE
    TYPE (maillage), INTENT(IN) :: Mesh  
    INTEGER, INTENT(IN) :: choicetest
    TYPE(Proximity), DIMENSION(:), INTENT(IN) :: VecProx !*Vector of proximity (for the gradient reconstruction)

    !*Loop variables
    INTEGER :: indi=0 !*Index for the loop on the DoFs
    INTEGER :: indj=0 !*Index for the loop on the DoFs
    INTEGER :: indt=0 !*Index for the loop on the triangles
    INTEGER :: indd=0 !*Index for the loop on the dimension
    INTEGER :: indc=0 !*Index for the loop on the components

    !*Support variables for the computation of the gradient
    REAL(DP), DIMENSION(:), ALLOCATABLE :: valloc, coeffloc
    !*valloc local values of the potential
    !*coeffloc local coefficients of the potential
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: MaBaDoF
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: invMaBaDoF
    !*For many reasons that I do not want to explain I generate the matrices to pass from coefficients to values here (essentially because e%base_at_DoFs is defined differently in the 2d and in the 1d code) 

    TYPE (element) :: e !*For the quick reference to the generic Mesh%e(indt) in the loops
 
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*STRUCTURE FOR SAFETY CHECK
    REAL(DP), DIMENSION(:), ALLOCATABLE :: supploc
    INTEGER :: numbtest=40
    INTEGER :: indl !*For a loop
    TYPE(PVAR), DIMENSION(:), ALLOCATABLE :: coeffsupp !*Support structure for the coefficients
    TYPE(PVAR) :: dersupp !*Support structure for the derivative
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    IF (ALLOCATED(PotentialAtDoFsUtils) .OR. ALLOCATED(GradPotentialAtDoFsUtils)) THEN !*Then there is sth wrong, it means that we have already initialized these structures and the sub is being recalled a second time. We stop the execution
       PRINT*, "You are trying to initialize again the potential and its gradient. Not possible." 
       PRINT*, "(I'm in utils_euler, InitializePotAndGrad)"
       STOP
    END IF
    !*Otherwise continue

    ALLOCATE(PotentialAtDoFsUtils(Mesh%Ns),GradPotentialAtDoFsUtils(Mesh%Ns,N_dim))

    !*PotentialAtDoFsUtils has 1 index 
    !*   -> DoF (1:Mesh%Ns)
    !*GradPotentialAtDoFsUtils has 2 indices
    !*   -> DoF (1:Mesh%Ns)
    !*   -> dim (1:N_Dim)

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*PRINT*, "Shape", SHAPE(GradPotentialAtDoFsUtils)
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !*Initialization
    PotentialAtDoFsUtils=0._DP
    GradPotentialAtDoFsUtils=0._DP


    !*Now we fill the vector PotentialAtDoFsUtils values of the potential at the DoFs
    DO indt=1,Mesh%Nt !*Loop on the elements
        e=Mesh%e(indt) !*Quick referencing
        DO indi=1,e%Nsommets !*Loop on the DoFs
           PotentialAtDoFsUtils(e%nu(indi))=potential(e%coor(indi),choicetest)
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*SAFETY CHECK
           !*PRINT*, indi, e%nu(indi), e%coor(indi)-potential(e%coor(indi),choicetest)
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        END DO
    END DO
    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indi=1,Mesh%Ns
    !*   PRINT*, indi, PotentialAtDoFsUtils(indi)
    !*END DO
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    !*Now we have the values at the DoFs
    !*To have the gradient at the DoFs we need to 
    !*1) Compute the gradient in the DoFs
    !*2) Averaging it thanks to VecProx

    !*In order to compute the gradient we need to pass to the coefficients
    !*To this end let's compute the matrix to pass from values to coefficients
   ALLOCATE(MaBaDoF(Mesh%e(1)%nsommets,Mesh%e(1)%nsommets),invMaBaDoF(Mesh%e(1)%nsommets,Mesh%e(1)%nsommets))
    MaBaDoF=0._DP
    invMaBaDoF=0._DP
    e=Mesh%e(1)
    DO indi=1,e%Nsommets
       DO indj=1,e%Nsommets
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, indj, e%x(:,indj)
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          MaBaDoF(indi,indj)=e%base(indj,e%x(:,indi)) !*row=point, column=basis
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, MaBaDoF(indi,indj)-e%base0(indi,indj), indi, indj
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO
    END DO
    invMaBaDoF=Inverse(MaBaDoF)
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK: 
    !*DO indi=1,e%Nsommets
    !*   DO indj=1,e%Nsommets
    !*      PRINT*, invMaBaDoF(indi,indj)-e%base1(indi,indj), indi, indj
    !*   END DO
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !*STRATEGY: We will first fill GradPotentialAtDoFsUtils adding all the values of the gradient in the DoFs
    !*NB: If a DoF is shared by more elements we will meet it as many times as the number of elements containing it
    !*In the end we divide the value in each DoF by the number of elements containing the DoF thanks to VecProx(indi)%numbofel

    ALLOCATE(valloc(Mesh%e(1)%nsommets),coeffloc(Mesh%e(1)%nsommets))
    !*SAFE INITIALIZATION
    valloc=0._DP !*Local values
    coeffloc=0._DP !*Local coefficients

    DO indt=1,Mesh%nt !*Loop on the elements
       e=Mesh%e(indt) !*Quick reference
       !*SAFE INITIALIZATION
       valloc=0._DP
       coeffloc=0._DP
       DO indi=1,e%nsommets !*Loop on the DoFs to extract the local values
          valloc(indi)=PotentialAtDoFsUtils(e%nu(indi))
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, indi, e%nu(indi)
          !*PRINT*, valloc(indi)
          !*PRINT*, valloc(indi)-Potential(e%coor(indi),choicetest)
          !*PRINT*, valloc(indi)-e%coor(indi)
          !*PRINT*
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO
       !*Now we need to pass to the local coefficients
       coeffloc=MATMUL(invMaBaDoF,valloc)
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*PRINT*, coeffloc-valloc
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*ALLOCATE(supploc(Mesh%e(1)%nsommets))
       !*supploc=0._DP
       !*DO indi=1,e%nsommets
       !*   DO indj=1,e%nsommets
       !*      supploc(indi)=supploc(indi)+invMaBaDoF(indi,indj)*valloc(indj)   
       !*   END DO
       !*   PRINT*, coeffloc(indi)-supploc(indi)!*, valloc(indi) 
       !*END DO
       !*DEALLOCATE(supploc)
       !*STOP
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !*Now let's go from the coefficients to the gradient in the DoFs
       DO indi=1,e%nsommets !*Loop on the DoF where we compute the gradient
          DO indj=1,e%nsommets !*Loop on the basis functions
             GradPotentialAtDoFsUtils(e%nu(indi),:)=GradPotentialAtDoFsUtils(e%nu(indi),:)+&
             e%gradient(indj,e%x(:,indi))*coeffloc(indj)
          END DO
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
          !*SAFETY CHECK
          !*IF (e%nu(indi)==numbtest) THEN
          !*   PRINT*, "Coor", e%coor(indi)
          !*   PRINT*, "Exact Grad", grad_phi(e%coor(indi),choicetest)
          !*   ALLOCATE(coeffsupp(Mesh%e(1)%Nsommets))
          !*   DO indl=1, Mesh%e(1)%Nsommets
          !*      coeffsupp(indl)=coeffloc(indl)
          !*   END DO
          !*   dersupp=0._DP
          !*   dersupp=e%eval_der(coeffsupp,e%x(:,indi))
          !*   PRINT*, "New values"
          !*   PRINT*, GradPotentialAtDoFsUtils(e%nu(indi),1)
          !*   PRINT*, "Because added"
          !*   PRINT*, dersupp
          !*   PRINT*
          !*   DEALLOCATE(coeffsupp)
          !*END IF
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO
    END DO
    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*PRINT*, "Elements", VecProx(numbtest)%numbofel
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !*Now let's deal with the averaging
    DO indi=1,Mesh%Ns !*Loop on the DoFs
       DO indd=1,n_dim !*Loop on the dimensions
          GradPotentialAtDoFsUtils(indi,indd)=GradPotentialAtDoFsUtils(indi,indd)/&
          & REAL(VecProx(indi)%numbofel,DP)
       END DO
    END DO

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indi=1,Mesh%Ns
    !*   PRINT*, indi, GradPotentialAtDoFsUtils(indi,:)
    !*   !*IF(indi==260) STOP
    !*END DO  
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indt=1,Mesh%nt
    !*   e=Mesh%e(indt)
    !*   DO indi=1,e%Nsommets
    !*      PRINT*, e%coor(indi), GradPotentialAtDoFsUtils(e%nu(indi),1)-grad_phi(e%coor(indi),choicetest),grad_phi(e%coor(indi),choicetest)
    !*   END DO
    !*END DO
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

    DEALLOCATE(valloc,coeffloc)
    DEALLOCATE(MaBaDoF,invMaBaDoF)

    
 END SUBROUTINE InitializePotAndGrad









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
