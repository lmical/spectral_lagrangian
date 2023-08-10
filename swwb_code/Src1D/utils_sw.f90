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
  USE PRECISION
  USE param2d
  USE algebra
  IMPLICIT NONE

  !*--------------------------------------------------------------------------------------
  !*We may want to use the bathymetry and its gradient also "outside" and not as analytical functions.
  !*Maybe we may want to have the values at the DoFs of the bathymetry and of the gradient (NB: with a gradient reconstruction). This happens for example when the bathymetry is not analytical and the gradient can just be reconstructed from the values
  !*This is why I SAVE here these PUBLIC structures
  !*They're filled/initialized down here in InitializeBathAndGrad and are at our disposal
  !*--------------------------------------------------------------------------------------
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: BathymetryAtDoFsUtils 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE, PUBLIC :: GradBathymetryAtDoFsUtils
  REAL(DP), PARAMETER :: gravity=9.81_DP

  !*Coefficients of the bathymetry
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: CoeffBathymetryAtDoFsUtils 

  !*Coefficients of the initial condition taken from scratch
  TYPE(Pvar), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: CoeffICUtils 
  !*Values of the initial condition taken from scratch
  TYPE(Pvar), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: ValICUtils 


CONTAINS


 !*--------------------------------------------------------------------------------------------------
 !*sourcetermfunction contains the different source terms depending on a parameter which is specified in the DATA
 !*-> 0 test for the elliptic problem
 !*-> 1 SW+bathymetry
 !*--------------------------------------------------------------------------------------------------
FUNCTION sourcetermfunction(x,U,GlobalIndexDoF,choicetest) RESULT(S)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x !*NB: No dimension because 1d and we have just a single real as cooridnate
    TYPE(Pvar), INTENT(IN) :: U !*Conservative
    INTEGER, INTENT(IN) :: GlobalIndexDoF
    INTEGER, INTENT(IN) :: choicetest
    INTEGER :: choicesourceterm
    TYPE(Pvar)             :: S !*Vector of the source term
    REAL(DP), DIMENSION(1) :: gradientbath

    SELECT CASE(choicetest)
        CASE(0)
           choicesourceterm=0

		 !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
		 !*3.2.2 Short channel: 100 m
		 !*150 Supercritical case   

       !*415 SWASHES Supercritical with friction to get the steady state
       !*425 SWASHES Subcritical with friction to get the steady state
       !*435 SWASHES Transcritical with friction to get the steady state


       !*416 SWASHES Supercritical with friction from steady state
       !*426 SWASHES Subcritical with friction from steady state
        CASE(150,440,441,415,416,425,426,435)
           choicesourceterm=2

        CASE DEFAULT !*A bit dangerous but there are other stops in case :)
           choicesourceterm=1
    END SELECT


    SELECT CASE(choicesourceterm) 
       CASE(1)            
          gradientbath=0._DP
#if(1==1)
          !*Analytical gradient
          gradientbath=deriv_bathymetry(choicetest,x)
#else
          !*Reconstructed gradient
          gradientbath=GradBathymetryAtDoFsUtils(GlobalIndexDoF,:)
#endif

          S%u(1) = 0._DP ! 0
          S%u(2) = -gravity*U%u(1)*gradientbath(1) ! -g*rho*b_x
       CASE(2)
          gradientbath=0._DP
#if(1==1)
          !*Analytical gradient
          gradientbath=deriv_bathymetry(choicetest,x)
#else
          !*Reconstructed gradient
          gradientbath=GradBathymetryAtDoFsUtils(GlobalIndexDoF,:)
#endif

          S%u(1) = 0._DP ! 0
          S%u(2) = -gravity*U%u(1)*gradientbath(1)-gravity*U%u(1)*n_manning**2*U%u(2)*ABS(U%u(2))/( U%u(1)**(10._DP/3._DP) ) ! -g*rho*b_x -g*h*S_f --- S_f=n^2*hu*ABS(hu)/h^(10/3)

       CASE DEFAULT
          PRINT*, "Wrong source term"
          STOP
    END SELECT
 END FUNCTION sourcetermfunction

  !*--------------------------------------------------------------------------------------
  !*Analytical bathymetry
  !*--------------------------------------------------------------------------------------
  REAL(DP) FUNCTION bathymetry(test,x)
    INTEGER, INTENT(IN):: test
    REAL(DP), INTENT(IN):: x
    REAL(DP) ::x0,r,b

   SELECT CASE(test) 
        CASE(-2)
            b=0._DP
        CASE(21,20,22,23,24,25,19,14,26)  !immersed smooth bumb
            x0=10._dp
            r=5._dp
            IF (x>x0-r .AND. x<x0+r ) THEN
                b = 0.2_dp*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
            ELSE
                b= 0._dp
            END IF
        CASE(0,2,3,4,5,6,10,11,12,13,40) !immersed 2 smooth bumbs
            b=0.1_dp*exp(-50._dp*(x-0.4_dp)**2) +0.04_dp*exp(-300._dp*(x-0.55_dp)**2) 
        CASE(301)
            b=0.1_dp*exp(-100._dp*(x-0.2_dp)**2) 
        CASE(30,29,34) !parabola
            b=0.5_dp*((x-2._dp)**2 -1._dp)
        CASE(31,33)
            b=0.17_dp*(15._dp-x)/10._dp+0.48_dp*(x-5._dp)/10._dp
        CASE(1,32,41,42,70,71,76,77)
            b=0._dp
        CASE(50,51)
            b = 0.5_dp+0.35_dp*SIN(6._dp*x*ACOS(-1._dp))
        CASE(62,63,67,68,97,98)
            b = 3._dp-x*0.002_dp !0.0026750049164762152
        CASE(60,61,65,66,95,96)
            b = 3._dp -x*0.01_dp
        CASE(27)
            IF (x>10._dp .AND. x<15._dp) THEN
              b=0.02_dp*(x-10._dp)
            ELSE
              IF (x>20._dp .AND. x<21._dp) THEN
                b=-0.1_DP*(x-21._dp)
              ELSE
                IF (x .GE. 15._dp .AND. x .LE. 20._dp) THEN
                  b= 0.1_dp
                ELSE 
                  b=0._dp
                END IF
              END IF
            END IF
        CASE(500)
           b=0.1_DP*EXP(-(x-0.5_DP)**2) 
        CASE(501)
           IF((x>-1._DP) .AND. (x<1._DP)) THEN
              b=EXP(-1._DP/(1._DP-x**2)) 
           ELSE
              b=0._DP
           ENDIF
       !*A New Approach for Designing Moving-Water Equilibria Preserving Schemes for the Shallow Water Equations - Yuanzhen Cheng, Alina Chertock, Michael Herty, Alexander Kurganov, Tong Wu, section 3 example 2
       !*100 supercritical transitory
       !*101 subcritical transitory
       !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 3.1.3 subcritical
       !*111 subcritical but with IC taken from scratch
       !*120 supercritical with IC taken from scratch but smooth
       !*121 subcritical with IC taken from scratch but smooth
       !*122 transcritical with IC taken from scratch but smooth
       CASE(100,101,111) 
           b=0._DP
           IF((x > 8._DP) .AND. (x < 12._DP)) THEN
              b=0.2_DP-0.05_DP*(x-10._DP)**2
           END IF      
        CASE(120,121,122)     
            x0=10._dp
            r=5._dp
            IF (x>x0-r .AND. x<x0+r ) THEN
                b = 0.2_dp*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
            ELSE
                b= 0._dp
            END IF



		 !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
		 !*3.2.2 Short channel: 100 m
		 !*150 Supercritical case 
       CASE(150,440,441)
			  PRINT*, "Analytical bathymetry not available, taken in input after the solution of an ODE"
           STOP




        !*Lake at rest perturbed
        CASE(200,201,210,211)     
            x0=10._dp
            r=5._dp
            IF (x>x0-r .AND. x<x0+r ) THEN
                b = 0.2_dp*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))
            ELSE
                b= 0._dp
            END IF
       CASE(230,231) 
           b=0._DP
           IF((x > 8._DP) .AND. (x < 12._DP)) THEN
              b=0.2_DP-0.05_DP*(x-10._DP)**2
           END IF      
		 !*Well-Balancing Via Flux Globalization: Applications to Shallow Water Equations with Wet/Dry Fronts
		 !*Alina Chertock, Alexander Kurganov, Xin Liu, Yongle Liu, and Tong Wu
		 CASE(300)
           b=0._DP
           IF((x > 0.1_DP) .AND. (x < 0.3_DP)) THEN
              b=1.25_DP*(COS(10._DP*Pi*(x-0.2_DP))+1._DP)
           END IF        
       !*---------------------------------------------------------------------------------
       !*PRINT*, Pi
       !*PRINT*, ACOS(-1.0)
       !*PRINT*, DACOS(-1._DP)
       !*PRINT*, COS(Pi)
       !*PRINT*, SIN(Pi)
       !*STOP
       !*---------------------------------------------------------------------------------        





		 !*SWASHES WITH BATHYMETRY NON SMOOTH
		 !*400 Lake at rest with an immersed bump 3.1.1
		 !*410 Supercritical
		 !*420 Subcritical
		 !*430 Transcritical
		 CASE(400,401,410,411,415,416,420,421,425,426,430,431,435)
           b=0._DP
           IF((x > 8._DP) .AND. (x < 12._DP)) THEN
              b=0.2_DP-0.05_DP*(x-10._DP)**2
           END IF      


		 !*IC taken from scratch
		 !*401 Lake at rest with an immersed bump 3.1.1 
		 !*411 Supercritical
		 !*421 Subcritical
		 !*431 Transcritical

       !*415 SWASHES Supercritical with friction to get the steady state
       !*425 SWASHES Supercritical with friction to get the steady state
       !*435 SWASHES Transcritical with friction to get the steady state

       !*416 SWASHES Supercritical with friction from steady state
       !*426 SWASHES Subcritical with friction from steady state


       CASE default
           PRINT*, "ERROR CASE NOT DEFINED IN bathymetry"
           STOP
    END SELECT
    bathymetry=b
  END FUNCTION bathymetry

  !*--------------------------------------------------------------------------------------
  !*Analytical gradient of the bathymetry
  !*--------------------------------------------------------------------------------------
  REAL(dp) FUNCTION deriv_bathymetry(test,x)
    INTEGER, INTENT(in):: test
    REAL(dp), INTENT(in):: x
    REAL(dp) ::x0,r,b
    REAL(dp) ::q, h_ex, dx_h_ex

    SELECT CASE(test) 
        CASE(-2)
            b=0._DP
        CASE(21,20,22,23,24,25,19,14,26)  !immersed smooth bumb
            x0=10._dp
            r=5._dp
            IF (x>x0-r .AND. x<x0+r ) THEN
                b = -0.2_dp*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))*&
                  &2._dp*(x-x0)/r**2/((1._DP-((x-x0)/r)**2)**2)
            ELSE
                b= 0._dp
            END IF

        CASE(0,2,3,4,5,6,10,11,12,13,40) !immersed 2 smooth bumbs
            b=-50._dp*(x-0.4_dp)*2._dp*0.1_dp*exp(-50._dp*(x-0.4_dp)**2) &
              &-300._dp*(x-0.55_dp)*2._dp*0.04_dp*exp(-300._dp*(x-0.55_dp)**2) 
        CASE(301)
            b=-100._dp*(x-0.2_dp)*2._dp*0.1_dp*exp(-100._dp*(x-0.2_dp)**2) 
        CASE(30,29,34) !parabola
            b=(x-2._dp)

        CASE(31,33)
            b=-0.17_dp/10._dp+0.48_dp/10._dp

        CASE(1,32,41,42,70,71,76,77)
            b=0._dp

        CASE(50,51)
            b = 6._dp*ACOS(-1._dp)*0.35_dp*COS(6._dp*x*ACOS(-1._dp))
        CASE(62,63,67,68,97,98)
            b = -0.002_dp !0.0026750049164762152
        CASE(60,61,65,66,95,96)
            b = -0.01_dp
        CASE(27)
            IF (x>10._dp .AND. x<15._dp) THEN
              b=0.02_dp
            ELSE
              IF (x>20._dp .AND. x<21._dp) THEN
                b=-0.1_dp
              ELSE
                IF (x .GE. 15._dp .AND. x .LE. 20._dp) THEN
                  b= 0._dp
                ELSE 
                  b=0._dp
                END IF
              END IF
            END IF
        CASE(500)
           b=0.1_DP*EXP(-(x-0.5_DP)**2)*(-2._DP*(x-0.5_DP))
        CASE(501)
           IF((x>-1._DP) .AND. (x<1._DP)) THEN
              b=EXP(-1._DP/(1._DP-x**2))*(1._DP/(1._DP-x**2)**2)*(-2._DP*x) 
           ELSE
              b=0._DP
           ENDIF
       !*A New Approach for Designing Moving-Water Equilibria Preserving Schemes for the Shallow Water Equations - Yuanzhen Cheng, Alina Chertock, Michael Herty, Alexander Kurganov, Tong Wu, section 3 example 2
       !*100 supercritical transitory
       !*101 subcritical transitory
       !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 3.1.3 subcritical
       !*111 subcritical but with IC taken from scratch
       !*120 supercritical with IC taken from scratch but smooth
       !*121 subcritical with IC taken from scratch but smooth
       !*122 transcritical with IC taken from scratch but smooth
        CASE(100,101,111) 
           b=0._DP
           IF((x > 8._DP) .AND. (x < 12._DP)) THEN
              b=-0.05_DP*2._DP*(x-10._DP)
           END IF           
        CASE(120,121,122)
           x0=10._dp
           r=5._dp
           IF (x>x0-r .AND. x<x0+r ) THEN
               b = -0.2_dp*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))*&
                 &2._dp*(x-x0)/r**2/((1._DP-((x-x0)/r)**2)**2)
           ELSE
               b= 0._dp
           END IF




		 !*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
		 !*3.2.2 Short channel: 100 m
		 !*150 Supercritical case 
       CASE(150,440,441)
			 q=2._DP
          h_ex=(4._DP/grav)**(1._DP/3._DP) * (1._DP-0.25_DP*EXP( -4._DP*(x/100._DP-0.5_DP)**2 ) )
          dx_h_ex=(4._DP/grav)**(1._DP/3._DP)*(-0.25_DP*EXP(-4._DP*(x/100._DP-0.5_DP)**2))* &
          & (-4._DP)*( 2._DP*(x/100._DP-0.5_DP) )/100._DP
		    b=(q**2/(grav*h_ex**3) - 1._DP)*dx_h_ex - n_manning**2*q*abs(q)/h_ex**(10._DP/3._DP)




        !*Lake at rest perturbed
        CASE(200,201,210,211)
           x0=10._dp
           r=5._dp
           IF (x>x0-r .AND. x<x0+r ) THEN
               b = -0.2_dp*EXP(1._dp-1._dp/(1._DP-((x-x0)/r)**2))*&
                 &2._dp*(x-x0)/r**2/((1._DP-((x-x0)/r)**2)**2)
           ELSE
               b= 0._dp
           END IF
       CASE(230,231) 
           b=0._DP
           IF((x > 8._DP) .AND. (x < 12._DP)) THEN
              b=-0.05_DP*2._DP*(x-10._DP)
           END IF
		 !*Well-Balancing Via Flux Globalization: Applications to Shallow Water Equations with Wet/Dry Fronts
		 !*Alina Chertock, Alexander Kurganov, Xin Liu, Yongle Liu, and Tong Wu
		 CASE(300)
           b=0._DP
           IF((x > 0.1_DP) .AND. (x < 0.3_DP)) THEN
              b=1.25_DP*(-SIN(10._DP*Pi*(x-0.2_DP)))*10._DP*Pi
           END IF   




		 !*SWASHES WITH BATHYMETRY NON SMOOTH
		 !*400 Lake at rest with an immersed bump 3.1.1
		 !*410 Supercritical
		 !*420 Subcritical
		 !*430 Transcritical
		 CASE(400,401,410,411,415,416,420,421,425,426,430,431,435)
           b=0._DP
           IF((x > 8._DP) .AND. (x < 12._DP)) THEN
              b=-0.05_DP*2._DP*(x-10._DP)
           END IF           


		 !*IC taken from scratch
		 !*401 Lake at rest with an immersed bump 3.1.1 
		 !*411 Supercritical
		 !*421 Subcritical
		 !*431 Transcritical

		 !*Same tests with friction
       !*415 SWASHES Supercritical with friction to get the steady state
       !*425 SWASHES Supercritical with friction to get the steady state
       !*435 SWASHES Transcritical with friction to get the steady state

       !*416 SWASHES Supercritical with friction from steady state     
       !*426 SWASHES Subcritical with friction from steady state     
   
        CASE default
            PRINT*, "ERROR CASE NOT DEFINED IN bathymetry"
            STOP

    END SELECT
    deriv_bathymetry=b
  END FUNCTION deriv_bathymetry

SUBROUTINE InitializeSourceTerm(choicetest,Mesh,VecProx)
    IMPLICIT NONE
    TYPE (maillage), INTENT(IN) :: Mesh  
    INTEGER, INTENT(IN) :: choicetest
    TYPE(Proximity), DIMENSION(:), INTENT(IN) :: VecProx


   
    SELECT CASE(choicetest)
    CASE(0)
       !*No potential required
    CASE DEFAULT !*A bit risky but there are other stops
       CALL InitializeBathAndGrad(choicetest,Mesh,VecProx)
    END SELECT

 END SUBROUTINE InitializeSourceTerm



 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !*Bathymetry at the DoFs and reconstruction of its gradient at the DoFs 
 !*NB: High-Order
 !*- Compute the gradient in every DoF
 !*- In the DoFs shared by more elements make the average between the gradients (in that DoF) in the different elements  
 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE InitializeBathAndGrad(choicetest,Mesh,VecProx)
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
    !*valloc local values of the bathymetry
    !*coeffloc local coefficients of the bathymetry
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


SELECT CASE(choicetest)


!*SWASHES: a compilation of Shallow Water Analytic Solutions for Hydraulic and Environmental Studies - Olivier Delestre, Carine Lucas, Pierre-Antoine Ksinant, Frédéric Darboux, Christian Laguerre, Thi Ngoc Tuoi Vo, Francois James, Stephane Cordier 
!*3.2.2 Short channel: 100 m
!*150 Supercritical case 
CASE(150,440,441)

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*PRINT*, "This should be allocated", ALLOCATED(BathymetryAtDoFsUtils), ALLOCATED(CoeffBathymetryAtDoFsUtils) 
    !*PRINT*, "This no", ALLOCATED(GradBathymetryAtDoFsUtils)
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          

    IF (ALLOCATED(GradBathymetryAtDoFsUtils) ) THEN !*Then there is sth wrong, it means that we have already initialized these structures and the sub is being recalled a second time. We stop the execution
       PRINT*, "You are trying to initialize again the bathymetry and its gradient. Not possible." 
       PRINT*, "(I'm in utils_sw, InitializeBathAndGrad)"
       STOP
    END IF
    !*Otherwise continue
    ALLOCATE(GradBathymetryAtDoFsUtils(Mesh%Ns,N_dim))
    GradBathymetryAtDoFsUtils=0._DP

CASE DEFAULT

    IF (ALLOCATED(BathymetryAtDoFsUtils) .OR. ALLOCATED(GradBathymetryAtDoFsUtils) &
       .OR. ALLOCATED(CoeffBathymetryAtDoFsUtils) ) THEN !*Then there is sth wrong, it means that we have already initialized these structures and the sub is being recalled a second time. We stop the execution
       PRINT*, "You are trying to initialize again the bathymetry and its gradient. Not possible." 
       PRINT*, "(I'm in utils_sw, InitializeBathAndGrad)"
       STOP
    END IF
    !*Otherwise continue

    ALLOCATE(BathymetryAtDoFsUtils(Mesh%Ns),GradBathymetryAtDoFsUtils(Mesh%Ns,N_dim),CoeffBathymetryAtDoFsUtils(Mesh%Ns))

    !*BathymetryAtDoFsUtils has 1 index 
    !*   -> DoF (1:Mesh%Ns)
    !*GradBathymetryAtDoFsUtils has 2 indices
    !*   -> DoF (1:Mesh%Ns)
    !*   -> dim (1:N_Dim)

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*PRINT*, "Shape", SHAPE(GradBathymetryAtDoFsUtils)
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !*Initialization
    BathymetryAtDoFsUtils=0._DP
    GradBathymetryAtDoFsUtils=0._DP
    CoeffBathymetryAtDoFsUtils=0._DP

    !*Now we fill the vector BathymetryAtDoFsUtils values of the bathymetry at the DoFs
    DO indt=1,Mesh%Nt !*Loop on the elements
        e=Mesh%e(indt) !*Quick referencing
        DO indi=1,e%Nsommets !*Loop on the DoFs
           BathymetryAtDoFsUtils(e%nu(indi))=bathymetry(choicetest,e%coor(indi))
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*SAFETY CHECK
           !*PRINT*, indi, e%nu(indi), e%coor(indi), bathymetry(choicetest,e%coor(indi))
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        END DO
    END DO
    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indi=1,Mesh%Ns
    !*   PRINT*, indi, BathymetryAtDoFsUtils(indi)
    !*END DO
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END SELECT


    
    !*Now we have the values at the DoFs
    !*To have the gradient at the DoFs we need to 
    !*1) Compute the gradient in the DoFs
    !*2) Averaging it thanks to VecProx

    !*In order to compute the gradient we need to pass to the coefficients (and we will store them)
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

    !*STRATEGY: We will first fill GradBathymetryAtDoFsUtils adding all the values of the gradient in the DoFs
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
          valloc(indi)=BathymetryAtDoFsUtils(e%nu(indi))
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, indi, e%nu(indi)
          !*PRINT*, valloc(indi)
          !*PRINT*, valloc(indi)-bathymetry(choicetest,e%coor(indi))
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
       DO indi=1,e%nsommets
          CoeffBathymetryAtDoFsUtils(e%nu(indi))=coeffloc(indi)
       END DO

       !*Now let's go from the coefficients to the gradient in the DoFs
       DO indi=1,e%nsommets !*Loop on the DoF where we compute the gradient
          DO indj=1,e%nsommets !*Loop on the basis functions
             GradBathymetryAtDoFsUtils(e%nu(indi),:)=GradBathymetryAtDoFsUtils(e%nu(indi),:)+&
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
          !*   PRINT*, GradBathymetryAtDoFsUtils(e%nu(indi),1)
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
          GradBathymetryAtDoFsUtils(indi,indd)=GradBathymetryAtDoFsUtils(indi,indd)/&
          & REAL(VecProx(indi)%numbofel,DP)
       END DO
    END DO

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indi=1,Mesh%Ns
    !*   PRINT*, indi, GradBathymetryAtDoFsUtils(indi,:)
    !*   !*IF(indi==260) STOP
    !*END DO  
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indt=1,Mesh%nt
    !*   e=Mesh%e(indt)
    !*   DO indi=1,e%Nsommets
    !*      !*PRINT*, e%coor(indi), BathymetryAtDoFsUtils(e%nu(indi))-bathymetry(choicetest,e%coor(indi))
    !*      PRINT*, e%coor(indi), GradBathymetryAtDoFsUtils(e%nu(indi),1),deriv_bathymetry(choicetest,e%coor(indi)),GradBathymetryAtDoFsUtils(e%nu(indi),1)-deriv_bathymetry(choicetest,e%coor(indi))
    !*   END DO
    !*END DO
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

    DEALLOCATE(valloc,coeffloc)
    DEALLOCATE(MaBaDoF,invMaBaDoF)

    
 END SUBROUTINE InitializeBathAndGrad


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

END MODULE  utils
