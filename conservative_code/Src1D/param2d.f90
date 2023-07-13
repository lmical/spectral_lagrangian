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
MODULE param2d
  USE element_class
  USE variable_def
  USE arete_class
  use precision
  IMPLICIT NONE
  REAL(8), PARAMETER:: Pi=ACOS(-1._DP) !*NB: IMPORTANT, it cannot be deleted and it is correct for it to be here because it is a parameter. It will be always recalled from this module !*ALERT ._DP was missing
  
  !*--------------------------------------------
  !*The type of the mesh
  !*--------------------------------------------
  TYPE maillage
     INTEGER:: ndofs, ns !*Total number of DoFs in the mesh
     INTEGER:: nt
     !*nt is the number of elements in the mesh
     TYPE(element), DIMENSION(:),ALLOCATABLE:: e !*vector of elements (look in elements_1D)
     REAL(dp),DIMENSION(:),ALLOCATABLE:: aires !*vector of the inverse of the areas of the dual cells for the DeC
     INTEGER:: nsegmt !*Very likely to be the boundaries of the elements i.e. points in 1d (segments in 2d) since they are 201 for a Mesh of 200 elements
     TYPE(arete), DIMENSION(:), ALLOCATABLE:: edge !*vector of edges (look in aretes)
  END TYPE maillage

  !*---------------------------------------------
  !*The type containing the variables (PVar) in all the DoFs and in all the subtimesteps 
  !*Type of Var
  !*I'd say our "main" variable, the solution we look for
  !*---------------------------------------------
  TYPE variables
     REAL(dp):: dt
     !*INTEGER:: Ncells !*Var%Ncells=Mesh%ns !*PROBABLY NOT NEEDED !*Var%Ncells=Mesh%ns is always 4096, probably an old structure not used anymore
     TYPE(PVar), DIMENSION(:,:),ALLOCATABLE:: ua, up
     !*Both are vectors of PVar (look in variable_def) for the definition. Essentially in Pvar we store the conservative (or primitive) variables (values or coefficients) in a single node (of a single subtimestep)
     !*In both ua and up we have 2 indices that will be allocated in the following way (in the main) 
     !*(Mesh%ndofs,0:DATA%iordret-1)
     !*1st->DoFs
     !*2nd->subtimesteps RMK: The subtimesteps are M+1=order i.e. 0,1,...,order-1
     !*up->PREVIOUS VALUES
     !*ua->ACTUAL VALUES (we want to pass from up to ua by adding L^2 divided by the measure of the dual cell)
     TYPE(Pvar), DIMENSION(:),ALLOCATABLE:: un
     !*In un is stored the operator L^2 for the single subtimestep 
     !*un has a single index referred to the DoFs (Mesh%Ndofs)
     
     !*SO THAT in the end we impose
     !*var%ua(is,k)=Var%up(is,k)-Var%un(is)*Mesh%aires(is)
     !*RMK: In Mesh%aires we have the inverse of the measures of the dual cells
  END TYPE variables

  !*---------------------------------------------------------------------
  !*All that is passed in input in the data file and that is user dependent for example the order, the choice of basis functions, ecc...
  !*The type of DATA 
  !*---------------------------------------------------------------------
  TYPE donnees
     INTEGER:: iordret ! also defines the number of levels in the Dec !*order of accuracy
     !*RMK: In the Dec iterations=P=subtimesteps=M+1 (0:M)= order
     REAL(dp):: cfl !*cfl with respect to the meximum dt allowed
     INTEGER:: ktmax !*maximal number of iterations carried for each execution to prevent running forever if dt is too small
     REAL(dp):: tmax !*final time
     INTEGER:: ifre !*frequency of storing of the results of the iterations (every ifre iterations we store the partial result)
     INTEGER:: ischema !*scheme chosen -> 1=supg, 2=psi, 3=mix, 4: galerkin+jump, 5: psi+jump 6 blend+jump 7: psi+galerkin2??
     INTEGER :: ijump !*jump chosen -> 1=conserved variables, 2=eta
     INTEGER:: iter !*number of iteration of the DeC !*RMK: iterations=P=subtimesteps=M+1 (0:M)= order
     !*So at least iterations=order (they may be more even if it would be useless)
     INTEGER:: nt, itype 
     !*nt number of elements
     !*itype type of elements 1: P1, 2: B2, 3: P2, 4:P3, 5: B3
     REAL(dp):: Length, domain_left !*Not read as an input but initialized in init_bc depending on the test chosen
     REAL(dp):: temps !*????????????????????????????????????????????????????? 
     LOGICAL:: restart=.FALSE. !*Possibility to restart from a partial result of a previous simulation in reality here it is always disbled but in 2d (where the simulations are longer it is essential)
     REAL(dp):: alpha_jump !*Stabilization parameter on the first derivative (used whenever we use the jump stabilization)
     REAL(dp):: alpha_jump2 !*Stabilization parameter on the second derivative (used whenever we use the jump stabilization)
     INTEGER:: test !*test chosen
     LOGICAL:: mood !*special (extra) technique to stabilize more
     LOGICAL:: periodicBC !*if periodic bc or nor

  END TYPE donnees


  TYPE Proximity !*Element of the vector referred to a single DoF i
     INTEGER  :: numbofel !*number of elements containing i
     INTEGER, DIMENSION(:), ALLOCATABLE :: vecEl !*Elements containing i
     INTEGER  :: numbofDoFs !*number of DoFs in the elements containing i (without redoundances) 
     INTEGER, DIMENSION(:), ALLOCATABLE :: vecDoFs !*Global indices of those DoFs
  END TYPE Proximity
END MODULE param2d
