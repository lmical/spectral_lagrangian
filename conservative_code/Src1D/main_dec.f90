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
PROGRAM main
  USE param2d
  USE scheme
  USE overloading
  USE Model
  USE geometry
  USE utils
  USE postprocessing
  USE init_bc
  USE timestepping
  use precision
  USE preprocessing
  IMPLICIT NONE

  !*Not used anywhere
  !*INTEGER, DIMENSION(5)::loc_ndofs=(/2,3,3,4,4/)

  !*Coefficients for the Dec
  !*P = iterations = M+1 = subtimesteps (counting 0) = basis functions=accuracy (order)
  !*So subtimesteps 0,1,...,M=accuracy-1
  !* In the DeC we have (for equispaced nodes)
  !* L^2 
  !*    theta_p^l=\int_0^{l/M} phi_p(t)dt 
  !*       where p=0,1,...,M=accuracy-1 basis functions & l=1,2,...,M=accuracy-1 subtimesteps
  !* L^1
  !*    beta^l=l/M !*NB: Never used in the last version of L^1, they cancel
  !*
  !*In this code they are stored in the following way (before let's present the notation)
  !*NOTATION:
  !*k->order
  !*l->subtimestep->1,2,...,k-1 (NB: The first subtimestep at 0 is in practice never used)
  !*p->basis function->0,1,...,k-1 (NB: The one in 0 is used obviously)
  
  !*NB: ONLY THE THETA COEFFICIENTS WILL BE USED

  REAL(dp),DIMENSION(4,2:5):: alpha
  !*The coefficients alpha are what we usually refer as beta in the DeC
  !*The structure alpha(:,:) has 2 indices
  !*alpha(l,k)=l/(k-1) subtimestep
  !*First index=l->subtimestep l=1,2,...,k-1
  !*Second index=k->order
  !*alpha(0,k) not present, if present it would be 0
  REAL(dp),DIMENSION(4,4,3:5):: beta, gamma !*?????????? gamma, anyway not needed
  !*Beta are coefficients that we would find in the "ugly" L^1 i.e. the interval between two consecutive subtimesteps, you use it if you make an Euler for each subtimesteps (but NB: we do not use this L^1, we use the good one)
  !*beta(l,k-1,k)
  !*l=1,2,...,k-1
  !*beta(1,k-1,k)=alpha(1,k)
  !*beta(2,k-1,k)=alpha(2,k)-alpha(1,k) ecc...   !*beta(l,k-1,k)=alpha(l,k)-alpha(l-1,k)
  !*They are not even initialized for order 2 (and 1 which is down here even if the DeC doesn't make sense for order 1)
  !*I have no clue about what the gammas are
  !*They are defined as the beta coefficients and basically he assigns
  !*gamma(2:,k-1,k)=beta(2:,k-1,k)
  !*So he copies the betas from the second value of the first index on
  INTEGER, DIMENSION(5):: n_theta !*Number of basis functions for each order, basically the order because the basis functions are M+1=accuracy=k so 0,1,...,M
  !*n_theta(k)=k
  REAL(dp),DIMENSION(0:4,1:4,1:5):: theta
  !*theta(p,l,k)=int_0^alpha(l,k) phi_p(s)ds
  !*p->0,1,...,k-1 basis functions
  !*l->1,2,...,k-1 subtimesteps
  !*k order

  !*REAL(dp), PARAMETER:: s=SQRT(3.0_dp) !*I commented and nothing changes

  TYPE(maillage):: mesh !*Type maillage in param2d (because it is the same as 2d)
  TYPE(variables):: var, debug !*Type variable in param2d
  TYPE(donnees):: DATA !*Type donnees in param2d
  TYPE(element):: e, e1, e2 !*Type element in elements_1D
  TYPE(arete):: ed !*Type aretes in arete
  TYPE(PVar):: error_L1, error_L2, error_Linf !*Type PVar in variable_def_euler
  !*You will get to know the structures when you'll see how they are initialized and used
  TYPE(Pvar),DIMENSION(:,:), ALLOCATABLE:: temp_debug, temp_var, temp
  TYPE(Pvar),DIMENSION(:,:), ALLOCATABLE:: u, up,ua1,ua2,up1,up2, up_p, u_p
  TYPE(Pvar), DIMENSION(:), ALLOCATABLE :: res,residu, uu, difference
  TYPE(Pvar), DIMENSION(:), ALLOCATABLE :: u1, u2
  !*TYPE(Pvar), DIMENSION(:,:), ALLOCATABLE:: resJ !*NOT USED
  !*REAL(dp), DIMENSION(3):: x=0._dp !*NOT USED
  INTEGER:: nt,itype, jt, i, kt, is, l, k_inter, p1, p2, k, iseg, jt1, jt2, ll, kt0
  REAL(dp):: dx, dt0, dt, temps
  CHARACTER(LEN = 1024) :: maille
  !*INTEGER:: nb_args, n, lp, liter !*NOT USED
  !*INTEGER:: Impre !*AS IF IT WAS NOT USED
  !*TYPE(Pvar), DIMENSION(:),ALLOCATABLE:: u_b, u_c !*NOT NEEDED
  !*TYPE(Pvar),DIMENSION(1):: FL,FR !*NOT USED
  INTEGER :: iarg, nargs
  REAL(dp)    :: num, tn
  CHARACTER(len=32) :: arg
  CHARACTER(len=100) :: folder
  !*INTEGER, DIMENSION(:), ALLOCATABLE :: ShIndArray!*NOT USED


  INTEGER, DIMENSION(:), ALLOCATABLE:: fluxes_mood 
  INTEGER :: nflux

  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*STRUCTURES FOR MY TESTS
  INTEGER :: indi, indj, indt, indq
  REAL(DP), DIMENSION(2) :: x_test
  TYPE(Pvar), DIMENSION(:), ALLOCATABLE :: u_test
  TYPE(arete) :: artest
  TYPE(Pvar) :: utest
  REAL(dp), DIMENSION(n_dim) :: ntest
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !*Preprocessing structure: a vector of elements of dimension=number of nodes whose generic component i is referred to the node i and gives the elements containing that node and the nodes in those elements
  !*The definition is in preprocessing but is recalled here for clarity
  !*TYPE Proximity !*Element of the vector referred to a single DoF i
  !*   INTEGER  :: numbofel !*number of elements containing i
  !*   INTEGER, DIMENSION(:), ALLOCATABLE :: vecEl !*Elements containing i
  !*   INTEGER  :: numbofDoFs !*number of DoFs in the elements containing i (without redoundances) 
  !*   INTEGER, DIMENSION(:), ALLOCATABLE :: vecDoFs !*Global indices of those DoFs
  !*END TYPE Proximity
  TYPE(Proximity), DIMENSION(:), ALLOCATABLE :: VecProx !*Vector calculated in precomputation in elliptic module  


  REAL(DP), DIMENSION(:), ALLOCATABLE :: error1,error2,errorinf

  !*I guess this will be the name of the folder where the results will be stored 
  folder='TEST'


!!$-------------------------------------------------------------------
!!$-------------------------------------------------------------------
!!$-------------------------------------------------------------------

  !*Subroutine down here in the same program to compute the coefficients needed for the DeC
  CALL theta_alpha_beta()
  
  !--------------------------- READ INITIAL DATA --------------------------------------------------
  !*Now we open the dile DATA and we read the input data and store them in DATA, we remark that the structure of DATA which is a type donnees is in param2d 
  OPEN(1, file="Data/don1d")
  READ(1,*) DATA%nt !*number of elements (segments)
  ! pass mesh size from command line
  nargs = command_argument_COUNT()
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*COMMAND_ARGUMENT_COUNT returns the number of arguments passed on the command line when the containing program was invoked.
  !*PRINT*, nargs
  !*STOP
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*It is used to understand if we specified the number of elements in the command window when running the program
  !*If yes he sets the number of elements to be the input given
  IF (nargs > 0) THEN
     CALL get_command_ARGUMENT(1, arg)
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*get_command_ARGUMENT(n, value)
     !*!*It stores in value the n-th argument that was passed on the command line when the containing program was invoked.
     !*PRINT*, arg
     !*STOP
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     PRINT*
     READ(arg, *) iarg !*It stores arg in iarg 
     DATA%nt = iarg !*It stores finally iarg in DATA%nt, the new number of triangles
     !*!*NB: 
     !*!*->arg is a character array
     !*!*->iarg and DATA%nt are integers

     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*PRINT*, arg, iarg, DATA%nt
     !*STOP
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END IF
  PRINT*, 'N = ', DATA%nt !*It prints the number of triangles (either default or in input if specified)
  READ(1,*) DATA%itype !*Choice of the basis functions
  !  DATA%itype_b=DATA%itype
  READ(1,*)DATA%iordret, DATA%iter
  !* DATA%iordret is the order
  !* DATA%iter is the number of iterations of the DeC
  !*RMK: P=iterations=M+1=subtimesteps counting 0=order
  !*So in the optimal case we have that the number of iterations is the same as the order 
  !*In any case we must have 
  !*iterations>=order
  !*If this condition is not fulfilled then... print some warnings but do not stop
  IF (DATA%iter.LT.DATA%iordret) THEN
     PRINT*, "wrong Data%iter,Data%iordret"
     PRINT*,"Data%iter should be >= than Data%iordret"
     PRINT*, "Data%iordret", DATA%iordret 
     PRINT*, "Data%iter", DATA%iter
  ENDIF
  READ(1,*) DATA%ischema, DATA%ijump !*Choice of the scheme and of the jump
  !*alpha_jump are parameters in the stabilization with jump (for example in CIP i.e. jump stabilization i.e. Burman)
  READ(1,*) DATA%alpha_jump
  READ(1,*) DATA%alpha_jump2
  READ(1,*) DATA%cfl !*You should know what the cfl is XD
  READ(1,*) DATA%ktmax !*Maximal number of iterations
  READ(1,*) DATA%tmax !*Final time but NOT USED, it will be reinitialized in init_bc
  READ(1,*) DATA%ifre !*Frequency of storing of the results
  !*Every DATA%ifre timesteps we store the results
  READ(1,*) DATA%test !*Choice of the test
  READ(1,*) DATA%mood !*It is a logical input 
  !*If false it uses the standard things
  !*If true it uses also mood which is an a posteriori technique to stabilize more
  !*I will skip the debugging of this part
  IF(DATA%mood) THEN
     READ(1,*)jt
     ALLOCATE(fluxes_mood(jt))
     DO k=1, jt
        READ(1,*) fluxes_mood(k)
     ENDDO
  ELSE !*If non mood
    ALLOCATE(fluxes_mood(1)) !*allocatable integer array, it is allocated with only one element and this is set to be the scheme integer
    fluxes_mood(1) = DATA%ischema
  ENDIF
  CLOSE(1) !*Close the file just read

  !----------- SET geometry of the problem test --------------------------
  CALL isPeriodic(DATA) !*Subroutine in init_bc that, depending on the problem, set periodic boundary conditions or not in DATA%periodicBC which is a logical field
  CALL geom(Mesh,DATA) !*Subroutine in geometry 

  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*SOME EXPERIMENTS TO UNDERSTAND WHAT THE VARIABLES EXACTLY ARE
  !*PRINT*, "Ndofs", Mesh%ndofs !*Total number of DoFs in the mesh
  !*For example if we run B3 with 200 cells -> ../bin1D/main_dec.out 200
  !*For each cell we have 4 DoFs so 3*numberofcells+1 DoFs=3*200+1=601 DoFs
  !*PRINT*, "ns, nt", Mesh%ns, Mesh%nt
  !*Mesh%ns is always 4096 for any number of nodes, probably an old structure not used anymore !*PROBABLY NOT NEEDED, REMOVED
  !*Mesh%nt number of elements
  !*PRINT*, "SIZE aires", SIZE(Mesh%aires)
  !*Just like in the 2D code, Mesh%aires are the inverses of the areas of the dual cells for the mass lumping performed in the DeC. They are obviously as many as the DoFs.
  !*PRINT*, "nsegment", Mesh%nsegmt !*Very likely to be the boundaries of the elements i.e. points in 1d (segments in 2d) since they are 201 for a Mesh of 200 elements
  !*PRINT*, "SIZE elements", SIZE(Mesh%e) !*Elements, pretty obvious
  !*PRINT*, "SIZE edges", SIZE(Mesh%edge) !*boudaries of the elements, 201 for 200 elements
  !*STOP
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO jt=1,Mesh%nt
     mesh%e(jt)%type_flux=DATA%ischema !*Initialize the field type_flux with the schema integer
  ENDDO
  !*Impre=6 !*NOT USED 
  jt=Mesh%nt !*Number of elements

  call descriptor(folder,data) !*In postprocessing. Creation of the folder TEST where to store the results

  ! geometrie
  !initialisation
  !*Var%Ncells=Mesh%ns !*PROBABLY NOT NEEDED !*Var%Ncells=Mesh%ns is always 4096, probably an old structure not used anymore
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*PRINT*, Var%Ncells, Mesh%ns
  !*STOP
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  kt0=0 !*Counter of the timesteps initialized to 0
  DATA%temps=0._dp !*Actual time initialized to 0

  ! initialisation Initialization of many structures
  ALLOCATE(Var%ua(Mesh%ndofs,0:DATA%iordret-1), &
       & Var%up(Mesh%ndofs,0:DATA%iordret-1)&
       & , var%un(Mesh%ndofs),&
      & debug%ua(Mesh%ndofs,0:DATA%iordret-1) &
       &,debug%up(Mesh%ndofs,0:DATA%iordret-1)&
       & , debug%un(Mesh%ndofs), &
       & temp_var(Mesh%ndofs,0:DATA%iordret-1), &
       & temp(Mesh%ndofs,0:DATA%iordret-1), &
       & temp_debug(Mesh%ndofs,0:DATA%iordret-1))

  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!*INVESTIGATION OF THE STRUCTURES IN param2d
  !*!*1) maillage
  !*!*2) variables
  !*!*3) donnees
  !*!*
  !*!*
  !*!*1) maillage -> Mesh is of type maillage
  !*!*We already investigated it
  !*!*TYPE maillage
  !*!*   INTEGER:: ndofs !*Total number of DoFs in the mesh
  !*!*   INTEGER:: nt!*, ns !*Mesh%ns is always 4096, probably an old structure not used anymore !*NOT NEEDED
  !*!*   !*nt is the number of elements in the mesh 
  !*!*   TYPE(element), DIMENSION(:),ALLOCATABLE:: e !*vector of elements (look in elements_1D)
  !*!*   REAL(dp),DIMENSION(:),ALLOCATABLE:: aires !*vector of the inverse of the areas of the dual cells for the DeC
  !*!*   INTEGER:: nsegmt !*Very likely to be the boundaries of the elements i.e. points in 1d (segments in 2d) since they are 201 for a Mesh of 200 elements
  !*!*   TYPE(arete), DIMENSION(:), ALLOCATABLE:: edge !*vector of edges (look in aretes)
  !*!*END TYPE maillage
  !*!*
  !*!*2) variables -> Var is of type variables
  !*!*TYPE variables
  !*!*   REAL(dp):: dt
  !*!*   !*INTEGER:: Ncells !*Var%Ncells=Mesh%ns !*PROBABLY NOT NEEDED !*Var%Ncells=Mesh%ns is always 4096, probably an old structure not used anymore
  !*!*   TYPE(PVar), DIMENSION(:,:),ALLOCATABLE:: ua, up
  !*!*   !*Both are vectors of PVar (look in variable_def) for the definition. Essentially in Pvar we store the conservative (or primitive) variables (values or coefficients) in a single node (of a single subtimestep)
  !*!*   !*In both ua and up we have 2 indices that will be allocated in the following way (in the main) 
  !*!*   !*(Mesh%ndofs,0:DATA%iordret-1)
  !*!*   !*1st->DoFs
  !*!*   !*2nd->subtimesteps RMK: The subtimesteps are M+1=order i.e. 0,1,...,order-1
  !*!*   !*up->PREVIOUS VALUES
  !*!*   !*ua->ACTUAL VALUES (we want to pass from up to ua by adding L^2 divided by the measure of the dual cell (at each iteration))
  !*!*   TYPE(Pvar), DIMENSION(:),ALLOCATABLE:: un
  !*!*   !*In un is stored the operator L^2 for the single subtimestep 
  !*!*   !*un has a single index referred to the DoFs (Mesh%Ndofs)
  !*!*   
  !*!*   !*SO THAT in the end we impose
  !*!*   !*var%ua(is,k)=Var%up(is,k)-Var%un(is)*Mesh%aires(is)
  !*!*   !*RMK: In Mesh%aires we have the inverse of the measures of the dual cells
  !*!*END TYPE variables
  !*!*
  !*!*TYPE donnees
  !*!*   INTEGER:: iordret ! also defines the number of levels in the Dec !*order of accuracy
  !*!*   !*RMK: In the Dec iterations=P=subtimesteps=M+1 (0:M)= order
  !*!*   REAL(dp):: cfl !*cfl with respect to the meximum dt allowed
  !*!*   INTEGER:: ktmax !*maximal number of iterations carried for each execution to prevent running forever if dt is too small
  !*!*   REAL(dp):: tmax !*final time
  !*!*   INTEGER:: ifre !*frequency of storing of the results of the iterations (every ifre iterations we store the partial result)
  !*!*   INTEGER:: ischema !*scheme chosen -> 1=supg, 2=psi, 3=mix, 4: galerkin+jump, 5: psi+jump 6 blend+jump 7: psi+galerkin2??
  !*!*   INTEGER:: iter !*number of iteration of the DeC !*RMK: iterations=P=subtimesteps=M+1 (0:M)= order
  !*!*   !*So at least iterations=order (they may be more even if it would be useless)
  !*!*   INTEGER:: nt, itype 
  !*!*   !*nt number of elements
  !*!*   !*itype type of elements 1: P1, 2: B2, 3: P2, 4:P3, 5: B3
  !*!*   REAL(dp):: Length, domain_left !*Not read as an input but initialized in init_bc depending on the test chosen
  !*!*   REAL(dp):: temps !*????????????????????????????????????????????????????? 
  !*!*   LOGICAL:: restart=.FALSE. !*Possibility to restart from a partial result of a previous simulation in reality here it is always disbled but in 2d (where the simulations are longer it is essential)
  !*!*   REAL(dp):: alpha_jump !*Stabilization parameter on the first derivative (used whenever we use the jump stabilization)
  !*!*   REAL(dp):: alpha_jump2 !*Stabilization parameter on the second derivative (used whenever we use the jump stabilization)
  !*!*   INTEGER:: test !*test chosen
  !*!*   LOGICAL:: mood !*special (extra) technique to stabilize more
  !*!*   LOGICAL:: periodicBC !*if periodic bc or nor
  !*!*END TYPE donnees
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !*!*IN INVESTIGATING THE STRUCTURES IN param2d WE HAVE ENCOUNTERED SOME SUBSTRUCTURES THAT WERE IN OTHER MODULES AND WE DIDN'T INVESTIGATE FOR EXAMPLE
  !*!*In maillage
  !*!*1)element->in elements_1D
  !*!*2)edge->in aretes
  !*!*In variables we have many vectors of
  !*!*3)PVar->in variable_def
  !*!*OUR AIM NOW IS TO TEST THESE STRUCTURES


  !*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !*!*1)element
  !*!*Let's recall the fields of the structures which are in elements_1D
  !*!*INTEGER:: type_flux=-10
  !*!*INTEGER:: diag=-1, diag2=-1
  !*!*INTEGER:: nsommets, itype, nvertex ! nbre de dofs, type element:
  !*!*   !1->P1
  !*!*   !2->B2
  !*!*   !3->P2
  !*!*   !4->P3
  !*!*   !5->B3
  !*!*   !6->B4
  !*!*   !7->P4
  !*!*   !11->PGL1
  !*!*   !12->PGL2
  !*!*   !13->PGL3
  !*!*   !14->PGL4
  !*!*   REAL(dp),  DIMENSION(:), POINTER :: coor =>NULL() ! 
  !*!*   !   *-----*-----* for quadratic, *---*---*---* for cubic elements
  !*!*   !   1     3     2                1   3   4   2
  !*!*   REAL(dp), DIMENSION(:,:), POINTER:: x=> NULL()
  !*!*   INTEGER, DIMENSION(:), POINTER   :: nu =>NULL() ! local connectivity table, see above for location
  !*!*   ! For Bezier, this corresponds to the Greville points
  !*!*   REAL(dp)                           :: volume =0.0_dp  ! volume
  !*!*   REAL(dp),  DIMENSION(:), POINTER   :: n =>NULL()     ! external normals 
  !*!*   INTEGER                        :: log    ! logic element  : this for boundary conditions, if needed
  !*!*!!!!   quadrature de surface
  !*!*   REAL(dp),   DIMENSION(:,:),POINTER :: quad =>NULL()  ! point de quadrature 
  !*!*   REAL(dp),   DIMENSION(:)  ,POINTER :: weight=>NULL() ! poids 
  !*!*   INTEGER                            :: nquad=0  ! nbre de points de quadrature
  !*!*   REAL(dp),DIMENSION(:,:),POINTER:: base0=>NULL(),base1=>NULL()
  !*!*   INTEGER, DIMENSION(:), POINTER :: dof2ind
  !*!*!!! int nabla phi_sigma * phi_sigma'
  !*!*REAL(dp), DIMENSION(:,:), POINTER:: coeff,coeff_b=> NULL(), mass=> NULL()
  !*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!*START TESTS
  !*e=Mesh%e(8) <-Quick referencing
  !*PRINT*, "e%type_flux", e%type_flux !*-10 at the beginning then set to be equal to the scheme number
  !*PRINT*, "e%diag, e%diag2", e%diag, e%diag2 !*-1,-1
  !*PRINT*, "nsommets, itype, nvertex", e%nsommets, e%itype, e%nvertex
  !*e%nsommets DoFs per element !*4 for B3
  !*e%itype type of elements (from the input) !*5 B3 (just because of our definition)
  !*e%nvertex number of vertices !*Clearly always 2 in 1dimension
  !*PRINT*, "coor", e%coor !*Coordinates of the DoFs of the element. Clearly only an abscissa for each DoF in 1d
  !*NB: As in 2D the natural order is: at the beginning we have the boundary vertices (2 in  this case) and then the internal DoFs, all with increasing abscissa
  !*!*Example: left vertex, right vertex, internal from lowest to highest abscissa
  !*PRINT*, "x", e%x !*Barycentric coordinates of the DoFs
  !*PRINT*, "size(x)", SHAPE(e%x)
  !*!*B3 
  !*PRINT*, "x(:,1)", e%x(:,1), "Associated to ->", e%coor(1) 
  !*PRINT*, "x(:,2)", e%x(:,2), "Associated to ->", e%coor(2) 
  !*PRINT*, "x(:,3)", e%x(:,3), "Associated to ->", e%coor(3)
  !*PRINT*, "x(:,4)", e%x(:,4), "Associated to ->", e%coor(4) 
  !*!*Annoying notation. 2 indices x(:,:)
  !*!*Second one referred to the DoF of which we consider the barycentric coordinates 1:nsommets
  !*!*First one referred to the 2 braycentric coordinates 1:2
  !*!*x(:,1) -> barycentric coordinates of the first DoF (left vertex)
  !*!*BUT NB: 
  !*!*Let's focus on the first index
  !*!*1->right vertex
  !*!*2->left vertex
  !*!*So
  !*!*x(:,1)=(0,1) which is, up to my opinion, a bit annoying because the order selected for the DoFs is with increasing abscissa. It is not wrong but a bit nonlogical.
  !*!*x(:,2)=(1,0)
  !*!*x(:,3)=(0,3;0,6) !*Closest to the left vertex associated to 2
  !*!*x(:,4)=(0,6;0,3) !*Closest to the right vertex associated to 1
  !*!*ACTUALLY if you want a logical way to imagine it just imagine that with this notation, in the reference interval [0,1], we have
  !*!*x(1)->x because you would start from x even if it is associated to the second node
  !*!*x(2)->1-x associated to the first node
  !*PRINT*, "e%nu", e%nu !*Global indices of the DoFs 
  !*!*The global numeration is
  !*!* vertices 1,2,...,N+1 (N number of cells)
  !*!* internal DoFs with increasing abscissa
  !*!*So for B3 cell K we have 
  !*!*k, k+1, N+2k, N+2k+1
  !*PRINT*, "Mesh%e(1)%nu", Mesh%e(1)%nu !*Global indices of the DoFs 
  !*PRINT*, "e%volume", e%volume !*Area of the cell (length in 1d) 
  !*PRINT*, "e%n", e%n !*INTERNAL normals !*I'd say internal since it is (1,-1) and, RMK, the first vertex is the left, the second is the right one??????????
  !*!*e%log NOT NEEDED
  !*!*PRINT*, "e%log", e%log !*Boundary tag associated to the element !*NOT NEEDED, removed
  !*!*PRINT*, "Mesh%e(1)%log", Mesh%e(1)%log !*Boundary tag associated to the element
  !*!*PRINT*, "Mesh%e(Mesh%nt)%log", Mesh%e(Mesh%nt)%log !*Boundary tag associated to the element
  !*PRINT*, "e%quad", e%quad !*Quadrature points
  !*PRINT*, "SHAPE(e%quad)", SHAPE(e%quad) !*Quadrature points in barycentric coordinates
  !*PRINT*, "e%quad(:,1)", e%quad(:,1) 
  !*PRINT*, "e%quad(:,2)", e%quad(:,2) 
  !*PRINT*, "e%quad(:,3)", e%quad(:,3)
  !*PRINT*, "e%quad(:,4)", e%quad(:,4)
  !*!*Again the we have 2 indices e%quad(:,:)
  !*!*Second index referred to the quadrature point
  !*!*First index referred to the two barycentric coordinates
  !*!*And also in this case (I guess) if we focus on the first index
  !*!*1->right vertex
  !*!*0->left vertex
  !*PRINT*, "e%weight", e%weight !*Quadrature weights
  !*PRINT*, "e%nquad", e%nquad !*Number of quadrature points
  !*PRINT*, "e%base0", e%base0 
  !*PRINT*, "e%base1", e%base1 
  !*!*base0 matrix of the values of the basis functions at the physical DoFs
  !*!*base1 inverse of the transposed of base0
  !*PRINT*, "SHAPE(e%base0)", SHAPE(e%base0) 
  !*PRINT*, "SHAPE(e%base1)", SHAPE(e%base1)
  !*PRINT*, "e%dof2ind", e%dof2ind !*local numeration 
  !*!*1 left vertex
  !*!*2 right vertex
  !*!*3,4,... internal DoFs with increasing abscissa
  !*PRINT*, "e%coeff", e%coeff
  !*!*NOT USED, in fact not initialized and they can give segmentation fault !*PRINT*, "e%coeff_b", e%coeff_b
  !*PRINT*, "e%mass", e%mass 
  !*!*mass and coeff are the matrices phi*phi and phi*gradphi

  !*!*The class element contains also some private procedures
  !*!*PROCEDURE, PUBLIC:: aire=>aire_element !*Area of the element, length in 1D
  !*!*PROCEDURE, PUBLIC:: gradient=>gradient_element 
  !*!*PROCEDURE, PUBLIC:: gradient2=>gradient2_element
  !*!*PROCEDURE, PUBLIC:: average => average_element
  !*!*PROCEDURE, PUBLIC:: base=>base_element
  !*!*PROCEDURE, PUBLIC:: normale=>normale_element !*Inward normals
  !*!*PROCEDURE, PUBLIC:: quadrature=>quadrature_element
  !*!*PROCEDURE, PUBLIC:: base_ref=>base_ref_element
  !*!*PROCEDURE, PUBLIC:: eval_func=>eval_func_element
  !*!*PROCEDURE, PUBLIC:: eval_der=>eval_der_element
  !*!*PROCEDURE, PUBLIC:: eval_der2=>eval_der2_element
  !*!*PROCEDURE, PUBLIC:: der_sec=>der_sec_element
  !*!*PROCEDURE, PUBLIC:: eval_coeff=>eval_coeff_element
  !*!*PROCEDURE, PUBLIC:: der2_poly=>der2_poly_element

  !*PRINT*, "e%aire", e%aire() !*Area of the element, length in 1D
  !*PRINT*, "e%normale", e%normale() !*Inward normals
  !*!*e%normale(1)=1 !*Left vertex, positive .->
  !*!*e%normale(2)=-1 !*Right vertex, negative <-.
  

  !*!*base_element(k,x)
  !*!*k->basis function 
  !*!*x->barycentric coordinates
  !*!*RKM: they are numerated in the following way
  !*!*1 left vertex
  !*!*2 right vertex
  !*!*other internal nodes with increasing abscissa
  !*!*x->barycentric coordinate REAL(DP), DIMENSION(2)
  !*!*RMK:
  !*!*Barycentric coordinates of the left vertex (0,1)
  !*!*Barycentric coordinates of the right vertex (1,0)
  !*!*So the second barycentric coordinate is associated to the left/first vertex, the first one is associated to the right/second vertex 
  !*!*RMK: To visualize it in a logical way in the reference interval
  !*!*x(1)=x even if associated to the second vertex
  !*!*x(2)=1-x
  !*x_test(1)=0.6_DP !*Closest to the second extremum
  !*x_test(2)=0.4_DP
  !*PRINT*, "e%base", e%base(2,x_test) 
  !*x_test(1)=0.5_DP 
  !*x_test(2)=0.5_DP
  !*PRINT*, "e%base", e%base(1,x_test) !*First one, left vertex
  !*PRINT*, "e%base", e%base(2,x_test) !*Second one, right vertex
  !*PRINT*, "e%base", e%base(3,x_test) !*Internal
  
  !*!*eval_func_element(u,y)
  !*!*eval_func_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the solution in y 
  !*ALLOCATE(u_test(e%nsommets))
  !*u_test=3.0_DP
  !*PRINT*, "u_test", u_test !*NB: u_test cannot be printed
  !*!*P3
  !*u_test(1)=1.0_DP !*Left vertex
  !*u_test(2)=2.0_DP !*Right vertex
  !*u_test(3)=3.0_DP !*First internal DoF
  !*u_test(4)=4.0_DP !*Second internal DoF
  !*PRINT*, "u_test(1)", u_test(1)
  !*PRINT*, "u_test(2)", u_test(2)
  !*PRINT*, "u_test(3)", u_test(3)
  !*PRINT*, "u_test(4)", u_test(4)
  !*!*Second internal DoF (4)
  !*x_test(1)=2._DP/3._DP
  !*x_test(2)=1._DP/3._DP
  !*PRINT*, "In the DoF 4", e%eval_func(u_test, x_test) !*4
  !*!*First DoF (1)
  !*x_test(1)=0._DP
  !*x_test(2)=1._DP
  !*PRINT*, "In the DoF 1", e%eval_func(u_test, x_test) !*1
  !*DEALLOCATE(u_test)

  !*!*eval_der_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the derivative of the solution in y 
  !*!*NB: We specify the barycentric coordinate but the output gradient is in the "real element"
  !*!*We keep into account the passage from the reference element to the real one dividing by the length the gradient in the reference element inside gradient_element (grad=fx/e%volume)
  !*!*P1
  !*ALLOCATE(u_test(e%nsommets))
  !*u_test(1)=1.0_DP !*Left vertex
  !*u_test(2)=2.0_DP !*Right vertex
  !*PRINT*, "u_test(1)", u_test(1)
  !*PRINT*, "u_test(2)", u_test(2)
  !*!*First DoF (1)
  !*x_test(1)=0._DP
  !*x_test(2)=1._DP
  !*PRINT*, "x_test", x_test
  !*PRINT*, "e%eval_der", e%eval_der(u_test,x_test) !*It should be 1/e%volume
  !*PRINT*, "1/e%volume", 1._DP/e%volume
  !*x_test(1)=0.2_DP
  !*x_test(2)=0.8_DP
  !*PRINT*, "e%eval_der", e%eval_der(u_test,x_test) !*It should be 1/e%volume
  !*x_test(1)=82913._DP
  !*x_test(2)=1.8_DP
  !*PRINT*, "e%eval_der", e%eval_der(u_test,x_test) !*It should be 1/e%volume
  !*DEALLOCATE(u_test)
  !*!*P2
  !*ALLOCATE(u_test(e%nsommets))
  !*u_test(1)=1.0_DP !*Left vertex
  !*u_test(2)=1.0_DP !*Right vertex
  !*u_test(3)=13.0_DP !*Internal DoF
  !*x_test(1)=0.5_DP
  !*x_test(2)=0.5_DP
  !*PRINT*, "e%eval_der", e%eval_der(u_test,x_test) !*It should be 0
  !*DEALLOCATE(u_test)

  !*!*eval_der2_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the SECOND derivative of the solution in y 
  !*NB: We specify the barycentric coordinate but the output gradient is in the "real element"
  !*We keep into account the passage from the reference element to the real one dividing by the length SQUARED the gradient in the reference element inside gradient2_element (grad2=fxx/(e%volume**2)) 
  !*!*P1
  !*ALLOCATE(u_test(e%nsommets))
  !*u_test(1)=1.0_DP !*Left vertex
  !*u_test(2)=1.0_DP !*Right vertex
  !*x_test(1)=3021._DP
  !*x_test(2)=321._DP
  !*PRINT*, "e%eval_der2", e%eval_der2(u_test,x_test) !*It should be 0
  !*DEALLOCATE(u_test)
  !*!*P2 
  !*ALLOCATE(u_test(e%nsommets))
  !*u_test(1)=2.0_DP !*Left vertex
  !*u_test(2)=2.0_DP !*Right vertex
  !*u_test(3)=13.0_DP !*Internal DoF
  !*x_test(1)=0.8_DP
  !*x_test(2)=0.2_DP
  !*PRINT*, "e%eval_der", e%eval_der2(u_test,x_test) 
  !*x_test(1)=0.2_DP
  !*x_test(2)=0.8_DP
  !*PRINT*, "e%eval_der", e%eval_der2(u_test,x_test) !*It should be symmetric 
  !*u_test(1)=1.0_DP !*Left vertex
  !*u_test(2)=1.0_DP !*Right vertex
  !*u_test(3)=1.0_DP !*Internal DoF
  !*x_test(1)=0.9_DP
  !*x_test(2)=0.1_DP
  !*PRINT*, "e%eval_der", e%eval_der2(u_test,x_test) !*It should be 0 
  !*DEALLOCATE(u_test)
  
  !*!*gradient_element(k,x) computes the first derivative of the basis function k in the point of baricentric coordinates x
  !*!*NB: The first derivative is already in the "real" element, not in the reference one, because of the final division by e%volume
  
  !*!*gradient2_element(k,x) computes the SECOND derivative of the basis function k in the point of baricentric coordinates x
  !*!*NB: The second derivative is already in the "real" element, not in the reference one, because of the final division by e%volume**2

  !*!*average(u) takes in input the coefficients of the solution (or a general function) in the DoFs and computes the average in the cell
  !*!*u->coefficients in the different nodes 
  
  !*!*quadrature_element fills the weights and the nodes of the quadrature formulas, it is a subroutine in fact
  !*!*Here I trust and do not check, to be checked and tested eventually
  !*CALL e%quadrature()
  !*PRINT*, "e%nquad", e%nquad !*4, order 7 for B3



  !*!*base0 matrix of the values of the basis functions at the physical DoFs
  !*!*base1 inverse of base0
  !*!*base0: 
  !*!*First index DoF
  !*!*Second index basis function
  
  !*!*e%base0=[ phi_1(x_1) phi_2(x_1) ... phi_N(x_1)
  !*!*          phi_1(x_2) phi_2(x_2) ... phi_N(x_2)
  !*!*              .           .
  !*!*              .           .
  !*!*              .           .
  !*!*          phi_1(x_N) phi_2(x_N) ... phi_N(x_n) ]
  !*!*
  !*!*NB: e%base0*vectorofthecoefficient=vectorofthevalues
  !*!*    vectorofthevalues=inv(e%base0)*vectorofthecoefficient=e%base1*vectorofthecoefficient

  !*!*NB: For Pn (any Lagrange) these matrices are identity matrices
  !*PRINT*, "e%base(2,2)", e%base0(2,2), "e%base(2,3)", e%base0(2,3)
  !*PRINT*, "e%base0", e%base0
  !*PRINT*, "e%base1", e%base1

  !*!*e%dof2ind(:) local indices of the nodes with increasing abscissa
  !*!*We have in the natural numeration
  !*!*At first the indices
  !*!*Then the internal DoF with increasing abscissa
  !*!*So e%dof2ind(:) is for example
  !*!*P1-> 1, 2
  !*!*B2-> 1, 3, 2
  !*!*B3,P3-> 1, 3, 4, 2
  !*PRINT*, "e%dof2ind", e%dof2ind

  !*!*eval_coeff_element evaluates some coefficients needed for the DeC
  !*!*eval_coeff_element evaluates some coefficients needed for the DeC
  !*!*e%mass(i,j)-> local mass matrix int_K phi_i phi_j
  !*!*e%coeff(i,j) -> another matrix int_K phi_i  d_x phi_j !*NB:DERIVATIVE ON THE SECOND INDEX

  !*PRINT*, "e%volume", e%volume 
  !*PRINT*, "e%mass", e%mass !*The quadrature formula chosen for P1 determines a lumping on the local mass matrix which is diagonal
  !*PRINT*, "e%coeff", e%coeff
  !*STOP
  !*!*END OF THE TESTS ON THE STRUCTURE ELEMENT
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!*2)arete 
  !*!*A little note: Arete is the class of the boundaries of the elements which are faces in 3d, segments in 2d, points in 1d, so it is clear that in this case (1d) it is almost useless
  !*!*We will have a lot of useless fields which just make sense for higher dimension and will not be used/initialized here !*THEY WILL BE MARKED BY USELESS IF IT IS POSSIBLE TO SAFELY ELIMINATE THEM !*(NB: If they are used I write GOOD)
  !*!*
  !*!*SPOILER: THREE FIELDS ONLY ARE IMPORTANT, THE OTHERS CAN BE DELETED
  !*!*I)bord !*GOOD !*It tells whether the arete is on the boundary of the domain or not
  !*!*NB: 
  !*!*Mesh%edge(1)%bord=T
  !*!*Mesh%edge(Mesh%nsegmt)%bord=Mesh%edge(Mesh%nt+1)%bord=T
  !*!*The others are F
  !*!*II) jt1=-1, jt2 =-1 ! les deux elements de par et d'autre de l'arete. !*GOOD
  !*The two (in 1d) elements sharing the arete 
  !*NB: The aretes at the boundary of the domain will have just one element !*IN THIS CASE IT IS SET jt1=jt2=that element !*Clearly this happens for the first and the last one
  !*Mesh%edge(1)%jt1=Mesh%edge(1)%jt2=1
  !*Mesh%edge(nsegmt)%jt1=Mesh%edge(nsegmt)%jt2=Mesh%nt
  !*For the rest we have
  !*Mesh%edge(indi)%jt1=indi-1
  !*Mesh%edge(indi)%jt2=indi
  !*!*Let's recall the fields of the structures which are in aretes
  !*!*TYPE, PUBLIC:: arete
  !*!*   INTEGER:: nsommets, itype, nvertex! nbre de dofs, type element: 1-> P1 !*USELESS
  !*!*   ! 2-> B2
  !*!*   !3->P2
  !*!*   ! nombre de sommet dans cet element arete
  !*!*   LOGICAL:: bord !*GOOD
  !*!*   INTEGER:: jt1=-1, jt2 =-1 ! les deux elements de par et d'autre de l'arete. !*GOOD
  !*!*   INTEGER, DIMENSION(:,:), POINTER :: nu =>Null() !*USELESS ! nu( indice des elements, indice des voisins): on cherche a connaitre le numero local des 
  !*!*   ! points communs (puisque le maillage est confome)  aux deux elements sur cette face dans chacun des elements jt1 et jt2
  !*!*   REAL(dp), DIMENSION(:,:), POINTER   :: coor=>Null() !*USELESS ! il s'agit des coordonnees physique des dofs (communs) sur la face (commune)
  !*!*   REAL(dp)                           :: volume=0  ! !*USELESS
  !*!*   REAL(dp)                           :: jump_flag=1.0 !*USELESS!(between 0 and 1 which gives the weight of the edge term
  !*!*   ! for this one check if there is a discontinuity around and take the maximum value)
  !*!*   REAL(dp),  DIMENSION(:), POINTER :: n   =>Null() !*USELESS  ! normales exterieures
  !*!*   INTEGER                        :: log    ! logique !*USELESS
!!!!   quadrature de surface !
  !*!*   REAL(dp),   DIMENSION(:,:),POINTER :: quad =>Null() !*USELESS ! point de quadrature 
  !*!*   REAL(dp),     DIMENSION(:),POINTER :: weight=>Null() !*USELESS ! poids 
  !*!*   INTEGER                        :: nquad !*USELESS ! nbre de points de quadrature 
!!! quadrature bord (dimension -1)
  !*!*   REAL(dp),   DIMENSION(:,:),POINTER :: quad_1 =>Null() !*USELESS ! point de quadrature 
  !*!*   REAL(dp),     DIMENSION(:),POINTER :: weight_1=>Null() !*USELESS ! poids 
  !*!*   INTEGER                        :: nquad_1 !*USELESS ! nbre de points de quadrature 
!!!
  !*!* CONTAINS
  !*!*   PROCEDURE, PUBLIC:: aire=>aire_arete !*USELESS
  !*!*   PROCEDURE, PUBLIC:: quadrature=>quadrature_arete !*USELESS
  !*!*   PROCEDURE, PUBLIC:: normale=>normale_arete !*USELESS
  !*!*   !FINAL:: clean_arete !*USELESS
  !*!*END TYPE arete
  !*!*
  !*!*TESTS
  !*artest=Mesh%edge(10)
  !*!*nsommets, itype, nvertex NOT USED ---SUCCEEDED IN DELETING---
  !*!*bord GOOD !*It tells whether the arete is on the boundary of the domain or not
  !*!*NB: 
  !*!*Mesh%edge(1)%bord=T
  !*!*Mesh%edge(Mesh%nsegmt)%bord=Mesh%edge(Mesh%nt+1)%bord=T
  !*!*The others are F
  !*PRINT*, "artest%bord", artest%bord
  !*!*First arete T
  !*PRINT*, "artest%bord first arete", Mesh%edge(1)%bord !*T
  !*!*Last arete T
  !*PRINT*, Mesh%nt+1, Mesh%nsegmt
  !*PRINT*, "artest%bord last arete", Mesh%edge(Mesh%nsegmt)%bord !*T
  !*!*Second arete F
  !*PRINT*, "artest%bord second arete", Mesh%edge(2)%bord !*F
  !*DO indi=1,Mesh%nsegmt
  !*   PRINT*, "artest%bord", indi, Mesh%edge(indi)%bord
  !*END DO
  !*!*jt1,jt2 GOOD !*The two (in 1d) elements sharing the arete 
  !*!*NB: The aretes at the boundary of the domain will have just one element !*IN THIS CASE IT IS SET jt1=jt2=that element !*Clearly this happens for the first and the last one
  !*!*Mesh%edge(1)%jt1=Mesh%edge(1)%jt2=1
  !*!*Mesh%edge(nsegmt)%jt1=Mesh%edge(nsegmt)%jt2=Mesh%nt
  !*!*For the rest we have
  !*!*Mesh%edge(indi)%jt1=indi-1
  !*!*Mesh%edge(indi)%jt2=indi
  !*DO indi=1,Mesh%nsegmt
  !*   PRINT*, "artest%jt1, artest%jt2 of", indi, "->", Mesh%edge(indi)%jt1, Mesh%edge(indi)%jt2
  !*END DO

  !*!*nu NOT USED ---SUCCEEDED IN DELETING---
  
  !*!*coor NOT USED ---SUCCEEDED IN DELETING---
  !*!*NOT EVEN INITIALIZED -> If you try to print you get segmentation fault
  !*!*CLEARLY ALSO THE FUNCTIONS WHICH INVOLVE coor ARE USELESS (moreover they are just for 2d)
  !*!*aire & normale

  !*!*volume and jump_flag USELESS ---SUCCEEDED IN DELETING---
  
  !*!*n USELESS ---SUCCEEDED IN DELETING---
  
  !*!*log USELESS ---SUCCEEDED IN DELETING---
  
  !*!*ALL THE QUADRATURE STRUCTURES USELESS ---SUCCEEDED IN DELETING---
  
  !*STOP
  !*!*END OF THE TESTS FOR THE MOMENT, WE CAN GO AHEAD
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-------------------------------------------------------------------------------------------------------------------------------
  ! ----------Aquire old data if test has been interrupted and does not need to restart from scratch --------
  IF (.NOT.DATA%restart) THEN ! we start with a new solution !*DATA%restart is not taken from the file read in input, it is in the defintion of the DATA structure and it is set to be .FALSE. by default. Probably bacause the simulations in 1d take always a reasonable time.
  !*In 2D it is taken from the input data file and if .TRUE. it restarts from the last result, from scratch
    
     !*Anyway, in this case we start from 0
     kt=0
 
     DO jt=1, Mesh%nt !*Loop over the elements
        e=Mesh%e(jt) !*Quick reference to the processed element through e
        DO l=1, e%nsommets !*Loop on the DoFs of th element
           Var%ua(e%nu(l),0)=IC(e%coor(l),DATA) !*IC is a function in init_bc which sets the initial condition, it takes in input the coordinates of the node and DATA to select through DATA%test the correct initial condition !*NB: The initialization is performed in all the nodes but only in the first sutimestep (the one in 0) !*RMK: subtimesteps 0,1,...,M where M+1=subtimesteps=P=iterations=oreder
           !*NB: The initial condition is given consists in the values (of the conserved variables)
           
           Var%un(e%nu(l))=Var%ua(e%nu(l),0) !*He stores the initial condition also in Var%un (he will work on this variable to make the transformation of the values of the initial condition into coefficients BUT NB: Apart from this initial use Var%un will be later used for the DeC to store the L^2 operator in every node needed for the updating var%ua(is,k-1)=Var%up(is,k-1)-Var%un(is)*Mesh%aires(is))

        ENDDO
        !*Now he works on Var%un (which at the moment contains the values) to get the coefficients
        SELECT CASE(e%itype) 
        CASE(1,3,4,7,11,12,13,14) !Lagrange-> nothing to do
           !*Lagrange: coefficients=values so do not modify Var%un
        CASE(2,5,6) ! B3 ! cubic Bezier: modif all dof except vertices
           !*Bezier, transformation needed
           !*NOT USED!*ALLOCATE(u_b(e%nsommets), u_c(e%nsommets))
           DO k=1, n_vars !*Loop on the components (clearly the transformation is needed for every component)
              DO l=1, e%nsommets !*Loop on the DoFs
                 Var%un(e%nu(l))%u(k)=SUM(e%base1(l,:)*Var%ua(e%nu(:),0)%u(k)) !*Passage from values to coefficients, look down here for the explaination
                 !*NB: Let's recall base1 from elements 
                 !*e%base0 matrix of the values of the basis functions at the physical DoFs
                 !*e%base1 inverse of e%base0
                 !*base0: 
                 !*First index DoF
                 !*Second index basis function
                 !*e%base0=[ phi_1(x_1) phi_2(x_1) ... phi_N(x_1)
                 !*          phi_1(x_2) phi_2(x_2) ... phi_N(x_2)
                 !*              .           .
                 !*              .           .
                 !*              .           .
                 !*          phi_1(x_N) phi_2(x_N) ... phi_N(x_n) ]
                 !*
                 !*NB: e%base0*vectorofthecoefficients=vectorofthevalues
                 !*    vectorofthecoefficients=inv(e%base0)*vectorofthevalues=e%base1*vectorofthevalues
                 !*With the expression SUM(e%base1(l,:)*Var%ua(e%nu(:),0)%u(k)) he's making the product of the row l (referred to DoF l) of the inverse by the vector of the values thus the coefficient of the DoF l
              ENDDO
           ENDDO
           !*NOT USED!*DEALLOCATE(u_b,u_c) NOT USED
           !*NOW in Var%un there are the coefficients either equal to the values (for Lagrange) or different (Bezier) 
        CASE default
           PRINT*, "correct initialisation not yet implemented for these elements,e%itype=",e%itype
           STOP
        END SELECT

     ENDDO

     Var%ua(:,0)=Var%un(:) !*Set the coefficients equal to the values !*NB The initial condition is saved in the first subtimestep
     !*!*RMK: Var is type variables
     !*!*TYPE variables
     !*!*   REAL(dp):: dt
     !*!*   !*INTEGER:: Ncells !*Var%Ncells=Mesh%ns !*PROBABLY NOT NEEDED !*Var%Ncells=Mesh%ns is always 4096, probably an old structure not used anymore
     !*!*   TYPE(PVar), DIMENSION(:,:),ALLOCATABLE:: ua, up
     !*!*   !*Both are vectors of PVar (look in variable_def) for the definition. Essentially in Pvar we store the conservative (or primitive) variables (values or coefficients) in a single node (of a single subtimestep)
     !*!*   !*In both ua and up we have 2 indices that are allocated in the following way (in the main) 
     !*!*   !*(Mesh%ndofs,0:DATA%iordret-1)
     !*!*   !*1st->DoFs
     !*!*   !*2nd->subtimesteps RMK: The subtimesteps are M+1=order i.e. 0,1,...,order-1
     !*!*   !*up->PREVIOUS VALUES
     !*!*   !*ua->ACTUAL VALUES (we want to pass from up to ua by adding L^2 divided by the measure of the dual cell (at each iteration))
     !*!*   TYPE(Pvar), DIMENSION(:),ALLOCATABLE:: un
     !*!*   !*In un is stored the operator L^2 for the single subtimestep 
     !*!*   !*un has a single index referred to the DoFs (Mesh%Ndofs)
     !*!*   
     !*!*   !*SO THAT in the end we impose
     !*!*   !*var%ua(is,k)=Var%up(is,k)-Var%un(is)*Mesh%aires(is)
     !*!*   !*RMK: In Mesh%aires we have the inverse of the measures of the dual cells
     !*!*END TYPE variables

     !*Kind of safety check, after the initialization of the IC he prints the maximum and the minimum values of the coefficients for each component
     DO l=1, n_vars 
        PRINT*, "min variable ",l,  MINVAL(Var%ua(:,0)%u(l))
        PRINT*, "max variable ",l,  MAXVAL(Var%ua(:,0)%u(l))
        REWIND(10)

        PRINT*
     ENDDO

  ELSE !*If DATA%restart is .TRUE. (RMK: it is always set to be .FALSE.) restart from scratch
     !*NB: It is not possible in this case because it's not just a matter to reinitialize the counter (which is actually useless also here) but we need also sth (which is present in 2d and absent here) to read the last result (in a certain folder eventually)
     kt=kt0 !*In reality kt should be read from the last result. (.sol in 2d)
     !*NB: NOT ONLY writesol is not present but also readsol
  ENDIF
  !-------------------------------------------------------------------------------------------------------------------------------

  !*Here it repeates the same safety check (the one that we had if we started from 0) printing the maximum and the minimum values of the coefficients for each component
  !*NB: If we started from scratch we see it for the first time  
  PRINT*, "Initialisation"
  DO l=1, n_vars
     PRINT*, "min variable ",l,  MINVAL(Var%ua(:,0)%u(l))
     PRINT*, "max variable ",l,  MAXVAL(Var%ua(:,0)%u(l))
     PRINT*
  ENDDO


  CALL visu(DATA, kt, mesh, Var,folder) !*In postprocessing, he stores in a file the result at this point, at the time 0 if we started from 0 or it resaves the last result if we started from the beginning

  !*-------------------------------------------------------------------------------------------
  !*HERE I CARRY MY PRECOMPUTATION
  Mesh%ns=Mesh%nDoFs
  CALL precomputation(Mesh, VecProx)
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*SAFETY CHECK
  !*PRINT*, SIZE(VecProx, DIM=1) !*OK
  !*STOP
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !*Last final Check
  !*indi=17
  !*PRINT*, "Node ", indi
  !*PRINT*, "It has ", VecProx(indi)%numbofDoFs, "close DoFs " 
  !*PRINT*, VecProx(indi)%vecDoFs(:)
  !*PRINT*, "It is contained by ", VecProx(indi)%numbofEl, "elements"
  !*PRINT*, VecProx(indi)%vecEl(:)
  !*DO indt=1,VecProx(indi)%numbofEl
  !*   PRINT*, "Element: ", VecProx(indi)%vecEl(indt)
  !*   PRINT*, Mesh%e(VecProx(indi)%vecEl(indt))%nu(:)
  !*END DO
  !*STOP
  !*REMARK: i is close to itself
  !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !*-------------------------------------------------------------------------------------------

  !*Initialization of the source term 
  !*Bathymetry and Grad of the Bathymetry for SW
  !*Potential and Grad of the Potential for Eul
  CALL InitializeSourceTerm(DATA%test,Mesh,VecProx)


  !-------------------------------------------------------------------------------------------------------------------------------
  ! -------------------- Time- stepping loop -------------------------

  DO kt=kt0+1, DATA%ktmax  ! loop over time steps 

     !*!*RMK: The initial condition is now in Var%ua(:,0)
     !*!*Var%ua has 2 indices
     !*!*first index referred to the DoFs 1:Mesh%NDoFs
     !*!*second index referred to the subtimesteps 0:order-1 i.e. 0:Data%iordret-1

     ! switch

     !*!*Copying the initial condition in other variables
     !*!*RMK on the structures
     !*!*Var is a variables type and contains
     !*!*Var%ua(Mesh%ndofs,0:DATA%iordret-1) vector of PVar
     !*!*Var%up(Mesh%ndofs,0:DATA%iordret-1) vector of PVar
     !*!*Var%un(Mesh%ndofs) vector of PVar
     !*!*debug is also a variables type
     !*!*debug%ua(Mesh%ndofs,0:DATA%iordret-1)
     !*!*debug%up(Mesh%ndofs,0:DATA%iordret-1)
     !*!*debug%un(Mesh%ndofs)
     !*!*temp are instead vectors of PVar, they are consistent with Var%ua and Var%up
     !*!*temp_var(Mesh%ndofs,0:DATA%iordret-1)
     !*!*temp(Mesh%ndofs,0:DATA%iordret-1)
     !*!*temp_debug(Mesh%ndofs,0:DATA%iordret-1)

     DO is=1, Mesh%ndofs !*Loop over the DoFs
        temp(is,:)=Var%ua(is,0) !*All the subtimesteps (second index) of temp set equal to the initial condition OR solution detected in the previous time stepping loop thanks to Var%ua(:,0)=Var%ua(:,DATA%iordret-1) at the end of the timestepping loop
        Var%ua(is,1:)=Var%ua(is,0) !*All the subtimesteps (second index) of Var%ua equal to initial condition OR solution detected in the previous time stepping loop thanks to Var%ua(:,0)=Var%ua(:,DATA%iordret-1) at the end of the timestepping loop
        Var%un(is)=0._dp !*Var%un initialized to 0. We rmk that it will contain the residuals
        debug%ua(is,:)=Var%ua(is,0) !*All the subtimesteps (second index) of debug%ua set equal to the initial condition OR solution detected in the previous time stepping loop thanks to Var%ua(:,0)=Var%ua(:,DATA%iordret-1) at the end of the timestepping loop
        temp_debug(is,:)=Debug%ua(is,:) !*All the subtimesteps (second index) of temp_debug set equal to the initial condition OR solution detected in the previous time stepping loop thanks to Var%ua(:,0)=Var%ua(:,DATA%iordret-1) at the end of the timestepping loop
        temp_var(is,:)=Var%ua(is,:) !*All the subtimesteps (second index) of temp_var set equal to the initial condition OR solution detected in the previous time stepping loop thanks to Var%ua(:,0)=Var%ua(:,DATA%iordret-1) at the end of the timestepping loop
     ENDDO

     dt0=  Pas_de_temps(Mesh,Var, DATA) !*Pas_de_temps is a function down here in the same program. In practice it determines the maximum dt allowed and multiplies it by the CFL
     dt=MIN(dt0,DATA%tmax-DATA%temps) !*Then we make sure not to overcome the final time by taking the minimum between the computed dt and the dt needed to arrive to the final time
     Var%dt=dt !*dt is stored in Var
     IF (dt.LE.1.e-13_dp) EXIT !*If dt is too small it stops to prevent the program to take ininite time to arrive at the final time
     tn=DATA%temps !*Actual time stored in tn
     !*DATA%temps is 0 at the beginning then it is updated in this loop
     !*So if it is the first time that we enter this cycle it is 0 (or the previous final value of the last simulation), if we here from the second iteration on we store in tn the updated final time of the previous iteration (the actual time)


     ! Euler time stepping
     ! initilalise the time step with the right scheme


     !*!!!!!!!!!!!!!!!!!!!!!!
     !*I WILL NOT DEBUG/COMMENT IN DEPTH MOOD SO GO FROM HERE--->
     !*!!!!!!!!!!!!!!!!!!!!!!
     ! loop over the cascade of schemes
     !*===========================================================================
     !*THE WHOLE MOOD PROCEDURE IS CONFINED HERE 
     !*It is a stabilization technique that I will skip
     !*I just report the following ALERT (in test which is a subroutine in timestepping)
     !*->test is basically all in SP, the real are declared as REAL and not REAL(DP)
     IF (DATA%mood) then 
         Mesh%e(:)%type_flux=DATA%ischema
         Mesh%e(:)%diag2=0

         DO nflux = 1,size(fluxes_mood)


            if (nflux.gt.1) then
               DO is=1, Mesh%ndofs
                  Debug%ua(is,:)= temp_debug(is,:)
                  Var%ua(is,:)=temp(is,:)
               END DO
            endif

            DO k_inter=1,DATA%iter
               DO is=1, Mesh%ndofs
                  debug%up(is,:)=debug%ua(is,:)
               ENDDO
               DO k=1, DATA%iordret-1 !*<-CHANGED
                  CALL debut()
               ENDDO ! k
            END DO

            ! do the test at the end of timestep, not at every iteration!!!!
            ! maybe positivity checks inside every iteration

            CALL test(Data%iordret-1,Debug,Var, mesh, DATA,fluxes_mood(nflux)) !*NB: test is in timestepping

       
         END DO ! nflux

    !!$     
         DO is=1, Mesh%ndofs
            Var%ua(is,:)=temp(is,:)
         END DO
     END IF
     !*===========================================================================

     !*!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*---->TO HERE. Basically I skipped the last if
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!

     !*Here we have the main of the main program. The DeC procedure
     !*NB: Up to now we have 
     !*->Var%ua(is,1:)=Var%ua(is,0) !*All the subtimesteps (second index) of Var%ua equal to initial condition OR solution detected in the previous time stepping loop thanks to Var%ua(:,0)=Var%ua(:,DATA%iordret-1) at the end of the timestepping loop
     !*In any case Var%ua(1:Mesh%NDoFs,:) is "constant" in the second component and we have there the starting value
     DO k_inter=1,DATA%iter !*Loop on the iterations !*RMK the iterations should be iterations=P>=M+1=subtimesteps=accuracy. In the optimal case P=M+1=accuracy
        !*SO we must have at least DATA%iter=DATA%iorder (they could be more but it would be useless)
        

        DO is=1, Mesh%ndofs !*Loop on the DoFs
           Var%up(is,:)=Var%ua(is,:) !*He's copying the initial value from which performing the iteration in Var%up (all the subtimesteps)
           !*!*In practice it is a (multi)vector "constant" in the second component with the starting value in each component at the first iteration, in the other iterations i 
        ENDDO
        !*!*At this point we have
        !*!*Var%up=Var%ua=u^(p-1) associated to the iteration p-1
        !*!*u^(p-1) is made of all the subtimesteps u^m,(p-1) m=0,1,...,M=order-1 and each of the subtimesteps contains all the DoFs
        !*!*u^(p-1)=Var%up(DoFs,m)=Var%ua(DoFs,m)
        !*!*First index->DoFs
        !*!*Second index->Subtimesteps

        !*!*OUR AIM IS TO COMPUTE u^(p) thorugh the DeC

        !*!*The updating formula of the DeC reads
        !*!*L^1(u^(p))=L^1(u^(p-1))-L^2(u^(p-1))
        !*!*This formula involves all the subtimesteps (and obviously all the nodes)
        !*!*We have in particular we have for each subtimestep m (and DoF i)
        !*!*u_i^m,(p)=u_i^m,(p-1)-1/|C_i|*{\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))+\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}]} (#)
        !*!*Formula (#) is performed in fin() inside the loop on the subtimesteps (inside the iteration loop of the DeC inside the timestepping loop)

        DO k=1, DATA%iordret-1 !*Loop on the subtimesteps
             !*NB: The subtimesteps are in total M+1=accuracy=DATA%iordret
             !*->subtimesteps 0,1,...,M=DATA%iordret-1
             !*BUT according to the DeC the solution in the first subtimestep (0) must not be updated
             !*SO this loop must cover all the subtimestep apart from 0
             !*SO the index goes from 1 to M=DATA%iordret-1
             

             DATA%temps = tn + alpha(k,DATA%iordret)*dt !*Set the time to be the one of the treated subtimestep

             CALL fin()!*<-THIS IS THE REAL TIMESTEPPING, THE DEC, AND IT IS A SUBROUTINE DOWN HERE
        ENDDO ! k

        !*We have updated Var%ua, we can go to the next iteration (or finish the iterations if we did all of them)

     ENDDO ! n_inter !*End loop on the iterations

     !*When we finish all the DeC iterations we are here

     DATA%temps= tn +dt !*We update the time
     Var%ua(:,0)=Var%ua(:,DATA%iordret-1) !*We update the value of the solution copying the last timestep of the last iteration in the first timestep.
     !*NB:OUR REFERENCE IS ALWAYS Var%ua(:,0)
     !*At the beginning of each loop of the timestepping we start with Var%ua(:,0) and we copy it in the other subtimesteps to get u^(0) initial value from which carrying the first iteration

     IF (MOD(kt,DATA%ifre)==0) THEN !*Every DATA%ifre we store the solution 
        PRINT*
        PRINT*, kt, DATA%temps, DATA%tmax, dt !*Print these general information
        DO l=1, n_vars !*Safety check to understand if sth is exploding
           PRINT*, "min ua, variable ",l, MINVAL(Var%ua(:,0)%u(l))
           PRINT*, "max ua, variable ",l, MAXVAL(Var%ua(:,0)%u(l))
           PRINT*
        ENDDO

        !100     FORMAT(1x,i5,3(1x,e10.4))
        CALL visu(DATA, kt, mesh, Var,folder) !*Storing of the actual solution

        PRINT*
     ENDIF
  ENDDO ! end  time stepping

  !*We reached now the final time

  PRINT*, "final time=", DATA%temps, DATA%tmax  !*Print the information

#if(1==1)
  PRINT*
  PRINT*, kt, DATA%temps, DATA%tmax, dt !*Print these general information
  DO l=1, n_vars !*Safety check to understand if sth is exploding
     PRINT*, "min ua, variable ",l, MINVAL(Var%ua(:,0)%u(l))
     PRINT*, "max ua, variable ",l, MAXVAL(Var%ua(:,0)%u(l))
     PRINT*
  ENDDO
#endif


  ALLOCATE(error1(n_vars),error2(n_vars),errorinf(n_vars))
  error1=0._DP
  error2=0._DP
  errorinf=0._DP
  CALL errorTestProjection(error1,error2,errorinf,Var%ua(:,0),Mesh,DATA)





  !post processing
  CALL visu(DATA, kt, mesh, Var,folder) !*Save the final solution
 
  !*Deallocate the used structures
  DEALLOCATE(Var%ua, Var%up , var%un, Mesh%e, Mesh%edge, Mesh%aires )
  DEALLOCATE( temp_var, temp_debug,temp, fluxes_mood)

CONTAINS

  !*-----------------------------------------------------
  !*Pas_de_temps determines the maximum dt allowed and multiplies it by the CFL
  !*NB: We take as an input the whole Var but the actual interesting reference value is in Var%ua(:,0) (it is also true that it has been copied in all the other subtimesteps)
  !*-----------------------------------------------------
  REAL(dp) FUNCTION pas_de_temps(Mesh,Var,DATA)
    IMPLICIT NONE
    TYPE(variables), INTENT(in):: Var 
    TYPE(Maillage), INTENT(in):: Mesh
    TYPE(donnees), INTENT(in):: DATA
    INTEGER:: jt, l
    TYPE(element):: e
    TYPE(Pvar), ALLOCATABLE, DIMENSION(:):: u
    REAL(dp):: dt, dx, umax
    REAL(dp), DIMENSION(n_dim):: x=0.0_dp, n

    dt=HUGE(1.00_dp) !*dt set to be the maximum possible for a real of DP precision. It will be then reduced to the maximum allowed 
    DO jt=1, Mesh%nt !*Loop on the elements
       e=Mesh%e(jt)
       ALLOCATE(u(e%nsommets))
        DO l=1, e%nsommets !*Loop on the DoFs
          dx=e%volume
          n=1.0_DP !*ALERT: IT WAS DECLARED WITHOUT _DP
          x=e%coor(l)
          u(l) = e%eval_func(Var%ua(e%nu,0),e%x(:,l)) !*It evaluates the solution in the processed DoF
          umax=u(l)%spectral_radius(x,n) !*It considers the spectral radius in that DoF of the (normal) Jacobian (actually in 1d "normal" is useless)
          dt=MIN(dt,e%volume/( ABS(umax)+1.e-6_dp)) !*Computes the maximum allowed dt by taking the minimum among all the DoFs of dx/|f'(u)| !*+1.e-6_dp is added
       ENDDO
       DEALLOCATE(u)
    ENDDO
    !*Then, once we have dt, we impose the CFL
    pas_de_temps=dt*DATA%cfl 
  END FUNCTION pas_de_temps


SUBROUTINE debut ()

    DO is=1, Mesh%Ndofs
       debug%un(is)=0.0_dp
    ENDDO

    !!    CALL shock_indicator(Mesh,Var%ua(k-1,:),ShIndArray)
    DO jt=1, Mesh%nt
       CALL main_update(k,dt,Mesh%e(jt),debug,DATA,alpha,beta,gamma,n_theta,theta,jt,Mesh)
    ENDDO
   !--------------------------------------------------------------------------------------
    ! EDGE UPDATE DONE ONLY IN CASE WE DO NOT HAVE LXF scheme
    IF (DATA%ischema==5 .OR.DATA%ischema==4) THEN !(with jumps--> Burman'stuff)
        CALL edge_main_update(k,DATA,Mesh,debug,dt,alpha,beta,gamma,n_theta,theta)
    END IF  
    CALL BC_Residuals(debug, Mesh%e(1)%nu(1), Mesh%e(Mesh%nt)%nu(2), DATA, k, n_theta,theta,alpha,dt,tn)

    ! Compute the 'trial' solution to make MOOD tests
    DO is=1, Mesh%ndofs
       debug%ua(is,k)=debug%up(is,k)-debug%un(is)*Mesh%aires(is) !*<-CHANGED
    ENDDO

    CALL Strong_BC(debug, Mesh%e(1)%nu(1), Mesh%e(Mesh%nt)%nu(2), DATA, k, n_theta,theta,alpha,dt,tn)
    

  END SUBROUTINE debut

  !*------------------------------------------------------------------------------------
  !*fin() performs the updating of the single subtimestep inside the loop on the subtimesteps (inside the iteration loop of the DeC inside the timestepping loop)
  !*
  !*Denoting by m the index of the subtimestep (and by p the index of the iteration), the updating formula performed in fin() reads
  !*u_i^m,(p)=u_i^m,(p-1)-1/|C_i|*{\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))+\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}]} (#)
  !*NB: In reality since we are in 1d we do not have the boundary residual
  !*\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}
  !*They are not needed here, there are escamotages to avoid them. So we just compute
  !*So we have:
  !*u_i^m,(p)=u_i^m,(p-1)-1/|C_i|*{\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))]}
  !*------------------------------------------------------------------------------------
  SUBROUTINE fin()
   !*!*Var%un will contain the updating part i.e.
   !*!*\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))+\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}]
   !*!*This must be multiplied by 1/|C_i| and subtracted to u_i^m,(p-1) to get the update
   !*!*We'll have in particular
   !*!*Var%ua(i,k-1)=Var%up(i,k-1)-Var%un(i)*Mesh%aires(i)
   !*!*RMK:
   !*!*i is the DoF index which goes through all the sutimesteps (we have a loop inside this sub)
   !*!*k-1 is the index of the subtimestep NB: the loop on the subtimesteps is outside. We do not have a loop on the subtimesteps inside fin()
   !*!*Mesh%aires(i) has a very misleading name since it contains the inverse of the areas (length because 1d) of the dual cells

    !*!*Loop to initialize Var%un which will contain the updating part
    DO is=1, Mesh%Ndofs
       var%un(is)=0.0_dp
    ENDDO

    !*!*Loop on the elements to compute the updating part
    !*!*NB: The updating part lacks the jump stabilization over the edges which is computer later if the scheme requires it (the schemes requing it are Galerkin+jump=4 and Psi+jump=5)
    !*!*Here we add to Var%un
    !*!*\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + 
    !*!*+ dt*\sum_{l=0}^M theta_l^m*\Sum_{K \in K_i} RES_i^K(u^l,(p-1))
    !*!*So just the mass matrix part 
    !*!*And the node-element (not boundary) residuals up to the edge stabilization
    DO jt=1, Mesh%nt 
       CALL main_update(k,dt,Mesh%e(jt),Var,DATA,alpha,beta,gamma,n_theta,theta,jt,Mesh) !*main_update in timestepping
    ENDDO

    !--------------------------------------------------------------------------------------
    ! EDGE UPDATE DONE ONLY IN CASE WE DO NOT HAVE LXF scheme
    !*!*
    IF (DATA%ischema==5 .OR.DATA%ischema==4) THEN !(with jumps--> Burman'stuff)
       CALL edge_main_update(k,DATA,Mesh,Var,dt,alpha,beta,gamma,n_theta,theta	) !*edge_main_update in timestepping
    ENDIF
    !--------------------------------------------------------------------------------------

    !*To keep into account the BCs
    !*BC in init_BC
    CALL BC_Residuals(Var, Mesh%e(1)%nu(1), Mesh%e(Mesh%nt)%nu(2), DATA, k, n_theta,theta,alpha,dt,tn) 
    !*IMPORTANT RMK: Also in the 1D case the boundary residuals of the different subtimesteps MUST be combined through the thetas.
    !*
    !*NB: the time DATA%time=tn+alpha(k,DATA%iordret)*dt so it is better to use tn and alpha to refer to the times of the subtimesteps if needed (or at least don't forget this)
    !*
    !*In practice we can often "avoid" to perform a combination inside BC_Residuals because we have: 
    !*-> PERIODIC BC: Var%un(1) and Var%un(Mesh%nt) (node residuals of the boundary nodes (first and last) and so ALREADY combined) are given to the first and the last node. We just need to keep into account the fact that Var%un(Mesh%nt) goes also to the first node and Var%un(1) goes also to the last node
    !*-> STRONG IMPOSITION: nothing to combine at the boundary. We just impose the values of Var%ua without caring for the updating.
    !*-> OUTFLOW: no contribution of the boundary residuals to the residuals at the nodes
    !*
    !*In 2d it is more critical: the weak boundary conditions generate boundary residuals that must be combined through the thetas
    !--------- Update of the solution -----------------------------------------------
    DO is=1, Mesh%ndofs
       var%ua(is,k)=Var%up(is,k)-Var%un(is)*Mesh%aires(is) !*Updating
    ENDDO

    CALL Strong_BC(Var, Mesh%e(1)%nu(1), Mesh%e(Mesh%nt)%nu(2), DATA, k, n_theta,theta,alpha,dt,tn)
    !*REPEATED FOR STRONG BCs BADLY IMPLEMENTED, essentially they update Var%ua which is overwritten once out of BC

   !*!*ALERT
   !*!*I report a problem on the strong imposition of the boundary conditions in 1d.Basically BC is called before updating the solution with in fin() 
   !*!*DO is=1, Mesh%ndofs
   !*!* var%ua(is,k-1)=Var%up(is,k-1)-Var%un(is)*Mesh%aires(is)
   !*!*ENDDO
   !*!*
   !*!*Whenever we want to impose strong Dirichlet (see for example case(4) in init_bc_euler), BC modifies Var%ua through
   !*!*Var%ua(k-1,i1)%u=whatever we want
   !*!*But then when we update we lose this information (because we set  var%ua(is,k-1)=Var%up(is,k-1)-Var%un(is)*Mesh%aires(is))
   !*!*I've talked a bit with Davide about it:One possibility is to call BC before and after the updating. This works but I think we lose a bit in "linearity" and "clarity" of the code. So I suggest to think a bit about it and to find a solution.
   !*!*NB:->Setting Var%un=0 doesn't work because later we would still have
   !*!*var%ua(is,k-1)=Var%up(is,k-1)-Var%un(is)*Mesh%aires(is)=Var%up(is,k-1) 
   !*!*->so you may think that we could modify Var%up and set Var%un=0This doesn't work either because at the next subtimestep you would have a different value for Var%up
   !*!*For the moment I suggest to apply BCs before and after the updating but as I wrote, this is not a proper solution if we want to fix the code forever.

   !*!* FIXED WITH Strong_BC (Which anyway recalls BC_Residuals, so it should be done better isolating the strong imposition anyway it is not bad for the moment) :)


  END SUBROUTINE fin

  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------

  SUBROUTINE theta_alpha_beta()

    ! loading the weights

    ! alpha: this defines the fractions of Delta t for the Euler method (operator L1)
    ! alpha( subtimestep, order in time)
    ! theta: \int_0^alpha(L,k) f = \sum_{p=0}^{order in time) theta(p,L, order in time)* f(u[p])
    ! beta contains the (time) weight for the Euler step of Version 1


    !*Coefficients for the Dec
    !*P = iterations = M+1 = subtimesteps (counting 0) = basis functions=accuracy (order)
    !*So subtimesteps 0,1,...,M=accuracy-1
    !* In the DeC we have (for equispaced nodes)
    !* L^2 
    !*    theta_p^l=\int_0^{l/M} phi_p(t)dt 
    !*       where p=0,1,...,M=accuracy-1 basis functions & l=1,2,...,M=accuracy-1 subtimesteps
    !* L^1
    !*    beta^l=l/M !*NB: Never used in the last version of L^1, they cancel
    !*
    !*In this code they are stored in the following way (before let's present the notation)
    !*NOTATION:
    !*k->order
    !*l->subtimestep->1,2,...,k-1 (NB: The first subtimestep at 0 is in practice never used)
    !*p->basis function->0,1,...,k-1 (NB: The one in 0 is used obviously)

    !*NB: ONLY THE THETA COEFFICIENTS WILL BE USED

    !*The coefficients alpha are what we usually refer as beta in the DeC ,in L^1
    !*The structure alpha(:,:) has 2 indices
    !*alpha(l,k)=l/(k-1) subtimestep
    !*First index=l->subtimestep l=1,2,...,k-1
    !*Second index=k->order
    !*alpha(0,k) not present, if present it would be 0

    !*Beta are coefficients that we would find in the "ugly" L^1 i.e. the interval between two consecutive subtimesteps, you use it if you make an Euler for each subtimesteps (but NB: we do not use this L^1, we use the good one)
    !*beta(l,k-1,k)
    !*l=1,2,...,k-1
    !*beta(1,k-1,k)=alpha(1,k)
    !*beta(2,k-1,k)=alpha(2,k)-alpha(1,k) ecc...   !*beta(l,k-1,k)=alpha(l,k)-alpha(l-1,k)
    !*They are not even initialized for order 2 (and 1 which is down here even if the DeC doesn't make sense for order 1)
    !*I have no clue about what the gammas are
    !*They are defined as the beta coefficients and basically he assigns
    !*gamma(2:,k-1,k)=beta(2:,k-1,k)
    !*So he copies the betas from the second value of the first index on
  
    !*n_theta=Number of basis functions for each order, basically the order because the basis functions are M+1=accuracy=k so 0,1,...,M
    !*n_theta(k)=k
  
    !*theta coefficients for L^2 in the DeC
    !*theta(p,l,k)=int_0^alpha(l,k) phi_p(s)ds
    !*p->0,1,...,k-1 basis functions
    !*l->1,2,...,k-1 subtimesteps
    !*k order



    alpha=0._dp
    theta=0._dp
    beta=0._dp
    gamma=0._dp
    n_theta=0 !*<- It is a vector of integers, no need to declare it as DP

    ! for order 1 !*The DeC doesn't make sense for order 1 pathlogic
    n_theta(1)=1
    theta(:,:,1)=0._dp
    theta(0,1,1)=1._dp
!!$
    ! for order 2 !*Betas not even defined/initialized in this case
    alpha(1,2)=1._dp !*beta in the normal L^1 of the DeC
    n_theta(2)=2 !*number of basis functions
    theta(:,:,2)=0._dp 
    theta(0:1,1,2)=0.5_dp !*Integral of the basis functions

    ! quadratic interpolation over [0,1]
    ! alpha, beta gamma coeffcients for order 3
    !
    alpha(1,3)=0.5_dp !*beta in the normal L^1 of the DeC
    alpha(2,3)=1._dp !*beta in the normal L^1 of the DeC
    beta(1,2,3)=alpha(1,3) !*beta in the "ugly" L^1 of the DeC
    beta(2,2,3)=alpha(2,3)-alpha(1,3) !*beta in the "ugly" L^1 of the DeC i.e. distance between consecutive subtimesteps
    GAMMA(2,2,3)=beta(2,2,3) !*Mystery XD, anyway only the theta are really useful
    n_theta(3)=3 !*Number of basis functions
    !
    ! integration coefficients for order 3
    theta(0:2,1:2,3)=RESHAPE( (/5.0_dp/24._dp,1._dp/3._dp,-1._dp/24._dp&
         &,1._dp/6._dp,2._dp/3._dp,1._dp/6._dp/),(/3,2/)) !*Theta coefficients

!!$  ! cubic interpolation
!!$#if (-2==0)
!!$  ! Here I choose the points 0.5*(1+cos(k*pi/N)), N=3, k=3,2,1,0
!!$!!!!!!!!!!!!!!
!!$  ! these are
!!$
!!$  alpha(1,4)=1./4.
!!$  alpha(2,4)=3./4.
!!$  alpha(3,4)=1.
!!$  n_theta(4)=4
!!$
!!$  theta(0:3,1:3,4)=RESHAPE( (/59./576., & !(0,1)
!!$       & 47./288.,  &                      !(1,1)
!!$       & -7./288., &                      !(2,1)
!!$       & 5./576.,  &                      !(3,1)
!!$       & 3./64.,   &                      !(0,2)
!!$       & 15./32.,  &                      !(1,2)
!!$       & 9./32.,   &                      !(2,2)
!!$       & -3./64.,  &                      !(3,2)
!!$       &1./18.,    &                      !(0,3)
!!$       &4./9.,     &                      !(1,3)
!!$       &4./9.,     &                      !(2,3)
!!$       &1./18./)   &                      !(3,3)
!!$       &,(/4,3/))
!!$#endif
#if (-2==-2) 
    ! Here I choose the points k/3, k=3,2,1,0
!!!!!!!!!!!!!!
    ! these are

    !*Beta in the normal L^1 of the DeC, location of the subtimesteps
    alpha(1,4)=1._dp/3._dp 
    alpha(2,4)=2._dp/3._dp
    alpha(3,4)=1._dp
    n_theta(4)=4

    !*Thetas
    theta(0:3,1:3,4)=RESHAPE( (/1._dp/8._dp, & !(0,1) 
         & 19._dp/72._dp,  &                      !(1,1)
         & -5._dp/72._dp, &                      !(2,1)
         & 1._dp/72._dp,  &                      !(3,1)
         & 1._dp/9._dp,   &                      !(0,2)
         & 4._dp/9._dp,  &                      !(1,2)
         & 1._dp/9._dp,   &                      !(2,2)
         & 0._dp,      &                      !(3,2)
         &1._dp/8._dp,    &                      !(0,3)
         &3._dp/8._dp,     &                      !(1,3)
         &3._dp/8._dp,     &                      !(2,3)
         &1._dp/8._dp/)   &                      !(3,3)
         &,(/4,3/))
#endif
    !! The following lines are good only for the version 1
    beta(1,3,4)=alpha(1,4) ! this gives the fraction of timestep for the Euler method, in the ugly L^1 not used
    beta(2,3,4)=alpha(2,4)-alpha(1,4)
    beta(3,3,4)=alpha(3,4)-alpha(2,4)
    GAMMA(2:3,3,4)=beta(2:3,3,4) !*Again mistery

!!$-------------------------------------------------------------------
!!$-------------------------------------------------------------------
!!$-------------------------------------------------------------------

!!$ 4th degree interpolation, 5th order

    !*Waste of time to comment the following, it's just like the previous ones
#if (-1==-1)
    !EQUISPACED POINTS
    alpha(1,5)=1._dp/4._dp
    alpha(2,5)=2._dp/4._dp
    alpha(3,5)=3._dp/4._dp
    alpha(4,5)=1._dp
    n_theta(5)=5

    theta(0:4,1:4,5)=RESHAPE( (/ 251._dp/2880._dp, & !(0,1)
         & 323._dp/1440._dp,  &                      !(1,1)
         & -11._dp/120._dp, &                      !(2,1)
         & 53._dp/1440._dp,  &                      !(3,1)
         & -19._dp/2880._dp,  &                      !(4,1)
         & 29._dp/360._dp,   &                      !(0,2)
         & 31._dp/90._dp,  &                      !(1,2)
         & 1._dp/15._dp,   &                      !(2,2)
         & 1._dp/90._dp,      &                      !(3,2)
         & -1._dp/360._dp,      &                    !(4,2)
         &27._dp/320._dp,    &                      !(0,3)
         &51._dp/160._dp,     &                      !(1,3)
         &9._dp/40._dp,     &                      !(2,3)
         &21._dp/160._dp,   &                      !(3,3)
         &-3._dp/320._dp,   &                      !(4,3)
         &7._dp/90._dp,    &                      !(0,4)
         &16._dp/45._dp,     &                      !(1,4)
         &2._dp/15._dp,     &                      !(2,4)
         &16._dp/45._dp,   &                      !(3,4)
         &7._dp/90._dp /)   &                      !(4,4)
         &,(/5,4/))
#else
    !Gauss Lobatto Chebishev POINTS 
    alpha(1,5)=1._dp/2._dp-1._dp/SQRT(2._dp)/2._dp
    alpha(2,5)=1._dp/2._dp
    alpha(3,5)=1._dp/2._dp+1._dp/SQRT(2._dp)/2._dp
    alpha(4,5)=1._dp
    n_theta(5)=5

    theta(0:4,1:4,5)=RESHAPE( (/SQRT(2._dp)/120._dp+23._dp/480._dp, & !(0,1)
         & -13._dp/480._dp*SQRT(2._dp)+2._dp/15._dp,  & !(1,1)
         & -3._dp/20._dp*SQRT(2._dp)+1._dp/5._dp, &     !(2,1)
         & -43._dp/480._dp*SQRT(2._dp)+2._dp/15._dp,  & !(3,1)
         & -1._dp/120._dp*SQRT(2._dp)-7._dp/480._dp,  &  !(4,1)
         & 1._dp/60._dp,   &                      !(0,2)
         &  1._dp/8._dp*SQRT(2._dp)+2._dp/15._dp,  &    !(1,2)
         & 1._dp/5._dp,   &                       !(2,2)
         & -1._dp/8._dp*SQRT(2._dp)+2._dp/15._dp,  &    !(3,2)
         & 1._dp/60._dp,     &                    !(4,2)
         & -SQRT(2._dp)/120._dp+23._dp/480._dp, &       !(0,3)
         &  43._dp/480._dp*SQRT(2._dp)+2._dp/15._dp,  & !(1,3)
         &  3._dp/20._dp*SQRT(2._dp)+1._dp/5._dp, &     !(2,3)
         &  13._dp/480._dp*SQRT(2._dp)+2._dp/15._dp,  & !(3,3)
         & -1._dp/120._dp*SQRT(2._dp)-7._dp/480._dp,  &  !(4,3)
         &1._dp/30._dp,    &                      !(0,4)
         &4._dp/15._dp,     &                      !(1,4)
         &2._dp/5._dp,     &                      !(2,4)
         &4._dp/15._dp,   &                      !(3,4)
         &1._dp/30._dp/)   &                      !(4,4)
         &,(/5,4/)) 
#endif

    !! The following lines are good only for the version 1
    beta(1,4,5)=alpha(1,5) ! this gives the fraction of timestep for the Euler method !*In the ugly L^1, which we do not use
    beta(2,4,5)=alpha(2,5)-alpha(1,5)
    beta(3,4,5)=alpha(3,5)-alpha(2,5)
    beta(4,4,5)=alpha(4,5)-alpha(3,5)
    GAMMA(2:4,4,5)=beta(2:4,4,5)


  END SUBROUTINE theta_alpha_beta

END PROGRAM
