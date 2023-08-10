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
   REAL(DP), DIMENSION(:), ALLOCATABLE :: error1EndPoints,error2EndPoints,errorinfEndPoints

   !*-----------------------------------------------------------
   !*GLOBAL FLUX STRUCTURES
   TYPE(Pvar), DIMENSION(:,:), ALLOCATABLE :: Rvector
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Rmatrix
   !*(Rmatrix)_ij=int_{x_L}^{x_i} phi_j
   !*row index -> x_i, DoF, end integration point
   !*column index -> phi_j, basis function 
   !*-----------------------------------------------------------
   !*STOP IF STEADY
   REAL(DP), DIMENSION(2) :: ConvergenceSteady !*It is ||Var%ua(:,DATA%iordret-1)-Var%ua(:,0)||_\infty*dx/dt COMPONENT BY COMPONENT, if it is around machine precision we stop
   !*-----------------------------------------------------------
   !*---------------------------------------------------------------------------------------------------------
   !*For the DeC staggered in time
   REAL(DP), DIMENSION(3,2) :: M1 !*Inteprpolation matrix from 2 subtimesteps to 3
   REAL(DP), DIMENSION(4,3) :: M2 !*Inteprpolation matrix from 3 subtimesteps to 4
   REAL(DP), DIMENSION(5,4) :: M3 !*Inteprpolation matrix from 4 subtimesteps to 5
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: M_interp !*Matrix actually used among the previous ones
   TYPE(variables):: VarSupp      !*Suppport structure to pass to Var with more subtimesteps
   REAL(DP), DIMENSION(:), ALLOCATABLE :: CoeffInterpOld, CoeffInterpNew !*Support structure to pass to more subtimesteps in each DoF
   INTEGER :: nSubOld,nSubNew !*Number of subtimesteps old and new in the context of an interpolation 
   INTEGER :: k_opt !*Optimal iterations  
   INTEGER :: indc !*Loop on the components
   INTEGER :: inds !*Loop on the subtimesteps
   !*---------------------------------------------------------------------------------------------------------  
   !*For the computational time
   REAL(DP) :: time_start, time_finish, time_total
   
   INTEGER :: indstart




   
   

   M1=RESHAPE( (/ 1._DP,     & 
               &   0.5_DP,    & 
               &   0._DP,     &
               &   0._DP,     &
               &   0.5_DP,    & 
               &   1._DP /)   &                      
               &,(/3,2/))

      M2=RESHAPE( (/ 1._DP,                    &
                  &   0.2222222222222223_DP,    &
                  &   -0.1111111111111111_DP,   &
                  &   0._DP,                    &
                  &   0._DP,                    &
                  &   0.8888888888888890_DP,    &
                  &   0.8888888888888890_DP,    &
                  &   0._DP,                    &
                  &   0._DP,                    &
                  &   -0.1111111111111111_DP,    & !*- inserted later
                  &   0.2222222222222222_DP,    &
                  &   1._DP   /)   &                      
                  &,(/4,3/))

      M3=RESHAPE( (/ 1._DP,          &  
                  &   0.1171875_DP,   &
                  &   -0.0625_DP,     &
                  &   0.0390625_DP,   &
                  &   0._DP,          &
                  &   0._DP,          &
                  &   1.0546875_DP,   &
                  &   0.5625_DP,      &
                  &   -0.2109375_DP,  &
                  &   0._DP,          &
                  &   0._DP,          &
                  &   -0.2109375_DP,  &
                  &   0.5625_DP,      &
                  &   1.0546875_DP,   &
                  &   0._DP,          &
                  &   0._DP,          &
                  &   0.0390625_DP,   &
                  &   -0.0625_DP,     &
                  &   0.1171875_DP,   &
                  &   1._DP   /)      &                      
                  &,(/5,4/))            





   !*I guess this will be the name of the folder where the results will be stored 
   folder='TEST'

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
   READ(1,*)DATA%iordret, DATA%iter, DATA%staggering
   PRINT*, "Staggering is: ", DATA%staggering
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
   READ(1,*) DATA%ischema, DATA%ijump, DATA%iGlobalFlux !*Choice of the scheme and of the jump
   !*alpha_jump are parameters in the stabilization with jump (for example in CIP i.e. jump stabilization i.e. Burman)
   READ(1,*) DATA%alpha_jump
   READ(1,*) DATA%alpha_jump2
   READ(1,*) DATA%cfl !*You should know what the cfl is XD
   READ(1,*) DATA%ktmax !*Maximal number of iterations
   READ(1,*) DATA%tmax !*Final time but NOT USED, it will be reinitialized in init_bc
   READ(1,*) DATA%ifre !*Frequency of storing of the results
   !*Every DATA%ifre timesteps we store the results
   READ(1,*) DATA%test, DATA%perturbation !*Choice of the test and of the perturbation
   CALL InitializeParameters(DATA%test) !*Right after the acquisition of the test case we initialize the parameters

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

   SELECT CASE(DATA%ijump)
   CASE(3,4,5,6,7,13,14,15,16,17)
      DATA%alpha_jump2=0._DP
   END SELECT

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

      !*----------------------------------------
      !*Initialize from scratch. It is used only when needed, there is an if inside.
      CALL InitializeFromScratch(DATA,Mesh,Var%un(:))
      !*Introduce a perturbation in the IC if needed
      CALL Perturbation(DATA,Mesh,Var%un(:))
      !*----------------------------------------

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


   !*NOT HERE BECAUSE YOU HAVE NOT INITIALIZED THE BATHYMETRY
   !*CALL visu(DATA, kt, mesh, Var,folder) !*In postprocessing, he stores in a file the result at this point, at the time 0 if we started from 0 or it resaves the last result if we started from the beginning

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


   !*----------------------------------------------------
   !*Depending on the order we generate the needed low-order meshes
	SELECT CASE(DATA%itype) 
		 CASE(1,11)
			 k_opt=2
	    CASE(2,3,12) !*B2, P2, PGL2
		 	 k_opt=3
		 CASE(4,5,13) !*P3, B3, PGL3
		 	 k_opt=4
		 CASE(6,7,14) !*B4, P4, PGL4
			 k_opt=5
		 CASE DEFAULT 
			 PRINT*, "Not implemented basis function. Problems in the main."
			 STOP
	END SELECT
   !*----------------------------------------------------


   CALL cpu_time(time_start)



   CALL visu(DATA, kt, mesh, Var,folder) !*In postprocessing, he stores in a file the result at this point, at the time 0 if we started from 0 or it resaves the last result if we started from the beginning


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

                  IF (DATA%ischema==24 .OR. DATA%ischema==25 .OR. DATA%ischema==-21 &
                     & .OR. DATA%ijump==6 .OR. DATA%ijump==7 &
                     & .OR. DATA%ijump==16 .OR. DATA%ijump==17) THEN
                     Rvector=0._DP
                     !*Loop over the scalar components of the global flux associated to each subtimestep
                     DO inds=0,DATA%iordret-1
                        CALL IntegrationSourceGlobalFLuxScalar(Rvector(:,inds),debug%up(:,inds),Rmatrix,Mesh,DATA)
                     END DO
                  END IF



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


         !*--------------------------------------------------------------------------------
         !*HERE IT IS MY TURN I HAVE TO DEAL WITH THE CHANGE OF SUBTIMESTEPS (INTERPOLATION)/CHANGE OF DATA BEFORE GOING AHEAD
         IF (DATA%staggering .EQV. .TRUE.) THEN 
               VarSupp=Var  !*Support var for the change of mesh at each iteration
               IF (k_inter==1) THEN
                  DEALLOCATE(Var%ua,Var%up) !*These are now in VarSupp
                  !*Now we need only two subtimesteps
                  ALLOCATE(Var%ua(Mesh%ndofs,0:1), Var%up(Mesh%ndofs,0:1)) !*NB: For the DoFs in space we are still using the same number as at the beginning, we do not touch the space in this project
                  Var%ua=0._DP
                  Var%up=0._DP
                  DO indi=1,Mesh%nDoFs 
                     Var%up(indi,0)=VarSupp%up(indi,0) !*Copying old values in VarSupp
                     Var%ua(indi,0)=VarSupp%up(indi,0) !*Copying also in Var%ua even if it is actually not needed because we'll work on Var%up to update Var%ua
                     Var%up(indi,1)=VarSupp%up(indi,DATA%iordret-1)
                     Var%ua(indi,1)=VarSupp%up(indi,DATA%iordret-1)
                  END DO
                  !*NB: Var%un not touched because it does not change
                  !*Now I also change DATA%iordret
                  DATA%iordret=2
               ELSE IF ( (k_inter.GE.2) .AND. (k_inter.LE. (k_opt-1)) ) THEN !*ALERT
                  DEALLOCATE(Var%ua,Var%up) !*These are now in VarSupp
                  !*Now we need one subtimestep more
                  ALLOCATE(Var%ua(Mesh%ndofs,0:k_inter), Var%up(Mesh%ndofs,0:k_inter)) !*NB: For the DoFs in space we are still using the same number as at the beginning, we do not touch the space in this project
                  Var%ua=0._DP
                  Var%up=0._DP
                  !*Interpolation
                  nSubOld=k_inter   !*From 0 to k_inter-1
                  nSubNew=k_inter+1 !*From 0 to k_inter
                  ALLOCATE(M_interp(nSubOld,nSubNew))
                  M_interp=0._DP
                  SELECT CASE(k_inter)
                     CASE(2)
                        M_interp=M1   
                     CASE(3)
                        M_interp=M2
                     CASE(4)
                        M_interp=M3            
                     CASE DEFAULT 
                        PRINT*, "Not possible interpolation in time staggered"
                        STOP
                  END SELECT
                  ALLOCATE(CoeffInterpOld(nSubOld), CoeffInterpNew(nSubNew))
                  CoeffInterpOld=0._DP
                  CoeffInterpNew=0._DP
                  
                  DO indi=1,Mesh%nDoFs !*Loop over all the DoFs to make the interpolation in each DoF
                     DO indc=1,n_vars !*Loop on the components so to interpolate each component
                        CoeffInterpOld=0._DP
                        CoeffInterpNew=0._DP
                        !*We need to fetch the old values from VarSupp
                        DO inds=0,nSubOld-1
                           CoeffInterpOld(inds+1)=VarSupp%up(indi,inds)%u(indc)
                        END DO 
                        CoeffInterpNew=MATMUL(M_interp,CoeffInterpOld) !*Interpolated coefficients for that component
                        !*We need to put them now in their place in Var%up
                        DO inds=0,nSubNew-1
                           Var%up(indi,inds)%u(indc)=CoeffInterpNew(inds+1)
                           Var%ua(indi,inds)%u(indc)=Var%up(indi,inds)%u(indc)
                        END DO 
                     END DO
                  END DO
                                 
                  DEALLOCATE(CoeffInterpOld, CoeffInterpNew,M_interp)
                  !*NB: Var%un not touched because it does not change
                  !*Now I also change DATA%iordret
                  DATA%iordret=k_inter+1
               ELSE
                  !*Do nothing, continue with the last number of subtimesteps
               END IF


            


         END IF




         !*-------------------------------------------------------------------------------------------
         !*Initialize GlobalFluxStuff
         IF (DATA%ischema==24 .OR. DATA%ischema==25 .OR. DATA%ischema==-21 &
            & .OR. DATA%ijump==6 .OR. DATA%ijump==7 &
            & .OR. DATA%ijump==16 .OR. DATA%ijump==17) THEN
            CALL GlobalFluxIntegrationMatrix()
            ALLOCATE(Rvector(Mesh%ndofs,0:DATA%iordret-1))
            Rvector=0._DP
         ELSE !*NONE
            !*In any case I initialize these structures with the right size but I leave them empty
            !*THIS IS DONE NOT TO HAVE MISMATCHING PROBLEMS
            ALLOCATE(Rmatrix(Mesh%e(1)%nsommets,Mesh%e(1)%nsommets))
            Rmatrix=0._DP
            ALLOCATE(Rvector(Mesh%ndofs,0:DATA%iordret-1))
            Rvector=0._DP
         END IF      
         !*-------------------------------------------------------------------------------------------







         !*----------------------------------------------------------------
         !*----------------------------------------------------------------
         !*----------------------------------------------------------------
         !*----------------------------------------------------------------
         !*----------------------------------------------------------------
         !*----------------------------------------------------------------
         !*INTEGRATE S FOR THE GLOBAL FLUX, INTRODUCE A LOOP ON THE SUBTIMESTEPS SO TO AVOID THE 0 BASED NUMERATION PROBLEM!
         !*R=-int^x S
         IF (DATA%ischema==24 .OR. DATA%ischema==25 .OR. DATA%ischema==-21 &
            &.OR. DATA%ijump==6 .OR. DATA%ijump==7 &
            &.OR. DATA%ijump==16 .OR. DATA%ijump==17) THEN
            Rvector=0._DP
            !*Loop over the scalar components of the global flux associated to each subtimestep
            DO inds=0,DATA%iordret-1
               CALL IntegrationSourceGlobalFLuxScalar(Rvector(:,inds),Var%up(:,inds),Rmatrix,Mesh,DATA)
            END DO
            !*------------------------------
            !*SAFETY CHECK
            !*DO inds=0,DATA%iordret-1
            !*   DO indi=1,Mesh%ndofs !*Loop on the elements
            !*      !*PRINT*, indi, Rvector(indi,inds)%u(1)
            !*      PRINT*, indi, Rvector(indi,inds)%u(2)+0.5_DP*gravity*Var%up(indi,inds)%u(1)**2
            !*   END DO
            !*END DO
            !*STOP
            !*------------------------------
         END IF
         !*----------------------------------------------------------------
         !*----------------------------------------------------------------
         !*----------------------------------------------------------------
         !*----------------------------------------------------------------
         !*----------------------------------------------------------------
         !*----------------------------------------------------------------

         !*!*OUR AIM IS TO COMPUTE u^(p) thorugh the DeC

         !*!*The updating formula of the DeC reads
         !*!*L^1(u^(p))=L^1(u^(p-1))-L^2(u^(p-1))
         !*!*This formula involves all the subtimesteps (and obviously all the nodes)
         !*!*We have in particular we have for each subtimestep m (and DoF i)
         !*!*u_i^m,(p)=u_i^m,(p-1)-1/|C_i|*{\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))+\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}]} (#)
         !*!*Formula (#) is performed in fin() inside the loop on the subtimesteps (inside the iteration loop of the DeC inside the timestepping loop)

         !*Not to compute all the subtimesteps in the final iteration of the DeC
         IF (k_inter /= DATA%iter) THEN !*NB: It cannot be k_opt otherwise you lose the possibility of making more iterations than necessary
            indstart=1
         ELSE
            indstart=DATA%iordret-1
         END IF


         DO k=indstart, DATA%iordret-1 !*Loop on the subtimesteps

            !*NB: The subtimesteps are in total M+1=accuracy=DATA%iordret
            !*->subtimesteps 0,1,...,M=DATA%iordret-1
            !*BUT according to the DeC the solution in the first subtimestep (0) must not be updated
            !*SO this loop must cover all the subtimestep apart from 0
            !*SO the index goes from 1 to M=DATA%iordret-1
               

            DATA%temps = tn + alpha(k,DATA%iordret)*dt !*Set the time to be the one of the treated subtimestep

            CALL fin()!*<-THIS IS THE REAL TIMESTEPPING, THE DEC, AND IT IS A SUBROUTINE DOWN HERE
         ENDDO ! k


         !*We have updated Var%ua, we can go to the next iteration (or finish the iterations if we did all of them)

         DEALLOCATE(Rmatrix,Rvector)

      ENDDO ! n_inter !*End loop on the iterations

      !*When we finish all the DeC iterations we are here

      DATA%temps= tn +dt !*We update the time

      !*BEFORE GOING AHEAD WE CHECK IF WE REACHED THE STEADY STATE
      !*We check ||Var%ua(:,DATA%iordret-1)-Var%ua(:,0)||_\infty*dx/dt, if it is around machine precision we stop
#if(1==0)
      IF (MOD(kt,DATA%ifre)==0) THEN !*Every DATA%ifre
        ConvergenceSteady=0._DP !*It is ||Var%ua(:,DATA%iordret-1)-Var%ua(:,0)||_\infty*dx/dt COMPONENT BY COMPONENT
         CALL CheckIfSteady(ConvergenceSteady,Var,Mesh,dt,DATA) !*StopIfSteady in timestepping
         PRINT*, DATA%temps, ConvergenceSteady
         IF (MAXVAL(ConvergenceSteady(:)) .LE.1.e-15_dp) EXIT
      END IF
#endif



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
         !*NO, JUST PRINT IN THE END AFTER THE TIMESTEPPING
         !*CALL MyPrintFinalSolution(Var%ua(:,0),Mesh,DATA,kt)

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


   !*We reached now the final time

   CALL cpu_time(time_finish)
   time_total=time_finish-time_start
   PRINT*, "Time of the computation= ",time_total," seconds." 
   PRINT*



   ALLOCATE(error1(n_vars),error2(n_vars),errorinf(n_vars))
   ALLOCATE(error1EndPoints(n_vars),error2EndPoints(n_vars),errorinfEndPoints(n_vars))

   error1=0._DP
   error2=0._DP
   errorinf=0._DP
   CALL errorTestProjection(error1,error2,errorinf,Var%ua(:,0),Mesh,DATA)
   error1EndPoints=0._DP
   error2EndPoints=0._DP
   errorinfEndPoints=0._DP
   CALL errorTestProjectionEndPoints(error1EndPoints,error2EndPoints,errorinfEndPoints,Var%ua(:,0),Mesh,DATA)

   !*--------------------------------------------------------------------------------------
   !*Print the numerical solution in a .dat file
   CALL MyPrintFinalSolution(Var%ua(:,0),Mesh,DATA,kt)
   !*--------------------------------------------------------------------------------------


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
       CALL main_update(k,dt,Mesh%e(jt),debug,DATA,alpha,beta,gamma,n_theta,theta,jt,Mesh,Rvector)
    ENDDO
   !--------------------------------------------------------------------------------------
    ! EDGE UPDATE DONE ONLY IN CASE WE DO NOT HAVE LXF scheme
    IF (DATA%ischema==5 .OR. DATA%ischema==4 .OR. DATA%ischema==14 .OR. DATA%ischema==15) THEN !(with jumps--> Burman'stuff)
        CALL edge_main_update(k,DATA,Mesh,debug,dt,alpha,beta,gamma,n_theta,theta,Rvector)
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
       CALL main_update(k,dt,Mesh%e(jt),Var,DATA,alpha,beta,gamma,n_theta,theta,jt,Mesh,Rvector) !*main_update in timestepping
    ENDDO

    !--------------------------------------------------------------------------------------
    ! EDGE UPDATE DONE ONLY IN CASE WE DO NOT HAVE LXF scheme
    !*!*
    IF (DATA%ischema==4 .OR.DATA%ischema==5 &
      & .OR. DATA%ischema==14 .OR. DATA%ischema==15 &
      & .OR. DATA%ischema==24 .OR. DATA%ischema==25 ) THEN !(with jumps--> Burman'stuff)
       CALL edge_main_update(k,DATA,Mesh,Var,dt,alpha,beta,gamma,n_theta,theta,Rvector) !*edge_main_update in timestepping
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
       !*PRINT*, "before", var%ua(is,k)
       var%ua(is,k)=Var%up(is,k)-Var%un(is)*Mesh%aires(is) !*Updating
       !*PRINT*, "after", var%ua(is,k)
       !*print*
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

   !*------------------------------------------------------
   !*Interpolation matrix for the global flux approach
   !*NB: Specific for PGL, if it is not PGL (or P) it will stop here
   !*------------------------------------------------------
   SUBROUTINE GlobalFluxIntegrationMatrix()
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: MYR
      REAL(DP) :: sqrt21, sqrt5
      REAL(DP), DIMENSION(:), ALLOCATABLE :: SWAPVECT
      INTEGER :: indr, indc


      SELECT CASE(DATA%itype)
         CASE(1,11) !*P1, PGL1
            ALLOCATE(Rmatrix(2,2))
            Rmatrix=0._DP
            Rmatrix= RESHAPE( (/ &
               !*First column
            &  0._DP, &
            &  0.5_DP, &
               !*Second column
            &  0._DP, &
            &  0.5_DP/), &
            &  (/2,2/) )
         CASE(3,12) !*P2, PGL2
            ALLOCATE(Rmatrix(3,3), MYR(3,3))
            Rmatrix=0._DP
            Rmatrix=RESHAPE(  (/ &
               !*First column
            &   0._DP, &
            &   0.16666666666666663_DP, &
            &   0.20833333333333329_DP, &
               !*Second column
            &   0._DP, &
            &   0.16666666666666657_DP, &
            &   -0.04166666666666667_DP, &
               !*Third column
            &   0._DP, &
            &   0.66666666666666674_DP, &
            &   0.33333333333333337_DP/),&      
            &   (/3,3/) )


            !*STILL NOT ORDERED IN MY LOCAL ORDERING
            MYR=0._DP
            MYR=RESHAPE(  (/ &
               !*First column
            &   0._DP, &
            &   5._DP/24._DP, &
            &   1._DP/6._DP, &
               !*Second column
            &   0._DP, &
            &   1._DP/3._DP, &
            &   2._DP/3._DP, &
               !*Third column
            &   0._DP, &
            &   -1._DP/24._DP, &
            &   1._DP/6._DP/),&      
            &   (/3,3/) )



         CASE(13) !*PGL3
            ALLOCATE(Rmatrix(4,4), MYR(4,4))
            Rmatrix=0._DP
            Rmatrix=RESHAPE(  (/ &
               !*First column
            &   0._DP, &
            &   0.08333333333333370_DP, &
            &   0.11030056647916493_DP, &
            &   0.07303276685416846_DP, &
               !*Second column
            &   0._DP, &
            &   0.08333333333333309_DP, &
            &   0.01030056647916491_DP, &
            &   -0.02696723314583164_DP, &
               !*Third column
            &   0._DP, &
            &   0.41666666666666741_DP, &
            &   0.18969943352083515_DP, &
            &   0.45057403089581088_DP, &
               !*Fourth column
            &   0._DP, &
            &   0.41666666666666696_DP,  &
            &   -0.03390736422914388_DP, &
            &   0.22696723314583167_DP/),&      
            &   (/4,4/) )



            !*STILL NOT ORDERED IN MY LOCAL ORDERING
            MYR=0._DP
            sqrt5=SQRT(5._DP)
            MYR=RESHAPE(  (/ &
               !*First column
            &   0._DP, &
            &   (11._DP+sqrt5)/120._DP, &
            &   (11._DP-sqrt5)/120._DP, &
            &   1._DP/12._DP, &
               !*Second column
            &   0._DP, &
            &   (25._DP-sqrt5)/120._DP, &
            &   (25._DP+13._DP*sqrt5)/120._DP, &
            &   5._DP/12._DP, &
               !*Third column
            &   0._DP, &
            &   (25._DP-13._DP*sqrt5)/120._DP, &
            &   (25._DP+sqrt5)/120._DP, &
            &   5._DP/12._DP, &
               !*Fourth column
            &   0._DP, &
            &   (-1._DP+sqrt5)/120._DP, &
            &   (-1._DP-sqrt5)/120._DP, &
            &   1._DP/12._DP/), &     
            &   (/4,4/) )


         CASE(14) !*PGL4
            ALLOCATE(Rmatrix(5,5), MYR(5,5))
            Rmatrix=0._DP
            Rmatrix=RESHAPE(  (/ &
                  !*First column
               &   0._DP, &
               &   0.05_DP, &
               &   0.06772843218615690_DP, &
               &   0.04062499999999993_DP, &
               &   0.05370013924241427_DP, &
                  !*Second column
               &   0._DP,     &
               &   0.05_DP,     &
               &   -0.00370013924241453_DP,    &
               &   0.00937500000000001_DP,     &
               &   -0.01772843218615693_DP,    &
                  !*Third column
               &   0._DP,     &
               &   0.27222222222222336_DP,     &
               &   0.11974476934341174_DP,     &
               &   0.30318418332304309_DP,     &
               &   0.26158639799680766_DP,     &
                  !*Fourth column
               &   0._DP,     &
               &   0.35555555555555579_DP,     &
               &   -0.02173572186655814_DP,    &
               &   0.17777777777777776_DP,     &
               &   0.37729127742211377_DP,     &
                  !*Fifth column
               &   0._DP,     &
               &   0.27222222222222203_DP,     &
               &   0.01063582422541549_DP,     &
               &   -0.03096196110082059_DP,    &
               &   0.15247745287881020_DP /),   &      
               &   (/5,5/) )
         
            !*STILL NOT ORDERED IN MY LOCAL ORDERING
            MYR=0._DP
            sqrt21=SQRT(21._DP)
            MYR=RESHAPE(  (/ &
                  !*First column
               &   0._DP, &
               &   (119._DP+3._DP*sqrt21)/1960._DP, &
               &   13._DP/320._DP, &
               &   (119._DP-3._DP*sqrt21)/1960._DP, &
               &   1._DP/20._DP, &
                  !*Second column
               &   0._DP,     &
               &   (343._DP-9._DP*sqrt21)/2520._DP,     &
               &   (392._DP+105._DP*sqrt21)/2880._DP,    &
               &   (343._DP+69._DP*sqrt21)/2520._DP,     &
               &   49._DP/180._DP,    &
                  !*Third column
               &   0._DP,     &
               &   (392._DP-96._DP*sqrt21)/2205._DP,     &
               &   8._DP/45._DP,     &
               &   (392._DP+96._DP*sqrt21)/2205._DP,     &
               &   16._DP/45._DP,     &
                  !*Fourth column
               &   0._DP,     &
               &   (343._DP-69*sqrt21)/2520._DP,     &
               &   (392._DP-105._DP*sqrt21)/2880._DP,    &
               &   (343._DP+9*sqrt21)/2520._DP,     &
               &   49._DP/180._DP,     &
                  !*Fifth column
               &   0._DP,     &
               &   (-21._DP+3._DP*sqrt21)/1960._DP,     &
               &   3._DP/320._DP,     &
               &   (-21._DP-3._DP*sqrt21)/1960._DP,    &
               &   1._DP/20._DP /),   &      
               &   (/5,5/) )

         CASE(4) !*P3
               ALLOCATE(Rmatrix(4,4))
               Rmatrix=0._DP
               Rmatrix=RESHAPE(  (/ &
                  !*First column
               &   0._DP, &
               &   0.125_DP, &
               &   0.125_DP, &
               &   0.1111111111111111_DP, &
                  !*Second column
               &   0._DP, &
               &   0.125_DP, &
               &   0.01388888888888889_DP, &
               &   0._DP, &
                  !*Third column
               &   0._DP, &
               &   0.375_DP, &
               &   0.26388888888888889_DP, &
               &   0.44444444444444444_DP, &
                  !*Fourth column
               &   0._DP, &
               &   0.375_DP,  &
               &   -0.06944444444444444_DP, &
               &   0.1111111111111111_DP/), &      
               &   (/4,4/) )

         CASE(7) !*P4
            ALLOCATE(Rmatrix(5,5))
            Rmatrix=0._DP
            Rmatrix=RESHAPE(  (/ &
                  !*First column
               &   0._DP, &
               &   0.07777777777777778_DP, &
               &   0.08715277777777778_DP, &
               &   0.08055555555555556_DP, &
               &   0.084375_DP, &
                  !*Second column
               &   0._DP, &
               &   0.07777777777777778_DP, &
               &   -0.00659722222222222_DP, &
               &   -0.00277777777777781_DP, &
               &   -0.009375_DP, &
                  !*Third column
               &   0._DP, &
               &   0.35555555555555556_DP, &
               &   0.22430555555555556_DP, &
               &   0.34444444444444444_DP, &
               &   0.31875_DP, &
                  !*Fourth column
               &   0._DP, &
               &   0.13333333333333333_DP, &
               &   -0.09166666666666667_DP, &
               &   0.06666666666666667_DP, &
               &   0.225_DP, &
                  !*Fifth column
               &   0._DP, &
               &   0.35555555555555556_DP, &
               &   0.03680555555555556_DP, &
               &   0.01111111111111111_DP, &
               &   0.13125_DP /),   &      
               &   (/5,5/) )

            CASE(2) !*B2
               ALLOCATE(Rmatrix(3,3))
               Rmatrix=0._DP
               Rmatrix=RESHAPE(  (/ &
                  !*First column
               &   0._DP, &
               &   0.33333333333333333_DP, &
               &   0.29166666666666667_DP, &
                  !*Second column
               &   0._DP, &
               &   0.33333333333333333_DP, &
               &   0.04166666666666667_DP, &
                  !*Third column
               &   0._DP, &
               &   0.33333333333333333_DP, &
               &   0.16666666666666669_DP/),&      
               &   (/3,3/) )

            CASE(5) !*B3
               ALLOCATE(Rmatrix(4,4))
               Rmatrix=0._DP
               Rmatrix=RESHAPE(  (/ &
                  !*First column
               &   0._DP, &
               &   0.25_DP, &
               &   0.20061728395061731_DP, &
               &   0.24691358024691357_DP, &
                  !*Second column
               &   0._DP, &
               &   0.25_DP, &
               &   0.00308641975308642_DP, &
               &   0.04938271604938269_DP, &
                  !*Third column
               &  0._DP, &
               &  0.25_DP, &
               &  0.10185185185185185_DP, &
               &  0.22222222222222222_DP, &
                  !*Fourth column
               &   0._DP, &
               &   0.25_DP,  &
               &   0.02777777777777778_DP, &
               &   0.14814814814814814_DP/), &      
               &   (/4,4/) )

         CASE(6) !*B4
            ALLOCATE(Rmatrix(5,5))
            Rmatrix=0._DP
            Rmatrix=RESHAPE(  (/ &
                  !*First column
               &   0._DP, &
               &   0.2_DP, &
               &   0.1525390625_DP, &
               &   0.19375_DP, &
               &   0.1998046875_DP, &
                  !*Second column
               &   0._DP, &
               &   0.2_DP, &
               &   0.0001953125_DP, &
               &   0.00625_DP, &
               &   0.0474609375_DP, &
                  !*Third column
               &   0._DP, &
               &   0.2_DP, &
               &   0.0734375_DP, &
               &   0.1625_DP, &
               &   0.196875_DP, &
                  !*Fourth column
               &  0._DP, &
               &  0.2_DP, &
               &  0.020703125_DP, &
               &  0.1_DP, &
               &  0.179296875_DP, &
                  !*Fifth column
               &  0._DP, &
               &  0.2_DP, &
               &  0.003125_DP, &
               &  0.0375_DP, &
               &  0.1265625_DP /),   &      
               &   (/5,5/) )

         CASE DEFAULT         
            PRINT*, "Stop in GlobalFluxIntegrationMatrix, Global flux not provided in this case"
            STOP
         END SELECT


         !*Now we need to swap a couple of things to respect our numeration
         !*In particular, both for rows and for columns
         !*FOR PGL4
         !*swap 4-5
         !*swap 3-5
         !*swap 2-3
         SELECT CASE(DATA%itype)         
            CASE(3,12,13,14)
               ALLOCATE(SWAPVECT(SIZE(MYR,DIM=1)))
               !*ROWS
               DO indr=1,SIZE(MYR,DIM=1)-2
                  SWAPVECT=0._DP
                  SWAPVECT=MYR(SIZE(MYR,DIM=1)+1-indr,:)
                  MYR(SIZE(MYR,DIM=1)+1-indr,:)=MYR(SIZE(MYR,DIM=1)+1-indr-1,:)
                  MYR(SIZE(MYR,DIM=1)+1-indr-1,:)=SWAPVECT
               END DO

               !*COLUMNS
               DO indc=1,SIZE(MYR,DIM=1)-2
                  SWAPVECT=0._DP
                  SWAPVECT=MYR(:,SIZE(MYR,DIM=1)+1-indc)
                  MYR(:,SIZE(MYR,DIM=1)+1-indc)=MYR(:,SIZE(MYR,DIM=1)+1-indc-1)
                  MYR(:,SIZE(MYR,DIM=1)+1-indc-1)=SWAPVECT
               END DO

               Rmatrix=MYR
               DEALLOCATE(MYR,SWAPVECT)
         CASE DEFAULT         
            !*Not necessary swap
         END SELECT
         
   END SUBROUTINE GlobalFluxIntegrationMatrix

   !*------------------------------------------------------
   !*Integration of the source for the global flux approach
   !*R=-int^x S
   !*NB: Routine working on the single subtimestep
   !*NB: Specific for P or PGL (coefficients are values)
   !*------------------------------------------------------
   SUBROUTINE IntegrationSourceGlobalFLuxScalar(RvectorSingleSub,UpSingleSubCoeff,Rmatrix,Mesh,DATA)
      IMPLICIT NONE
      TYPE(Pvar), DIMENSION(:), INTENT(INOUT)   :: RvectorSingleSub    !*Vector of R=-int^x S for the single subtimestep
      TYPE(Pvar), DIMENSION(:), INTENT(IN)      :: UpSingleSubCoeff    !*Solution Coefficients in all the DOFs in a single subtimestep
      REAL(DP),   DIMENSION(:,:),   INTENT(IN)  :: Rmatrix          !*Interpolation matrix
      TYPE(maillage),   INTENT(IN)              :: Mesh       
      TYPE(donnees),    INTENT(IN)              :: DATA
      !*--------------------
      !*Local variables
      TYPE(Pvar), DIMENSION(:), ALLOCATABLE      :: UpSingleSubVal     !*Solution Values in all the DOFs in a single subtimestep
      INTEGER :: indi, indj !*Loop on the DoFs
      INTEGER :: indt !*Loop on the elements
      TYPE(element) :: e !*For the quick referencing to the processed element
      !*--------------------
      !*GF1 WB
      REAL(DP),   DIMENSION(:), ALLOCATABLE          :: eta_el     !*Local eta
      REAL(DP),   DIMENSION(:), ALLOCATABLE          :: b_el       !*Local b
      REAL(DP),   DIMENSION(:), ALLOCATABLE          :: coeff_b_el !*Local coeff of b
      REAL(DP),   DIMENSION(:,:), ALLOCATABLE        :: db_el      !*Local d_x b 
      REAL(DP),   DIMENSION(:), ALLOCATABLE          :: getadb_el  !*Local g*eta*db
      REAL(DP),   DIMENSION(:), ALLOCATABLE          :: gn2qabsq_h73 !*Local g*n^2*q |q|/h^{7/3}
      !*--------------------
      !*GF2 NOT WB
      TYPE(Pvar),   DIMENSION(:), ALLOCATABLE             :: S_el
      REAL(DP),     DIMENSION(:), ALLOCATABLE             :: ms2_el !* -  Local second component of the flux !  -(-g*rho*b_x -g*h*S_f) --- S_f=n^2*hu*ABS(hu)/h^(10/3)
      !*--------------------
      REAL(DP),   DIMENSION(:), ALLOCATABLE          :: valloc, coeffloc
      !*---------------------
      !*SAFETY CHECK
      !*PRINT*, Rmatrix
      !*PRINT*
      !*PRINT*, Rmatrix(1,:)
      !*PRINT*, Rmatrix(2,:)
      !*stop
      !*---------------------

      RvectorSingleSub=0._DP
      
      !*NB: 
      !*FIRST COMPONENT 0: I need to care about the second component only, the first one must be zero because the source is 0 in that component
      !*RvectorSingleSub 0 ON THE LEFT BOUNDARY

      !*Support vector for eta, b and db in the element
      ALLOCATE(UpSingleSubVal(Mesh%Ns), &
               & eta_el(Mesh%e(1)%nsommets), b_el(Mesh%e(1)%nsommets), db_el(Mesh%e(1)%nsommets,N_dim), getadb_el(Mesh%e(1)%nsommets), &
               & gn2qabsq_h73(Mesh%e(1)%nsommets), S_el(Mesh%e(1)%nsommets), ms2_el(Mesh%e(1)%nsommets), &
               & coeff_b_el(Mesh%e(1)%nsommets), valloc(Mesh%e(1)%nsommets), coeffloc(Mesh%e(1)%nsommets) )

      
      !*Safe intialization 
      UpSingleSubVal=0._DP
      eta_el=0._DP 
      b_el=0._DP 
      db_el=0._DP
      getadb_el=0._DP
      gn2qabsq_h73=0._DP
      ms2_el=0._DP
      S_el=0._DP
      coeff_b_el=0._DP
      valloc=0._DP
      coeffloc=0._DP

      !*From coefficients to values
      UpSingleSubVal=global_Control_to_Cons(UpSingleSubCoeff,Mesh)

      DO indt=1,Mesh%nt !*Loop on the elements
         e=Mesh%e(indt)
         
         IF(DATA%iGlobalFlux==1) THEN
            !*Safe intialization before processing the element
            eta_el=0._DP 
            b_el=0._DP 
            db_el=0._DP
            getadb_el=0._DP
            gn2qabsq_h73=0._DP
            coeff_b_el=0._DP
            valloc=0._DP
            coeffloc=0._DP


            !*Local b and eta
            DO indi=1,e%nsommets
               b_el(indi)=BathymetryAtDoFsUtils(e%nu(indi))
               eta_el(indi)=b_el(indi)+UpSingleSubVal(e%nu(indi))%u(1)
               coeff_b_el(indi)=CoeffBathymetryAtDoFsUtils(e%nu(indi))
               !*--------------------------------------------------------
               !*SAFETY CHECK
               !*PRINT*, e%coor(indi), bathymetry(DATA%test,e%coor(indi))-BathymetryAtDoFsUtils(e%nu(indi))
               !*PRINT*, eta_el(indi)
               !*--------------------------------------------------------
            END DO
            !*Local grad b
            DO indi=1,e%nsommets !*Loop on the DoF where we compute the gradient
               DO indj=1,e%nsommets !*Loop on the basis functions
                  db_el(indi,:)=db_el(indi,:)+&
                  e%gradient(indj,e%x(:,indi))*coeff_b_el(indj)
               END DO
               !*--------------------------------------------------------
               !*SAFETY CHECK
               !*PRINT*, indi, db_el(indi,:)-(b_el(2)-b_el(1))/(e%coor(2)-e%coor(1))
               !*--------------------------------------------------------
            END DO

            !*Local g*eta*d_x_b
            DO indi=1,e%nsommets
               getadb_el(indi)=gravity*eta_el(indi)*db_el(indi,1)
            END DO

            !*--------------------------
            !*SAFETY CHECK
            !*DO indi=1,e%nsommets
            !*   PRINT*, indi, db_el(indi,1)-GradBathymetryAtDoFsUtils(e%nu(indi),1)
            !*END DO
            !*--------------------------
            

            !*Local g*n^2*q |q|/h^{7/3}
            DO indi=1,e%nsommets
               gn2qabsq_h73(indi)=gravity*n_manning**2 &
                                 & *UpSingleSubVal(e%nu(indi))%u(2)*ABS(UpSingleSubVal(e%nu(indi))%u(2)) &
                                 & /(UpSingleSubVal(e%nu(indi))%u(1)**(7._DP/3._DP))
            END DO

            !*-----------------------------------------
            !*-----------------------------------------
            !*In order to make the wb for Bernstein one should pass to the coefficients of getadb_el+gn2qabsq_h73
            !*-----------------------------------------
            !*-----------------------------------------

            valloc=getadb_el+gn2qabsq_h73
            SELECT CASE(e%itype)
            CASE(1,3,4,7,11,12,13,14) ! Lagrange
               !*In this case there is nothing to do because coeff=values
               coeffloc=valloc
            CASE(2,5,6)
               !*Projection needed
               coeffloc = MATMUL(e%base1,valloc)
            CASE default
               PRINT*, "erreur dans Model/Control_to_Cons"
               STOP
            END SELECT

            !*Going on with the integration of R
            !*NB: It is not an error to start from 2
            DO indi=2,e%nsommets
               !*----------------------------------
               !*SAFETY CHECK
               !*PRINT*
               !*PRINT*, Rmatrix(indi,:)
               !*PRINT*, getadb_el
               !*PRINT*, Rmatrix(indi,:)*getadb_el
               !*----------------------------------
               RvectorSingleSub(e%nu(indi))%u(2)=RvectorSingleSub(e%nu(1))%u(2) &
                                                & +SUM(Rmatrix(indi,:)* coeffloc)*(e%coor(2)-e%coor(1)) & !*ALERT TO FRICTION
                                                & -(0.5_DP*gravity*b_el(indi)**2-0.5_DP*gravity*b_el(1)**2) 

            END DO         

            

            !*-------------------------------------
            !*SAFETY CHECK
            !*DO indi=1,e%nsommets
            !*   PRINT*, indi, RvectorSingleSub(e%nu(indi))%u(2)+0.5_DP*gravity*UpSingleSub(e%nu(indi))%u(1)**2
            !*END DO         
            !*-------------------------------------

            !*------------------------------------
            !*SAFETY CHECK
            !*PRINT*, getadb_el
            !*PRINT*, getadb_el+getadb_el
            !*STOP
            !*------------------------------------
         
         ELSE IF (DATA%iGlobalFlux==2)  THEN
            !*Safe intialization before processing the element
            S_el=0._DP
            ms2_el=0._DP 
            !*Local second component of the flux
            DO indi=1,e%nsommets
               S_el(indi)=sourcetermfunction(e%coor(indi),UpSingleSubVal(e%nu(indi)),e%nu(indi),DATA%test)
               ms2_el(indi)=-S_el(indi)%u(2)
            END DO

            !*Going on with the integration of R
            !*NB: It is not an error to start from 2
            DO indi=2,e%nsommets
               RvectorSingleSub(e%nu(indi))%u(2)=RvectorSingleSub(e%nu(1))%u(2) &
                                                & +SUM(Rmatrix(indi,:)*ms2_el)*(e%coor(2)-e%coor(1)) 
            END DO         
         ELSE
            PRINT*, "No global flux provided in IntegrationSourceGlobalFLuxScalar in main"
            STOP
         END IF



      END DO !*End loop on the elements      

      !*-------------------------------------------------------------
      !*SAFETY CHECK
      !*DO indi=1,Mesh%ndofs !*Loop on the elements
      !*   !*PRINT*, indi, RvectorSingleSub(indi)%u(2)
      !*   PRINT*, indi, RvectorSingleSub(indi)%u(2)+0.5_DP*gravity*UpSingleSub(indi)%u(1)**2
      !*END DO
      !*stop
      !*DO indi=1,Mesh%ndofs !*Loop on the elements
      !*   PRINT*, indi, UpSingleSub(indi)%u(1)+BathymetryAtDoFsUtils(indi)
      !*   PRINT*, indi, UpSingleSub(indi)%u(2)
      !*   PRINT*
      !*END DO
      !*-------------------------------------------------------------
      

      DEALLOCATE(UpSingleSubVal, eta_el, b_el, db_el, getadb_el, gn2qabsq_h73, S_el, ms2_el,coeff_b_el,&
                 & valloc, coeffloc)

   END SUBROUTINE IntegrationSourceGlobalFLuxScalar



END PROGRAM
