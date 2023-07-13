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

!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*THIS MODULE CONTAINS
!*-> SUBROUTINE precomputation(Mesh, VecProx)
!*   It takes the mesh and gives as output VecProx which is a vector. Each element corresponds to a DoF and has
!*      - number of the elements containing that DoF
!*      - global indices of those elements
!*      - number of DoFs in those elements (without repetitions if node is shared)
!*      - global indices of those elements (without repetitions)
!*      NB:The order is always min->max
!*-> SUBROUTINE internalandboundarydofs(DATA,bDoFs,iDoFs,NewNumeration,OldNumeration)
!*   It takes in inputs the DATA and reads the mesh like in geom.f90 and gives
!*      - bDoFs global indices of the boundary DoFs
!*      - iDoFs global indices of the internal DoFs
!*      - OldNumeration, Vector containing first the global indices of the bDoFs and then the golobal indices of the iDoFs. OldNumeration(48) is the global index of the node which in the new numeration is number 48. So it is the old index in the new numeration. It is meant to pass from the new numeration to the old one.
!*      - NewNumeration, Vector containing the position of the (old) DoF in the new numeration. For example NewNumeration(24) is the place in the new numeration of the DoF with global index 24.
!*      NB: I need it for the BCs but it is actually a bit useless since gmesh generates first the boundary. So in practice there is never reordering. Anyway I implemented everything in the general context in which the boundary DoFs are not mandatorily the first ones.
!*-> SUBROUTINE stiffnessmatrix(CRSStMaSc1,Mesh,VecProx,bDoFs,iDoFs,NewNumeration,OldNumeration,switchBCs)
!*      It takes the Mesh and the previous structures and computes the CRS of the stiffness matrix in 1-based-index-numeration. I take care of the boundary conditions in this way: I compute the nonzero entries of the matrix in the original numeration and I store them in an intermediate structure tempNnzSc which is a vector. Each element is referred to a row (so to a DoF i) and contains the number of nonzero entries, their position (the columns i.e. the close DoFs to i taken from VecProx) and the nonzero entries. Then I pass to NnzSc in the new numeration: I move the rows of the boundary nodes up and change the numeration of the columns (not the entries obviously) and I reorder everything because in principle in this moving one loses the order of the columns (but I remark, in practice nothing happens because gmesh generates from the boundary and the boundary nodes are already the first ones). Then if I have Neumann in switchbs(=1) I do nothing, if I have Dirichlet(=0) I change the boundary rows and put the identity and remove the other elements of those columns. 
!*      Then there is another thing that I have implemented before passing to the CRS: Remotion of the zero elements -> For example in P2 the entries 1,5 2,6 3,4 of any local stiffness matrix are 0. Obviously the calculator does not know so he calculates and gets something close to zero machine in those positions and treats them as nonzero elements. So I remove the elements whose ABS() is less than a parameter eps set to 10^(-15). This part anyway can be skipped setting #if(1==0). It just makes the stiffness matrix smaller.
!*      I call then a subroutine contained in this subroutine that passes from the nonzero entries (NnzSc) to the CRS.
!*      NB: It is very "clean" from the point of view of the memory. I have for every structure that I use a cleaning routine which avoid the wasting of memory (in every SUBROUTINE, not just in this one). So basically everything that is not needed anymore is cleaned in the end. The matrix is the scalar one. But we need that one :)
!*      NB: If you will ever need the mass matrix the code is exactly the same provided that you substitute the calculation of the local stiffness matrix with the calculation of the local mass matrix.
!*-> FUNCTION global_Control_to_Cons(Ucoeff,Mesh) RESULT (Uval)
!*-> FUNCTION global_Cons_to_Control(Uval,Mesh) RESULT (Ucoeff)
!*      They are functions that pass from the vector of PVar of dimension equal to the number of DoFs containing the coefficients in the nodes to the values in the nodes and viceversa. We had just the local version in model so I thought It was nice to avoid the loop on the elements that we do every time. This is how I discovered the bug in model for P2 :)
!*-> SUBROUTINE CleaningVecProx(VecProx)
!*      Because VecProx is not used anymore after the construction of the CRS structure of the stiffness matrix
!*-> SUBROUTINE CleaningCRS(CRSstr)
!*      
!*-> SUBROUTINE errortest(error1,error2,errorinf,CoeffNumSol,Mesh)
!*      It calculates the difference between my solution and the exact solution and its L^1,L^2 and L^\infty norm
!*      The inputs are
!*         1) CoeffNumSol, coefficients of the solution to the elliptic problem in the original numeration
!*         2) Mesh
!*      The (in)outputs are
!*         error1,error2,errorinf, respectively the L^1,L^2 and L^\infty norm of the error w.r.t. the exact solution. 
!*-> SUBROUTINE gradient_reconstruction(GradAtDoFs,CoeffNumSol,VecProx,Mesh)
!*      It takes in input the coeffients CoeffNumSol, VecProx and Mesh and computes (for each component) the gradient of the solution in each DoF by averaging between in the DoFs contained by more than one element
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




MODULE preprocessing
  USE PRECISION
  USE param2d
  USE OVERLOADING
  USE MODEL
  USE UTILS

  
  IMPLICIT NONE
   CHARACTER(LEN = *), PARAMETER :: mod_name = "ellipticmodule"
  
  REAL(DP), PARAMETER, PUBLIC :: smallepsilon = 1.e-15_DP !*It will be the parameter for the remotion of the zero entries

  !*Needed types

  !*MOVED IN Param2D
  !*TYPE Proximity !*Element of the vector referred to a single DoF i
  !*   INTEGER  :: numbofel !*number of elements containing i
  !*   INTEGER, DIMENSION(:), ALLOCATABLE :: vecEl !*Elements containing i
  !*   INTEGER  :: numbofDoFs !*number of DoFs in the elements containing i (without redoundances) 
  !*   INTEGER, DIMENSION(:), ALLOCATABLE :: vecDoFs !*Global indices of those DoFs
  !*END TYPE Proximity
  
  !*Temporary, to fill the list of triangles containing a node and of nodes close to one node
  TYPE cellule !*cellule=integer(val)+pointer to cellule(suiv)
     INTEGER :: val
     TYPE(cellule), POINTER :: next
  END TYPE cellule
  
  !*NB: I'll use 1-based indices.
  TYPE CRS !*Compressed row storage
     REAL(DP), DIMENSION(:), ALLOCATABLE :: values !*Nonzero entries
     INTEGER, DIMENSION(:), ALLOCATABLE :: columns !*columns of the nonzero entries
     INTEGER, DIMENSION(:), ALLOCATABLE :: rowstart !*Index of elements in val that is at the start of the given row
  END TYPE



!*GMSH_ELE(type_ele=8 or 9,1 or 2)
!*1->dim of nu
!*2->dim of vv

!*vertex_ele(type_ele=8 or 9)
!*number of vertices

CONTAINS
  !*--------------------------------------------------------------------------------------------
  !*Fill VecProx
  !*VecProx is a vector of elements of dimension=number of nodes whose generic component i is referred to the node i and gives the elements containing that node and the nodes in those elements
  !*--------------------------------------------------------------------------------------------
  SUBROUTINE precomputation(Mesh, VecProx)
     IMPLICIT NONE
     TYPE (maillage), INTENT(IN) :: Mesh
     TYPE(Proximity), DIMENSION(:), ALLOCATABLE :: VecProx !*Vector 
     INTEGER :: indt=0 !*Index for the loop on the triangles
     INTEGER :: indi=0 !*Index for the loop on the DoFs i
     INTEGER :: indj=0 !*Index for the loop on the DoFs j
          
     TYPE (element) :: e !*For the quick reference to the generic Mesh%e(indt) in the loop

     !*Temporary pointer structures
     TYPE(cellule), DIMENSION(:), ALLOCATABLE, TARGET :: tempEl, tempDoFs
     !*tempEL is a temporary vector of cellules, the component i will carry the list of elements containing i
     !*tempDoFs is a temporary vector of cellules, the component i will carry the list of DoFs close to i (in at least one element containing i)
     !*(tempDoFs gives us the nonzero entries of the mass matrix: colptr and rowptx
      
     TYPE(cellule), POINTER :: Cell, NextCell, AddCell=>NULL() !*Support structures for the inserting 
      
     
     ALLOCATE(VecProx(Mesh%Ns))
     DO indi=1,Mesh%Ns
        VecProx(indi)%numbofel=0
        VecProx(indi)%numbofDoFs=0 
     END DO
     
     !*We need some preprecomputations involving temporary structures

     !*1) We make a loop on all the elements indt
     !*2) We make a loop on all the nodes indi to add indt to the list of triangles containing indi
     !*3) We make a double loop on all the the DoFs indi and on all the DoFs indj to add the nodes indj(/=indi) to the list of nodes close to indi
     
     !*NB: The addition of new triangles or new nodes must be in order min->max and must avoid redoundances 
     !*We work on some temporary pointer structures, then later we will fill VecProx

     !*Initialization of the temporary structures
     ALLOCATE(tempEl(Mesh%Ns),tempDoFs(Mesh%Ns))
     DO indi=1,Mesh%Ns
        tempEl(indi)%val=indi !*First value set equal to the index of the node
        tempDoFs(indi)%val=indi !*First value set equal to the index of the node
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !*SAFETY CHECK
        !*PRINT*, indi, tempEl(indi)%val, tempDoFs(indi)%val
        !*PRINT*, ASSOCIATED(tempEl(indi)%next), ASSOCIATED(tempDoFs(indi)%next)
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(ASSOCIATED(tempEl(indi)%next)) NULLIFY(tempEl(indi)%next)
        IF(ASSOCIATED(tempDoFs(indi)%next)) NULLIFY(tempDoFs(indi)%next) !*Fixed ;)
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !*SAFETY CHECK
        !*PRINT*, ASSOCIATED(tempEl(indi)%next), ASSOCIATED(tempDoFs(indi)%next) !*F,F
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     END DO

     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*TEST
     !*ALLOCATE(AddCell)
     !*AddCell%val=1122334455
     !*IF(ASSOCIATED(AddCell%next)) NULLIFY(AddCell%next)
     !*PRINT*, AddCell%val, ASSOCIATED(AddCell%next)
     !*tempEl(1)%next=>AddCell
     !*PRINT*, AddCell%val, ASSOCIATED(AddCell%next)
     !*PRINT*, tempEl(1)%val, ASSOCIATED(tempEl(1)%next), tempEl(1)%next%val, ASSOCIATED(tempEl(1)%next%next)
     !*ALLOCATE(AddCell)
     !*AddCell%val=987654321
     !*IF(ASSOCIATED(AddCell%next)) NULLIFY(AddCell%next)
     !*PRINT*, tempEl(1)%val, ASSOCIATED(tempEl(1)%next), tempEl(1)%next%val, ASSOCIATED(tempEl(1)%next%next)
     !*tempEl(1)%next%next=>AddCell
     !*PRINT*, tempEl(1)%val, ASSOCIATED(tempEl(1)%next), tempEl(1)%next%val, &
     !*& ASSOCIATED(tempEl(1)%next%next), tempEl(1)%next%next%val, ASSOCIATED(tempEl(1)%next%next%next)
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     DO indt=1,Mesh%nt !*Loop on the triangles
        e=Mesh%e(indt)
        !*Let's save that the triangle indt contains its DoFs in e%nu(:)
        DO indi=1,e%nsommets !*Loop on the DoFs i of the triangle
           !*We now scroll the list of the triangles of that DoF in order to find the right position where to insert that triangle (NB:We'll also check if it is already there)
           Cell=>tempEl(e%nu(indi)) 
           NextCell=>Cell%next
           DO WHILE (ASSOCIATED(NextCell))
              IF(NextCell%val .GE. indt) THEN !*We found the position, between Cell and NextCell
                 EXIT
              ELSE
                 Cell=>NextCell
                 NextCell=>Cell%next
              END IF
           END DO
           !*Now we can have two cases: 
           !*1) We are at the end of the list 
           !*NextCell is not associated
           !*The triangle is not already there otherwise we would have exited before
           IF (.NOT. ASSOCIATED(NextCell)) THEN
              VecProx(e%nu(indi))%numbofel=VecProx(e%nu(indi))%numbofel+1
              ALLOCATE(AddCell)
              IF(ASSOCIATED(AddCell%next)) NULLIFY(AddCell%next)
              AddCell%val=indt
              !*Create the bridge with Cell
              Cell%next=>AddCell
           ELSE !*We are in the middle of the list
              !*We have to memorize indt only if it is not already present
              IF(indt<NextCell%val) THEN
                 VecProx(e%nu(indi))%numbofel=VecProx(e%nu(indi))%numbofel+1
                 ALLOCATE(AddCell)
                 AddCell%val=indt
                 !*Create the bridges with Cell and NextCell
                 Cell%next=>AddCell
                 AddCell%next=>NextCell
              END IF
           END IF
        END DO !*End loop on the DoFs
     END DO !*End loop on the triangles

     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*SAFETY CHECK
     !*indi=2000
     !*PRINT*, "DoF ", indi, "has ", VecProx(indi)%numbofel
     !*Cell=>tempEl(indi) 
     !*DO indt=1,VecProx(indi)%numbofel
     !*   Cell=>Cell%next
     !*   PRINT*, "Element", indt, "containing the node", indi, "is", Cell%val
     !*   PRINT*, Mesh%e(Cell%val)%nu(:)
     !*END DO
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !*Now before passing to the nodes let's store the information on the triangles in the structure VecProx
     DO indi=1,Mesh%Ns
        ALLOCATE(VecProx(indi)%vecEl(VecProx(indi)%numbofel))
        Cell=>tempEl(indi)
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !*SAFETY CHECK
        !*PRINT*, indi, "Numb of el: ", VecProx(indi)%numbofel
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO indt=1,VecProx(indi)%numbofel
           Cell=>Cell%next
           VecProx(indi)%vecEl(indt)=Cell%val
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*SAFETY CHECK
           !*PRINT*, indi, VecProx(indi)%vecEL(indt)
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        END DO
     END DO

     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*SAFETY CHECK
     !*DO indi=1,Mesh%Ns
     !*   PRINT*, "Node ", indi, "->", VecProx(indi)%vecEl(:)
     !*END DO 
     !*STOP
     !*indi=300
     !*PRINT*, "Node ", indi, "->", VecProx(indi)%vecEl(:)
     !*DO indt=1,SIZE(VecProx(indi)%vecEl(:),DIM=1)
     !*   PRINT*,"Element", VecProx(indi)%vecEl(indt), "has the nodes", &
     !*   & Mesh%e(VecProx(indi)%vecEl(indt))%nu(:)
     !*END DO
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     CALL CleaningTemp(tempEl) 
     DEALLOCATE(tempEl) !*Not needed anymore   

     !*Now let's go to the nodes
     !*We want to do the following. We have now the number of triangles containing each node and their gloabl indices. We want to store, for each DoF i, the number of the DoFs in the elements containing i and their global indices.
     !*NB: They are the nonzero entries of the Stiffness-matrix of the elliptic problem 
     !*Just like we did for the elements, we work with a temporary structure of pointers
     !*tempDoFs which has already been initialized with tempEl
     
     !*We first make a loop on the DoFs, indi
     !*Then we make a loop on the elements containing that DoF indi thanks to VecProx(indi)%vecEl, indt
     !*We make a loop over the DoFs indj storing indj in the list of the nodes close to indi

     DO indi=1,Mesh%Ns !*Loop on the DoFs, indi
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !*SAFETY CHECK
        !*PRINT*, indi, "Elements: ", VecProx(indi)%vecEl(:)
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO indt=1,VecProx(indi)%numbofel !*Loop on the elements containing indi, indt 
           e=Mesh%e(VecProx(indi)%vecEl(indt))
           DO indj=1,e%nsommets !*Loop on the DOFs of the element containing the node indi, indj
           !*We want to store these nodes in the list tempDoFs as usual from min to max and without repetitions
              !*NB: We refer to indi and not to e%nu(indi) but we also refer to e%nu(indj)
              
              !*We now scroll the list of the triangles of that DoF in order to find the right position where to insert that triangle (NB:We'll also check if it is already there)
              Cell=>tempDoFs(indi) 
              NextCell=>Cell%next
              DO WHILE (ASSOCIATED(NextCell))
                 IF(NextCell%val .GE. e%nu(indj)) THEN !*We found the position, between Cell and NextCell
                    EXIT
                 ELSE
                    Cell=>NextCell
                    NextCell=>Cell%next
                 END IF
              END DO
             !*Now we can have two cases: 
             !*1) We are at the end of the list 
             !*NextCell is not associated
             !*The DoF is not already there otherwise we would have exited before
             IF (.NOT. ASSOCIATED(NextCell)) THEN
                 VecProx(indi)%numbofDoFs=VecProx(indi)%numbofDoFs+1
                 ALLOCATE(AddCell)
                 IF(ASSOCIATED(AddCell%next)) NULLIFY(AddCell%next)
                 AddCell%val=e%nu(indj)
                 !*Create the bridge with Cell
                 Cell%next=>AddCell
              ELSE !*We are in the middle of the list
                 !*We have to memorize e%nu(indj) only if it is not already present
                 IF(e%nu(indj)<NextCell%val) THEN
                    VecProx(indi)%numbofDoFs=VecProx(indi)%numbofDoFs+1
                    ALLOCATE(AddCell)
                    AddCell%val=e%nu(indj)
                    !*Create the bridges with Cell and NextCell
                    Cell%next=>AddCell
                    AddCell%next=>NextCell
                 END IF
              END IF


           END DO !*End loop on the DoFs of the element, indj
        END DO !*End loop on the elements, indt
     END DO !*End first loop on the nodes, indi

     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*SAFETY CHECK
     !*indi=4320 
     !*indi=1315 
     !*PRINT*, "DoF ", indi, "has ", VecProx(indi)%numbofel 
     !*PRINT*, "The node", indi, "is contained by the elements: ", VecProx(indi)%vecEl(:)    
     !*DO indt=1, VecProx(indi)%numbofel
     !*   PRINT*, VecProx(indi)%vecEl(indt), "Nodes: ", Mesh%e(VecProx(indi)%vecEl(indt))%nu(:)
     !*END DO
     !*Close nodes
     !*PRINT*, "DoF ", indi, "has ", VecProx(indi)%numbofDoFs, "voisins" 
     !*Cell=>tempDoFs(indi) 
     !*DO indj=1,VecProx(indi)%numbofDoFs
     !*     Cell=>Cell%next
     !*     PRINT*, "DoF", indi, "close to the node", Cell%val
     !*END DO
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !*Let's store the information on the DoFs in the structure VecProx
     DO indi=1,Mesh%Ns
        ALLOCATE(VecProx(indi)%vecDoFs(VecProx(indi)%numbofDoFs))
        Cell=>tempDoFs(indi)
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !*SAFETY CHECK
        !*PRINT*, indi, "Numb of DoFs: ", VecProx(indi)%numbofDoFs
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO indj=1,VecProx(indi)%numbofDoFs
           Cell=>Cell%next
           VecProx(indi)%vecDoFs(indj)=Cell%val
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*SAFETY CHECK
           !*PRINT*, indi, VecProx(indi)%vecDoFs(indj)
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        END DO
     END DO

     CALL CleaningTemp(tempDoFs)      
     DEALLOCATE(tempDoFs) !*Not needed anymore


     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*SAFETY CHECK
     !*DO indi=1,Mesh%Ns
     !*   PRINT*, "Node ", indi, "->", VecProx(indi)%vecDoFs(:)
     !*   IF(indi==400) STOP
     !*END DO      
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
     !*Last final Check
     !*indi=807
     !*PRINT*, "Node ", indi
     !*PRINT*, "It has ", VecProx(indi)%numbofDoFs, "close DoFs " 
     !*PRINT*, VecProx(indi)%vecDoFs(:)
     !*PRINT*, "It is contained by ", VecProx(indi)%numbofEl, "elements"
     !*PRINT*, VecProx(indi)%vecEl(:)
     !*DO indt=1,VecProx(indi)%numbofEl
     !*   PRINT*, "Triangle: ", VecProx(indi)%vecEl(indt)
     !*   PRINT*, Mesh%e(VecProx(indi)%vecEl(indt))%nu(:)
     !*END DO
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
     !*REMARK: i is close to itself


     CONTAINS




        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !*Cleaning of a temporary structure
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !*A SMALL THEORY RECAP :)

        !!POINTERS are locations in memory not variables
        
        !!- They can be associated to existing variables through '=>' and deassociated NULLIFY()
        !!TIP: When declared set => NULL() to avoid the uncertainity of the state UNDEFINED
        
        !!- They can be allocated i.e. associated to an unnamed part of the memory through ALLOCATE

        !!AN IMPORTANT NOTE ON THE DEALLOCATION
        !!ALLOCATE not only creates the desired object in type and shape but it PRESERVES A GIVEN (UNNAMED) PART of the MEMORY 
        !!That part of memory is not a variable and it is only accessible through the pointer
        !!If we nullify the pointer (all the pointers to that part of memory) or associate it through another TARGET without deallocating first we cannot access anymore to that part of memory and moreover it is still there preserved until the end of the program
        !!So we have to DEALLOCATE a given part of memory when we don't need that information anymore
        !!NB: If we have more pointers to a given part of memory and we DEALLOCATE one of them the memory is free but the other pointers are still linked to it and it is not safe since the computer is now free to use it as it likes
        !!We have to NULLIFY the others or better their links to that part of memory witht the command NULLIFY.
        !!AND NB: We cannot DEALLOCATE because even if the link is still to be destroyed, the memory is already free and the command would result in a segmentation fault.
        !!DEALLOCATE in practice both frees the memory and NULLIFY the link to it.
        !!In this case it cannot free the memory because it is already free
        !!So one must be DEALLOCATED and the other NULLIFIED or assigned to other variables
        
        !!The DEALLOCATION only can be used on arguments created by an ALLOCATE statement, not with pointers associated to existing variables



        SUBROUTINE CleaningTemp(temp)
           TYPE(cellule), DIMENSION(:), ALLOCATABLE, TARGET :: temp
           !*Support structures for the dellocation loop
           TYPE(cellule), POINTER :: Cell, NextCell
           INTEGER :: indtemp, j
           INTEGER :: counter
           

           IF (ALLOCATED(temp)) THEN
              DO indtemp=1,SIZE(temp,DIM=1)
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !*SAFETY CHECK
                 !*PRINT*, "I'm here", indtemp
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !*We make a deallocaction loop
                 Cell=>temp(indtemp)%next 
                 NextCell=>Cell%next
                 counter=1
                 DO WHILE (ASSOCIATED(NextCell))
                    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !*SAFETY CHECK
                    !*PRINT*, 1, ASSOCIATED(Cell), ASSOCIATED(NextCell)
                    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    counter=counter+1
                    DEALLOCATE(Cell)
                    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !*SAFETY CHECK
                    !*PRINT*, 2, ASSOCIATED(Cell), ASSOCIATED(NextCell)
                    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    Cell=>NextCell !*Which still exists
                    NextCell=>Cell%next !*Which now exists
                    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !*SAFETY CHECK
                    !*PRINT*, 3, ASSOCIATED(Cell), ASSOCIATED(NextCell)
                    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 END DO
                 !*Last nullify
                 DEALLOCATE(Cell)
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !*SAFETY CHECK
                 !*PRINT*, "Number of deallocations in ", indtemp, "is ", counter
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              END DO
           

              DO indtemp=1, SIZE(temp, DIM=1)
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !*SAFETY CHECK
                 !*PRINT*, "Association: ", ASSOCIATED(temp(indtemp)%next) !*Of course still associated but to a free memory location in fact if I print the value I have a random value
                 !*PRINT*, temp(indtemp)%next%val
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 NULLIFY(temp(indtemp)%next)
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !*SAFETY CHECK
                 !*PRINT*, ASSOCIATED(temp(indtemp)%next)
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              END DO
           END IF

        END SUBROUTINE CleaningTemp

  END SUBROUTINE precomputation



  !*--------------------------------------------------------------------------------------------
  !*Research for the boundary and internal nodes
  !*NB:
  !*In 2d I read them from the mesh
  !*In 1d I directly know that the boundary DoFs are 
  !*1 and nt+1
  !*So I can directly fill bDoFs, iDoFs, NewNumeration and OldNumeration
  !*In this case I do not need to search, so I do not need temporary structures
  !*--------------------------------------------------------------------------------------------
  SUBROUTINE internalandboundarydofs(Mesh,bDoFs,iDoFs,NewNumeration,OldNumeration) !*Read the mesh and return the vector of the global indices of the DoFs on the boundary
     IMPLICIT NONE
     TYPE (maillage), INTENT(IN) :: Mesh
     INTEGER, DIMENSION(:), ALLOCATABLE :: bDoFs !*Vector of the boundary DoFs after processing i.e. ordered min->max and without redoundaces
     INTEGER, DIMENSION(:), ALLOCATABLE :: iDoFs !*Vector of the internal DoFs ordered min->max
     INTEGER, DIMENSION(:), ALLOCATABLE :: NewNumeration, OldNumeration 
     !*We make a new numeration to solve the elliptic problem
     !*We'll move the bDoFs up and the then we'll have the iDoFs
     !*NewNumeration->Vector containing the position of the DoF in the new numeration
     !*For example NewNumeration(24) is the place in the new numeration of the DoF with global index 24
     !*OldNumeration->Vector containing first the global indices of the bDoFs and then the golobal indices of the iDoFs OldNumeration(48) is the global index of the node which in the new numeration is number 48
     !*NB: n 
     !*NewNumeration(n)=new position
     !*OldNumeration(NewNumeration(n))=n

     !*Loop variable
     INTEGER :: indi, i, j
     !*Useful variable used to fill iDoFs
     INTEGER :: counter, check 

     ALLOCATE(bDoFs(2), iDoFs(Mesh%nDoFs-2), NewNumeration(Mesh%nDoFs), OldNumeration(Mesh%nDoFs))   
     !*Initialization of the allocated structures
     bDoFs=0
     iDoFs=0
     NewNumeration=0
     OldNumeration=0

     !*!!!!!!!!!!!!!!!!!!!!!!!!
     !*SAFETY CHECK to assess that they are all 0
     !*PRINT*, bDoFs
     !*PRINT*, iDoFs
     !*PRINT*, NewNumeration
     !*PRINT*, OldNumeration
     !*!!!!!!!!!!!!!!!!!!!!!!!!

     !*Boundary DoFs: 1 and nt+1
     bDoFs(1)=1
     bDoFs(2)=Mesh%nt+1



     !*Now we isolate the internal DoFs    
     check=0
     counter=0
     !*Now we store the internal DoFs
     DO i=1,Mesh%nDoFs !*Loop on all the nodes
        DO j=1,SIZE(bDoFs,DIM=1) !*Loop on the boundary nodes
           IF (i==bDoFs(j)) THEN
              check=1 !*The Node is a boundary node and we must not store it
              EXIT
           END IF 
        END DO
        IF (check==0) THEN !*The internal Node must be stored
           counter=counter+1
           iDoFs(counter)=i
        END IF
        check=0
     END DO
     
     !*Now we have the vector of the internal nodes: iDoFs
     !*and the vector of the boundary nodes:bDoFs
 
     !*!!!!!!!!!!!!!!!!!!!!!!!!
     !*SAFETY CHECK 
     !*PRINT*, "bDoFs", bDoFs
     !*PRINT*, "iDoFs", iDoFs
     !*!!!!!!!!!!!!!!!!!!!!!!!!
 
     OldNumeration(1:SIZE(bDoFs,DIM=1))=bDoFs(:)
     OldNumeration(SIZE(bDoFs,DIM=1)+1:SIZE(bDoFs,DIM=1)+SIZE(iDoFs,DIM=1))=iDoFs(:)
 
     DO i=1,SIZE(NewNumeration,DIM=1)
        NewNumeration(OldNumeration(i))=i
     END DO

     !*!!!!!!!!!!!!!!!!!!!!!!!!
     !*SAFETY CHECK 
     !*PRINT*, "OldNumeration", OldNumeration
     !*PRINT*, "NewNumeration", NewNumeration
     !*!!!!!!!!!!!!!!!!!!!!!!!!

    
  END SUBROUTINE internalandboundarydofs







 !*--------------------------------------------------------------------------------------------
 !*Construction of the CRS FORMAT of the Stiffness-matrix
 !*--------------------------------------------------------------------------------------------
 SUBROUTINE stiffnessmatrix(CRSStMaSc1,Mesh,VecProx,bDoFs,iDoFs,NewNumeration,OldNumeration,switchBCs)
    IMPLICIT NONE
    TYPE(CRS), INTENT(INOUT) :: CRSStMaSc1 !*REMARK: boundary nodes moved at the top and columns switched coherently with the new numeration 
    !*Identity submatrix considered in place of the first boundary rows if Dirichlet (swithcBCs=0) or classical matrix considered if Neumann (swithcBCs=1)
    !*REMARK: 1-based indices will be used.
    
    !*TYPE(CRS) :: CRSStMaSc0 !*Remark: CRSStMaSc0 -> CRS structure of the scalar stiffness matrix in 0 based system !*Can be useful but not needed its calculation is in a commented part

    !*TYPE(CRS) :: CRSStMaSc1supp !*CRS structure of the scalar stiffness matrix in 1 based system !*CRSStMaSc1supp actually is not used, I used it in a part which is now commented
        
    !*NB: We do not compute the CRS structure of the matrix of the whole vectorial problem. We solve separately each component.

    TYPE (maillage), INTENT(IN) :: Mesh
    TYPE(Proximity), DIMENSION(:), INTENT(IN) :: VecProx !*Vector of proximity
    INTEGER, DIMENSION(:), INTENT(IN) :: bDoFs,iDoFs
    INTEGER, DIMENSION(:), INTENT(IN) :: NewNumeration, OldNumeration
    !*NewNumeration->Vector containing the position of the DoF in the new numeration
    !*OldNumeration->Vector containing the original index of the DoF in the new numeration vector
    INTEGER, INTENT(IN) :: switchBCs
    !*0=Dirichlet
    !*1=Neumann
  
    !*We refer now to the scalar Stiffness-Matrix
    TYPE NzeScStMaDoF !*Nonzero entries of the Scalar Stiffness matrix for a certain DoF i !*It will be the element of a vector
       INTEGER :: numbofDoFs !*Number of DoFs "close" to i !*Basically VecProx(i)%numbofDoFs 
       INTEGER, DIMENSION(:), ALLOCATABLE :: vecDoFs !*Global indices of those DoFs !*Basically VecProx(i)%vecDoFs
       REAL(DP), DIMENSION(:), ALLOCATABLE :: entries !*Nonzero entries
    END TYPE NzeScStMaDoF

    TYPE(NzeScStMaDoF), DIMENSION(:), ALLOCATABLE :: tempNzeSc !*temporary vector of the nonzero entries row by row
    TYPE(NzeScStMaDoF), DIMENSION(:), ALLOCATABLE :: NzeSc !*final vector of the nonzero entries row by row
    !*We have a final vector because tempNzeSc is in the original order. We then have to move the boundary DoFs up according to the new numeration and switch the columns coherently
    !*and trick a bit those DoFs in the Dirichlet case

    !*Loop variables
    INTEGER :: i,j,k

    !*Loop variables for the evaluation the integrals (for the nonzero elements of the Stiffness-Matrix)
    INTEGER :: indt=0 !*Index for the loop on the triangles
    INTEGER :: indi=0 !*Index for the first loop on the DoFs i
    INTEGER :: indq=0 !*Index for the loop on the quadrature points to make the integral over the element
    INTEGER :: indj=0 !*Index for the second loop on the DoFs j
    
    !*Support variables for the integration
    TYPE (element) :: e !*For the quick reference to the generic Mesh%e(indt) in the loop
     
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dummyStMa !*Local Stiffness-matrix in the loop 
    !*We fill this one at first and then we send the components to the right places in tempNzeSc thanks to e%nu
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: grad !*Gradient of the basis functions in a specific quadrature point
    !*First index=dimension
    !*Second index=basis function
    INTEGER, PARAMETER :: ndim=1 !*number of the physical dimensions needed to allocate grad

    !*Support variable for scrolling tempNzeSc
    INTEGER :: scroll
    !*Support variable for check
    INTEGER :: check

    !*Support variables for reordering after the change of numeration
    INTEGER :: suppint
    REAL(DP) :: suppreal
     
    !*Support variables to impose Dirichlet 
    INTEGER :: counter !*Used to count the elements to delete to impose the Dirichelt BCs
    REAL(DP), DIMENSION(:), ALLOCATABLE :: suppvalues !*Used as support vector in the imposition of Dirichlet
    INTEGER, DIMENSION(:), ALLOCATABLE :: suppcolumns !*Used as support vector in the imposition of Dirichlet
    INTEGER :: intdummy !*Used in the imposition of the Dirichlet BCs

    !*eps very small to remove the "zero" elements
    REAL(DP), PARAMETER :: eps = smallepsilon
    !*For example in P2 the local entries 1,5 2,6 3,4 are 0. Obviously the calculator does not know so he calculates and get something close to zero machine in those positions and treats them as nonzero elements. This parameter is used to remove those elements if we set 1==1 in the #if


    !*Allocation of tempNzeSc
    ALLOCATE(tempNzeSc(Mesh%Ns))
    !*We first initialize the vector tempNzeSc
    !*Some structures have to be copied from VecProx
    DO i=1,Mesh%Ns
       tempNzeSc(i)%numbofDoFs=VecProx(i)%numbofDoFs 
       ALLOCATE(tempNzeSc(i)%vecDoFs(SIZE(VecProx(i)%vecDoFs,DIM=1)))
       tempNzeSc(i)%vecDoFs=VecProx(i)%vecDoFs 
       ALLOCATE(tempNzeSc(i)%entries(SIZE(VecProx(i)%vecDoFs,DIM=1)))
       tempNzeSc(i)%entries(:)=0._DP !*Initialization before the calculation
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*IF(i==200) THEN
       !*   PRINT*, tempNzeSc(i)%numbofDoFs
       !*   PRINT*, tempNzeSc(i)%vecDoFs(:)
       !*   PRINT*, tempNzeSc(i)%entries(:)
       !*   STOP
       !*END IF
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END DO

    !*Filling of the nonzero entries
    !*To fill this mass-matrix I'm going to make
    !*1)Loop over the elements K
    !*2)Loop over the DoFs i\inK
    !*3)Loop over the quadrature points for the integration
    !*4)Loop over the DoFs j\inK
 
    !*For each element we fill a dummy local stiffness-matrix and we send its components to the right place in the end 
    ALLOCATE(dummyStMa(Mesh%e(1)%Nsommets,Mesh%e(1)%Nsommets))
    dummyStMa=0._DP !*Initialization
    ALLOCATE(grad(ndim,Mesh%e(1)%Nsommets))
    grad=0._DP !*Initialization
    

    DO indt=1,Mesh%nt !*Loop on the triangles
       dummyStMa=0.0_DP
       grad=0.0_DP
       e=Mesh%e(indt) !*Quick reference to the processed triangle
       DO indi=1,e%nsommets !*Loop on the row DoF i
          DO indq=1,e%nquad !*Loop on the quadrature points
             DO indj=1,e%nsommets !*Loop on the column DoF j->To collect the values of the gradient of the basis functions in grad(:,indj)
                grad(:,indj)=e%gradient(indj, e%quad(:,indq))
             ENDDO !*End loop on the column DoF j
             DO indj=1,e%nsommets !*Integration loop
                dummyStMa(indi,indj)=dummyStMa(indi,indj)+ &
                & DOT_PRODUCT(grad(:,indi),grad(:,indj))*e%weight(indq)
             END DO
          END DO !*End loop on the quadrature points 
       END DO !*End loop on the row DoF i
       !*Now dummyStMa is filled
       !*Multiplication for the Jacobian
       dummyStMa=dummyStMa*e%volume
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*PRINT*, "Local Stiffness-matrix: ", dummyStMa
       !*IF(indi==10) STOP
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !*Now before going to the next element let's map this local Stiffness matrix to its place in tempNzeSc
       DO indi=1,e%nsommets !*Loop on the row DoF i
          DO indj=1,e%nsommets !*Loop on the column DoF j
             !*We need to check in which position e%nu(indj) is in tempNzeSc(e%nu(indi))%vecDoFs
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !*SAFETY CHECK
             check=0
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !*SAFETY CHECK
             !*PRINT*, e%nu(indj)
             !*PRINT*, tempNzeSc(e%nu(indi))%vecDoFs(:)
             !*IF(indt==10) STOP
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             DO scroll=1,tempNzeSc(e%nu(indi))%numbofDoFs !*Scrolling the DoFs close to e%nu(indi)
                IF(e%nu(indj)==tempNzeSc(e%nu(indi))%vecDoFs(scroll)) THEN !*We found the position
                   tempNzeSc(e%nu(indi))%entries(scroll)=tempNzeSc(e%nu(indi))%entries(scroll)+&
                   & dummyStMa(indi,indj)
                   !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !*SAFETY CHECK
                   check=1
                   !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                END IF
             END DO
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !*SAFETY CHECK
             IF(check==0) THEN
                PRINT*, "It seems we have an error because e%nu(indj) is not in the list of the DoFs close to e%nu(indi)"
                STOP
             END IF
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          END DO !*End loop on the column DoF j   
       END DO !*End loop on the row DoF i
    END DO !*End loop on the triangles 
    !*Now the Stiffness-matrix is filled    
    DEALLOCATE(dummyStMa) !*Not needed anymore
    DEALLOCATE(grad) !*Not needed anymore



    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*Let'check if there are zero elements in the nonzero entries
    !*DO i=1,Mesh%Ns
    !*   DO j=1,tempNzeSc(i)%numbofDoFs
    !*      IF(tempNzeSc(i)%entries(j)==0._DP) PRINT*, "Zero element", i, j 
    !*   END DO
    !*END DO
    !*STOP
    !*No zero elements present
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !*Now we can define NzeSc and move the boundary DoFs up 
    !*NB: After we move the bDoFs up ALSO THE COLUMNS have to be switched coherently with the new numeration
    !*Allocation of NzeSc
    ALLOCATE(NzeSc(Mesh%Ns))

    !*First we move the rows accoding to the new numeration
    !*NB: OldNumeration has the original indices of the DoFs
    !*First the global indices of the bDoFs, then the global indices of the iDoFs
    DO i=1,Mesh%Ns
       NzeSc(i)%numbofDoFs=VecProx(OldNumeration(i))%numbofDoFs 
       ALLOCATE(NzeSc(i)%vecDoFs(SIZE(VecProx(OldNumeration(i))%vecDoFs,DIM=1))) !*<- ALERT: OldNumeration(i)
       NzeSc(i)%vecDoFs=VecProx(OldNumeration(i))%vecDoFs !*<- ALERT: OldNumeration(i)
       ALLOCATE(NzeSc(i)%entries(SIZE(VecProx(OldNumeration(i))%vecDoFs,DIM=1)))
       NzeSc(i)%entries(:)=0._DP !*Initialization before copying from tempNzeSc
       NzeSc(i)%entries=tempNzeSc(OldNumeration(i))%entries
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*PRINT*, i, Oldnumeration(i)
       !*PRINT*, NzeSc(i)%entries
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END DO  

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO i=1,Mesh%Ns
    !*   PRINT*, NzeSc(i)%numbofDoFs!*-tempNzeSc(Oldnumeration(i))%numbofDoFs
    !*   PRINT*, NzeSc(i)%vecDoFs!*-tempNzeSc(Oldnumeration(i))%vecDoFs
    !*   PRINT*, NzeSc(i)%entries!*-tempNzeSc(Oldnumeration(i))%entries
    !*END DO
    !*PRINT*, NzeSc(2)%numbofDoFs-tempNzeSc(Mesh%nt+1)%numbofDoFs
    !*PRINT*, NzeSc(2)%vecDoFs-tempNzeSc(Mesh%nt+1)%vecDoFs
    !*PRINT*, NzeSc(2)%entries-tempNzeSc(Mesh%nt+1)%entries
    !*PRINT*, NzeSc(7)%numbofDoFs-tempNzeSc(6)%numbofDoFs
    !*PRINT*, NzeSc(7)%vecDoFs-tempNzeSc(6)%vecDoFs
    !*PRINT*, NzeSc(7)%entries-tempNzeSc(6)%entries
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL CleaningNze(tempNzeSc)
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO i=1,SIZE(tempNzeSc)
    !*    IF(ALLOCATED(tempNzeSc(i)%vecDoFs) .OR. ALLOCATED(tempNzeSc(i)%entries)) THEN
    !*       PRINT*, "Error in stiffnessmatrix. I can't deallocate tempNzeSc(i)%vecDoFs or tempNzeSc(i)%entries"
    !*       STOP
    !*    END IF
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DEALLOCATE(tempNzeSc) !*Not needed anymore
    
   
    !*NB:Now we have to switch also the columns because 
    !*NzeSc(i)%vecDoFs=VecProx(OldNumeration(i))%vecDoFs and NzeSc(i)%entries=tempNzeSc(OldNumeration(i))%entries are referred to the previous numeration
    !*We use Newnumeration which has the position in the new numeration of the DoF
    !*For example NewNumeration(24) is the place in the new numeration of the DoF with global index 24
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*PRINT*, "This number should be 1: ", NewNumeration(bDoFs(1))
    !*PRINT*, "This number should be 14: ", NewNumeration(bDoFs(14))
    !*NB: In principle bDoFs(1) and bDoFs(14) shouldn't be 1 and 14 but they are because gmesh generates starting from the boudary
    !*PRINT*, bDoFs(1), bDoFs(14) 
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !*So let's switch the columns
    DO i=1,Mesh%Ns !*Loop on the rows
       DO j=1,NzeSc(i)%numbofDoFs !*Loop on the columns to change the indices
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, NzeSc(i)%vecDoFs(j), NewNumeration(NzeSc(i)%vecDoFs(j)) 
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          NzeSc(i)%vecDoFs(j)=NewNumeration(NzeSc(i)%vecDoFs(j))
          !*REMARK: NewNumeration has the new position
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*IF(i==300) PRINT*, NzeSc(i)%vecDoFs(j), NewNumeration(NzeSc(i)%vecDoFs(j))
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO
    END DO

    
    
    !*Now let's reorder the columns from min->max
    !*In general the order is wasted if we apply NewNumeration
    !*REORDERING
    !*NB: This is just a reordering, no indices are changed. It's just to have also in the new numeration min->max but actually the indices of the columns and their entries will not be touched
    !*Moreover it is useless because gmesh generates the Mesh starting by the boundary
    DO i=1,Mesh%Ns !*Loop on the rows
       DO j=1,NzeSc(i)%numbofDoFs-1 !*Loop on the columns to order the indices
          !*We make a loop to reorder NzeSc(i)%vecDoFs(:) AND we move coherently also NzeSc(i)%entries(:)
          DO k=j+1,NzeSc(i)%numbofDoFs
             IF(NzeSc(i)%vecDoFs(k)<NzeSc(i)%vecDoFs(j)) THEN !*We need to change the order
                 !*We change the order of the columns
                 suppint=NzeSc(i)%vecDoFs(j)
                 NzeSc(i)%vecDoFs(j)=NzeSc(i)%vecDoFs(k)
                 NzeSc(i)%vecDoFs(k)=suppint
                 !*Coherently we need to change also the entries
                 suppreal=NzeSc(i)%entries(j)
                 NzeSc(i)%entries(j)=NzeSc(i)%entries(k)
                 NzeSc(i)%entries(k)=suppreal
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !*SAFETY CHECK
                 !*PRINT*, "I made a change!", i, NzeSc(i)%vecDoFs(j), NzeSc(i)%vecDoFs(k)
                 !*PRINT*, "Entry moved", NzeSc(i)%entries(j)
                 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             END IF
          END DO
       END DO
    END DO
    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK to assess the order
    !*DO i=1,Mesh%Ns
    !*   PRINT*, NzeSc(i)%vecDoFs
    !*END DO
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

    !*Now we have the matrix in the new numeration
    !*Let's insert the switch on the BCs
    SELECT CASE(switchBCs)
       CASE(0) !*Dirichlet: remove the boundary nodes
          !*The first rows corresponding to the boudary have to be turned into the identity matrix
          !*Only one entry in position ii and that entry has to be 1
          DO i=1,SIZE(bDoFs, DIM=1)
             NzeSc(i)%numbofDoFs=1 !*Only one entry
             DEALLOCATE(NzeSc(i)%vecDoFs,NzeSc(i)%entries) !*NB:It was allocated for sure
             ALLOCATE(NzeSc(i)%vecDoFs(1),NzeSc(i)%entries(1))
             NzeSc(i)%vecDoFs(1)=i
             NzeSc(i)%entries(1)=1.0_DP
          END DO

          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, NzeSc(1)%vecDoFs
          !*PRINT*, NzeSc(1)%entries
          !*PRINT*, NzeSc(2)%vecDoFs
          !*PRINT*, NzeSc(2)%entries
          !*PRINT*, NzeSc(3)%vecDoFs
          !*PRINT*, NzeSc(3)%entries
          !*PRINT*, NzeSc(Mesh%nt)%vecDoFs
          !*PRINT*, NzeSc(Mesh%nt)%entries
          !*PRINT*, NzeSc(Mesh%nt+1)%vecDoFs
          !*PRINT*, NzeSc(Mesh%nt+1)%entries
          !*STOP
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !*Not only, we also have to remove those columns (it's not mandatory actually)
          !*We make now a loop on the iDoFs rows to remove the columns from 1 to SIZE(bDoFs,DIM=1) in each row
          DO i=SIZE(bDoFs,DIM=1)+1,SIZE(NzeSc,DIM=1) !*Loop on the DoFs corresponding to the internal nodes
             !*Let's count the elements to delete in the row i
             counter=0 !*Initialize the counter 
             DO j=1,NzeSc(i)%numbofDoFs !*Loop on the columns
                IF (NzeSc(i)%vecDoFs(j)<=SIZE(bDoFs,DIM=1)) THEN !*If the element belongs to a boundary column we need to remove it
                   NzeSc(i)%entries(j)=0._DP !*For security
                   counter=counter+1
                   !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !*SAFETY CHECK
                   !*PRINT*, NzeSc(i)%vecDoFs(j), NzeSc(i)%entries(j), counter
                   !*PRINT*, "Before removing"
                   !*PRINT*, "internal DoF", i, "counter", counter
                   !*PRINT*, NzeSc(i)%vecDoFs
                   !*PRINT*, NzeSc(i)%entries
                   !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                END IF
             END DO

             !*Now counter is the number of elements to delete in that row
             !*The correction is performed ONLY if counter>0
             IF(counter>0) THEN
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !*SAFETY CHECK
                !*Before the correction
                !*PRINT*, i, "Before"
                !*PRINT*, NzeSc(i)%numbofDoFs, "among which remove", counter
                !*PRINT*, NzeSc(i)%VecDoFs
                !*PRINT*, NzeSc(i)%entries
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                intdummy=0
                ALLOCATE(suppvalues(NzeSc(i)%numbofDoFs-counter),suppcolumns(NzeSc(i)%numbofDoFs-counter))
                suppvalues=0._DP
                suppcolumns=0
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !*SAFETY CHECK
                !*PRINT*, suppvalues
                !*PRINT*, suppcolumns
                !*PRINT*, "Number of bDoFs", SIZE(bDoFs,DIM=1)
                !*PRINT*, "Elements to remove", counter
                !*PRINT*, "Still to clean. The columns are now: ", NzeSc(i)%vecDoFs
                !*STOP
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                DO j=1,NzeSc(i)%numbofDoFs !*Loop on the columns to save the right elements
                   IF (NzeSc(i)%vecDoFs(j)>SIZE(bDoFs,DIM=1)) THEN !*We have to keep the element
                      intdummy=intdummy+1
                      suppcolumns(intdummy)=NzeSc(i)%vecDoFs(j)
                      suppvalues(intdummy)=NzeSc(i)%entries(j)
                   END IF
                END DO
                !*Now suppvalues(NzeSc(i)%numbofDoFs-counter==intdummy
                !*!!!!!!!!!!!!!!!!!!!!!!!!!
                !*SAFETY CHECK
                !*IF (counter/=0) THEN
                !*   PRINT*, "Node", i, "-> remove", counter 
                !*   PRINT*, NzeSc(i)%vecDoFs
                !*   PRINT*, NzeSc(i)%entries
                !*   PRINT*, suppcolumns
                !*   PRINT*, suppvalues
                !*END IF
                !*!!!!!!!!!!!!!!!!!!!!!!!!!
                !*SAFETY CHECK
                !*IF (NzeSc(i)%numbofDoFs-counter/=intdummy) THEN
                !*   PRINT*, "ERROR in the imposition of BCs Dirichlet"
                !*END IF
                !*!!!!!!!!!!!!!!!!!!!!!!!!!
                !*Now we correct the row
                NzeSc(i)%numbofDoFs=NzeSc(i)%numbofDoFs-counter !*Correction of the number of entries
                DEALLOCATE(NzeSc(i)%vecDoFs,NzeSc(i)%entries) !*Deallocation
                ALLOCATE(NzeSc(i)%vecDoFs(NzeSc(i)%numbofDoFs),NzeSc(i)%entries(NzeSc(i)%numbofDoFs)) !*Allocation with the right size
                NzeSc(i)%vecDoFs=0
                NzeSc(i)%entries=0._DP
                NzeSc(i)%vecDoFs(:)=suppcolumns(:) !*Correction of the columns
                NzeSc(i)%entries(:)=suppvalues(:) !*Correction of the values
                DEALLOCATE(suppvalues,suppcolumns) !*Deallocation of the support structures
             END IF
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !*SAFETY CHECK
             !*IF (counter/=0) THEN !*It means that we corrected something
             !*   !*After the correction
             !*   PRINT*, i, "After"
             !*   PRINT*, NzeSc(i)%numbofDoFs
             !*   PRINT*, NzeSc(i)%VecDoFs
             !*   PRINT*, NzeSc(i)%entries
             !*END IF
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          END DO
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*DO i=SIZE(bDoFs,DIM=1)+1,SIZE(NzeSc,DIM=1)
          !*   DO j=1,NzeSc(i)%numbofDoFs
          !*      IF(NzeSc(i)%vecDoFs(j)<=SIZE(bDoFs,DIM=1)) THEN
          !*         PRINT*, "Very big problem. Not removed correctly the columns of the boundary DoFs"
          !*         STOP
          !*      END IF
          !*   END DO
          !*END DO
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CASE(1) !*Neumann: do nothing
          !*Really nothing to do
       CASE DEFAULT
          PRINT*, "Not defined BCs. I'm stopping in stiffnessmatrix. switchBCs is: ", switchBCs
    END SELECT

#if(1==0)
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*Remove the zero elements
    !*For example in P2 in 2D the local entries 1,5 2,6 3,4 are 0. Obviously the calculator does not know so he calculates and get something close to zero machine in those positions and treats them as nonzero elements. Here I remove them if they are less than eps
    PRINT*, "Remotion of the zero terms activated in the computation of the Stiffness Matrix"
    DO i=1,SIZE(NzeSc,DIM=1) !*Loop on all the DoFs
       !*Let's count the elements to delete in the row i
       counter=0 !*Initialize the counter 
       DO j=1,NzeSc(i)%numbofDoFs !*Loop on the columns
          IF (ABS(NzeSc(i)%entries(j))<=eps) THEN !*If the element is "very small" we need to remove it
             NzeSc(i)%entries(j)=0._DP !*For security
             counter=counter+1
          END IF
       END DO
       !*Now counter is the number of elements to delete in that row
       !*The correction is performed ONLY if counter>0
       IF(counter>0) THEN
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*Before the correction
          !*PRINT*, i, "Before"
          !*PRINT*, NzeSc(i)%numbofDoFs, "among which remove", counter
          !*PRINT*, NzeSc(i)%VecDoFs
          !*PRINT*, NzeSc(i)%entries
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
          intdummy=0
          ALLOCATE(suppvalues(NzeSc(i)%numbofDoFs-counter),suppcolumns(NzeSc(i)%numbofDoFs-counter))
          suppvalues=0._DP
          suppcolumns=0
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, suppvalues
          !*PRINT*, suppcolumns
          !*PRINT*, "Number of bDoFs", SIZE(bDoFs,DIM=1)
          !*PRINT*, "Elements to remove", counter
          !*PRINT*, "Still to clean. The columns are now: ", NzeSc(i)%vecDoFs
          !*STOP
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO j=1,NzeSc(i)%numbofDoFs !*Loop on the columns to save the right elements
             IF (ABS(NzeSc(i)%entries(j))>eps) THEN !*We have to keep the element
                intdummy=intdummy+1
                suppcolumns(intdummy)=NzeSc(i)%vecDoFs(j)
                suppvalues(intdummy)=NzeSc(i)%entries(j)
             END IF
          END DO
          !*Now suppvalues(NzeSc(i)%numbofDoFs-counter==intdummy
          !*!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*IF (counter/=0) THEN
          !*   PRINT*, "Node", i, "-> remove", counter 
          !*   PRINT*, NzeSc(i)%vecDoFs
          !*   PRINT*, NzeSc(i)%entries
          !*   PRINT*, suppcolumns
          !*   PRINT*, suppvalues
          !*END IF
          !*!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*IF (NzeSc(i)%numbofDoFs-counter/=intdummy) THEN
          !*   PRINT*, "ERROR in the imposition of BCs Dirichlet"
          !*END IF
          !*!!!!!!!!!!!!!!!!!!!!!!!!!
          !*Now we correct the row
          NzeSc(i)%numbofDoFs=NzeSc(i)%numbofDoFs-counter !*Correction of the number of entries
          DEALLOCATE(NzeSc(i)%vecDoFs,NzeSc(i)%entries) !*Deallocation
          ALLOCATE(NzeSc(i)%vecDoFs(NzeSc(i)%numbofDoFs),NzeSc(i)%entries(NzeSc(i)%numbofDoFs)) !*Allocation with the right size
          NzeSc(i)%vecDoFs=0
          NzeSc(i)%entries=0._DP
          NzeSc(i)%vecDoFs(:)=suppcolumns(:) !*Correction of the columns
          NzeSc(i)%entries(:)=suppvalues(:) !*Correction of the values
          DEALLOCATE(suppvalues,suppcolumns) !*Deallocation of the support structures
       END IF
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*IF (counter/=0) THEN !*It means that we corrected something
       !*   !*After the correction
       !*   PRINT*, i, "After"
       !*   PRINT*, NzeSc(i)%numbofDoFs
       !*   PRINT*, NzeSc(i)%VecDoFs
       !*   PRINT*, NzeSc(i)%entries
       !*END IF
       !*IF (counter>6) THEN 
       !*   PRINT*, "We are a looooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooot"
       !*   STOP
       !*END IF
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif   




    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO i=1,Mesh%Ns 
    !*   PRINT*, NzeSc(i)%numbofDoFs
    !*   PRINT*, NzeSc(i)%vecDoFs(:)
    !*   PRINT*, NzeSc(i)%entries(:)
    !*   IF(i==1000) STOP
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*Checking simmetry
    !*Thanks to the sky it's symmetric, fiuuuuuuuuuuuuuuu ;)SS
    !*i=1000
    !*PRINT*, i
    !*PRINT*, "Number of entries", NzeSc(i)%numbofDoFs
    !*PRINT*, NzeSc(i)%vecDoFs(:)
    !*PRINT*, NzeSc(i)%entries(:)
    !*j=342
    !*PRINT*, j
    !*PRINT*, "Number of entries", NzeSc(j)%numbofDoFs
    !*PRINT*, NzeSc(j)%vecDoFs(:)
    !*PRINT*, NzeSc(j)%entries(:)
    !*PRINT*, Mesh%e(2000)%nu(6)
    !*PRINT*, Mesh%e(2000)%coorL(:,6)
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !*Now I can consider the CRS structures
    !*I will first construct the structures in 1-based numeration system i.e.
    !*vectors indices (1:end)
    !*first column index=1
    !*first element in rowstart=1
    !*Then I will construct the structures in 0-based numeration system thanks to a SUBROUTINE CRS1to0 contained in this SUBROUTINE
    !*vectors indices (0:end-1)
    !*first column index=0
    !*first element in rowstart=0


    !*First I construct the CRS structure of the scalar stiffness matrix in 1 based system
    !*CRSStMaSc1
    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*Recalling the CRS structure
    !*TYPE CRS !*Compressed row storage
    !*    REAL(DP), DIMENSION(:), ALLOCATABLE :: values !*Nonzero entries
    !*    INTEGER, DIMENSION(:), ALLOCATABLE :: columns !*columns of the nonzero entries
    !*    INTEGER, DIMENSION(:), ALLOCATABLE :: rowstart !*Index of elements in val that is at the start of the given row
     !*END TYPE
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     !*The following is commented because I moved it into the subroutine FromNzeToCRS so that it is easier to debug
     !*This sub gives from the nonzero entries the CRS structure
     !*I made some safety check to make sure that the result is the same 

     !*SO GO FROM HERE TO------>

     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!
     !*!!!!*Let's count the nonzero vectors so that we have the dimension of values and columns, we already know the dimension of rowstart which is number of rows+1
     !*!!!    
     !*!!!
     !*!!!
     !*!!!counter=0
     !*!!!DO i=1,Mesh%Ns
     !*!!!   !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!   !*SAFETY CHECK
     !*!!!   !*IF(NzeSc(i)%numbofDoFs==0) THEN
     !*!!!   !*   PRINT*, "Something terribly wrong, singular matrix. An empty row."
     !*!!!   !*   STOP
     !*!!!   !*END IF
     !*!!!   counter=counter+NzeSc(i)%numbofDoFs
     !*!!!   !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!   !*PRINT*, counter
     !*!!!   !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!END DO
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!*PRINT*, counter
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!
     !*!!!!*Now we can allocate the structures of CRSStMaSc1
     !*!!!ALLOCATE(CRSStMaSc1%values(1:counter))
     !*!!!ALLOCATE(CRSStMaSc1%columns(1:counter))
     !*!!!ALLOCATE(CRSStMaSc1%rowstart(1:Mesh%Ns+1))
     !*!!!
     !*!!!!*Let's store the values row by row
     !*!!!k=0 !*k will be the index in values and columns
     !*!!!DO i=1,Mesh%Ns
     !*!!!   DO j=1,NzeSc(i)%numbofDoFs
     !*!!!      k=k+1
     !*!!!      CRSStMaSc1%values(k)=NzeSc(i)%entries(j)
     !*!!!      CRSStMaSc1%columns(k)=NzeSc(i)%vecDoFs(j)
     !*!!!   END DO
     !*!!!END DO
     !*!!!
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!*SAFETY CHECK
     !*!!!!*PRINT*, CRSStMaSc1%values(77000:)
     !*!!!!*PRINT*, CRSStMaSc1%columns(77000:)
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!
     !*!!!!*Let's now store the starting row index
     !*!!!CRSStMaSc1%rowstart(1)=1
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!*SAFETY CHECK
     !*!!!!*PRINT*, "Start row ", 1, "with", CRSStMaSc1%rowstart(1)
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!DO i=1,Mesh%Ns
     !*!!!   CRSStMaSc1%rowstart(i+1)=CRSStMaSc1%rowstart(i)+NzeSc(i)%numbofDoFs
     !*!!!   !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!   !*SAFETY CHECK
     !*!!!   !*PRINT*, "Adding ", NzeSc(i)%numbofDoFs 
     !*!!!   !*PRINT*, "Start row ", i+1, "with", CRSStMaSc1%rowstart(i+1)
     !*!!!   !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!END DO
     !*!!!
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!*SAFETY CHECK
     !*!!!!*PRINT*, CRSStMaSc1%rowstart(1)
     !*!!!!*PRINT*, NzeSc(1)%numbofDoFs
     !*!!!!*PRINT*, CRSStMaSc1%rowstart(2)
     !*!!!!*PRINT*, CRSStMaSc1%rowstart(256)
     !*!!!!*PRINT*, CRSStMaSc1%rowstart(257)
     !*!!!!*PRINT*, NzeSc(257)%numbofDoFs
     !*!!!!*PRINT*, CRSStMaSc1%rowstart(258)
     !*!!!!*i=Mesh%Ns
     !*!!!!*PRINT*, CRSStMaSc1%rowstart(i)
     !*!!!!*PRINT*, NzeSc(i)%numbofDoFs
     !*!!!!*PRINT*, CRSStMaSc1%rowstart(i+1)
     !*!!!!*PRINT*, "End", counter
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!
     !*!!!
     !*!!!CALL FromNzeToCRS(CRSStMaSc1supp,NzeSc)
     !*!!!
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!*SAFETY CHECK
     !*!!!DO i=1,SIZE(CRSStMaSc1%values,DIM=1)
     !*!!!   IF(CRSStMaSc1supp%values(i)/=CRSStMaSc1%values(i)) THEN
     !*!!!      PRINT*, "Brutta storia"
     !*!!!      STOP
     !*!!!   END IF
     !*!!!END DO
     !*!!!DO i=1,SIZE(CRSStMaSc1supp%values,DIM=1)
     !*!!!   IF(CRSStMaSc1supp%values(i)/=CRSStMaSc1%values(i)) THEN
     !*!!!      PRINT*, "Brutta storia"
     !*!!!      STOP
     !*!!!   END IF
     !*!!!END DO
     !*!!!DO i=1,SIZE(CRSStMaSc1%columns,DIM=1)
     !*!!!   IF(CRSStMaSc1supp%columns(i)/=CRSStMaSc1%columns(i)) THEN
     !*!!!      PRINT*, "Brutta storia"
     !*!!!      STOP
     !*!!!   END IF
     !*!!!END DO
     !*!!!DO i=1,SIZE(CRSStMaSc1supp%columns,DIM=1)
     !*!!!   IF(CRSStMaSc1supp%columns(i)/=CRSStMaSc1%columns(i)) THEN
     !*!!!      PRINT*, "Brutta storia"
     !*!!!      STOP
     !*!!!   END IF
     !*!!!END DO
     !*!!!DO i=1,SIZE(CRSStMaSc1%rowstart,DIM=1)
     !*!!!   IF(CRSStMaSc1supp%rowstart(i)/=CRSStMaSc1%rowstart(i)) THEN
     !*!!!      PRINT*, "Brutta storia"
     !*!!!      STOP
     !*!!!   END IF
     !*!!!END DO
     !*!!!DO i=1,SIZE(CRSStMaSc1supp%rowstart,DIM=1)
     !*!!!   IF(CRSStMaSc1supp%rowstart(i)/=CRSStMaSc1%rowstart(i)) THEN
     !*!!!      PRINT*, "Brutta storia"
     !*!!!      STOP
     !*!!!   END IF
     !*!!!END DO
     !*!!!!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
     !*TO------> HERE

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*If I modify a random element of the matrix
    !*NzeSc(SIZE(bDoFs, DIM=1)+3)%entries(1)=0.005_DP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
     CALL FromNzeToCRS(CRSStMaSc1,NzeSc)
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !"SAFETY CHECK
     !*DO i=1,SIZE(CRSStMaSc1%values,DIM=1)
     !*   IF(CRSStMaSc1%values(i)==0._DP) PRINT*, "Problem"
     !*END DO
     !*PRINT*, CRSStMaSc1%values
     !*PRINT*, CRSStMaSc1%columns
     !*PRINT*, CRSStMaSc1%rowstart
     !*PRINT*, SIZE(CRSStMaSc1%values,DIM=1)
     !*PRINT*, SIZE(CRSStMaSc1%columns,DIM=1)
     !*PRINT*, SIZE(CRSStMaSc1%rowstart,DIM=1)
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     CALL CleaningNze(NzeSc)
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*SAFETY CHECK
     !*DO i=1,SIZE(NzeSc)
     !*    IF(ALLOCATED(NzeSc(i)%vecDoFs) .OR. ALLOCATED(NzeSc(i)%entries)) THEN
     !*       PRINT*, "Error in stiffnessmatrix. I can't deallocate NzeSc(i)%vecDoFs or NzeSc(i)%entries"
     !*       STOP
     !*    END IF
     !*END DO
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DEALLOCATE(NzeSc) !*Not needed anymore


     !*Now the CRS structure of the scalar matrix is ready in system of numeration 1-based
     !*Let's compute the CRS structure 0-based
     !*We actually don't need it but can be useful
     !*CALL CRS1to0(CRSStMaSc0,CRSStMaSc1)

     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !*SAFETY CHECK
     !*PRINT*, SIZE(CRSStMaSc1%values), SIZE(CRSStMaSc0%values)
     !*PRINT*, SIZE(CRSStMaSc1%columns), SIZE(CRSStMaSc0%columns) 
     !*PRINT*, SIZE(CRSStMaSc1%rowstart), SIZE(CRSStMaSc0%rowstart) 
     !* 
     !*DO i=0,SIZE(CRSStMaSc1%values,DIM=1)-1
     !*   IF(CRSStMaSc0%values(i)/=CRSStMaSc1%values(i+1)) THEN
     !*      PRINT*, "Error in the passage from 1-based to 0-based - values", i
     !*      STOP
     !*   END IF 
     !*END DO
     !*DO i=0,SIZE(CRSStMaSc1%columns,DIM=1)-1
     !*   IF(CRSStMaSc0%columns(i)/=CRSStMaSc1%columns(i+1)-1) THEN
     !*      PRINT*, "Error in the passage from 1-based to 0-based - columns", i
     !*      STOP
     !*   END IF 
     !*END DO
     !*DO i=0,SIZE(CRSStMaSc1%rowstart,DIM=1)-1
     !*   IF(CRSStMaSc0%rowstart(i)/=CRSStMaSc1%rowstart(i+1)-1) THEN
     !*      PRINT*, "Error in the passage from 1-based to 0-based - rowstart", i
     !*      STOP
     !*   END IF 
     !*END DO
     !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !*We're done. It is not necessary to evaluate the Stiffness matrix associated to the vectorial problem. We solve separately numbofcomp scalar problems

     CONTAINS
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !*Compute the CRS structure starting from the nonzero entries
        !*Let's recall the structure of the nonzero entries
        !*
        !*Single element referred to the single DoF (row i)
        !*TYPE NzeScStMaDoF 
        !*   INTEGER :: numbofDoFs !*Number nonzerocolumns  
        !*   INTEGER, DIMENSION(:), ALLOCATABLE :: vecDoFs !*Global indices of those columns
        !*   REAL(DP), DIMENSION(:), ALLOCATABLE :: entries !*Nonzero entries
        !*END TYPE NzeScStMaDoF
        !*
        !*final vector of the nonzero entries row by row
        !*TYPE(NzeScStMaDoF), DIMENSION(:), ALLOCATABLE :: NzeSc 
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE FromNzeToCRS(CRSstructure,Nze)
           IMPLICIT NONE
           TYPE(CRS), INTENT(INOUT) :: CRSstructure
           TYPE(NzeScStMaDoF), DIMENSION(:), INTENT(IN) :: Nze

           !*Support variables
           INTEGER :: i,j,k !*Loops
           INTEGER :: counter !*To count the nonzero entries


           counter=0

           DO i=1,SIZE(Nze,DIM=1)
              !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !*SAFETY CHECK
              !*IF(Nze(i)%numbofDoFs==0) THEN
              !*   PRINT*, "Something terribly wrong, singular matrix. An empty row."
              !*   STOP
              !*END IF
              counter=counter+Nze(i)%numbofDoFs
              !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !*PRINT*, counter
              !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           END DO
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*PRINT*, counter
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !*Now we can allocate the structures of CRS
           IF(ALLOCATED(CRSstructure%values)) DEALLOCATE(CRSstructure%values)
           IF(ALLOCATED(CRSstructure%columns)) DEALLOCATE(CRSstructure%columns)
           IF(ALLOCATED(CRSstructure%rowstart)) DEALLOCATE(CRSstructure%rowstart)

           ALLOCATE(CRSstructure%values(1:counter))
           ALLOCATE(CRSstructure%columns(1:counter))
           ALLOCATE(CRSstructure%rowstart(1:SIZE(Nze,DIM=1)+1))
     
           !*Let's store the values row by row
           k=0 !*k will be the index in values and columns
           DO i=1,SIZE(Nze,DIM=1)
              DO j=1,Nze(i)%numbofDoFs
                 k=k+1
                 CRSstructure%values(k)=Nze(i)%entries(j)
                 CRSstructure%columns(k)=Nze(i)%vecDoFs(j)
              END DO
           END DO

          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, CRSstructure%values(138:)
          !*PRINT*, CRSstructure%columns(138:)
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !*Let's now store the starting row index
          CRSstructure%rowstart(1)=1
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, "Start row ", 1, "with", CRSstructure%rowstart(1)
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO i=1,SIZE(Nze,DIM=1)
             CRSstructure%rowstart(i+1)=CRSstructure%rowstart(i)+Nze(i)%numbofDoFs
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !*SAFETY CHECK
             !*PRINT*, "Adding ", Nze(i)%numbofDoFs 
             !*PRINT*, "Start row ", i+1, "with", CRSstructure%rowstart(i+1)
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          END DO

          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, "values", CRSstructure%values
          !*PRINT*, "columns", CRSstructure%columns
          !*PRINT*, "rowstart", CRSstructure%rowstart
          !*STOP          
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        END SUBROUTINE FromNzeToCRS


        
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !*CRS1to0 evaluates the CRS in 0-based system starting from the one in 1-based system
        !*We actually don't need it but it can be useful
        !*1-based system
        !*vectors indices (1:end)
        !*first column index=1
        !*first element in rowstart=1
        !*0-based system
        !*vectors indices (0:end-1)
        !*first column index=0
        !*first element in rowstart=0
        !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE CRS1to0(CRS0,CRS1)
           IMPLICIT NONE
           TYPE(CRS), INTENT(INOUT) :: CRS0
           TYPE(CRS), INTENT(IN) :: CRS1
           INTEGER :: i !*For the loop

           IF(ALLOCATED(CRS0%values)) DEALLOCATE(CRS0%values)
           IF(ALLOCATED(CRS0%columns)) DEALLOCATE(CRS0%columns)
           IF(ALLOCATED(CRS0%rowstart)) DEALLOCATE(CRS0%rowstart)

           ALLOCATE(CRS0%values(0:SIZE(CRS1%values,DIM=1)-1))
           ALLOCATE(CRS0%columns(0:SIZE(CRS1%columns,DIM=1)-1))
           ALLOCATE(CRS0%rowstart(0:SIZE(CRS1%rowstart,DIM=1)-1))

           !*The values have just to be copied with -1 index
           !*The columns have to be copied with -1 index but their values has to be decreased by 1
           DO i=1,SIZE(CRS1%values,DIM=1)
              CRS0%values(i-1)=CRS1%values(i)
              CRS0%columns(i-1)=CRS1%columns(i)-1
           END DO

           !*The startrows have to be copied with -1 index and their values has to be decreased by 1
           DO i=1,SIZE(CRS1%rowstart,DIM=1)
              CRS0%rowstart(i-1)=CRS1%rowstart(i)-1
           END DO
        END SUBROUTINE CRS1to0

        SUBROUTINE CleaningNze(Nze)
           TYPE(NzeScStMaDoF), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Nze
           INTEGER :: i
           
           IF (ALLOCATED(Nze)) THEN
              DO i=1,SIZE(Nze,DIM=1) 
                 IF (ALLOCATED(Nze(i)%vecDoFs)) DEALLOCATE(Nze(i)%vecDoFs)
                 IF (ALLOCATED(Nze(i)%entries)) DEALLOCATE(Nze(i)%entries)
              END DO
           END IF

        END SUBROUTINE CleaningNze
        
 END SUBROUTINE stiffnessmatrix


 !*--------------------------------------------------------------------------------------------
 !*global_Control_to_Cons 
 !*coefficients->values
 !*NB: It is impossible to make a check on the dimension of the output so be careful with matching before entering the function
 !*Up to now we had just a local transformation in model i.e. the function Control_to_Cons(u,e) RESULT (u_cons) taking as input the local coefficients in an element and giving as output the local values
 !*But often the projection to do is global and we make a loop on the elements recalling Control_to_Cons.
 !*I make now a global projection function that takes the mesh and the coefficients stored in a vector of PVar of dimension equal to the number of DoFs and performs the projection to get the values
 !*--------------------------------------------------------------------------------------------
 FUNCTION global_Control_to_Cons(Ucoeff,Mesh) RESULT (Uval) 
    IMPLICIT NONE
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: Ucoeff !*Coefficients of the solution in all the nodes. INPUT
    TYPE (maillage), INTENT(IN) :: Mesh !*INPUT
    TYPE(PVar), DIMENSION(SIZE(Ucoeff,DIM=1)) :: Uval !*Values of the solution in all the nodes. OUTPUT
    
    !*Support structures for the projection
    TYPE(PVar), DIMENSION(:), ALLOCATABLE :: Ucoeffloc,Uvalloc
    !*Ucoeffloc local coefficients extracted from Ucoeff
    !*Uvalloc local values got from Ucoeffloc thanks to Control_to_Cons

    !*Loop variables used in the loop on the elements and on the DoFs for the projection
    INTEGER :: indt=0 !*Index for the loop on the triangles
    INTEGER :: indi=0 !*Index for the loop on the DoFs i

    TYPE (element) :: e !*For the quick reference to the generic Mesh%e(indt) in the loops
    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*PRINT*, "Sizes", SIZE(Ucoeff,DIM=1), SIZE(Uval,DIM=1), Mesh%Ns 
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(SIZE(Ucoeff,DIM=1)/=Mesh%Ns) THEN
       PRINT*, "Error in global_Control_to_Cons. Mismatch of the dimensions"
       STOP
    END IF
    !*If it didn't stop we can continue

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*DO indi=1,Mesh%Ns
    !*   PRINT*, indi, Ucoeff(indi)%u
    !*END DO
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ALLOCATE(Ucoeffloc(Mesh%e(1)%nsommets),Uvalloc(Mesh%e(1)%nsommets))
    Ucoeffloc=0._DP
    Uvalloc=0._DP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*PRINT*, Ucoeffloc
    !*PRINT*, Uvalloc    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO indt=1,Mesh%Nt !*Loop on the elements
       e=Mesh%e(indt) !*Quick reference to the processed element
       Ucoeffloc=0._DP
       Uvalloc=0._DP
       DO indi=1,e%nsommets !*Loop on the DoFs of the processed element to extract the local coefficients 
          Ucoeffloc(indi)=Ucoeff(e%nu(indi))
       END DO !*End loop of the acquisition at the DoFs
       Uvalloc=Control_to_Cons(Ucoeffloc,e) !*Local values !*It is a function in the module model
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*PRINT*, "Triangle", indt
       !*PRINT*, "Local coefficients"
       !*DO indi=1,e%nsommets
       !*   PRINT*, Ucoeffloc(indi)
       !*END DO
       !*PRINT*, "Local values"
       !*DO indi=1,e%nsommets
       !*   PRINT*, Uvalloc(indi)
       !*END DO
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO indi=1,e%nsommets !*We put everything at its place
          Uval(e%nu(indi))=Uvalloc(indi)
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, indt, e%nu(indi), indi, Uvalloc(indi), Uval(e%nu(indi))
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO       
    END DO !*End loop on the elements

    !*End of the projection
    DEALLOCATE(Ucoeffloc,Uvalloc)
    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indi=1,Mesh%Ns
    !*   PRINT*,indi, Uval(indi)
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

 END FUNCTION global_Control_to_Cons


 !*--------------------------------------------------------------------------------------------
 !*global_Cons_to_Control 
 !*values->coefficients
 !*NB: It is impossible to make a check on the dimension of the output so be careful with matching before entering the function
 !*Up to now we had just a local transformation in model i.e. the function Cons_to_Control(u_cons,e) RESULT (u) taking as input the local values in an element and giving as output the local coefficients
 !*But often the projection to do is global and we make a loop on the elements recalling Cons_to_Control.
 !*I make now a global projection function that takes the mesh and the values stored in a vector of PVar of dimension equal to the number of DoFs and performs the projection to have the coefficients
 !*--------------------------------------------------------------------------------------------
 FUNCTION global_Cons_to_Control(Uval,Mesh) RESULT (Ucoeff) 
    IMPLICIT NONE
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: Uval !*Values of the solution in all the nodes. INPUT
    TYPE (maillage), INTENT(IN) :: Mesh !*INPUT
    TYPE(PVar), DIMENSION(SIZE(Uval,DIM=1)) :: Ucoeff !*Coefficients of the solution in all the nodes. OUTPUT
    
    !*Support structures for the projection
    TYPE(PVar), DIMENSION(:), ALLOCATABLE :: Uvalloc, Ucoeffloc
    !*Uvalloc local values extracted from Uval
    !*Ucoeffloc local coefficients got from Uvalloc through Cons_to_Control
     
    !*Loop variables used in the loop on the elements and on the DoFs for the projection
    INTEGER :: indt=0 !*Index for the loop on the elements
    INTEGER :: indi=0 !*Index for the loop on the DoFs i

    TYPE (element) :: e !*For the quick reference to the generic Mesh%e(indt) in the loops
    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*PRINT*, "Sizes", SIZE(Ucoeff,DIM=1), SIZE(Uval,DIM=1), Mesh%Ns 
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IF(SIZE(Uval,DIM=1)/=Mesh%Ns) THEN
       PRINT*, "Error in global_Cons_to_Control. Mismatch of the dimensions"
       STOP
    END IF
    !*If it didn't stop we can go ahead
    ALLOCATE(Uvalloc(Mesh%e(1)%nsommets),Ucoeffloc(Mesh%e(1)%nsommets))
    Uvalloc=0._DP
    Ucoeffloc=0._DP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*PRINT*, Uvalloc
    !*PRINT*, Ucoeffloc    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO indt=1,Mesh%Nt !*Loop on the elements
       e=Mesh%e(indt) !*Quick reference to the processed element
       Uvalloc=0._DP
       Ucoeffloc=0._DP
       DO indi=1,e%nsommets !*Loop on the DoFs of the processed element to extract the local values
          Uvalloc(indi)=Uval(e%nu(indi))
       END DO !*End loop of the acquisition at the DoFs
       Ucoeffloc=Cons_to_Control(Uvalloc,e) !*Local coefficients !*It is a function in the module model
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*PRINT*, "Triangle", indt
       !*PRINT*, "Local values"
       !*DO indi=1,e%nsommets
       !*   PRINT*, Uvalloc(indi)
       !*END DO
       !*PRINT*, "Local coefficients"
       !*DO indi=1,e%nsommets
       !*   PRINT*, Ucoeffloc(indi)
       !*END DO
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO indi=1,e%nsommets !*We put everything at its place
          Ucoeff(e%nu(indi))=Ucoeffloc(indi)
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, indt, e%nu(indi), indi, Ucoeffloc(indi), Ucoeff(e%nu(indi))
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO       
    END DO !*End loop on the elements

    !*End of the projection
    DEALLOCATE(Uvalloc,Ucoeffloc)

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indi=1,Mesh%Ns
    !*   PRINT*,indi, Ucoeff(indi)
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 END FUNCTION global_Cons_to_Control

 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !*It cleans VecProx (which is not needed anymore after the construction of the stiffness matrix)
 !*REMARK: VecProx is an allocatable vector of dimension equal to the number of DoFs
 !*Each element is
 !*TYPE Proximity !*Element of the vector referred to a single DoF i
 !*   INTEGER  :: numbofel !*number of elements containing i
 !*   INTEGER, DIMENSION(:), ALLOCATABLE :: vecEl !*Elements containing i
 !*   INTEGER  :: numbofDoFs !*number of DoFs in the elements containing i (without redoundances) 
 !*   INTEGER, DIMENSION(:), ALLOCATABLE :: vecDoFs !*Global indices of those DoFs
 !*END TYPE Proximity
 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE CleaningVecProx(VecProx)
    IMPLICIT NONE
    TYPE(Proximity), DIMENSION(:), ALLOCATABLE :: VecProx
    INTEGER :: indi

    IF(ALLOCATED(VecProx)) THEN !*Clean only if it is allocated
       DO indi=1,SIZE(VecProx,DIM=1) !*Loop on the rows
          IF(ALLOCATED(VecProx(indi)%vecEl)) DEALLOCATE(VecProx(indi)%vecEl)
          IF(ALLOCATED(VecProx(indi)%vecDoFs)) DEALLOCATE(VecProx(indi)%vecDoFs)
       END DO
    END IF

 END SUBROUTINE CleaningVecProx

 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !*It cleans a CRS structure (which is not needed anymore after the construction of the stiffness matrix)
 !*REMARK: a CRS structure is
 !*TYPE CRS !*Compressed row storage
 !*   REAL(DP), DIMENSION(:), ALLOCATABLE :: values !*Nonzero entries
 !*   INTEGER, DIMENSION(:), ALLOCATABLE :: columns !*columns of the nonzero entries
 !*   INTEGER, DIMENSION(:), ALLOCATABLE :: rowstart !*Index of elements in val that is at the start of the given row
 !*END TYPE
 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE CleaningCRS(CRSstr)
    IMPLICIT NONE
    TYPE(CRS) :: CRSstr

    IF(ALLOCATED(CRSstr%values)) DEALLOCATE(CRSstr%values)
    IF(ALLOCATED(CRSstr%columns)) DEALLOCATE(CRSstr%columns)
    IF(ALLOCATED(CRSstr%rowstart)) DEALLOCATE(CRSstr%rowstart)    
    
 END SUBROUTINE CleaningCRS

 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !*errorTest calculates the difference between my solution and the exact solution and its L^1,L^2 and L^\infty norm
 !*- In particular we start from the coefficients of our solution as an input FinalSolutionElliptic NB:They are already in the initial numeration -> CoeffNumSol=FinalSolutionElliptic
 !*- We calculate the values of the analytical exact solution -> ValExactSol
 !*- We pass to the coefficients of the analytical exact solution -> CoeffExactSol
 !*- We make the difference between the coefficients of our solution and the analytical one and we have the difference -> CoeffDiff
 !*We properly integrate it to get the L^1,L^2 norms
 !*The L^\infty norm is made directly on the coefficients CoeffDiff of the difference
 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE errorTest(error1,error2,errorinf,CoeffNumSol,Mesh)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(OUT) :: error1,error2,errorinf !*errors L^1,L^2,L^infty
    TYPE(Pvar), DIMENSION(:), INTENT(IN) :: CoeffNumSol !*Coefficients of our solution to the elliptic problem in the original numeration
    TYPE (maillage), INTENT(IN) :: Mesh 

    TYPE(PVar), DIMENSION(:), ALLOCATABLE :: ValExactSol !*Values of the exact solution
    TYPE(PVar), DIMENSION(:), ALLOCATABLE :: CoeffExactSol !*Coeff of the exact solution
    TYPE(PVar), DIMENSION(:), ALLOCATABLE :: CoeffDiff !*Coeff of the difference
    TYPE(PVar), DIMENSION(:), ALLOCATABLE :: ValDiff !*Values of the difference

    !*Loops variables used in the integration of the difference 
    INTEGER :: indt=0 !*Index for the loop on the triangles
    INTEGER :: indi=0 !*Index for the first loop on the DoFs i
    INTEGER :: indq=0 !*Index for the loop on the quadrature points

    REAL(DP), DIMENSION(:), ALLOCATABLE :: base !*Vector of the basis functions in a specific quadrature point for the integrations required for the calculation of the norms of the error
    
    TYPE (element) :: e !*For the quick reference to the generic Mesh%e(indt) in the loops
    REAL(DP) :: suppreal=0._DP !*Support real variable in the calculation of the values of the exact solution
    REAL(DP), DIMENSION(:), ALLOCATABLE :: eloc1!*local l1 error in an element
    REAL(DP), DIMENSION(:), ALLOCATABLE :: eloc2 !*local l2 error SQUARED in an element
    !*NB:
    !*error in norm L^1=sum eloc1
    !*error in norm L^2=sqrt(sum eloc2)
    REAL(DP), DIMENSION(:), ALLOCATABLE :: valueofthefunction !*Value of the function in the quadrature point=sum basis in the quadrature points * coefficients

    INTEGER :: component, numbofcomp


    ALLOCATE(ValExactSol(Mesh%Ns)) !*Allocation and initialization of the values of the exact solution
    ValExactSol=0._DP

    !*Calculation of the values of the exact solution
    DO indt=1,Mesh%Nt !*Loop on the elements
       e=Mesh%e(indt) !*Quick reference to the processed element
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*PRINT*, "Processed element: ", indt
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO indi=1,e%nsommets !*Loop on the DoFs of the element
          suppreal=0._DP
          suppreal=exactsolution(e%coor(indi))
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, indi, e%nu(indi), e%coor(indi), exactsolution(e%coor(indi)), suppreal
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ValExactSol(e%nu(indi))=suppreal !*PVar=real
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, e%nu(indi), e%coor(indi), ValExactSol(e%nu(indi))%u, suppreal
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
       END DO
    END DO

    !*Now that we have the values, let's pass to the coefficients
    ALLOCATE(CoeffExactSol(Mesh%Ns)) !*Allocation
    CoeffExactSol=0._DP !*Initialization
    CoeffExactSol=global_Cons_to_Control(ValExactSol,Mesh) !*Filling
    
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indi=1,Mesh%Ns
    !*   PRINT*, "DoF", indi
    !*   PRINT*, CoeffExactSol(indi)%u !*Coefficients
    !*   PRINT*, CoeffExactSol(indi)%u-ValExactSol(indi)%u !*For Lagrange they should be equal
    !*   PRINT*, CoeffExactSol(indi)%u-CoeffNumSol(indi)%u !*Error on the coefficients
    !*END DO
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DEALLOCATE(ValExactSol) !*Not needed anymore

    !*Now we can consider the coefficients of the error function
    ALLOCATE(CoeffDiff(Mesh%Ns)) !*Allocation
    CoeffDiff=0._DP !*Initialization
    !*Filling
    DO indi=1,Mesh%Ns
       CoeffDiff(indi)%u=CoeffExactSol(indi)%u-CoeffNumSol(indi)%u !*NB: without absolute value
    END DO
    DEALLOCATE(CoeffExactSol) !*Not needed anymore

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indi=1,Mesh%Ns
    !*   PRINT*, CoeffDiff(indi)%u
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
    !*We can now start the integration loop to get the norms of the error
    !*Initialization
    error1=0._DP !*L^1
    error2=0._DP !*L^2

    ALLOCATE(base(Mesh%e(1)%nsommets)) !*Basis functions in a specific quadrature point
    base=0._DP !*Initialization

    numbofcomp=SIZE(error1,DIM=1) !*Number of components
    ALLOCATE(eloc1(numbofcomp),eloc2(numbofcomp))

    eloc1=0._DP !*local l1 error over an element
    eloc2=0._DP !*local l2 error SQUARED over an element
    !*NB:
    !*error in norm L^1=sum eloc1
    !*error in norm L^2=sqrt(sum eloc2)

    ALLOCATE(valueofthefunction(numbofcomp)) !*Value of the function in the quadrature point=sum basis in the quadrature points * coefficients
    valueofthefunction=0._DP

    DO indt=1,Mesh%Nt !*Loop on the elements
        e=Mesh%e(indt) !*Quick reference
        !*For each element we initialize
        eloc1=0._DP 
        eloc2=0._DP
        DO indq=1,e%nquad !*Integration loop on the quadratures points
           valueofthefunction=0._DP
           base=0._DP
           !*We need the basis functions in the quadrature point to reconstruct the value in that point of the function to integrate 
           DO indi=1,e%nsommets !*Loop to evaluate the basis functions in that quadrature point
              base(indi)=e%base(indi,e%quad(:,indq))
              !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !*SAFETY CHECK
              !*PRINT*, indi, e%quad(:,indq), base(indi)
              !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           END DO
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*SAFETY CHECK
           !*PRINT*, SUM(base)
           !*STOP
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !*Reconstruction of the value in that point
           DO component=1,numbofcomp !*For each component
              DO indi=1,e%nsommets !*Loop on the basis functions
                 valueofthefunction(component)=valueofthefunction(component)+&
base(indi)*CoeffDiff(e%nu(indi))%u(component)
              END DO  
           END DO
           !*Now we pass to the absolute value
           DO component=1,numbofcomp
              valueofthefunction(component)=ABS(valueofthefunction(component))
           END DO
           !*We give the contribute of the quadrature point to the local integral
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*SAFETY CHECK
           !*PRINT*, indt, indq, "Old value of eloc1", eloc1
           !*PRINT*, indt, indq, "Old value of eloc2", eloc2
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           eloc1=eloc1+valueofthefunction*e%weight(indq)
           eloc2=eloc2+(valueofthefunction**2)*e%weight(indq)
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*SAFETY CHECK
           !*PRINT*, indt, indq, "New value of eloc1", eloc1
           !*PRINT*, indt, indq, "New value of eloc2", eloc2
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !*PRINT*, "value", valueofthefunction
           !*PRINT*, "square", valueofthefunction**2
           !*PRINT*, "numbofcomp", numbofcomp
           !*STOP
           !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        END DO !*End of the integration loop on the quadratures points
        eloc1=eloc1*e%volume !*Multiplication by the Jacobian
        eloc2=eloc2*e%volume !*Multiplication by the Jacobian
        error1=error1+eloc1 !*Adding the contribute of the element to the global error
        error2=error2+eloc2 !*Adding the contribute of the element to the global error
    END DO
    
    error2=SQRT(error2) !*sqrt on error 2 to make the L^2 norm
    !*PRINT*, error1
    !*PRINT*, error2

    ALLOCATE(ValDiff(Mesh%Ns)) !*Allocation
    ValDiff=0._DP !*Initialization
    ValDiff=global_Control_to_Cons(CoeffDiff,Mesh)
    DEALLOCATE(CoeffDiff)

    !*PRINT*, errorinf
    errorinf=0._DP
    !*PRINT*, errorinf

    !*Now let's compute the errorinf
    DO indi=1,Mesh%Ns
       DO component=1,numbofcomp
          errorinf(component)=MAX(errorinf(component),ABS(ValDiff(indi)%u(component)))
       END DO
    END DO



    CONTAINS

    FUNCTION exactsolution(x) RESULT(sol) 
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x !*Cartesian coordinates
       REAL(DP) :: sol
       REAL(DP) :: coeff, coeff1, coeff2

#if(1==0)       
       coeff=valueoftheknownterm/2._DP
       sol=coeff*(-x**2+1._DP) !*Correct not to put DP on 2, reference Chapman (but ALERT)
#endif
#if(1==0)       
       sol=-(-SIN(2._DP*PI*x)/(2._DP*PI)**2)
#endif       
#if(1==1)       
       coeff1=-1/2._DP*(EXP(1._DP)-EXP(-1._DP))
       coeff2=-1/2._DP*(EXP(1._DP)+EXP(-1._DP))
       sol=-(EXP(x)+coeff1*x+coeff2)
#endif
#if(1==0)       
       sol=-(x**3/6._DP-x/6._DP)
#endif       

    END FUNCTION exactsolution
    
 END SUBROUTINE errorTest

 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !*It performs the reconstruction at the DoFs of the gradient of the solution to the elliptic problem 
 !*- Compute the gradient in every DoF
 !*- In the DoFs shared by more elements make the average between the gradients (in that DoF) in the different elements
 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE gradient_reconstruction(GradAtDoFs,CoeffNumSol,VecProx,Mesh)
    IMPLICIT NONE
    TYPE(PVAR), DIMENSION(:,:), INTENT(OUT) :: GradAtDoFs !*Gradient reconstructed at the DoFs
    !*Vector of PVar with 2 indices
    !*First index->DoFs
    !*Second index-> spatial components
    !*Then since it is a PVar in %u(:) we have n_vars components for the components of the vectorial problem
    TYPE(PVAR), DIMENSION(:), INTENT(IN) :: CoeffNumSol !*Coefficients of the solution to the elliptic problem in the original numeration
    TYPE(Proximity), DIMENSION(:), INTENT(IN) :: VecProx !*VecProx is a vector of proximity elements. The definition of the generic element is recalled here
    !*TYPE Proximity !*Element of the vector referred to a single DoF i
    !*   INTEGER  :: numbofel !*number of elements containing i
    !*   INTEGER, DIMENSION(:), ALLOCATABLE :: vecEl !*Elements containing i
    !*   INTEGER  :: numbofDoFs !*number of DoFs in the elements containing i (without redoundances) 
    !*   INTEGER, DIMENSION(:), ALLOCATABLE :: vecDoFs !*Global indices of those DoFs
    !*END TYPE Proximity

    TYPE(maillage), INTENT(IN) :: Mesh

    TYPE(PVAR), DIMENSION(:), ALLOCATABLE :: CoeffLoc !*Support structure for the coefficients of the solution in a specific element to compute the gradients in the DoFs

    INTEGER :: indi=0 !*Index for the loop on the DoFs
    INTEGER :: indj=0 !*Index for the loop on the DoFs
    INTEGER :: indt=0 !*Index for the loop on the triangles
    INTEGER :: indd=0 !*Index for the loop on the dimension
    INTEGER :: indc=0 !*Index for the loop on the components

    TYPE (element) :: e !*For the quick reference to the generic Mesh%e(indt) in the loops

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*STRUCTURE FOR SAFETY CHECK
    TYPE(PVAR) :: dersupp !*Support structure for the derivative
    TYPE(PVAR), DIMENSION(:), ALLOCATABLE :: coeffsupp !*Support structure for the coeffcients
    INTEGER :: numbtest=102
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !*STRATEGY: We will first fill GradAtDoFs adding all the values of the gradient in the DoFs
    !*NB: If a DoF is shared by more elements we will meet it as many times as the number of elements containing it
    !*In the end we divide the value in each DoF by the number of elements containing the DoF thanks to VecProx(indi)%numbofel

    GradAtDoFs=0._DP !*Safe intialization of the output

    ALLOCATE(CoeffLoc(Mesh%e(1)%nsommets)) !*Allocation of the support structure for the local coefficients in an element to compute the gradient
    CoeffLoc=0._DP !*Safe initialization

    DO indt=1,Mesh%nt !*Loop on the elements
       e=Mesh%e(indt)
       CoeffLoc=0._DP !*Safe initialization
       DO indi=1,e%nsommets !*Loop on the DoFs to extract the local coefficients in the element
          CoeffLoc(indi)=CoeffNumSol(e%nu(indi))
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, CoeffLoc(indi)%u(1)-CoeffLoc(indi)%u(2)
          !*PRINT*, CoeffLoc(indi)%u(2)-CoeffLoc(indi)%u(3)
          !*PRINT*, CoeffLoc(indi)%u(3)-PotentialAtDoFsUtils(e%nu(indi))
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO
       !*Now I compute the gradient in each DoF and I add it in its position in GradAtDoFs
       !*RMK
       !*GradAtDoFs is a vector of PVar with 2 indices
       !*First index->DoFs
       !*Second index-> spatial components
       DO indi=1,e%nsommets
          DO indj=1,e%nsommets
             GradAtDoFs(e%nu(indi),:)=GradAtDoFs(e%nu(indi),:)+&
             e%gradient(indj,e%x(:,indi))*CoeffLoc(indj)
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !*SAFETY CHECK
             !*IF (indt==837 .AND. indi==1) THEN
             !*   PRINT*, indt 
             !*   PRINT*, "Basis function of barycentric coordinates", e%coorL(:,indi)
             !*   PRINT*, "Computed in", e%coorL(:,indj)
             !*   PRINT*, e%gradient(indi,e%coorL(:,indj))
             !*   PRINT*
             !*END IF
             !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
          END DO
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*ALLOCATE(coeffsupp(e%nsommets))
          !*dersupp=0._DP
          !*coeffsupp=0._DP
          !*DO indj=1,e%Nsommets
          !*   coeffsupp(indj)=CoeffLoc(indj)%u(2)
          !*END DO
          !*dersupp=e%eval_der(coeffsupp,e%coorL(:,indi))
          !*DEALLOCATE(coeffsupp)
          !*!*PRINT*, GradAtDoFs(e%nu(indi),1)%u(2), GradAtDoFs(e%nu(indi),2)%u(2)
          !*PRINT*, GradAtDoFs(e%nu(indi),1)%u(2)-dersupp(1)
          !*PRINT*, GradAtDoFs(e%nu(indi),2)%u(2)-dersupp(2)
          !*IF(indt==50) STOP
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
          !*SAFETY CHECK
          !*IF (e%nu(indi)==numbtest) THEN
          !*   PRINT*, e%coor(indi)
          !*   dersupp=0._DP
          !*   dersupp=e%eval_der(CoeffLoc,e%x(:,indi))
          !*   PRINT*, "If correct I added"
          !*   PRINT*, dersupp !*Computed here
          !*   PRINT*, "So the new values are"
          !*   PRINT*, GradAtDoFs(e%nu(indi),1) !*Computed before
          !*END IF
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO
    END DO

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*PRINT*, "Elements containing the DoF", VecProx(numbtest)%numbofel
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !*Now we need to average i.e. to divide every component (indd and indc) of a specific DoF indi by the number of elements containing it 
    DO indi=1,Mesh%Ns !*Loop on the DoFs
       DO indd=1,n_dim !*Loop on the dimensions
          DO indc=1,n_vars !*Loop on the component
             GradAtDoFs(indi,indd)%u(indc)=GradAtDoFs(indi,indd)%u(indc)/&
             & REAL(VecProx(indi)%numbofel,DP)
          END DO
       END DO
    END DO

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*DO indi=1,Mesh%Ns
    !*   !*PRINT*, GradPotentialAtDoFsUtils(indi,1)
    !*   PRINT*, GradAtDoFs(indi,1)%u-GradPotentialAtDoFsUtils(indi,1)
    !*END DO    
    !*STOP
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    DEALLOCATE(CoeffLoc)

 END SUBROUTINE gradient_reconstruction




END MODULE preprocessing
