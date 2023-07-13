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
MODULE geometry
  USE param2d
  USE init_bc
  USE PRECISION
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: geom

CONTAINS
  !*--------------------------------------------------
  !*geom builds the mesh
  !*--------------------------------------------------
  SUBROUTINE geom(Mesh, DATA)
    TYPE (maillage), INTENT(inout):: mesh
    TYPE(donnees),   INTENT(inout):: DATA
    REAL(dp), DIMENSION(:), ALLOCATABLE:: base

    INTEGER, PARAMETER,DIMENSION(14)::loc_ndofs=(/2,3,3,4,4,5,5,0,0,0,2,3,4,5/)
    !*loc_ndofs(itype)=number of DoFs in a single element for the chosen type
    !*For example 
    !*itype=5=B3 -> loc_ndofs(5)=4
    !*itype=12=P2 in Gauss-Lobatto loc_ndofs(12)=3
    REAL(dp):: dx, a, alpha
    TYPE(element):: e
    INTEGER:: nt, jt, itype, p1, p2, k, i, iq, l
    nt=DATA%nt; itype=DATA%itype
    !*Choice on the type of elements
    !*RMK:
    !*itype 
    !*1: P1, 
    !*2: B2, 
    !*3: P2, 
    !*4: P3, 
    !*5: B3
    !*11, 12, 13, 14: P1, P2, P3, P4 in Gauss-Lobatto
    SELECT CASE(DATA%itype)
    CASE(1,11) !*P1 (eventually in Gauss-Lobatto points)
       Mesh%nt=DATA%nt
       Mesh%ndofs=DATA%nt+1 !*2 DoFs per element so nt+1
    CASE(2,3,12) !*B2, P2 (eventually in Gauss-Lobatto points)
       Mesh%nt=DATA%nt 
       Mesh%ndofs=2*DATA%nt+1 !*3 DoFs per element so 2*nt+1 (essentially the ones of case(1) + nt internal)
    CASE(4,5,13) !*P3, B3 (P3 eventually in Gauss-Lobatto points)
       mesh%nt=DATA%nt
       Mesh%ndofs=nt+1+2*nt !*4 DoFs per element so 3*nt+1 (essentially the ones of case(1) + 2*nt internal)
    case(6,7,14)
       Mesh%nt=Data%nt
       Mesh%ndofs=nt+1+3*nt
    CASE default
       PRINT*, "Error: this element is not yet defined", DATA%itype
       STOP
    END SELECT

    CALL init_geom(DATA) !*the called subroutine is in init_bc
    !*Depending on the choice of the test init_geom sets
    !*DATA%domain_left = starting point (the left one) of the (one dimensional) domain
    !*DATA%Length =length of the domain      

    dx=DATA%Length/REAL(nt,DP) !*Length of the single element =length domain/numer of element
    a =DATA%domain_left !*Starting point (set in init_geom)

    ALLOCATE(Mesh%e(nt)) !*Allocate the vector of elements of the mesh

    DO jt=1, nt !*Loop on the elements
       !*Allocate the structures of the elements (type defined in elements_1d)
       ALLOCATE(Mesh%e(jt)%coor(loc_ndofs(itype)),Mesh%e(jt)%nu(loc_ndofs(itype)))
       !*loc_ndofs(itype)=number of DoFs in a single element for the chosen type
       !*For example 
       !*itype=5=B3 -> loc_ndofs(5)=4
       !*itype=12=P2 in Gauss-Lobatto loc_ndofs(12)=3
       !*So he's allocating coor and nu for each element with dimension equal to the number of DoFs in the element
       !*e%coor = Coordinates of the DoFs of the element. Clearly only an abscissa for each DoF in 1d
       !*NB: As in 2D the natural order is: at the beginning we have the boundary vertices (2 in  this case) and then the internal DoFs, all with increasing abscissa
       !*Example: left vertex, right vertex, internal from lowest to highest abscissa
       !*   *-----*-----* for quadratic, *---*---*---* for cubic elements
       !*   1     3     2                1   3   4   2
       !*e%nu=Global indices of the DoFs 
       !*The global numeration is
       !* vertices 1,2,...,N+1 (N number of cells)
       !* internal DoFs with increasing abscissa
       !*So for B3 cell K we have 
       !*k, k+1, N+2k, N+2k+1
       ALLOCATE(Mesh%e(jt)%x(2,loc_ndofs(itype)))
       !*e%x(:,:)=Barycentric coordinates of the DoFs
       !*Annoying notation. 2 indices x(:,:)
       !*Second one referred to the DoF of which we consider the barycentric coordinates 1:nsommets
       !*First one referred to the 2 barycentric coordinates 1:2
       !*x(:,1) -> barycentric coordinates of the first DoF (left vertex)
       !*BUT NB: 
       !*Let's focus on the first index
       !*1->right vertex
       !*2->left vertex
       !*So
       !*First vertex is x(:,1)=(0,1) which is, up to my opinion, a bit annoying because the order selected for the DoFs is with increasing abscissa. It is not wrong but a bit nonlogical
       !*x(:,2)=(1,0) !*Second (right) vertex
       !*x(:,3)=(0,3;0,6) !*Closer to the left vertex associated to (0,1)
       !*x(:,4)=(0,6;0,3) !*Closer to the right vertex associated to (1,0)
       !*ACTUALLY if you want a logical way to imagine it just imagine that with this notation, in the reference interval [0,1], we have
       !*x(1)->x because you would start from x even if it is associated to the second node
       !*x(2)->1-x associated to the first node

       mesh%e(jt)%itype=itype !*Element chosen
       Mesh%e(jt)%nvertex = 2 !*Number of vertices (2 in 1d)
       mesh%e(jt)%nsommets=loc_ndofs(itype) !*Number of DoFs !*loc_ndofs(itype)=number of DoFs in a single element for the chosen type


       !*AGAIN WE REMARK THAT:
       !*The global numeration is
       !* vertices 1,2,...,N+1 (N number of cells)
       !* internal DoFs with increasing abscissa
       !*So for B3 cell K we have 
       !*k, k+1, N+2k, N+2k+1
       !*Coordinates of the first DoF, left vertex
       mesh%e(jt)%coor(1)=(jt-1)*dx+a
       !*Coordinates of the second DoF, right vertex
       mesh%e(jt)%coor(2)=(jt  )*dx+a
       !*Global index of the first DoF, left vertex
       mesh%e(jt)%nu(1)=jt
       !*Global index of the second DoF, right vertex
       mesh%e(jt)%nu(2)=jt+1
       
       !*Now let's deal with the internal nodes
       SELECT CASE(itype) !*Checked, correct
       !*NB: Internal nodes present from thrird order on
       CASE(2,3,12) !*B2, P2 (eventually in Gauss-Lobatto points)
          mesh%e(jt)%coor(3)=mesh%e(jt)%coor(1)+dx/2._DP 
          mesh%e(jt)%nu(3)=nt+1+jt 
       CASE(4,5) !*P3, B3
          mesh%e(jt)%coor(3)=mesh%e(jt)%coor(1)+dx/3._DP 
          mesh%e(jt)%coor(4)=mesh%e(jt)%coor(1)+2._DP*dx/3._DP 
          mesh%e(jt)%nu(3)=nt+1+2*(jt-1)+1
          mesh%e(jt)%nu(4)=nt+1+2*(jt-1)+2
       CASE(13) !*P3 in Gauss-Lobatto points
          alpha= 0.5_DP-SQRT(5._DP)/10._DP
          mesh%e(jt)%coor(3)=mesh%e(jt)%coor(1)+alpha*dx
          mesh%e(jt)%coor(4)=mesh%e(jt)%coor(1)+(1._DP-alpha)*dx
          mesh%e(jt)%nu(3)=nt+1+2*(jt-1)+1
          mesh%e(jt)%nu(4)=nt+1+2*(jt-1)+2
       case(6,7)
          mesh%e(jt)%coor(3)=mesh%e(jt)%coor(1)+dx/4._DP
          mesh%e(jt)%coor(4)=mesh%e(jt)%coor(1)+dx/2._DP
          mesh%e(jt)%coor(5)=mesh%e(jt)%coor(1)+dx*0.75_DP
          mesh%e(jt)%nu(3)=nt+1+3*(jt-1)+1
          mesh%e(jt)%nu(4)=nt+1+3*(jt-1)+2
          mesh%e(jt)%nu(5)=nt+1+3*(jt-1)+3
       case(14)
          alpha = 0.5_DP-SQRT(21._DP)/14._DP
          mesh%e(jt)%coor(3)=mesh%e(jt)%coor(1)+alpha*dx
          mesh%e(jt)%coor(4)=mesh%e(jt)%coor(1)+dx/2._DP
          mesh%e(jt)%coor(5)=mesh%e(jt)%coor(1)+(1._DP-alpha)*dx
          mesh%e(jt)%nu(3)=nt+1+3*(jt-1)+1
          mesh%e(jt)%nu(4)=nt+1+3*(jt-1)+2
          mesh%e(jt)%nu(5)=nt+1+3*(jt-1)+3

       END SELECT

       !*He's filling some fields of the structure element with internal procedures, once the coordinates have been set (have a look at elements_1D)
       mesh%e(jt)%volume=mesh%e(jt)%aire() !*Length
       ALLOCATE(mesh%e(jt)%n(2))
       mesh%e(jt)%n     =mesh%e(jt)%normale() !*Inward normals
       CALL mesh%e(jt)%quadrature() !*Quadratures
       ALLOCATE(Mesh%e(jt)%base0(mesh%e(jt)%nsommets,mesh%e(jt)%nsommets))
       ALLOCATE(Mesh%e(jt)%base1(mesh%e(jt)%nsommets,mesh%e(jt)%nsommets))
       CALL mesh%e(jt)%base_ref() !*base_ref_element computes the values of the basis functions in the DoFs which are stored in e%base0(:,:)
       !*First index DoF
       !*Second index basis function
       !*e%base0=[ phi_1(x_1) phi_2(x_1) ... phi_N(x_1)
       !*          phi_1(x_2) phi_2(x_2) ... phi_N(x_2)
       !*              .           .              .
       !*              .           .              .
       !*              .           .              .
       !*          phi_1(x_N) phi_2(x_N) ... phi_N(x_n) ]
       !*
       !*NB: e%base0*vectorofthecoefficient=vectorofthevalues
       !*vectorofthecoefficients=inv(e%base0)*vectorofthevalues=e%base1*vectorofthevalues
       !*It also computes the inverse of this matrix stored in e%base1(:,:) 

       call mesh%e(jt)%eval_coeff() !*eval_coeff_element evaluates some coefficients needed for the DeC
       !*eval_coeff_element evaluates some coefficients needed for the DeC
       !*e%mass(i,j)-> local mass matrix int_K phi_i phi_j
       !*e%coeff(i,j) -> another matrix int_K phi_i  d_x phi_j !*NB:DERIVATIVE ON THE SECOND INDEX


    ENDDO

    !*NB: before starting with the edges we remark that they are points in 1d

    ! connectivite: for each edge, we give the element before and after
    IF (DATA%periodicBC) THEN !*IF MESH WITH PERIODIC BC
        Mesh%nsegmt=nt ! periodicite: the last edge is the first one, one should not count it twice
        !*Basically he doesn't count the last edge twice, the last edge i.e. nt+1 coincides with the first one

        ALLOCATE(Mesh%edge(Mesh%nsegmt)) !*Allocate the vector of aretes edge 

        DO jt=2,Mesh%nsegmt !*We start from the second edge so we count just the internal ones
           Mesh%edge(jt)%jt1=jt-1 !*Left element containing the edge
           Mesh%edge(jt)%jt2=jt !*Right element containing the edge
	        Mesh%edge(jt)%bord=.FALSE. !*Not on the boundary
	     ENDDO

        !*NOw let's deal with the "boundary" node i.e. the 1=nt+1
        ! for periodic BCs
        ! here I take into account the periodicity of the mesh
        Mesh%edge(1)%jt2= 1 !*The right element containing it is the first element 
        Mesh%edge(1)%jt1=Mesh%nt  !*The left element containing it is the last element
        Mesh%edge(1)%bord=.FALSE. !*Not on the boundary
        !*With this trick the edge is on the boundary and it has two well-defined neighbours
    
    
    ELSE !*IF MESH WITH NO PERIODIC BC

        !*We have to count separately first and last edge
        Mesh%nsegmt=nt+1 ! also first and last one

        ALLOCATE(Mesh%edge(Mesh%nsegmt))
        DO jt=2,Mesh%nsegmt-1 !*Loop on the internal edges
           Mesh%edge(jt)%jt1=jt-1 !*Left element containing the edge
           Mesh%edge(jt)%jt2=jt !*Right element containing the edge
           Mesh%edge(jt)%bord=.false. !*Not on the boundary
        ENDDO

        !*As a convenction, for what concerns the boundary edges which have only one element containing them we set as jt1 and jt2 that only neighbour
        !*Obviously we set the logical bord = .TRUE. because they are on the boudary

        Mesh%edge(1)%jt2= 1
        Mesh%edge(1)%jt1= 1
        Mesh%edge(1)%bord=.true.

        Mesh%edge(Mesh%nsegmt)%jt2=Mesh%nt !Mesh%nt !1
        Mesh%edge(Mesh%nsegmt)%jt1=Mesh%nt
        Mesh%edge(Mesh%nsegmt)%bord=.true.
  
    END IF


    !*We will now deal with the dual cells neeeded for the DeC

    ALLOCATE(Mesh%aires(Mesh%ndofs)) !*Vector of the inverses of the measures of the dual cells associated to the DoFs needed for the DeC
    !*RMK: Despite the misleading name it contains the inverse of the measures of the dual cells and not the measures
    !*We will compute the measures of the dual cells and then make the inverse in the end
    Mesh%aires=0.0_dp !*Initialization
      
    DO jt=1, nt !*Loop on the elements
       e=Mesh%e(jt)
#if (0==1)
	    ! VERSION GOOD ONLY FOR B^k
       !*C^K_i=|K|/number of DoFs (clearly to sum over the K containing i but this is embedded in the loop)
       Mesh%aires(e%nu)=Mesh%aires(e%nu)+e%volume/REAL(e%nsommets,DP)
#else
		 !!!!! MORE GENERAL VERSION
       !*C^K_i=\int_K phi_i (clearly to sum over the K containing i but this is embedded in the loop)
		 ALLOCATE( base(e%nsommets))
       DO iq=1, e%nquad !*Loop on the quadrature points
		    DO l=1, e%nsommets !*Loop on the DoFs/basis functions
		       base(l  )=e%base    (l, e%quad(:,iq)) !*All the basis functions in the processed quadrature point
		    ENDDO
		    DO i=1, e%nsommets !*Loop on the DoFs
		       Mesh%aires(e%nu(i)) = Mesh%aires(e%nu(i)) + (base(i))*e%weight(iq)* e%volume !*For each DoF we integrate over K the basis function
             !*NB: In Mesh%aires we have the global index, in base the local one clearly
		    ENDDO ! i
       ENDDO

		 DEALLOCATE(base)
#endif			
    ENDDO

    !*NB: In case of periodic BCs we are missing sth
    !*We only adjusted the edges in case of periodic BCs 
    !*To adjust the dual cells we must identify the first left vertx (local DoF 1) of the first element with the right vertex (local DoF 2) of the last cell

    !*SO...
    IF (DATA%periodicBC) THEN !*In case of periodic BCs
      !!ONLY IF PERIODIC BOUNDARY
      p1=Mesh%e(1)%nu(1) !*p1 global index of the left vertex of the first element
      p2=Mesh%e(Mesh%nt)%nu(2) !*p2 global index of the right vertex of the last element
      mesh%aires(p1)=Mesh%aires(p1)+Mesh%aires(p2) !*SUM of the two contribution in the dual cell of p1
      Mesh%aires(p2)=Mesh%aires(p1) !*Copying it in p2
		END IF

    ! this to avoid further divisions
    Mesh%aires=1.0_dp/Mesh%aires !*Inverse to avoid every time the division :)
    !*Even if I have to say that another name could have been chosen 

  END SUBROUTINE geom

END MODULE geometry
