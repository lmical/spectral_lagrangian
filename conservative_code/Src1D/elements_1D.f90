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

!*----------------------------------------------------------
!*Basically a segment in 1d
!*----------------------------------------------------------
MODULE element_class
  ! In this module are listed all the tools linked to the basis functions and quadrature formulas
  ! Structure of the module:
  ! ELEMENT CLASS
  ! - aire: area of the element
  ! - normale: normal of the element
  ! - base: basis function -- here we do have Lagrangian and Bernstein polynomials
  ! - eval_function: evaluation of the solution via basis functions SUM(base(:)*(u(:)%u(l))
  ! -  eval_der: evaluation of the derivative SUM( (u(:)%u(l) *grad(1,:)) 
  ! -  eval_der2: evaluation of the second derivative SUM( (u(:)%u(l) *grad2(1,:))
  ! - der_sec: alternative to eval_der2
  ! - gradient: corresponds to the gradient of the basis function
  ! - gradient2: corresponds to the second order gradient of the basis function
  ! - quadrature:  collects all the points and weights for each typology of element
  ! - base_ref:  computes the values of the basis functions at the physical dofs
  ! - bary: defines the barycentric coordinates of the points
  ! - clean: cleans the memory of the features defined within the ELEMENT

  USE algebra
  USE variable_def
  USE overloading
  USE PRECISION
  IMPLICIT NONE


  TYPE, PUBLIC:: element
     INTEGER:: type_flux=-10 !*Set to be (later) equal to the scheme number
     INTEGER:: diag=-1, diag2=-1 !*?????????? !*Used fo the mood procedure
     INTEGER:: nsommets, itype, nvertex ! nbre de dofs, type element:
     !1->P1
     !2->B2
     !3->P2
     !4->P3
     !5->B3
     !6->B4
     !7->P4
     !11->PGL1
     !12->PGL2
     !13->PGL3
     !14->PGL4
     !*e%nsommets DoFs per element !*4 for B3
     !*e%itype type of elements (from the input) !*5 B3 (just because of our definition)
     !*e%nvertex number of vertices !*Clearly always 2 in 1dimension

     REAL(dp),  DIMENSION(:), POINTER :: coor =>NULL() ! 
     !   *-----*-----* for quadratic, *---*---*---* for cubic elements
     !   1     3     2                1   3   4   2
     !*Coordinates of the DoFs of the element. Clearly only an abscissa for each DoF in 1d
     !*NB: As in 2D the natural order is: at the beginning we have the boundary vertices (2 in  this case) and then the internal DoFs, all with increasing abscissa
     !*Example: left vertex, right vertex, internal from lowest to highest abscissa

     REAL(dp), DIMENSION(:,:), POINTER:: x=> NULL() !*Barycentric coordinates of the DoFs
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
     INTEGER, DIMENSION(:), POINTER   :: nu =>NULL() ! local connectivity table, see above for location !*Global indices of the DoFs 
     !*The global numeration is
     !* vertices 1,2,...,N+1 (N number of cells)
     !* internal DoFs with increasing abscissa
     !*So for B3 cell K we have 
     !*k, k+1, N+2k, N+2k+1
     ! For Bezier, this corresponds to the Greville points
     REAL(dp)                          :: volume =0.0_dp  ! volume !*Area of the cell (length in 1d)
     REAL(dp),  DIMENSION(:), POINTER   :: n =>NULL()     ! INTERNAL normals !*I'd say internal since it is (1,-1) and, RMK, the first vertex is the left, the second is the right one??????????
     !*INTEGER                        :: log    ! logic element  : this for boundary conditions, if needed !*NOT NEEDED
!!!!   quadrature de surface
     REAL(dp),   DIMENSION(:,:),POINTER :: quad =>NULL()  ! point de quadrature !*Quadrature points in barycentric coordinates
     !*Again the we have 2 indices e%quad(:,:)
     !*Second index referred to the quadrature point
     !*First index referred to the two barycentric coordinates
     !*And also in this case (I guess) if we focus on the first index
     !*1->right vertex
     !*0->left vertex
     REAL(dp),   DIMENSION(:)  ,POINTER :: weight=>NULL() ! poids !*Quadrature weights
     INTEGER                            :: nquad=0  ! nbre de points de quadrature !*Number of quadrature points
     REAL(dp),DIMENSION(:,:),POINTER:: base0=>NULL(),base1=>NULL() 
     !*base0 matrix of the values of the basis functions at the physical DoFs
     !*base1 inverse of base0
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
     !*NB: e%base0*vectorofthecoefficient=vectorofthevalues
     !*vectorofthecoefficients=inv(e%base0)*vectorofthevalues=e%base1*vectorofthevalues
     INTEGER, DIMENSION(:), POINTER :: dof2ind !*e%dof2ind(:) local indices of the nodes with increasing abscissa
     !*We have in the natural numeration
     !*At first the indices
     !*Then the internal DoF with increasing abscissa
     !*So e%dof2ind(:) is for example
     !*P1-> 1, 2
     !*B2-> 1, 3, 2
     !*B3,P3-> 1, 3, 4, 2

!!! int nabla phi_sigma * phi_sigma'
     REAL(dp), DIMENSION(:,:), POINTER:: coeff=> NULL(), mass=> NULL()!*, coeff_b=> NULL() !*NOT USED, not initialized, they can give segmentation fault in fact
     !*mass and coeff are the matrices phi*phi and phi*gradphi
     !*eval_coeff_element evaluates some coefficients needed for the DeC
     !*e%mass(i,j)-> local mass matrix int_K phi_i phi_j
     !*e%coeff(i,j) -> another matrix int_K phi_i  d_x phi_j !*NB:DERIVATIVE ON THE SECOND INDEX

   CONTAINS
     PRIVATE
     PROCEDURE, PUBLIC:: aire=>aire_element !*Area of the element, length in 1D
     PROCEDURE, PUBLIC:: gradient=>gradient_element !*gradient_element(k,x) computes the first derivative of the basis function k in the point of baricentric coordinates x
     !*NB: The first derivative is already in the "real" element, not in the reference one, because of the final division by e%volume
     PROCEDURE, PUBLIC:: gradient2=>gradient2_element !*gradient2_element(k,x) computes the SECOND derivative of the basis function k in the point of baricentric coordinates x
     !*NB: The second derivative is already in the "real" element, not in the reference one, because of the final division by e%volume**2
     PROCEDURE, PUBLIC:: average => average_element !*average(u) takes in input the coefficients of the solution (or a general function) in the DoFs and computes the average in the cell
     !*u->coefficients in the different nodes 
     PROCEDURE, PUBLIC:: base=>base_element !*base_element(k,x) !*It calculates a specific basis function k in the point of barycentric coordinates x
     !*k->basis function 
     !*x->barycentric coordinates
     !*RKM: basis functions/nodes are numerated in the following way
     !*1 left vertex
     !*2 right vertex
     !*other internal nodes with increasing abscissa
     !*x->barycentric coordinate REAL(DP), DIMENSION(2)
     !*RMK:
     !*Barycentric coordinates of the left vertex (0,1)
     !*Barycentric coordinates of the right vertex (1,0)
     !*So the second barycentric coordinate is associated to the left/first vertex, the first one is associated to the right/second vertex 
     !*RMK: To visualize it in a logical way, in the reference interval
     !*x(1)=x even if associated to the second vertex
     !*x(2)=1-x
     PROCEDURE, PUBLIC:: normale=>normale_element !*Inward normals
     PROCEDURE, PUBLIC:: quadrature=>quadrature_element !*quadrature_element fills the weights and the nodes of the quadrature formulas, it is a subroutine in fact
    !*Here I trust and do not check, to be checked and tested eventually
     PROCEDURE, PUBLIC:: base_ref=>base_ref_element !*base_ref_element computes the values of the basis functions in the DoFs which are stored in e%base0(:,:)
     !*First index DoF
     !*Second index basis function
     !*It also computes the inverse of this matrix stored in e%base1(:,:) 


     PROCEDURE, PUBLIC:: eval_func=>eval_func_element !*eval_func(u,y) 
     !*eval_func_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the solution in y 
     !*u->coefficients in the different nodes 
     !*y->barycentric coordinates of the point 
     PROCEDURE, PUBLIC:: eval_der=>eval_der_element !*eval_der(u,y) 
     !*eval_der_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the derivative of the solution in y 
     !*NB: We specify the barycentric coordinates but the output gradient is in the "real element"
     !*We keep into account the passage from the reference element to the real one dividing the gradient in the reference element by the length, inside the function gradient_element (grad=fx/e%volume)
     !*u->coefficients in the different nodes 
     !*y->barycentric coordinates of the point 
     PROCEDURE, PUBLIC:: eval_der2=>eval_der2_element !*eval_der2_element(u,y)
     !*eval_der2_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the SECOND derivative of the solution in y 
     !*NB: We specify the barycentric coordinates but the output gradient is in the "real element"
     !*We keep into account the passage from the reference element to the real one dividing by the length SQUARED the gradient in the reference element inside gradient2_element (grad2=fxx/(e%volume**2))

     PROCEDURE, PUBLIC:: eval_coeff=>eval_coeff_element  !*eval_coeff_element evaluates some coefficients needed for the DeC
     !*eval_coeff_element evaluates some coefficients needed for the DeC
     !*e%mass(i,j)-> local mass matrix int_K phi_i phi_j
     !*e%coeff(i,j) -> another matrix int_K phi_i  d_x phi_j !*NB:DERIVATIVE ON THE SECOND INDEX


     !*NB: There is another procedure bary_element(e) which has not been made public and in fact is recalled here inside e%base_ref_element and not outside. 
     !*It is used to fill e%x(:,:) barycentric coordinates of the DoFs of the element and e%dof2ind(:) local indices of the nodes with increasing abscissa
     !*RMK: e%x(:,:) has two indices
     !*Second one referred to the DoF of which we consider the barycentric coordinates 1:nsommets
     !*First one referred to the 2 barycentric coordinates 1:2
     !*
     !*For what concerns e%dof2ind(:) local indices of the nodes with increasing abscissa
     !*We have in the natural numeration
     !*At first the indices
     !*Then the internal DoF with increasing abscissa
     !*So e%dof2ind(:) is for example
     !*P1-> 1, 2
     !*B2-> 1, 3, 2
     !*B3,P3-> 1, 3, 4, 2
  
     FINAL:: clean !*It cleans the element
  END TYPE element

CONTAINS
  !*--------------------
  !*Length in 1D
  !*--------------------
  REAL(dp) FUNCTION aire_element(e) ! area of a element 
    CLASS(element), INTENT(in):: e
    !*REAL(dp), DIMENSION(1):: a!*,b
    REAL(dp) :: a
    a= e%coor(2)-e%coor(1)
    aire_element=ABS ( a )
  END FUNCTION aire_element

  !*--------------------------------------------
  !*Inward normals
  !*e%normale(1)=1 !*Left vertex, positive .->
  !*e%normale(2)=-1 !*Right vertex, negative <-.
  !*--------------------------------------------
  FUNCTION normale_element(e) RESULT(n) ! inward normals
    CLASS(element), INTENT(in)::e
    REAL(dp), DIMENSION(2):: n
    !*INTEGER:: l !*USELESS
    !*DO l=1,2 !*USELESS
       n(1)= 1.0_dp !* at e%nu(1), i.e. at x_{i}
       n(2)=-1.0_dp !* at e%nu(2), i.e. at x_{i+1}
    !*ENDDO
  END FUNCTION normale_element


  !*----------------------------------------------------
  !*Basis functions in the barycentric coordinates
  !*It calculates a specific basis function k in the point of barycentric coordinates x
  !*base_element(k,x)
  !*k->basis function 
  !*x->barycentric coordinates
  !*RKM: they are numerated in the following way
  !*1 left vertex
  !*2 right vertex
  !*other internal nodes with increasing abscissa
  !*x->barycentric coordinate REAL(DP), DIMENSION(2)
  !*RMK:
  !*Barycentric coordinates of the left vertex (0,1)
  !*Barycentric coordinates of the right vertex (1,0)
  !*So the second barycentric coordinate is associated to the left/first vertex, the first one is associated to the right/second vertex 
  !*RMK: To visualize it in a logical way in the reference interval
  !*x(1)=x even if associated to the second vertex
  !*x(2)=1-x
  !*-----------------------------------------------------
  REAL(dp) FUNCTION base_element(e,k,x) ! basis functions
    CHARACTER(Len=*), PARAMETER :: mod_name="base_element"
    CLASS(element), INTENT(in):: e
    INTEGER, INTENT(in):: k ! index of basis function
    REAL(dp), DIMENSION(2), INTENT(in):: x ! barycentric coordinate
    REAL(dp) :: alpha, beta, s

    SELECT CASE(e%itype) !*Choice on the type of the element
    CASE(1,11) ! P1
       SELECT CASE(k) !*Choice on k, on the basis function
       CASE(1) !*First basis function, so the one associated to the left vertex BUT NB: the left vertex is associated to the second component of the barycentric coordinates so we have the following
          base_element=x(2) !*When we have x(2)=1 [i.e. x=(0,1)] we are talking about the left vertex and base_element(1,x)=base_element(1,(0,1))=1 while base_element(1,(1,0))=0
          !*x(2)->1-x associated to the first vertex
       CASE(2) !*k=2, right vertex, x(1) because the first component is associated to the right vertex
          base_element=x(1) 
          !*x(1)->x associated to the second node
       CASE default
          PRINT*,  mod_name

          PRINT*, "P1, numero base ", k
          STOP
       END SELECT

    CASE(2) ! B2 Bezier
       SELECT CASE(k) !*Choice on k
       CASE(1) !*Again we remark that the second barycentric coordinate is associated to the first node, the left vertex so base_element(1,x)=x(2)*x(2)
          base_element=x(2)*x(2)!(1.0-x(1))*(1.0-x(1)) (alternatively)
       CASE(2) !*Right vertex (so x(1))
          base_element=x(1)*x(1)
       CASE(3) !*Only one internal DoF
          base_element=2.0_dp*x(1)*x(2)!(1.0-x(1))
       CASE default
          PRINT*,  mod_name

          PRINT*, "B2, numero base ", k
          STOP
       END SELECT

    CASE(3,12)! P2
       SELECT CASE(k) !*Choice on k
       !*Here it will be even more immediate to understand the relations between nodes (basis funcitons) and barycentric coordinates
       CASE(1) !*Associated to the left vertex (x(2))
          base_element=x(2)*(x(2)-x(1) )!(1.0-x(1))*(1.00-2.0*x(1))
          !*RMK: 
          !*x(1)=x -> second node
          !*x(2)=1-x -> first node
          !*(0,1) first vertex -> base_element=1
          !*(1,0) other extremum ->base_element=0
       CASE(2)
          base_element=x(1)*(x(1)-x(2) )!x(1)*(2.0*x(1)-1.00)
          !*Specular w.r.t. the previous case
       CASE(3)
          base_element=4._dp*x(1)*x(2) !4.00*x(1)*(1.00-x(1))
          !*0 in (0,1) and (1,0)
          !*1 in (0.5,0.5)
       CASE default
          PRINT*,  mod_name

          PRINT*, "P2, numero base ", k
          STOP
       END SELECT

    CASE(4) ! P3
       SELECT CASE(k)
       CASE(1)
          base_element = -0.5_dp*(3._dp*x(1)-1._dp)*(3._dp*x(1)-2._dp)*(x(1)-1._dp)
          !*(0,1) -> -0.5*(-1)*(-2)*(-1)=1
          !*(1,0) -> -0.5*(3-1)*(3-2)*(1-1)=0
          !*(1/3,2/3) and (2/3,1/3) -> 0
       CASE(2)
          base_element = 0.5_dp*x(1)*(3._dp*x(1)-1._dp)*(3._dp*x(1)-2._dp)
          !*(1,0) -> 0.5*1*2*1
          !*(0,1) -> 0
          !*(1/3,2/3) and (2/3,1/3) -> 0
       CASE(3)
          base_element = 1.5_dp*x(1)*(3._dp*x(1)-2._dp)*(3._dp*x(1)-3._dp)
       CASE(4)
          base_element = -1.5_dp*x(1)*(3._dp*x(1)-1._dp)*(3._dp*x(1)-3._dp)
       CASE default
          PRINT*,  mod_name

          PRINT*, "P3, numero base ", k
          STOP
       END SELECT

    CASE(5) ! B3
       SELECT CASE(k)
       CASE(1) !*Left vertex -> x(2)
          base_element =x(2)**3
       CASE(2) !*Right vertex -> x(1)
          base_element = x(1)**3
       CASE(3)
          base_element = 3.0_dp*x(1)*x(2)*x(2)
       CASE(4)
          base_element = 3.0_dp*x(1)*x(1)*x(2)
       CASE default
          PRINT*,  mod_name

          PRINT*, "P3, numero base ", k
          STOP
       END SELECT
    CASE(6) ! B4
       SELECT CASE(k)
       CASE(1)
          base_element = x(2)**4! (1.-x(1))**3
       CASE(2)
          base_element = x(1)**4
       CASE(3)
          base_element = 4.0_dp*x(1)*x(2)*x(2)*x(2)!3.*x(1)*( (1-x(1))**2)
       CASE(4)
          base_element = 6.0_dp*x(1)*x(1)*x(2)*x(2)!3.*x(1)*x(1)*((1-x(1)))
       CASE(5)
          base_element = 4.0_dp*x(1)*x(1)*x(1)*x(2)!3.*x(1)*x(1)*((1-x(1)))
       CASE default
          PRINT*, "B4, numero base ", k
          STOP
       END SELECT
    CASE(7) ! P4 
       s= x(1)
       SELECT CASE(k)
       CASE(1)
          base_element = (32._dp/3._dp*(s - 1._dp)*(s - 1._dp/2._dp)*(s - 1._dp/4._dp)*(s - 3._dp/4._dp))
       CASE(2)
          base_element = (32._dp*s*(s - 1._dp/2._dp)*(s - 1._dp/4._dp)*(s - 3._dp/4._dp))/3._dp
       CASE(3)
          base_element = -(128._dp*s*(s - 1._dp)*(s - 1._dp/2._dp)*(s - 3._dp/4._dp))/3._dp
       CASE(4)
          base_element = 64._dp*s*(s - 1._dp)*(s - 1._dp/4._dp)*(s - 3._dp/4._dp)
       CASE(5)
          base_element = -(128._dp*s*(s - 1._dp)*(s - 1._dp/2._dp)*(s - 1._dp/4._dp))/3._dp
    CASE default
          PRINT*, "P4, numero base ", k
          STOP
       END SELECT

    CASE(13) ! P3 Gauss Lobatto
			 alpha=0.5_dp-SQRT(5._dp)/10._dp
			 beta = 1._dp -alpha
       SELECT CASE(k)
       CASE(1)
          base_element = (-x(1)+alpha)*(x(1)-beta)*(x(1)-1.0_dp)/(alpha*beta)
       CASE(2)
          base_element = x(1)*(x(1)-alpha)*(x(1)-beta)/alpha/beta
       CASE(3)
          base_element = x(1)*(x(1)-beta)*(x(1)-1.0_dp)/alpha/(-2._dp*alpha+1._dp)/beta
       CASE(4)
          base_element = x(1)*(x(1)-alpha)*(x(1)-1.0_dp)/alpha/(2._dp*alpha-1._dp)/beta
       CASE default
          PRINT*, "P3 Gauss Lobatto, numero base ", k
          STOP
       END SELECT

    CASE(14) ! P4 Gauss Lobatto
			 alpha=0.5_dp-SQRT(21._dp)/14._dp
			 beta = 1._dp -alpha
       SELECT CASE(k)
       CASE(1)
          base_element = (2._dp*(alpha-x(1))*(x(1)-1._dp)*(x(1)-0.5_dp)*(x(1)-beta))/(-alpha*beta)
       CASE(2)
          base_element = (2._dp*x(1)*(-alpha + x(1))*(x(1) - 0.5_dp)*(x(1) - beta))/(alpha*beta)
       CASE(3)
          base_element = (x(1)*(x(1)-1._dp)*(x(1)-0.5_dp)*(x(1)-beta))/(-alpha*(2._dp*alpha-1)*beta*(alpha-0.5_dp))
       CASE(4)
          base_element = -(4._dp*x(1)*(alpha -x(1))*(x(1)-1._dp)*(x(1)-beta))/(alpha-0.5_dp)**2
       CASE(5)
          base_element =-(x(1)*(alpha-x(1))*(x(1)-1._dp)*(x(1)-0.5_dp))/(-alpha*(2._dp*alpha-1._dp)*beta*(alpha-0.5_dp))
       CASE default
          PRINT*, "P4 Gauss Lobatto, numero base ", k
          STOP
       END SELECT

    CASE default
       PRINT*, "Type non existant", e%itype
       STOP
    END SELECT

  END FUNCTION base_element

  !*--------------------------------------------------------------------------------
  !*eval_func_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the solution in y 
  !*--------------------------------------------------------------------------------
  TYPE(PVar) FUNCTION eval_func_element(e, u, y)
    CLASS(element),                       INTENT(in):: e
    TYPE(PVar),    DIMENSION(e%nsommets), INTENT(in):: u
    REAL(dp),       DIMENSION(2),          INTENT(in):: y
    !*REAL(dp),       DIMENSION(2)                     :: x
    REAL(dp),       DIMENSION(e%nsommets)            :: base!*, alpha, beta 
    !*TYPE(PVar)                                      :: a,b,c, aa, cc
    INTEGER                                         :: l
    !*LOGICAL                                         :: flag

    DO l=1, e%nsommets !*Loop on the basis functions
       base(l)=e%base(l,y) !*All the basis functions in that point are stored in base
       !*NB e%base is a function
    ENDDO
    DO l=1,n_vars !*Loop on the components
       !*For each component the value in y is given by the SUM of the product between base and the coefficients. So we are just summing all the basis functions evaluated in that DoF multiplied by their coefficients 
       eval_func_element%u(l)=u(1)%u(l)+SUM(base(:)*(u(:)%u(l)-u(1)%u(l))) 
       !*NB: There is this extra u(1)%u(l) outside and inside the sum. Since the sum of the basis functions in any point of the element is = 1 we are just adding and subtracting u(1)%u(l) 
       !*It is a trick used to set exactly u=const and not up to machine precision (which would in any case be good). If u is constant SUM(base(:)*(u(:)%u(l)) is constant up to machine precision. If we use the trick we have (u(:)%u(l)-u(1)%u(l))=0._DP and so
       !*eval_func_element%u(l)=u(1)%u(l)+SUM(base(:)*(u(:)%u(l)-u(1)%u(l)))=u(1)%u(l)
       !*In the general case nothing changes up to machine precision 
    END DO


  END FUNCTION eval_func_element

  !*--------------------------------------------------------------------------------
  !*eval_der_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the derivative of the solution in y 
  !*NB: We specify the barycentric coordinate but the output gradient is in the "real element"
  !*We keep into account the passage from the reference element to the real one dividing by the length the gradient in the reference element inside gradient_element (grad=fx/e%volume) 
  !*--------------------------------------------------------------------------------
  TYPE(PVar) FUNCTION eval_der_element(e,u,y)
    CLASS(element), INTENT(in)::e
    TYPE(PVar), DIMENSION(e%nsommets),INTENT(in):: u
    REAL(dp), DIMENSION(2), INTENT(in):: y
    REAL(dp), DIMENSION(1,e%nsommets):: grad !*IMPORTANT: the first component in the shape is only 1
    INTEGER:: l

    DO l=1, e%nsommets !*Loop on the basis functions
       grad(:,l)=gradient_element(e,l,y) !*Collecting in grad the gradient of the basis functions in the point of barycentric coordinates y !*IMPORTANT: the first component in the shape is only 1 so even if there are the : we are talking about a scalar value. In 2D it would make more sense to write : because we'd have
    ENDDO
    DO l=1,n_vars !*Loop on the components
       eval_der_element%u(l)=SUM( ( u(:)%u(l)-u(1)%u(l) ) *grad(1,:) ) !*For each component the derivative is the sum of the derivatives of the basis functions in that point multiplied by the respective coefficients (the product u(:)%u(l)*grad(1,:) )
       !*Again we have the extra u(1)%u(l). It is equivalent up to machine precision because
       !*The sum of the basis fuctions is constant=1 => The gradient of the sum (or the sum of the gradients) is 0 => SUM(u(1)%u(l) *grad(1,:))=0
       !*The pro is that if the solution is constant in the cell we get exactly 0 because u(:)%u(l)-u(1)%u(l) would cancel otherwise it would be only 0 up to machine precision
    END DO

  END FUNCTION eval_der_element

  !*--------------------------------------------------------------------------------
  !*eval_der2_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the SECOND derivative of the solution in y 
  !*NB: We specify the barycentric coordinate but the output gradient is in the "real element"
  !*We keep into account the passage from the reference element to the real one dividing by the length SQUARED the gradient in the reference element inside gradient2_element (grad2=fxx/(e%volume**2)) 
  !*--------------------------------------------------------------------------------
  TYPE(PVar) FUNCTION eval_der2_element(e,u,y)
    CLASS(element), INTENT(in)::e
    TYPE(PVar), DIMENSION(e%nsommets),INTENT(in):: u
    REAL(dp), DIMENSION(2), INTENT(in):: y
    REAL(dp), DIMENSION(1,e%nsommets):: grad2
    INTEGER:: l

    DO l=1, e%nsommets
       grad2(:,l)=e%gradient2(l, y) !*Standard way: second derivative of all the basis functions in the point
    ENDDO
    DO l=1,n_vars
       eval_der2_element%u(l)=SUM( ( u(:)%u(l)-u(1)%u(l) ) *grad2(1,:)  ) !*Sum of the second derivative of the basis functions in the point y multiplied by the coefficients
       !*Standard trick with u(1)%u(l) to achieve exactly 0 in case of a constant solution
    END DO
  END FUNCTION eval_der2_element


  !*---------------------------------------------------------------------
  !*gradient_element(k,x) computes the first derivative of the basis function k in the point of baricentric coordinates x
  !*NB: The first derivative is already in the "real" element, not in the reference one, because of the final division by e%volume
  !*---------------------------------------------------------------------
  FUNCTION gradient_element(e,k,x) RESULT (grad) !*derivative in the real element (due to the final division by the length)
    CLASS(element), INTENT(in):: e
    INTEGER, INTENT(in):: k ! numero de la fonction de base
    REAL(dp), DIMENSION(2), INTENT(in):: x ! coordonnees barycentriques
    REAL(dp),DIMENSION(n_dim):: grad
    REAL(dp):: fx, alpha, s

    SELECT CASE(e%itype) !*Choice on the type of basis function
    CASE(1,11)! P1
       SELECT CASE(k) !*Choice on the basis function
       CASE(1) !*RMK: The first basis function is associated to the left node so the derivative is -1 in this case
          fx=-1.0_dp
       CASE(2) !*OK
          fx= 1.0_dp
       END SELECT
    CASE(3,12) !P2
       SELECT CASE(k)
       CASE(1)
          fx=-3.0_dp+4.0_dp*x(1)
       CASE(2)
          fx=4.0_dp*x(1)-1.0_dp
       CASE(3) !*Associated to the internal DoF, NB: It is 0 in x=(1/2,1/2)
          fx=4.0_dp-8.0_dp*x(1)
       END SELECT
    CASE(2) ! B2
       SELECT CASE(k)
       CASE(1)
          fx=-2.0_dp*(1.0_dp-x(1))
       CASE(2)
          fx=2.0_dp*x(1)
       CASE(3)
          fx=2.0_dp-4.0_dp*x(1)
       END SELECT
    CASE(4) ! P3
       SELECT CASE(k)
       CASE(1)
          fx=-11._dp/2._dp+(-27._dp/2._dp*x(1)+18._dp)*x(1)
       CASE(2)
          fx=1.0_dp+x(1)*(27._dp/2._dp*x(1)-9._dp)
       CASE(3)
          fx=9._dp+x(1)*(81._dp/2._dp*x(1)-45._dp)
       CASE(4)
          fx=-9._dp/2._dp+x(1)*(-81._dp/2._dp*x(1)+36._dp)
       END SELECT
    CASE(5)
       SELECT CASE(k)
       CASE(1)
          fx=-3.0_dp*x(2)*x(2)
          !fx=-3.0*x(2)*x(2)
       CASE(2)
          fx=3.0_dp*x(1)*x(1)
          !fx=3.0*x(1)*x(1)
       CASE(3)
          fx= 3.0_dp*x(2)*( x(2)-2._dp*x(1) )
          !fx=3.0*x(2)*( x(2)-2.*x(1) )
       CASE(4)
          fx= 3.0_dp*x(1)*(2._dp*x(2)-x(1))
          !fx=6*x(1)*x(2)-3*x(1)*x(1)
       END SELECT
    CASE(6) !B4
       SELECT CASE(k)
       CASE(1)
          fx=-4.0_dp*x(2)*x(2)*x(2)
       CASE(2)
          fx= 4.0_dp* x(1)*x(1)*x(1)
       CASE(3)
          fx= 4.0_dp*( x(2)-3.0_dp*x(1) )*x(2)*x(2)
       CASE(4) 
          fx=12.0_dp * x(1)*x(2) * ( x(2)-x(1) )
       CASE(5)
          fx=4.0_dp*( 3._dp*x(2)-x(1) )*x(1)*x(1)
       END SELECT

    CASE(7) !P4
       s=x(1)
       SELECT CASE(k)
       CASE(1)
          !fx=-4.0+12.0*x(1)-12.0*x(1)*x(1)+4.0*x(1)*x(1)*x(1)
          fx= (((128._dp*s)/3._dp - 80._dp)*s + 140._dp/3._dp)*s - 25._dp/3._dp
       CASE(2)
          !fx=4.0*x(1)*x(1)*x(1)
          fx=(((128._dp*s)/3._dp - 48._dp)*s + (44._dp)/3._dp)*s - 1._dp
       CASE(3)
          !fx=4.0-24.0*x(1)+36.0*x(1)*x(1)-16.0*x(1)*x(1)*x(1)
          fx= ((- (512._dp*s)/3._dp + 288._dp)*s - (416._dp)/3._dp)*s + 16._dp
       CASE(4) 
          !fx=12.0*x(1)-36.0*x(1)*x(1)+24.0*x(1)*x(1)*x(1)
          fx=((256._dp*s - 384._dp)*s + 152._dp)*s - 12._dp
       CASE(5)
          !fx=12.0*x(1)*x(1)-16.0*x(1)*x(1)*x(1)
          fx=  s*(s*(-(512._dp*s)/3._dp + 224._dp) - 224._dp/3._dp) + 16._dp/3._dp
       END SELECT

    CASE(13) !P3 Gauss Lobatto
			 alpha=0.5_dp-SQRT(5._dp)/10._dp
       SELECT CASE(k)
       CASE(1)
          fx=(- alpha**2 + alpha + 3._dp*x(1)**2._dp - 4._dp*x(1) + 1._dp)/(alpha*(alpha - 1._dp))
          !fx=-3.0*x(2)*x(2)
       CASE(2)
          fx=-(- alpha**2 + alpha + 3._dp*x(1)**2._dp - 2._dp*x(1))/(alpha*(alpha - 1._dp))
          !fx=3.0*x(1)*x(1)
       CASE(3)
          fx=(2._dp*alpha*x(1) - 4._dp*x(1) - alpha + 3._dp*x(1)**2._dp + 1._dp)/(alpha*(2._dp*alpha**2 - 3._dp*alpha + 1._dp))
          !fx=3.0*x(2)*( x(2)-2.*x(1) )
       CASE(4)
          fx= -(alpha - 2._dp*x(1) - 2._dp*alpha*x(1) + 3._dp*x(1)**2)/(alpha*(2._dp*alpha**2 - 3._dp*alpha + 1._dp))
          !fx=6*x(1)*x(2)-3*x(1)*x(1)
       END SELECT
    CASE(14) !P4 Gauss Lobatto
			 alpha=0.5_dp-SQRT(21._dp)/14._dp
			 s=x(1)
       SELECT CASE(k)
       CASE(1)
          fx=(4._dp*alpha**2*s - 3._dp*alpha**2 - 4._dp*alpha*s + 3._dp*alpha - 8._dp*s**3 + 15._dp*s**2 - 8._dp*s + 1._dp)/(alpha*(alpha - 1._dp))
          !fx=-3.0*x(2)*x(2)
       CASE(2)
          fx=-(- 4._dp*alpha**2*s + alpha**2 + 4._dp*alpha*s - alpha + 8._dp*s**3 - 9._dp*s**2 + 2*s)/(alpha*(alpha - 1._dp))
          !fx=3.0*x(1)*x(1)
       CASE(3)
          fx=(alpha + 8._dp*s - 6._dp*alpha*s + 6._dp*alpha*s**2 - 15._dp*s**2 + 8._dp*s**3 - 1._dp)/(alpha*(2._dp*alpha - 1._dp)**2*(alpha - 1._dp))
          !fx=3.0*x(2)*( x(2)-2.*x(1) )
       CASE(4)
          fx= (16._dp*(2._dp*s - 1._dp)*(- alpha**2 + alpha + 2._dp*s**2 - 2._dp*s))/(2._dp*alpha - 1._dp)**2
          !fx=6*x(1)*x(2)-3*x(1)*x(1)
	   CASE(5)
	      fx = -(alpha - 2._dp*s - 6._dp*alpha*s + 6._dp*alpha*s**2 + 9._dp*s**2 - 8._dp*s**3)/(alpha*(2._dp*alpha - 1._dp)**2*(alpha - 1._dp))
       END SELECT

    CASE default
       PRINT*, "Type non existant", e%itype
       STOP
    END SELECT
    grad=fx/e%volume !*Division by the length to pass from the reference element to the real element
  END FUNCTION gradient_element

  !*---------------------------------------------------------------------
  !*gradient2_element(k,x) computes the SECOND derivative of the basis function k in the point of baricentric coordinates x
  !*NB: The second derivative is already in the "real" element, not in the reference one, because of the final division by e%volume**2
  !*---------------------------------------------------------------------
  FUNCTION gradient2_element(e,k,x) RESULT (grad2) !*second derivative in the real element (due to the final division by the length^2)
    CLASS(element), INTENT(in):: e
    INTEGER, INTENT(in):: k ! numero de la fonction de base
    REAL(dp), DIMENSION(2), INTENT(in):: x ! coordonnees barycentriques
    REAL(dp),DIMENSION(n_dim):: grad2
    REAL(dp):: fxx, alpha,s

    SELECT CASE(e%itype)
    CASE(1,11)! P1 !*OK 0
       SELECT CASE(k) 
       CASE(1) 
          fxx=0.0_dp
       CASE(2)
          fxx=0.0_dp
       END SELECT
    CASE(3,12) !P2 !*OK constant
       SELECT CASE(k) 
       CASE(1) 
          fxx=4.0_dp
       CASE(2)
          fxx=4.0_dp
       CASE(3)
          fxx=-8.0_dp
       END SELECT
    CASE(2) ! B2 !*OK CONSTANT
       SELECT CASE(k)
       CASE(1)
          fxx=2.0_dp
       CASE(2)
          fxx=2.0_dp
       CASE(3)
          fxx=-4.0_dp
       END SELECT
    CASE(4) ! P3
			 alpha=1._dp/3._dp
       SELECT CASE(k)
       CASE(1)
          !fxx=-6.0*x(1)+6.0
          fxx=(6._dp*x(1) - 4._dp)/(alpha*(alpha - 1._dp))
       CASE(2)
          !fxx=6.0*x(1)
          fxx=-(6._dp*x(1) - 2._dp)/(alpha*(alpha - 1._dp))
       CASE(3)
          !fxx=18.0*x(1)-12.0
          fxx= (2._dp*alpha + 6._dp*x(1) - 4._dp)/(alpha*(2._dp*alpha**2._dp - 3._dp*alpha + 1._dp))
       CASE(4)
          !fxx=-18.0*x(1)+6.0
          fxx=(2._dp*alpha - 6._dp*x(1) + 2._dp)/(alpha*(2._dp*alpha**2._dp - 3._dp*alpha + 1._dp))
       END SELECT
    CASE(5)
       SELECT CASE(k)
       CASE(1)
          fxx=6.0_dp*x(2)
       CASE(2)
          fxx=6.0_dp*x(1)
       CASE(3)
          fxx= 6.0_dp*x(1)-12.0_dp*x(2)
       CASE(4)
          !fxx=-18.0*x(1)+6.0
          fxx=6.0_dp*x(2)-12.0_dp*x(1)
       END SELECT

    CASE(6)
       SELECT CASE(k)
       CASE(1)
          fxx=12.0_dp-24.0_dp*x(1)+12.0_dp*x(1)*x(1)
       CASE(2)
          fxx=12.0_dp*x(1)*x(1)
       CASE(3)
          fxx=-24.0_dp+72.0_dp*x(1)-48.0_dp*x(1)*x(1)
       CASE(4)
          fxx=12.0_dp-72.0_dp*x(1)+72.0_dp*x(1)*x(1)
       CASE(5)
          fxx=24.0_dp*x(1)-48.0_dp*x(1)*x(1)
       END SELECT

    CASE(7) !P4
       s=x(1)
       SELECT CASE(k)
       CASE(1)
          fxx=(128._dp*s - 160._dp)*s + 140._dp/3._dp
       CASE(2)
          fxx= (128._dp*s - 96._dp)*s + 44._dp/3._dp
       CASE(3)
          fxx= (-512._dp*s + 576._dp)*s - 416._dp/3._dp
       CASE(4)
          fxx= 768._dp*s*(s - 1._dp) + 152._dp
       CASE(5)
          fxx= (- 512._dp*s + 448._dp)*s - 224._dp/3._dp
       END SELECT

    CASE(13)
			 alpha=0.5_dp-SQRT(5._dp)/10._dp
       SELECT CASE(k)
       CASE(1)
          !fxx=-6.0*x(1)+6.0
          fxx=(6._dp*x(1) - 4._dp)/(alpha*(alpha - 1._dp))
       CASE(2)
          !fxx=6.0*x(1)
          fxx=-(6._dp*x(1) - 2._dp)/(alpha*(alpha - 1._dp))
       CASE(3)
          !fxx=18.0*x(1)-12.0
          fxx= (2._dp*alpha + 6._dp*x(1) - 4._dp)/(alpha*(2._dp*alpha**2 - 3._dp*alpha + 1._dp))
       CASE(4)
          !fxx=-18.0*x(1)+6.0
          fxx=(2._dp*alpha - 6._dp*x(1) + 2._dp)/(alpha*(2._dp*alpha**2 - 3._dp*alpha + 1._dp))
       END SELECT

    CASE(14)
			 alpha=0.5_dp-SQRT(21._dp)/14._dp
			 s=x(1)
       SELECT CASE(k)
       CASE(1)
          !fxx=-6.0*x(1)+6.0
          fxx=-(- 4._dp*alpha**2 + 4._dp*alpha + 24._dp*s**2 - 30._dp*s + 8._dp)/(alpha*(alpha - 1._dp))
       CASE(2)
          !fxx=6.0*x(1)
          fxx=-(- 4._dp*alpha**2 + 4._dp*alpha + 24._dp*s**2 - 18._dp*s + 2._dp)/(alpha*(alpha - 1._dp))
       CASE(3)
          !fxx=18.0*x(1)-12.0
          fxx= (12._dp*alpha*s - 30._dp*s - 6._dp*alpha + 24._dp*s**2 + 8._dp)/(alpha*(2._dp*alpha - 1._dp)**2*(alpha - 1._dp))
       CASE(4)
          !fxx=-18.0*x(1)+6.0
          fxx=(32._dp*(- alpha**2 + alpha + 6._dp*s**2 - 6._dp*s + 1._dp))/(2._dp*alpha - 1._dp)**2
       CASE(5)
          fxx = (6._dp*alpha - 18._dp*s - 12._dp*alpha*s + 24._dp*s**2 + 2._dp)/(alpha*(2._dp*alpha - 1._dp)**2*(alpha - 1._dp))
       END SELECT

    CASE default
       PRINT*, "Type non existant", e%itype
       STOP
    END SELECT
    grad2=fxx/(e%volume**2)
  END FUNCTION gradient2_element

  !*-----------------------------------------------------------
  !*average(u) takes in input the coefficients of the solution (or a general function) in the DoFs and computes the average in the cell
  !*u->coefficients in the different nodes 
  !*-----------------------------------------------------------
  FUNCTION average_element(e,u) RESULT(u_bar)
    CLASS(element), INTENT(in):: e		
    TYPE(PVar), DIMENSION(e%nsommets),INTENT(in):: u
    TYPE(PVar)  :: u_bar, u_loc
    INTEGER:: iq!*, i
    u_bar=0._dp

    DO iq = 1, e%nquad !*Loop on the quadrature point

       u_loc   =e%eval_func(u,e%quad(:,iq)) !*Evaluation of the function in the quadrature point
       u_bar = u_bar + (u_loc)*e%weight(iq) !*Add the point evaluation multiplied by the weight
       !*NB: This is the average because he should multiply by the area in the end to get the integral and divide by the area to get the average. The area simplifies in this process so it is not kept into account

    ENDDO ! iq
  END FUNCTION average_element

  !*--------------------------------------------------------------------
  !*quadrature_element fills the weights and the nodes of the quadrature formulas, it is a subroutine in fact
  !*Here I trust and do not check, to be checked and tested eventually
  !*--------------------------------------------------------------------
  SUBROUTINE quadrature_element(e)  ! quadrature points and weights
    ! for triangle and edges
    CLASS(element), INTENT(inout):: e
    REAL(dp):: w,xo,yo,zo,s
    INTEGER:: nquad, nquad_1
    ! degree 5
    REAL(dp), PARAMETER:: s1_3=SQRT(0.6_dp)
    REAL(dp), PARAMETER:: s2_3=0.5_dp*(1.0_dp - s1_3)
    REAL(dp), PARAMETER:: s3_3=0.5_dp*(1.0_dp + s1_3)
    REAL(dp), PARAMETER:: weight1_3=5._dp/18._dp
    REAL(dp), PARAMETER:: weight2_3=8._dp/18._dp
    ! degre 7
    REAL(dp), PARAMETER:: s1_7=(18._dp-SQRT(30._dp))/72._dp
    REAL(dp), PARAMETER:: s2_7=SQRT(525._dp+70._dp*SQRT(30._dp))/35._dp
    REAL(dp), PARAMETER:: s3_7=0.50_dp*(1.00_dp-s2_7)
    REAL(dp), PARAMETER:: s4_7=0.50_dp*(1.00_dp+s2_7)
    REAL(dp), PARAMETER:: s5_7=(18._dp+  SQRT(30._dp))/72._dp 
    REAL(dp), PARAMETER:: s6_7=SQRT(525._dp-70._dp*SQRT(30._dp))/35._dp
    REAL(dp), PARAMETER:: s7_7=0.50_dp*(1.00_dp-s6_7)
    REAL(dp), PARAMETER:: s8_7=0.50_dp*(1.00_dp+s6_7)
    ! degree 9
    REAL(dp),PARAMETER:: s1_9=1.0_dp/3._dp*SQRT(5._dp-2._dp*SQRT(10._dp/7._dp))
    REAL(dp),PARAMETER:: s2_9=0.50_dp*(1.00_dp - s1_9)
    REAL(dp),PARAMETER:: s3_9=0.50_dp*(1.00_dp + s1_9)
    REAL(dp),PARAMETER:: weight1_9=(322._dp+13._dp*SQRT(70.0_dp))/1800.0_dp


    REAL(dp),PARAMETER:: s5_9=1.0_dp/3._dp*SQRT(5._dp+2._dp*SQRT(10._dp/7._dp))
    REAL(dp),PARAMETER:: s6_9=0.50_dp*(1.00_dp - s5_9)
    REAL(dp),PARAMETER:: s7_9=0.50_dp*(1.00_dp + s5_9)
    REAL(dp),PARAMETER:: weight2_9=(322._dp-13._dp*SQRT(70.0_dp))/1800.0_dp
    REAL(dp),PARAMETER:: s4_9=1._dp/3._dp*SQRT(5._dp+2._dp*SQRT(10.0_dp/7.0_dp))

    SELECT CASE(e%itype) !*Depending on the type of element (in particular) on the degree of the polynomial chosen a suitable quadrature is chosen
    CASE(1,11)
       e%nquad=2
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad))

       ! quadrature for edges (3 point gauss formula)

       e%quad(1,1) = 1.0_dp
       e%quad(2,1) = 0._dp
       e%weight(1) = 0.5_dp

       e%quad(1,2) = 0.0_dp
       e%quad(2,2) = 1.0_dp
       e%weight(2) = 0.5_dp
       
       !*NB: For P1 the quadrature is
       !*int=1/2(f(x_1)+f(x_2)) [error O(h^3)]
       !*This results in a lumping of the local mass matrix which is diagonal 

    CASE(2,3)!! exact for degree 5
       e%nquad=3
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad))



       !       s=sqrt(0.6_dp) !0.7745966692414834_dp  !SQRT(0.6_dp)
       e%quad(1,1)=s2_3 !0.5_dp*(1.0_dp - s)
       e%quad(2,1)=s3_3 !0.5_dp*(1.0_dp + s)
       e%weight(1) =weight1_3 !5._dp/18._dp! 0.2777777777777778 !5.0_dp/18.0_dp

       e%quad(1,2)=s3_3 !0.5_dp*(1.0_dp + s)
       e%quad(2,2)=s2_3 !0.5_dp*(1.0_dp - s)
       e%weight(2) =weight1_3 !5._dp/18._dp! 0.2777777777777778 !5.0_dp/18.0_dp

       e%quad(1,3)=0.5_dp
       e%quad(2,3)=0.5_dp
       e%weight(3)= weight2_3 !8._dp/18._dp!0.4444444444444444 ! 8.0_dp/18.0_dp

    CASE(12)!! exact for degree 3, underintegration
       e%nquad = 3
       nquad = e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad))
       
       e%quad(1,1)=0._dp
       e%quad(2,1)=1.0_dp
       e%weight(1)=1._dp/6._dp


       e%quad(1,2)=0.5_dp
       e%quad(2,2)=0.5_dp
       e%weight(2)=2._dp/3._dp

       e%quad(1,3)=1._dp
       e%quad(2,3)=0._dp
       e%weight(3)=1._dp/6._dp

    CASE(4,5) !! exact for degree 7
       e%nquad=4 ! ordre 7
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad) )
       !      s=(18._dp-sqrt(30._dp))/72._dp !0.1739274225687269_dp  !(18.- SQRT(30.))/72. !SQRT(30._dp)*(-5.0_dp+3._dp*SQRT(30._dp))/360._dp
       e%weight(1:2)=s1_7
       !      s=SQRT(525._dp+70._dp*SQRT(30._dp))/35._dp !0.8611363115940526_dp!SQRT(525._dp+70._dp*SQRT(30._dp))/35._dp
       e%quad(1,1)=s3_7! 0.50_dp*(1.00_dp-s)
       e%quad(1,2)=s4_7 !0.50_dp*(1.00_dp+s)
       e%quad(2,1)=e%quad(1,2) !0.5_dp*(1.0_dp+s)
       e%quad(2,2)=e%quad(1,1) !0.5_dp*(1.0_dp-s)

       !      s=(18._dp+  SQRT(30._dp))/72._dp !0.3260725774312731_dp !(18.+  SQRT(30.))/72.! SQRT(30._dp)*(5.0_dp+3._dp*SQRT(30._dp))/360._dp
       e%weight(3:4)=s5_7
       !       s=SQRT(525._dp-70._dp*SQRT(30._dp))/35._dp !0.3399810435848563_dp !SQRT(525._dp-70._dp*SQRT(30._dp))/35._dp
       e%quad(1,3)=s7_7 !0.50_dp*(1.00_dp-s)
       e%quad(1,4)=s8_7 !0.50_dp*(1.00_dp+s)
       e%quad(2,3)=e%quad(1,4) !0.5_dp*(1.0_dp+s)
       e%quad(2,4)=e%quad(1,3) !0.5_dp*(1.0_dp-s)

       !       e%quad(2,:)=1.0_dp-e%quad(1,:)
    CASE(13)!! exact for order 4, underintegration
       e%nquad=4 
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad) )

       s=0.5_dp-SQRT(5._dp)/10._dp

       e%quad(1,1)=0._dp
       e%quad(2,1)=1.0_dp
       e%weight(1)=1._dp/12._dp


       e%quad(1,2)=s
       e%quad(2,2)=1._dp-s
       e%weight(2)=5._dp/12._dp

       e%quad(1,3)=1._dp -s
       e%quad(2,3)=s
       e%weight(3)=5._dp/12._dp

       e%quad(1,4)=1._dp
       e%quad(2,4)=0._dp
       e%weight(4)=1._dp/12._dp

    CASE(6,7) ! exact for degree 9

       e%nquad=5
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad))

       !s=1.0_dp/3._dp*SQRT(5._dp-2._dp*SQRT(10._dp/7._dp))
       e%quad(1,1)=s2_9 !0.50_dp*(1.00_dp - s)
       e%quad(2,1)=s3_9 !0.50_dp*(1.00_dp + s)
       e%weight(1) = weight1_9 !(322._dp+13._dp*SQRT(70.0_dp))/1800.0_dp
       e%quad(1,2)=s3_9 !0.5_dp*(1.0_dp + s)
       e%quad(2,2)=s2_9 !0.5_dp*(1.0_dp - s)
       e%weight(2) = weight1_9 !(322._dp+13._dp*SQRT(70.0_dp))/1800.0_dp

       !      s=1._dp/3._dp*SQRT(5._dp+2._dp*SQRT(10.0_dp/7.0_dp))
       e%quad(1,3)=s6_9 !0.5_dp*(1.0_dp - s)
       e%quad(2,3)=s7_9 !0.5_dp*(1.0_dp + s)
       e%weight(3) = weight2_9 !(322.0_dp-13.0_dp*SQRT(70.0_dp))/1800.0_dp
       e%quad(1,4)=s7_9 !0.5_dp*(1.0_dp + s)
       e%quad(2,4)=s6_9 !0.5_dp*(1.0_dp - s)
       e%weight(4) = weight2_9 ! (322.0_dp-13.0_dp*SQRT(70.0_dp))/1800.0_dp

       e%quad(1,5)=0.5_dp
       e%quad(2,5)=0.5_dp
       e%weight(5)=64.0_dp/225.0_dp


    CASE(14)!! exact for order 5, underintegration
       e%nquad=5
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad) )

       s=0.5_dp-SQRT(21._dp)/14._dp

       e%quad(1,1)=0._dp
       e%quad(2,1)=1.0_dp
       e%weight(1)=1._dp/20._dp


       e%quad(1,2)=s
       e%quad(2,2)=1._dp-s
       e%weight(2)=49._dp/180._dp

       e%quad(1,3)=0.5_dp
       e%quad(2,3)=0.5_dp
       e%weight(3)=16._dp/45._dp

       e%quad(1,4)=1._dp -s
       e%quad(2,4)=s
       e%weight(4)=49._dp/180._dp

       e%quad(1,5)=1._dp
       e%quad(2,5)=0._dp
       e%weight(5)=1._dp/20._dp

    END SELECT

    e%weight = e%weight/SUM(e%weight)
  END SUBROUTINE quadrature_element

  !*----------------------------------------------------------------------
  !*base_ref_element computes the values of the basis functions in the DoFs which are stored in e%base(:,:)
  !*First index DoF
  !*Second index basis function
  !*NB: THE COMPUTATION IS DONE IN THE OTHER ORDER (First index basis function, second index DoF) but in the end he TRANSPOSE everything
  !*----------------------------------------------------------------------
  SUBROUTINE base_ref_element(e)
    ! compute the values of the basis functions at the physical dofs
    REAL(dp), PARAMETER:: ut=1._DP/3._DP, det=1._DP-ut
    CLASS(element), INTENT(inout)::e
    INTEGER:: l, l_base, l_dof
    REAL(dp), DIMENSION(2,4):: x ! barycentric coordinates for B3
    REAL(dp), DIMENSION(2,5):: z ! barycentric coordinates for B4
    REAL(dp), DIMENSION(2,3):: ba
    ! ebase0( number of the basis function, point l)
    e%base0=0._DP
    SELECT CASE(e%itype)
    CASE(1,3,4,7,11,12,13,14) ! lagrangian functions !*Lagrange, easy, they are 1 in their DoF and 0 in the other DoFs
       DO l=1,e%nsommets
          e%base0(l,l)=1.0_DP
       ENDDO
    CASE(2) ! B2 Bezier !*Bezier
       ba=0.0_DP !*Barycentric coordinates 
       !*Two indices for ba(:,:)
       !*first index->1,2 barycentric coordinates
       !*second index-> DoF       

       !*RMK: Let's focus on the first index, on the barycentric coordinates of a given node
       !*first component associated to the second node x 
       !*second component associated to the first node 1-x

       !*He fills the first components 
       ba(1,1)=0.0_DP !*First vertex->Second node=0
       ba(1,2)=1.0_DP !*Second vertex-SSecond node=1
       ba(1,3)=0.5_DP !*Internal DoF

       !*Second components got by difference with 1
       ba(2,:)=1.0_DP-ba(1,:) 

       DO l_dof=1,3 !*Loop on the DoFs
          DO l_base=1,3 !*Loop on the basis functions
             e%base0(l_base,l_dof)=base_element(e,l_base,ba(:,l_dof)) 
             !*RMK: 
             !*base_element(k,x) !*It calculates a specific basis function k in the point of barycentric coordinates x
             !*k->basis function -> l_base
             !*x->barycentric coordinates -> l_dof
             !*So we'll have in e%base0(:,:) !*BUT NB: HE WILL TRANSPOSE IN THE END
             !*FOR THE MOMENT WE HAVE
             !*First index basis function
             !*Second index DoF
             !*AFTER THE TRANSPOSITION WE'LL HAVE
             !*First index DoF
             !*Second index basis function
          ENDDO
          !write(*,11)l_dof,e%base0(:,l_dof)
          !11        FORMAT(i2,3(1x,f10.4))
       ENDDO


    CASE(5) !*B3
       x=0._DP
       !*Filling the first barycentric coordinate
       x(1,1)=0.0_DP; !*First vertex->Second node=0
       x(1,2)=1.0_DP !*Second vertex->Second node=1
       x(1,3)=1.0_DP/3.0_DP !*First internal DoF->Closer to first vertex (1/3)
       x(1,4)=2.0_DP/3.0_DP !*Second internal DoF->Closer to second vertex (2/3)
       !*Second barycentric coordinate got by difference with 1
       x(2,:)=1.0_DP-x(1,:)
       DO l_dof=1,4 !*Loop on the DoFs
          DO l_base=1,4 !*Loop on the basis functions
             e%base0(l_base,l_dof)=base_element(e,l_base,x(:,l_dof)) !*As before
          ENDDO
       ENDDO

    CASE(6) !*B4
       z=0._DP
       !*Filling first barycentric coordinate
       z(1,1)=0.0_DP; !*First vertex->Second node=0
       z(1,2)=1.0_DP !*Second vertex->Second node=1
       z(1,3)=1.0_DP/4.0_DP !*Closer to first vertex
       z(1,4)=1.0_DP/2.0_DP !*Central
       z(1,5)=0.75_DP !*Closer to second vertex
       !*Second barycentric coordinate
       z(2,:)=1.0_DP-z(1,:) 
       DO l_dof=1,5 !*Loop on the DoFs
          DO l_base=1,5 !*Loop on the basis functions
             e%base0(l_base,l_dof)=base_element(e,l_base,z(:,l_dof)) !*As before
          ENDDO
       ENDDO

       !10     FORMAT(i2,10(1x,f10.4))
    CASE default
       PRINT*, "Element type not yet take into account in base_ref_element", e%itype
       STOP
    END SELECT

    e%base0 = TRANSPOSE(e%base0)
    !*AFTER THE TRANSPOSITION WE HAVE
    !*First index DoF
    !*Second index basis function
 
    e%base1=inverse(e%base0) !*Inverse
 
    !*e%base0=[ phi_1(x_1) phi_2(x_1) ... phi_N(x_1)
    !*          phi_1(x_2) phi_2(x_2) ... phi_N(x_2)
    !*              .           .
    !*              .           .
    !*              .           .
    !*          phi_1(x_N) phi_2(x_N) ... phi_N(x_n) ]
    !*
    !*NB: e%base0*vectorofthecoefficients=vectorofthevalues
    !*vectorofthecoefficients=inv(e%base0)*vectorofthevalues=e%base1*vectorofthevalues
     

    CALL bary_element(e) !*Down here in the same module

  END SUBROUTINE base_ref_element


  !*-----------------------------------------------------
  !*bary_element has not been made public and in fact is recalled here in this module inside e%base_ref_element and not outside. It is used to fill e%x(:,:) barycentric coordinates of the DoFs of the element and e%dof2ind(:) local indices of the nodes with increasing abscissa
  !*RMK: e%x(:,:) has two indices
  !*Second one referred to the DoF of which we consider the barycentric coordinates 1:nsommets
  !*First one referred to the 2 barycentric coordinates 1:2
  !*
  !*For what concerns e%dof2ind(:) local indices of the nodes with increasing abscissa
  !*We have in the natural numeration
  !*At first the indices
  !*Then the internal DoF with increasing abscissa
  !*So e%dof2ind(:) is for example
  !*P1-> 1, 2
  !*B2-> 1, 3, 2
  !*B3,P3-> 1, 3, 4, 2
  !*-----------------------------------------------------
  SUBROUTINE bary_element(e)
    ! defines barycentric coordinates of Lagrange points
    CLASS(element), INTENT(inout)::e
    ALLOCATE(e%x(2,e%nsommets)) !*e%x(:,:) has two indices
    !*Second one referred to the DoF of which we consider the barycentric coordinates 1:nsommets
    !*First one referred to the 2 barycentric coordinates 1:2
    ALLOCATE(e%dof2ind(e%nsommets)) 



    !*NB: HE WILL FILL THE FIRST BARYCENTRIC COORDINATE OF ALL THE DoFs (1,:) AND THEN COMPUTE THE OTHER ONE COMPLEMENTING WITH 1 -> e%x(2,:)=1._DP-e%x(1,:)
 
    !*(:,1) and (:,2) referred to the extrema, the first 2 DoFs, always present.
    e%x(1,1)=0.0_DP !*First coordinate (referred to the second node) of the first vertex -> 0
    e%x(1,2)=1.0_DP !*First coordinate of the second vertrex ->1
    !*NB: (:,1) and (:,2) always present, the others will be eventually filled but just the first coordinate, the other one wil be got in the end by complementing with 1


    SELECT CASE (e%itype)
    CASE(1,11) !*P1, no other x
       e%dof2ind = (/1, 2/) !*Clear
    CASE(2,3,12) !*B2,P2 another x
       e%x(1,3)=0.5_DP
       e%dof2ind = (/1, 3, 2/) !*Clear
    CASE(4,5) !*B3,P3 2 other x
       e%x(1,3)=1._DP/3._DP !*Closer to the first 
       e%x(1,4)=2._DP/3._DP !*Closer to the second
       e%dof2ind = (/1, 3, 4, 2/)
    CASE(13) !*Gauss Lobatto
       e%x(1,3)=0.5_dp-SQRT(5._dp)/10._dp
       e%x(1,4)=0.5_dp+SQRT(5._dp)/10._dp
       e%dof2ind = (/1, 3, 4, 2/)
    CASE(6,7)
       e%x(1,3)=0.25_DP !*Closer to first
       e%x(1,4)=0.5_DP !*Central
       e%x(1,5)=0.75_DP !*Closer to second
       e%dof2ind = (/1, 3, 4, 5, 2/) !*Clear
    CASE(14)
       e%x(1,3)=0.5_dp-SQRT(21._dp)/14._dp
       e%x(1,4)=0.5_DP
       e%x(1,5)=0.5_dp+SQRT(21._dp)/14._dp
       e%dof2ind = (/1, 3, 4, 5, 2/)
    CASE default
       PRINT*, "bary_element, cas not implemented"
       STOP
    END SELECT
    e%x(2,:)=1._DP-e%x(1,:) !*Complementing with 1 to get the other coordinate
  END SUBROUTINE bary_element



!*  !*------------------------------------------------------------------------
!*  !*NOT USED, IT CAN BE REMOVED. IT IS NOT USED HERE AND NOT EVEN A PUBLIC PROCEDURE
!*  !*ANYWAY IT IS NOT WRONG
!*  !*Evaluation of some coefficients needed for the DeC
!*  !*e%mass(i,j)-> local mass matrix int_K phi_i phi_j
!*  !*e%coeff(i,j) -> another matrix int_K d_x phi_i phi_j 
!*  !*------------------------------------------------------------------------
!*  SUBROUTINE eval_coeff_element_old(e)
!*    ! compute \int_K \nabla phi_sigma * phi_sigma'
!*    CLASS(element), INTENT(inout)::e
!*    REAL(dp), DIMENSION(e%nsommets):: base
!*    REAL(dp), DIMENSION(n_dim,e%nsommets):: grad
!*
!*    INTEGER:: i, iq, l
!*
!*    ALLOCATE(e%coeff(e%nsommets,e%nsommets)) !*int_K d_x phi_i phi_j
!*    ALLOCATE(e%mass (e%nsommets,e%nsommets)) !*int_K phi_i phi_j
!*    e%coeff=0._dp
!*    e%mass=0._dp
!*    DO iq=1,e%nquad !*Loop on the quadrature points, for each quadrature point we'll make point wise evaluations of the quantities to integrate and we will sum them after multiplying by the weight
!*       DO i=1, e%nsommets !*Loop on the basis functions i
!*          DO l=1, e%nsommets !*Loop on the basis funtions l
!*             base(l  )=e%base    (l, e%quad(:,iq)) !*Values of the basis functions in the quadrature point iq
!*             grad(:,l)=e%gradient(l, e%quad(:,iq)) !*Value of the derivative of the basis functions in the quadrature point iq
!*          ENDDO
!*          DO l=1, e%nsommets !*We can reuse the index l. Loop on the basis functions l
!*             e%coeff(l,i)=e%coeff(l,i)+ grad(1,l)*base(i)*e%weight(iq) !*NB: Derivative on the first index e%coeff(,)=d_x phi * phi
!*             e%mass (l,i)=e%mass (l,i)+ base(l)  *base(i)*e%weight(iq) !*e%mass(,)=phi * phi
!*          ENDDO
!*       ENDDO
!*    ENDDO
!*    !*To complete the integration with the quadrature we have to keep into account the real dimension of the segment so we multiply by the length (e%volume)
!*    e%coeff=e%coeff * e%volume 
!*    e%mass=e%mass * e%volume
!*
!*  END SUBROUTINE eval_coeff_element_old


  !*------------------------------------------------------------------------
  !*eval_coeff_element evaluates some coefficients needed for the DeC
  !*e%mass(i,j)-> local mass matrix int_K phi_i phi_j
  !*e%coeff(i,j) -> another matrix int_K phi_i  d_x phi_j !*NB:DERIVATIVE ON THE SECOND INDEX
  !*------------------------------------------------------------------------
  SUBROUTINE eval_coeff_element(e)
    ! compute \int_K \nabla phi_sigma * phi_sigma'
    CLASS(element), INTENT(inout)::e
    REAL(dp), DIMENSION(e%nsommets):: base
    REAL(dp), DIMENSION(n_dim,e%nsommets):: grad

    !*REAL(dp):: eps !*NOT NEEDED
    INTEGER:: i, iq, l, j, k
    !    REAL(dp),  DIMENSION(2,2),parameter:: xp=reshape( (/0._dp,1._dp, 1._dp,0._dp/),(/2,2/) )
    !*NOT NEEDED !*REAL(dp), DIMENSION(2,2):: xp
    !*NOT NEEDED !*xp(1,1)=0.0_dp; xp(2,1)=1.0_dp; xp(1,2)=1.0_dp; xp(2,2)=0.0_dp


    ALLOCATE(e%coeff  (e%nsommets,e%nsommets)) !*e%coeff(i,j) -> another matrix int_K d_x phi_i phi_j !*NB:DERIVATIVE ON THE SECOND INDEX
    ALLOCATE(e%mass   (e%nsommets,e%nsommets)) !*e%mass(i,j)-> local mass matrix int_K phi_i phi_j
!    ALLOCATE(e%coeff_b(e%nsommets,e%nsommets))

    e%coeff=0._dp
    e%mass=0._dp
    DO iq=1,e%nquad !*Loop on the quadrature points, for each quadrature point we'll make point wise evaluations of the quantities to integrate and we will sum them after multiplying by the weight
       DO i=1, e%nsommets !*Loop on the basis function i
          DO l=1, e%nsommets !*Loop on the basis functions to compute the value of the basis functions and of their derivative in the quadrature point
             base(l  )=e%base    (l, e%quad(:,iq))
             grad(:,l)=e%gradient(l, e%quad(:,iq))
          ENDDO
          DO l=1, e%nsommets !*Loop on the basis functions !*Actual quadrature
             e%coeff(l,i)=e%coeff(l,i)+ grad(1,i)*base(l)*e%weight(iq) !*NB: Derivative on the SECOND INDEX
             e%mass (l,i)=e%mass (l,i)+ base(l)  *base(i)*e%weight(iq)
          ENDDO
       ENDDO
    ENDDO
    !*To complete the integration with the quadrature we have to keep into account the real dimension of the segment so we multiply by the length (e%volume)
    e%coeff=e%coeff * e%volume
    e%mass=e%mass   * e%volume

!!$    ! coefficients for the boundary term
!!$    e%coeff_b=0._dp
!!$    eps=1._dp
!!$    
!!$    DO k=1,2
!!$       
!!$       DO l=1, e%nsommets
!!$          base(l)=e%base(l,xp(:,k))
!!$       ENDDO
!!$       
!!$       DO i=1,e%nsommets
!!$          DO j=1, e%nsommets
!!$             e%coeff_b(i,j)=e%coeff_b(i,j)+ eps*base(i)*base(j)
!!$          ENDDO
!!$       ENDDO
!!$
!!$       eps = -eps
!!$    ENDDO


  END SUBROUTINE eval_coeff_element

  SUBROUTINE clean(e) !*CLEANING SUB
    TYPE(element), INTENT(inout)::e
    IF (ASSOCIATED(e%coor))  NULLIFY(e%coor)
    IF (ASSOCIATED(e%n))  NULLIFY(e%nu)
    IF (ASSOCIATED(e%n))  NULLIFY(e%n)
    IF (ASSOCIATED(e%quad))  NULLIFY(e%quad)
    IF (ASSOCIATED(e%weight))  NULLIFY(e%weight)
    IF (ASSOCIATED(e%base0))  NULLIFY(e%base0)
    IF (ASSOCIATED(e%base1))  NULLIFY(e%base1)
  END SUBROUTINE clean

END MODULE element_class
