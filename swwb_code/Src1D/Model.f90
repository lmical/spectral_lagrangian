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

!*----------------------------------------------------------------------------------------
!*Module to pass from coefficients to values and viceversa
!*contr=coefficient
!*cons=values (cons stays for conserved variables)
!*----------------------------------------------------------------------------------------
MODULE Model

  ! This module allows to pass from control to physical variables and viceversa

  USE param2d
  USE PRECISION
  IMPLICIT NONE
CONTAINS

  !*--------------------------------------------------------------------------------------
  !*coeff->values
  !*--------------------------------------------------------------------------------------
  FUNCTION Control_to_Cons(u,e) RESULT (u_cons)
    ! crom control point, compute the values a physical dofs
    TYPE(Pvar), DIMENSION(:), INTENT(in):: u !*coeff in input
    TYPE(element), INTENT(in):: e
    TYPE(Pvar),DIMENSION(SIZE(u,dim=1)):: u_cons !*values as output
    INTEGER:: k,l

    u_cons=u !*Set u_cons=u

    SELECT CASE(e%itype)
    CASE(1,3,4,7,11,12,13,14) ! Lagrange
       !*In this case there is nothing to do because values=coeff
    CASE(2,5,6) ! Bezier
       !*In this case we need to do a projection
       DO l=1, e%nsommets !*Loop on the DoFs
          u_cons(l) = e%eval_func(u(:),e%x(:,l)) !*He's storing in u_cons the pointwise evaluations of the solution in the DoFs
          !*eval_func_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the solution in y 
          !*u->coefficients in the different nodes 
          !*y->barycentric coordinates of the point 

          !*e%x are the barycentric coordinates of the DoFs
          !*2 indices x(:,:)
          !*Second one referred to the DoF of which we consider the barycentric coordinates 1:nsommets
          !*First one referred to the 2 barycentric coordinates 1:2
       ENDDO

    CASE default
       PRINT*, "erreur dans Model/Control_to_Cons"
       STOP
    END SELECT
  END FUNCTION Control_to_cons


  !*--------------------------------------------------------------------------------------
  !*values->coeff
  !*--------------------------------------------------------------------------------------
  FUNCTION Cons_to_Control(u_cons,e) RESULT (u)
    ! from values at dofs, compute control points
    TYPE(Pvar), DIMENSION(:), INTENT(in):: u_cons !*values in input
    TYPE(element), INTENT(in):: e
    TYPE(Pvar),DIMENSION(SIZE(u_cons,dim=1)):: u !*coefficients as output
    INTEGER:: l, k

    u=u_cons !*He's copying the values in the output variable

    SELECT CASE(e%itype)
    CASE(1,3,4,7,11,12,13,14) ! Lagrange
       !*In this case there is nothing to do because coeff=values
    CASE(2,5,6)! cubic bezier
       !*Projection needed

       DO k=1, n_vars        !*Loop on the components of the vector u
          u(:)%u(k) = MATMUL(e%base1,u_cons%u(k))
          !*NB: u_cons%u(k) has the component k in all the DoFs, we mean u_cons(:)%u(k), it's the same

          !*RMK, from the fields of element
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
       ENDDO
    CASE default
       PRINT*, "erreur dans Model/Control_to_Cons"
       STOP
    END SELECT
    
  END FUNCTION Cons_to_Control
END MODULE Model
