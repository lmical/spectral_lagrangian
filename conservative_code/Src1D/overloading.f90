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
!*Overloading equiparates some operators/functions user defined involving data structures user defined (like PVars) to existing operators on standard data strucutres (like + on real variables)
!*(For example we can define a function to add two PVar (in some sense that we decide) and we refer to it as the operator +. Whenever we invoke Pvar+Pvar we call this function) 
!*----------------------------------------------------------------------------------------
MODULE overloading
  USE variable_def
  use precision
  IMPLICIT NONE

  !*Assignments
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE assign_var !*PVar=PVar
     MODULE PROCEDURE assign_var_vec !*Vector of PVar=Vector of PVar !*NB: There is no safety check for eventual mismatchings
     MODULE PROCEDURE assign_var_real !*PVar=Real !*All the components=to that real
  END INTERFACE ASSIGNMENT (=)

  !*Addition
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE addition_var !*PVar+PVar=PVar
     MODULE PROCEDURE addition_Vecvar !*Vector of PVar+Vector of PVar=Vector of PVar !*NB: There is no safety check for eventual mismatchings
     MODULE PROCEDURE addition_Vecvar_var !*Vector of PVar+ single PVar=Vector of PVar !*single PVar added to all the elements of the vector
     MODULE PROCEDURE addition_var_Vecvar !*single PVar+Vector of PVar=Vector of PVar !*single PVar added to all the elements of the vector
  
  END INTERFACE OPERATOR (+)

  !*Subtraction
  INTERFACE OPERATOR (-)
     MODULE PROCEDURE subtraction_var !*PVar-PVar=PVar
     MODULE PROCEDURE subtraction_Vecvar !*Vector of PVar-Vector of PVar=Vector of Pvar !*NB: There is no safety check for eventual mismatchings
     MODULE PROCEDURE subtraction_Vecvar_var !*Vector of PVar-single PVar=Vector of PVar !*single PVar subtracted to all the elements of the Vector
     MODULE PROCEDURE subtraction_var_Vecvar !*PVar-Vector of Pvar=Vector of PVar !*each component of the result is the difference between the single PVar and the corresponding element of the argument vector
  END INTERFACE OPERATOR (-)

  !*Multiplication
  INTERFACE OPERATOR (*)
     MODULE PROCEDURE multiple_var_var !*PVar*PVar=PVar !*every component is the product between the corresponding components
     MODULE PROCEDURE multiple_real_var !*real*PVar=PVar
     MODULE PROCEDURE multiple_var_real !*PVar*real=PVar
     MODULE PROCEDURE multiple_Vec_Vecvar !*Vector of PVar*Vector of real=Vector of PVar !*each component of the result is the corresponding PVar*the corresponding real !*No check for eventual mismatchings
     MODULE PROCEDURE multiple_Vecvar_Vec !*Vector of real*Vector of PVar=Vector of PVar !*each component of the result is the corresponding PVar*the corresponding real !*No check for eventual mismatchings
     MODULE PROCEDURE multiple_Vecvar_real !*Vector of PVar*real=Vector of PVar
     MODULE PROCEDURE multiple_real_Vecvar !*real*Vector of PVar=Vector of PVar
     MODULE PROCEDURE multiple_var_realvec !*PVar*Vector of real=Vector of Pvar !*every element of the result is the corresponding real*PVar
     MODULE PROCEDURE multiple_realvec_var !*Vector of real*PVar=Vector of Pvar !*every element of the result is the corresponding real*PVar
     MODULE PROCEDURE multiple_tensor_vec !*tensor(N_var,N_vars,N_dim)*vec(N_dim)=matrix(N_vars,N_vars)
     !*Typical application: normal Jacobian
     !*df/du(N_var,N_vars,N_dim)*n(N_dim)=df/du*n(N_vars,N_vars)
     MODULE PROCEDURE multiple_vec_tensor !*vec(N_dim)*tensor(N_var,N_vars,N_dim)=matrix(N_vars,N_vars)
     !*Typical application: normal Jacobian
     !*n(N_dim)*df/du(N_var,N_vars,N_dim)=df/du*n(N_vars,N_vars)
     MODULE PROCEDURE multiple_matrix_var !*matrix*PVar=PVar !*PVar is treated as a vector and the resulting PVar is the standard matrix*vector !*No check for eventual mismatching
  END INTERFACE OPERATOR (*)

  !*Division
  INTERFACE OPERATOR (/)
     MODULE PROCEDURE divide_var_real !*PVar/real=Pvar
     MODULE PROCEDURE divide_Vecvar_Vecreal !*Vector of Pvar/Vector of real=Vector of Pvar !*each element od the result is the corresponding element of the PVar vector divided by the corresponding element of the real vector !*No check for eventual mismatchings
  END INTERFACE OPERATOR (/)

  !*Sum
  INTERFACE SUM
     MODULE PROCEDURE pvar_sum   !*SUM(Vector of Pvar)=Pvar !*We sum component by component all the elements of the vector of Pvar
  END INTERFACE SUM



CONTAINS

  !*--------------------------------------------------------------------------------------
  !*PVar=PVar
  !*--------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE assign_var( var_result, var)
    TYPE(PVar), INTENT(in)  :: var
    TYPE(PVar), INTENT(out) :: var_result

    var_result%NVars = var%NVars
    var_result%u     = var%u

  END SUBROUTINE assign_var

  !*--------------------------------------------------------------------------------------
  !*Vector of PVar=Vector of PVar
  !*--------------------------------------------------------------------------------------
  SUBROUTINE assign_var_vec( var_result, var)
    TYPE(PVar), DIMENSION(:),INTENT(in)  :: var
    TYPE(PVar), DIMENSION(:),INTENT(out) :: var_result
    INTEGER:: i

    DO i=1, SIZE(Var,dim=1) !*Loop on the elements of the vector of PVars
       !*Single assignment
       var_result(i)%NVars = var(i)%NVars
       var_result(i)%u     = var(i)%u
    ENDDO
    !*NB: There is no safety check for eventual mismatchings
  END SUBROUTINE assign_var_vec

  !*--------------------------------------------------------------------------------------
  !*PVar=Real
  !*--------------------------------------------------------------------------------------  
  ELEMENTAL  SUBROUTINE assign_var_real( var_result, alpha)
    REAL(dp),       INTENT(IN)  :: alpha
    TYPE(PVar), INTENT(out) :: var_result

    var_result%u(:) = alpha !*All the components=to that real

  END SUBROUTINE assign_var_real

  !*--------------------------------------------------------------------------------------
  !*PVar+PVar=PVar
  !*--------------------------------------------------------------------------------------
  FUNCTION addition_var(var1, var2) RESULT(add_var)
    TYPE(PVar), INTENT(IN) :: var1
    TYPE(PVar), INTENT(IN) :: var2
    TYPE(PVar)             :: add_var

    add_var%u = var1%u + var2%u

  END FUNCTION addition_var

  !*--------------------------------------------------------------------------------------
  !*Vector of PVar+Vector of PVar=Vector of PVar
  !*--------------------------------------------------------------------------------------
  FUNCTION addition_Vecvar(var1, var2) RESULT(add_var)
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var1
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var2
    TYPE(PVar), DIMENSION(SIZE(var1))    :: add_var
    INTEGER:: i

    DO i=1, SIZE(var1) !*Loop on the elements of the vector of PVar 
       add_var(i)%u(:) = var1(i)%u(:) + var2(i)%u(:)
    ENDDO
    !*NB: There is no safety check for eventual mismatchings
  END FUNCTION addition_Vecvar
  
  !*--------------------------------------------------------------------------------------
  !*Vector of PVar+ single PVar=Vector of PVar (single PVar added to all the elements of the vector)
  !*--------------------------------------------------------------------------------------
  FUNCTION addition_Vecvar_var(var1, var2) RESULT(add_vecvar_var)
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var1
    TYPE(PVar),               INTENT(IN) :: var2
    TYPE(PVar), DIMENSION(SIZE(var1))    :: add_vecvar_var
    INTEGER:: l

    DO l=1, SIZE(var1) !*Loop on the elements of the vector
       add_vecvar_var(l)%u = var1(l)%u + var2%u
    ENDDO
  END FUNCTION addition_Vecvar_var

  !*--------------------------------------------------------------------------------------
  !*single PVar+Vector of PVar=Vector of PVar (single PVar added to all the elements of the vector)
  !*--------------------------------------------------------------------------------------
  FUNCTION addition_var_Vecvar(var1, var2) RESULT(add_vecvar_var)
    TYPE(PVar),               INTENT(IN) :: var1
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var2
    TYPE(PVar), DIMENSION(SIZE(var2))    :: add_vecvar_var
    INTEGER:: l

    DO l=1, SIZE(var2) !*Loop on the elements of the vector
       add_vecvar_var(l)%u = var1%u + var2(l)%u
    ENDDO
  END FUNCTION addition_var_Vecvar

  !*--------------------------------------------------------------------------------------
  !*PVar-PVar=PVar
  !*--------------------------------------------------------------------------------------
  FUNCTION subtraction_var(var1, var2) RESULT(subtract_var)
    TYPE(PVar), INTENT(IN) :: var1
    TYPE(PVar), INTENT(IN) :: var2
    TYPE(PVar)             :: subtract_var

    subtract_var%u = var1%u - var2%u

  END FUNCTION subtraction_var

  !*--------------------------------------------------------------------------------------
  !*Vector of PVar-Vector of PVar=Vector of Pvar
  !*--------------------------------------------------------------------------------------
  FUNCTION subtraction_Vecvar(var1, var2) RESULT(subtract_var)
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var1
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var2
    TYPE(PVar), DIMENSION(SIZE(var1))    :: subtract_var
    INTEGER:: l

    DO l=1, SIZE(var1)
       subtract_var(l)%u = var1(l)%u - var2(l)%u
    ENDDO
  END FUNCTION subtraction_Vecvar

  !*--------------------------------------------------------------------------------------
  !*Vector of PVar-single PVar=Vector of PVar (single PVar subtracted to all the elements of the Vector)
  !*--------------------------------------------------------------------------------------
  FUNCTION subtraction_Vecvar_var(var1, var2) RESULT(subtract_vecvar_var)
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var1
    TYPE(PVar),               INTENT(IN) :: var2
    TYPE(PVar), DIMENSION(SIZE(var1))    :: subtract_vecvar_var
    INTEGER:: l

    DO l=1, SIZE(var1)
       subtract_vecvar_var(l)%u = var1(l)%u - var2%u
    ENDDO
  END FUNCTION subtraction_Vecvar_var

  !*--------------------------------------------------------------------------------------
  !*PVar-Vector of Pvar=Vector of PVar (each component of the result is the difference between the single PVar and the corresponding element of the argument vector)
  !*--------------------------------------------------------------------------------------
  FUNCTION subtraction_var_Vecvar(var1, var2) RESULT(subtract_vecvar_var)
    TYPE(PVar),               INTENT(IN) :: var1
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var2
    TYPE(PVar), DIMENSION(SIZE(var2))    :: subtract_vecvar_var
    INTEGER:: l
    DO l=1, SIZE(var2)
       subtract_vecvar_var(l)%u = var1%u - var2(l)%u
    ENDDO
  END FUNCTION subtraction_var_Vecvar

  !*--------------------------------------------------------------------------------------
  !*PVar*PVar=PVar (every component is the product between the corresponding components)
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_var_var(var1, var2) RESULT(mult_var_var)
    TYPE(PVar), INTENT(IN) :: var1
    TYPE(PVar), INTENT(IN) :: var2
    TYPE(PVar)             :: mult_var_var

    mult_var_var%u = var1%u * var2%u

  END FUNCTION multiple_var_var

  !*--------------------------------------------------------------------------------------
  !*real*PVar=PVar
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_real_var(alpha, var) RESULT(mult_real_var)
    REAL(dp),       INTENT(IN) :: alpha
    TYPE(PVar), INTENT(IN) :: var
    TYPE(PVar)             :: mult_real_var

    mult_real_var%u(:) = alpha * var%u(:)

  END FUNCTION multiple_real_var

  !*--------------------------------------------------------------------------------------
  !*PVar*real=PVar
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_var_real(var,alpha) RESULT(mult_var_real)
    REAL(dp),       INTENT(IN) :: alpha
    TYPE(PVar), INTENT(IN) :: var
    TYPE(PVar)             :: mult_var_real

    mult_var_real%u(:) = alpha * var%u(:)

  END FUNCTION multiple_var_real

  !*--------------------------------------------------------------------------------------
  !*Vector of PVar*Vector of real=Vector of PVar (each component of the result is the corresponding PVar*the corresponding real)
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_Vec_Vecvar(var1, var2) RESULT(mult_var_var)
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var1
    REAL(dp),       DIMENSION(:), INTENT(IN) :: var2
    TYPE(PVar), DIMENSION(SIZE(var1))    :: mult_var_var
    INTEGER:: i

    DO i=1, SIZE(Var1)
       mult_var_var(i)%u = var1(i)%u * var2(i)
    ENDDO
    !*No check for eventual mismatchings
  END FUNCTION multiple_Vec_Vecvar

  !*--------------------------------------------------------------------------------------
  !*Vector of real*Vector of PVar=Vector of PVar (each component of the result is the corresponding PVar*the corresponding real)
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_Vecvar_Vec(var1, var2) RESULT(mult_var_var)
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var2
    REAL(dp),       DIMENSION(:), INTENT(IN) :: var1
    TYPE(PVar), DIMENSION(SIZE(var1))    :: mult_var_var
    INTEGER:: i

    DO i=1, SIZE(Var1)
       mult_var_var(i)%u = var2(i)%u * var1(i)
    ENDDO
    !*No check for eventual mismatchings
  END FUNCTION multiple_Vecvar_Vec

  !*--------------------------------------------------------------------------------------
  !*Vector of PVar*real=Vector of PVar
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_Vecvar_real(var,alpha) RESULT(mult_var_real)
    REAL(dp),       INTENT(IN) :: alpha
    TYPE(PVar), DIMENSION(:),INTENT(IN) :: var
    TYPE(PVar),DIMENSION(SIZE(var))             :: mult_var_real
    INTEGER:: l
    DO l=1, SIZE(var)
       mult_var_real(l)%u(:) = alpha * var(l)%u(:)
    ENDDO

  END FUNCTION multiple_Vecvar_real

  !*--------------------------------------------------------------------------------------
  !*real*Vector of PVar=Vector of PVar 
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_real_Vecvar(alpha,var) RESULT(mult_var_real)
    REAL(dp),       INTENT(IN) :: alpha
    TYPE(PVar), DIMENSION(:),INTENT(IN) :: var
    TYPE(PVar),DIMENSION(SIZE(var))             :: mult_var_real
    INTEGER:: l
    DO l=1, SIZE(var)
       mult_var_real(l)%u(:) = alpha * var(l)%u(:)
    ENDDO

  END FUNCTION multiple_real_Vecvar

  !*--------------------------------------------------------------------------------------
  !*PVar*Vector of real=Vector of Pvar (every element of the result is the corresponding real*PVar)
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_var_realvec(var,alpha) RESULT(mult_var_vecreal)
    REAL(dp),    DIMENSION(:), INTENT(IN) :: alpha
    TYPE(PVar),               INTENT(IN) :: var
    TYPE(PVar), DIMENSION(SIZE(alpha))   :: mult_var_vecreal
    INTEGER :: l
    DO l=1, SIZE(alpha)
       mult_var_vecreal(l)%u(:) = alpha(l) * var%u(:)
    ENDDO

  END FUNCTION multiple_var_realvec

  !*--------------------------------------------------------------------------------------
  !*Vector of real*PVar=Vector of Pvar (every element of the result is the corresponding real*PVar)
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_realvec_var(alpha,var) RESULT(mult_var_vecreal)
    REAL(dp),    DIMENSION(:), INTENT(IN) :: alpha
    TYPE(PVar),               INTENT(IN) :: var
    TYPE(PVar), DIMENSION(SIZE(alpha))     :: mult_var_vecreal
    INTEGER :: l
    DO l=1, SIZE(alpha)
       mult_var_vecreal(l)%u(:) = alpha(l) * var%u(:)
    ENDDO

  END FUNCTION multiple_realvec_var

  !*--------------------------------------------------------------------------------------
  !*tensor(N_var,N_vars,N_dim)*vec(N_dim)=matrix(N_vars,N_vars)
  !*Typical application: normal Jacobian
  !*df/du(N_var,N_vars,N_dim)*n(N_dim)=df/du*n(N_vars,N_vars)
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_tensor_vec(tensor, vec) RESULT (matrix)
    REAL, DIMENSION(N_vars, N_vars, N_dim), INTENT(IN) :: tensor
    REAL, DIMENSION(n_dim),                 INTENT(IN) :: vec
    REAL, DIMENSION(n_vars,n_vars)                     :: matrix
    INTEGER:: i

    matrix=0.0_dp
    DO i=1, N_dim
       matrix(:,:)= matrix(:,:) + vec(i) * tensor(:,:,i)
    ENDDO

  END FUNCTION multiple_tensor_vec

  !*--------------------------------------------------------------------------------------
  !*vec(N_dim)*tensor(N_var,N_vars,N_dim)=matrix(N_vars,N_vars)
  !*Typical application: normal Jacobian
  !*n(N_dim)*df/du(N_var,N_vars,N_dim)=df/du*n(N_vars,N_vars)
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_vec_tensor(vec, tensor) RESULT (matrix)
    REAL(dp), DIMENSION(N_vars, N_vars, N_dim), INTENT(IN) :: tensor
    REAL(dp), DIMENSION(n_dim),                 INTENT(IN) :: vec
    REAL(dp), DIMENSION(n_vars,n_vars)                     :: matrix
    INTEGER:: i

    matrix=0.0_dp
    DO i=1, N_dim
       matrix(:,:)= matrix(:,:)+vec(i) * tensor(:,:,i)
    ENDDO

  END FUNCTION multiple_vec_tensor



  !*--------------------------------------------------------------------------------------
  !*matrix*PVar=PVar (PVar is treated as a vector and the resulting PVar is the standard matrix*vector)
  !*--------------------------------------------------------------------------------------
  FUNCTION multiple_matrix_var(A,var1) RESULT(var2)
    REAL, INTENT(in), DIMENSION(:,:):: A
    TYPE(PVar), INTENT(in):: var1
    TYPE(PVar):: var2

    var2%u=MATMUL(A,var1%u)
    !*No check for eventual mismatching
  END FUNCTION multiple_matrix_var

  !*--------------------------------------------------------------------------------------
  !*PVar/real=Pvar
  !*--------------------------------------------------------------------------------------
  FUNCTION divide_var_real(var,alpha) RESULT(div_var_real)
    REAL(dp),       INTENT(IN) :: alpha
    TYPE(PVar), INTENT(IN) :: var
    TYPE(PVar)             :: div_var_real

    div_var_real%u(:) =  var%u(:)/alpha

  END FUNCTION divide_var_real

  !*--------------------------------------------------------------------------------------
  !*Vector of Pvar/Vector of real=Vector of Pvar (each element od the result is the corresponding element of the PVar vector divided by the corresponding element of the real vector)
  !*--------------------------------------------------------------------------------------
  FUNCTION divide_Vecvar_Vecreal(var,alpha) RESULT(div_var_real)
    REAL(dp),       DIMENSION(:), INTENT(IN) :: alpha
    TYPE(PVar), DIMENSION(:), INTENT(IN) :: var
    TYPE(PVar), DIMENSION(SIZE(var))     :: div_var_real
    INTEGER:: i

    DO i=1, SIZE(Var)
       div_var_real(i)%u(:) =  var(i)%u(:)/alpha(i)
    ENDDO
    !*No check for eventual mismatching
  END FUNCTION divide_Vecvar_Vecreal

  !*--------------------------------------------------------------------------------------
  !*SUM(Vector of Pvar)=Pvar (we sum component by component all the elements of the vector of Pvar)
  !*--------------------------------------------------------------------------------------  
  FUNCTION pvar_sum(pvar_array)
    TYPE(PVar), DIMENSION(:), INTENT(in) :: pvar_array
    TYPE(PVar)                           :: pvar_sum
    INTEGER :: i

    pvar_sum%u(:)=0.0_dp

    DO i=1,SIZE(pvar_array)
       pvar_sum%u(:) = pvar_sum%u(:) + pvar_array(i)%u(:)
    END DO

  END FUNCTION pvar_sum

END MODULE overloading
