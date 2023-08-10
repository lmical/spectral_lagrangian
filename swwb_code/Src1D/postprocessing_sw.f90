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
!*Module used to write the solution at a certain timestep so that it can be plotted
!*In case there is an exact solution the errors are computed here 
!*----------------------------------------------------------------------------------------
MODULE postprocessing

  ! This module is intended to output the postprocessing files with the data to generate plots.
  ! The data can in general easily read via gnuplot or xmgrace (as example)
  ! In the data files are printed ALL the DOFs of the domain

  USE param2d
  USE init_bc
  USE Model
  USE variable_def
  USE PRECISION
  USE preprocessing
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: visu, descriptor, errorTestProjection, errorTestProjectionEndPoints, MyPrintFinalSolution!, visudeb !,errors

CONTAINS

  !*--------------------------------------------------------------------------------------
  !*visu prints the numerical solution at a certain timestep
  !*--------------------------------------------------------------------------------------
  SUBROUTINE visu(DATA, kt, mesh, var,folder)
    TYPE(Donnees), INTENT(inout):: DATA
    INTEGER, INTENT(in):: kt
    TYPE(maillage), INTENT(in):: mesh
    TYPE(variables), INTENT(in):: var

    CHARACTER(len=100), INTENT(in)::  folder
    INTEGER::jt, l, k, ios, lp, is
    TYPE(element):: e
    TYPE(Pvar),DIMENSION(:),ALLOCATABLE::v,u,w
    REAL(dp), DIMENSION(:,:), ALLOCATABLE:: x
    TYPE(PVar), DIMENSION(:), ALLOCATABLE:: sol
    CHARACTER(LEN = 75) :: nombre, nt_str

    !*Print the actual iteration number and the number of elements in the shell
    WRITE(nombre, *) kt
    WRITE(nt_str, *) Mesh%nt
    !*The SYSTEM subroutine executes a system command as if from the command line.
    call system('mkdir -p ' //folder )

    !*We open the files (name.dat) on which we will write the solution
    OPEN( UNIT = 1, FILE = trim(FOLDER)//"/"//"/sol_cons_last.dat", IOSTAT = ios )
    OPEN( UNIT = 2, FILE = trim(FOLDER)//"/"//"/ex_prim_last.dat", IOSTAT = ios )
    OPEN( UNIT = 3, FILE = trim(FOLDER)//"/"//"/error_last.dat", IOSTAT = ios )
    OPEN( UNIT = 4, FILE = trim(FOLDER)//"/"//"/sol_prim_last.dat", IOSTAT = ios )
    OPEN( UNIT = 80, FILE = trim(FOLDER)//"/"//"/sol_prim_nodal_last.dat", IOSTAT = ios )
    OPEN( UNIT = 5, FILE = trim(FOLDER)//"/"//"/sol_prim_"//TRIM(ADJUSTL(nombre))//".dat", IOSTAT = ios )
    OPEN( UNIT = 50, FILE = trim(FOLDER)//"/"//"/ex_prim_"//TRIM(ADJUSTL(nombre))//".dat", IOSTAT = ios )
    OPEN( UNIT = 7, FILE = trim(FOLDER)//"/"//"/FluxIndicator_"//TRIM(ADJUSTL(nt_str))//"_"//TRIM(ADJUSTL(nombre))//".dat", IOSTAT = ios )
    OPEN( UNIT = 9, FILE = trim(FOLDER)//"/"//"/DiagIndicator_"//TRIM(ADJUSTL(nt_str))//"_"//TRIM(ADJUSTL(nombre))//".dat", IOSTAT = ios )
    OPEN( UNIT = 60, FILE = trim(FOLDER)//"/"//"/DiagIndicator2_"//TRIM(ADJUSTL(nt_str))//"_"//TRIM(ADJUSTL(nombre))//".dat", IOSTAT = ios )


    OPEN( UNIT = 70, FILE = trim(FOLDER)//"/"//"/FluxIndicator_last.dat", IOSTAT = ios )
    OPEN( UNIT = 90, FILE = trim(FOLDER)//"/"//"/DiagIndicator_last.dat", IOSTAT = ios )




    OPEN( UNIT = 8, FILE = trim(FOLDER)//"/"//"/Entropy_"//TRIM(ADJUSTL(nt_str))//".dat", IOSTAT = ios )
!!$
!!$    
!!$
    !*The REWIND statement positions the file associated with the specified unit to its initial point. It means that we will write from the first line
    REWIND(1); REWIND(2); REWIND(3); REWIND(4); REWIND(5);REWIND(7); REWIND(8); REWIND(9);REWIND(70);REWIND(90);REWIND(80);REWIND(50);REWIND(60)
!!$
    !*Loop on the elements
    DO jt=1, Mesh%nt
       e=Mesh%e(jt) !*Quick reference to the processed element
       ALLOCATE(sol(e%nsommets),u(e%nsommets),v(e%nsommets),w(e%nsommets))

       DO l=1, e%nsommets !*Loop on the DoFs of the element
          sol(l)=0._DP !*In case it is known, this is the exact solution which one can use to have feedbacks or to make convergence analyses. For example in the linear advection with speed a=1 we have sol(l)IC(e%coor(l)-DATA%temps,DATA)
          u(l)=Var%ua(e%nu(l),0) !*Copying the coefficients of the DoFs of the processed triangle in u
       ENDDO
       DO l=1, e%nsommets !*Loop on the DoFs
          v(l)=e%eval_func(u,e%x(:,l)) !*Value of the solution in that DoF (we only had the coefficients so we are now storing in v the values)
#if(1==0)
          w(l)=convert_cons2prim(v(l)) !*We are now passing from conserved to primitive variables
#else
          !*Conserved variables
          w(l)=v(l)
#if(2==2)
          w(l)%u(1)=w(l)%u(1)+BathymetryAtDoFsUtils(e%nu(l)) !*BETTER BECAUSE SOMETIMES YOU HAVE THE DISCRETE BATHYMETRY ONLY
#else
          w(l)%u(1)=w(l)%u(1)+bathymetry(DATA%test,e%coor(l)) 
#endif
#endif
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, BathymetryAtDoFsUtils(e%nu(l))-bathymetry(DATA%test,e%coor(l))
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       ENDDO

       !*Writing of the needed information
       !*NB: In some files we write the same information.
       !*This is due to the fact that some of them are overwritten every time (the "last" files), while some are written for the first time 
!!$
       WRITE(7,*) SUM(e%coor)/REAL(e%nsommets,DP), Mesh%e(jt)%type_flux
       WRITE(9,*) SUM(e%coor)/REAL(e%nsommets,DP), Mesh%e(jt)%diag
       WRITE(60,*) SUM(e%coor)/REAL(e%nsommets,DP), Mesh%e(jt)%diag2
  
       WRITE(70,*) SUM(e%coor)/REAL(e%nsommets,DP), Mesh%e(jt)%type_flux
       WRITE(90,*) SUM(e%coor)/REAL(e%nsommets,DP), Mesh%e(jt)%diag
!!$
       SELECT CASE(e%itype) !*Choice on the type of elements
       CASE(1,11)
          !*Order 2, 2 DoFs per element 
          !*Local DoFs with increasing abscissa: 1, 2

          !*With increasing abscissa of the DoF of the element we write 
          !*in file 1
          !*coordinate, value in the DoF, coefficients in the DoF
          WRITE(1,*) e%coor(1),v(1)%u, u(1)%u
          WRITE(1,*) e%coor(2),v(2)%u, u(2)%u
          !*in file 2
          !*coordinate, value of the exact solution in the DoF
          WRITE(2,*) e%coor(1),sol(1)%u
          WRITE(2,*) e%coor(2),sol(2)%u
          !*in file 3
          !*coordinate, error of the numerical solution w.r.t. the exact solution in the DoF
          WRITE(3,*) e%coor(1),sol(1)%u-v(1)%u
          WRITE(3,*) e%coor(2),sol(2)%u-v(2)%u
          !*in file 4
          !*coordinate, values of the primitive variables in the DoF
          WRITE(4,*) e%coor(1),w(1)%u
          WRITE(4,*) e%coor(2),w(2)%u
          !*in file 5
          !*coordinate, values of the primitive variables in the DoF
          WRITE(5,*) e%coor(1),w(1)%u
          WRITE(5,*) e%coor(2),w(2)%u
          !*in file 50
          !*coordinate, value of the exact solution in the DoF
          WRITE(50,*) e%coor(1),sol(1)%u
          WRITE(50,*) e%coor(2),sol(2)%u
          !*in file 80
          !*coordinate, values of the primitive variables in the DoF
          WRITE(80,*) e%coor(1),w(1)%u
          WRITE(80,*) e%coor(2),w(2)%u
          !*It's the same in the other cases   
       CASE(2,3,12)
          !*Order 3, 3 DoFs per element
          !*Local DoFs with increasing abscissa: 1, 3, 2
          WRITE(1,*) e%coor(1),v(1)%u, u(1)%u
          WRITE(1,*) e%coor(3),v(3)%u, u(3)%u
          WRITE(1,*) e%coor(2),v(2)%u, u(2)%u
          WRITE(2,*) e%coor(1),sol(1)%u
          WRITE(2,*) e%coor(3),sol(3)%u
          WRITE(2,*) e%coor(2),sol(2)%u
          WRITE(3,*) e%coor(1),sol(1)%u-v(1)%u
          WRITE(3,*) e%coor(3),sol(3)%u-v(3)%u
          WRITE(3,*) e%coor(2),sol(2)%u-v(2)%u
          WRITE(4,*) e%coor(1),w(1)%u
          WRITE(4,*) e%coor(3),w(3)%u
          WRITE(4,*) e%coor(2),w(2)%u
          WRITE(5,*) e%coor(1),w(1)%u
          WRITE(5,*) e%coor(3),w(3)%u
          WRITE(5,*) e%coor(2),w(2)%u
          WRITE(50,*) e%coor(1),sol(1)%u
          WRITE(50,*) e%coor(3),sol(3)%u
          WRITE(50,*) e%coor(2),sol(2)%u
          WRITE(80,*) e%coor(1),w(1)%u
          WRITE(80,*) e%coor(2),w(2)%u

       CASE(4,5,13)
          !*Order 4, 4 DoFs per element
          !*Local DoFs with increasing abscissa: 1, 3, 4, 2
          WRITE(1,*) e%coor(1),v(1)%u, u(1)%u
          WRITE(1,*) e%coor(3),v(3)%u, u(3)%u
          WRITE(1,*) e%coor(4),v(4)%u, u(4)%u
          WRITE(1,*) e%coor(2),v(2)%u, u(2)%u
          WRITE(2,*) e%coor(1),sol(1)%u
          WRITE(2,*) e%coor(3),sol(3)%u
          WRITE(2,*) e%coor(4),sol(4)%u
          WRITE(2,*) e%coor(2),sol(2)%u
          WRITE(3,*) e%coor(1),sol(1)%u-v(1)%u
          WRITE(3,*) e%coor(3),sol(3)%u-v(3)%u
          WRITE(3,*) e%coor(4),sol(4)%u-v(4)%u
          WRITE(3,*) e%coor(2),sol(2)%u-v(2)%u
          WRITE(4,*) e%coor(1),w(1)%u
          WRITE(4,*) e%coor(3),w(3)%u
          WRITE(4,*) e%coor(4),w(4)%u
          WRITE(4,*) e%coor(2),w(2)%u
          WRITE(5,*) e%coor(1),w(1)%u
          WRITE(5,*) e%coor(3),w(3)%u
          WRITE(5,*) e%coor(4),w(4)%u
          WRITE(5,*) e%coor(2),w(2)%u

          WRITE(50,*) e%coor(1),sol(1)%u
          WRITE(50,*) e%coor(3),sol(3)%u
          WRITE(50,*) e%coor(4),sol(4)%u
          WRITE(50,*) e%coor(2),sol(2)%u
          WRITE(80,*) e%coor(1),w(1)%u
          WRITE(80,*) e%coor(2),w(2)%u

       CASE(6,7,14)
          !*Order 5, 5 DoFs per element
          !*Local DoFs with increasing abscissa: 1, 3, 4, 5, 2
          WRITE(1,*) e%coor(1),v(1)%u, u(1)%u
          WRITE(1,*) e%coor(3),v(3)%u, u(3)%u
          WRITE(1,*) e%coor(4),v(4)%u, u(4)%u
          WRITE(1,*) e%coor(5),v(5)%u, u(5)%u
          WRITE(1,*) e%coor(2),v(2)%u, u(2)%u
          WRITE(2,*) e%coor(1),sol(1)%u
          WRITE(2,*) e%coor(3),sol(3)%u
          WRITE(2,*) e%coor(4),sol(4)%u
          WRITE(2,*) e%coor(5),sol(5)%u
          WRITE(2,*) e%coor(2),sol(2)%u
          WRITE(3,*) e%coor(1),sol(1)%u-v(1)%u
          WRITE(3,*) e%coor(3),sol(3)%u-v(3)%u
          WRITE(3,*) e%coor(4),sol(4)%u-v(4)%u
          WRITE(3,*) e%coor(5),sol(5)%u-v(5)%u
          WRITE(3,*) e%coor(2),sol(2)%u-v(2)%u
          WRITE(4,*) e%coor(1),w(1)%u
          WRITE(4,*) e%coor(3),w(3)%u
          WRITE(4,*) e%coor(4),w(4)%u
          WRITE(4,*) e%coor(5),w(5)%u
          WRITE(4,*) e%coor(2),w(2)%u
          WRITE(5,*) e%coor(1),w(1)%u
          WRITE(5,*) e%coor(3),w(3)%u
          WRITE(5,*) e%coor(4),w(4)%u
          WRITE(5,*) e%coor(5),w(5)%u
          WRITE(5,*) e%coor(2),w(2)%u

          WRITE(50,*) e%coor(1),sol(1)%u
          WRITE(50,*) e%coor(3),sol(3)%u
          WRITE(50,*) e%coor(4),sol(4)%u
          WRITE(50,*) e%coor(5),sol(5)%u
          WRITE(50,*) e%coor(2),sol(2)%u
          WRITE(80,*) e%coor(1),w(1)%u
          WRITE(80,*) e%coor(2),w(2)%u
       END SELECT
       DEALLOCATE(sol,v,u,w)
    ENDDO

    !*Close all the files where we wrote our stuff
    CLOSE(1);CLOSE(2);CLOSE(3); CLOSE(4);CLOSE(5);CLOSE(7);CLOSE(8);CLOSE(9);CLOSE(70); CLOSE(90);CLOSE(80); CLOSE(50); close(60);

  END SUBROUTINE visu

  !*--------------------------------------------------------------------------------------
  !*descriptor creates the folder TEST where to store the results
  !*--------------------------------------------------------------------------------------
  SUBROUTINE descriptor(folder,data)
      character(len=100)::folder
      character(len=100)::filename
      TYPE(donnees):: DATA
      call system('mkdir -p ' //folder )
      filename = '/Descriptor'
      open(unit=10,status="REPLACE",file= trim(folder)//trim(filename)//'.dat')
      write(10,*)  data%nt, data%itype,data%ischema,data%test

      !write(10,*) 'itype = ', data%itype
      !write(10,*) 'ischema = ', data%ischema
      !write(10,*) 'test = ', data%test
      close(10)
   END SUBROUTINE descriptor





 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !*errorTestProjection calculates the difference between my solution and the exact solution and its L^1,L^2 and L^\infty norm
 !*- In particular we start from the coefficients of our solution as an input FinalSolutionElliptic NB:They are already in the initial numeration -> CoeffNumSol=FinalSolutionElliptic
 !*- We calculate the values of the analytical exact solution -> ValExactSol
 !*- We pass to the coefficients of the analytical exact solution -> CoeffExactSol
 !*- We make the difference between the coefficients of our solution and the analytical one and we have the difference -> CoeffDiff
 !*We properly integrate it to get the L^1,L^2 norms
 !*The L^\infty norm is made directly on the coefficients CoeffDiff of the difference
 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE errorTestProjection(error1,error2,errorinf,CoeffNumSol,Mesh,DATA)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(OUT) :: error1,error2,errorinf !*errors L^1,L^2,L^infty
    TYPE(Pvar), DIMENSION(:), INTENT(IN) :: CoeffNumSol !*Coefficients of our solution to the elliptic problem in the original numeration
    TYPE (maillage), INTENT(IN) :: Mesh 
    TYPE(donnees), INTENT(INOUT) :: DATA

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

IF(DATA%test==12 .OR. DATA%test==500 .OR. DATA%test==501 &
  & .OR. DATA%test==100 .OR. DATA%test==101 .OR. DATA%test==111 &
  & .OR. DATA%test==120 .OR. DATA%test==121 .OR. DATA%test==122 &
  & .OR. (DATA%test==300 .AND. DATA%perturbation==0) .OR. DATA%test==210 &
  & .OR. DATA%test==150 .OR. (DATA%test==400 .AND. DATA%perturbation==0) ) THEN

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
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, indi, e%nu(indi), e%coor(indi), exactsolution(e%coor(indi)), suppreal
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ValExactSol(e%nu(indi))=ExactSolutionHyperbolic(e%coor(indi),DATA,e%nu(indi)) !*PVar=real
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

    PRINT*, "error1", error1
    PRINT*, "error2", error2
    PRINT*, "errorinf", errorinf

    CALL printerrors(DATA, Mesh, error1, error2, errorinf)

END IF

    CONTAINS

    FUNCTION ExactSolutionHyperbolic(x,DATA,GlobalIndex) RESULT(sol) 
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x !*Cartesian coordinates
       TYPE(PVar) :: sol
       TYPE(donnees), INTENT(INOUT) :: DATA
       INTEGER, INTENT(IN) :: GlobalIndex 
       
       SELECT CASE(DATA%test) 
          CASE(12,500,501,100,101,300,210,150,400)
             sol=IC(x,DATA)
          CASE(111,120,121,122)
             sol=ValICUtils(GlobalIndex)
          CASE default
             PRINT*, "ERROR CASE NOT DEFINED IN ExactSolutionHyperbolic in postprocessing sw"
             STOP
       END SELECT



    END FUNCTION ExactSolutionHyperbolic
    
 END SUBROUTINE errorTestProjection



 !*---------------------------------------------------------------------------------------
 !*This function is used to print the errors in a data file
 !*---------------------------------------------------------------------------------------
 SUBROUTINE printerrors(DATA, Mesh, error_L1, error_L2, error_linf) 
    IMPLICIT NONE
    TYPE(donnees),  INTENT(in)  :: DATA
    TYPE(maillage), INTENT(in)  :: Mesh
    REAL(DP), DIMENSION(:), INTENT(in) :: error_L1, error_L2, error_Linf
    CHARACTER(LEN=255)   :: FILENAME

    IF(DATA%test==12 .OR. DATA%test==500 .OR. DATA%test==501 &
      & .OR. DATA%test==100 .OR. DATA%test==101 .OR. DATA%test==111 &
      & .OR. DATA%test==120 .OR. DATA%test==121 .OR. DATA%test==122 &
      & .OR. (DATA%test==300 .AND. DATA%perturbation==0) .OR. DATA%test==210 &
      & .OR. DATA%test==150 .OR. (DATA%test==400 .AND. DATA%perturbation==0) ) THEN


       FILENAME = 'ErrorL1_XXXX.dat'
       WRITE(FILENAME(9:12),FMT="(I4.4)") Mesh%nt   
       OPEN(666,FILE=TRIM(FILENAME))
       WRITE(666,*) Mesh%nt,  error_L1(1:n_vars)
       WRITE(666,*) "number of elements, errors in l1 norm"
       CLOSE(666)

       FILENAME = 'ErrorL2_XXXX.dat'
       WRITE(FILENAME(9:12),FMT="(I4.4)") Mesh%nt   
       OPEN(667,FILE=TRIM(FILENAME))
       WRITE(667,*) Mesh%nt,  error_L2(1:n_vars)
       WRITE(667,*) "number of elements, errors in l2 norm"
       CLOSE(667)

       FILENAME = 'ErrorLinf_XXXX.dat'
       WRITE(FILENAME(11:14),FMT="(I4.4)") Mesh%nt   
       OPEN(668,FILE=TRIM(FILENAME))
       WRITE(668,*) Mesh%nt,  error_linf(1:n_vars)
       WRITE(668,*) "number of elements, errors in linf norm"
       CLOSE(668)


    END IF
 END SUBROUTINE printerrors



 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !*errorTestProjectionEndPoints calculates the L^1,L^2 and L^\infty norms of the error keeping into account the ending points of the segments only
 !*NB: the end points are already values for our basis functions. 
 !*So no need to pass from coefficients to values
 !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE errorTestProjectionEndPoints(error1EndPoints,error2EndPoints,errorinfEndPoints,CoeffNumSol,Mesh,DATA)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(OUT) :: error1EndPoints,error2EndPoints,errorinfEndPoints !*errors L^1,L^2,L^infty
    TYPE(Pvar), DIMENSION(:), INTENT(IN) :: CoeffNumSol !*Coefficients of our solution to the elliptic problem in the original numeration
    !*NB: The end points are already the values
    TYPE (maillage), INTENT(IN) :: Mesh 
    TYPE(donnees), INTENT(INOUT) :: DATA

    TYPE(PVar)  :: ValExactSol !*Value of the exact solution in an end point
    TYPE(PVar)  :: ValDiff     !*Value of the difference (Error) in an end point

    !*Loops variables used in the integration of the difference 
    INTEGER :: indt=0 !*Index for the loop on the triangles
    INTEGER :: indi=0 !*Index for the first loop on the DoFs i
    INTEGER :: indc=0 !*Index for the component

    TYPE (element) :: e !*For the quick reference to the generic Mesh%e(indt) in the loops


    INTEGER :: component, numbofcomp

IF(DATA%test==12 .OR. DATA%test==500 .OR. DATA%test==501 &
  & .OR. DATA%test==100 .OR. DATA%test==101 .OR. DATA%test==111 &
  & .OR. DATA%test==120 .OR. DATA%test==121 .OR. DATA%test==122 &
  & .OR. (DATA%test==300 .AND. DATA%perturbation==0) .OR. DATA%test==210 &
  & .OR. DATA%test==150 .OR. (DATA%test==400 .AND. DATA%perturbation==0) ) THEN

   error1EndPoints=0._DP !*L^1
   error2EndPoints=0._DP !*L^2
   errorinfEndPoints=0._DP !*L^inf

   DO indt=1,Mesh%nt !*Loop on the elements
      e=Mesh%e(indt)
      DO indi=1,2 !*Loop on the end points
         ValExactSol=ExactSolutionHyperbolic(e%coor(indi),DATA,e%nu(indi))
         DO indc=1,n_vars !*For each component
            error1EndPoints(indc)=error1EndPoints(indc)+ABS(CoeffNumSol(e%nu(indi))%u(indc)-ValExactSol%u(indc))*e%volume*0.5_DP
            error2EndPoints(indc)=error2EndPoints(indc)+ABS(CoeffNumSol(e%nu(indi))%u(indc)-ValExactSol%u(indc))**2*e%volume*0.5_DP
            errorinfEndPoints(indc)=MAX(errorinfEndPoints(indc),ABS(CoeffNumSol(e%nu(indi))%u(indc)-ValExactSol%u(indc)))         
         END DO
      END DO
   END DO
   error2EndPoints=SQRT(error2EndPoints) !*sqrt on error 2 to make the L^2 norm


   PRINT*, "error1EndPoints", error1EndPoints
   PRINT*, "error2EndPoints", error2EndPoints
   PRINT*, "errorinfEndPoints", errorinfEndPoints

    CALL printerrorsEndPoints(DATA, Mesh, error1EndPoints, error2EndPoints, errorinfEndPoints)

END IF

    CONTAINS

    FUNCTION ExactSolutionHyperbolic(x,DATA,GlobalIndex) RESULT(sol) 
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x !*Cartesian coordinates
       TYPE(PVar) :: sol
       TYPE(donnees), INTENT(INOUT) :: DATA
       INTEGER, INTENT(IN) :: GlobalIndex 
       
       SELECT CASE(DATA%test) 
          CASE(12,500,501,100,101,300,210,150,400)
             sol=IC(x,DATA)
          CASE(111,120,121,122)
             sol=ValICUtils(GlobalIndex)
          CASE default
             PRINT*, "ERROR CASE NOT DEFINED IN ExactSolutionHyperbolic in postprocessing sw"
             STOP
       END SELECT

    END FUNCTION ExactSolutionHyperbolic

 END SUBROUTINE errorTestProjectionEndPoints



 !*---------------------------------------------------------------------------------------
 !*This function is used to print the errors in a data file
 !*---------------------------------------------------------------------------------------
 SUBROUTINE printerrorsEndPoints(DATA, Mesh, error_L1, error_L2, error_linf) 
    IMPLICIT NONE
    TYPE(donnees),  INTENT(in)  :: DATA
    TYPE(maillage), INTENT(in)  :: Mesh
    REAL(DP), DIMENSION(:), INTENT(in) :: error_L1, error_L2, error_Linf
    CHARACTER(LEN=255)   :: FILENAME

    IF(DATA%test==12 .OR. DATA%test==500 .OR. DATA%test==501 &
      & .OR. DATA%test==100 .OR. DATA%test==101 .OR. DATA%test==111 &
      & .OR. DATA%test==120 .OR. DATA%test==121 .OR. DATA%test==122 &
      & .OR. (DATA%test==300 .AND. DATA%perturbation==0) .OR. DATA%test==210 &
      & .OR. DATA%test==150 .OR. (DATA%test==400 .AND. DATA%perturbation==0) ) THEN


       FILENAME = 'EndPointsErrorL1_XXXX.dat'
       WRITE(FILENAME(18:21),FMT="(I4.4)") Mesh%nt   
       OPEN(666,FILE=TRIM(FILENAME))
       WRITE(666,*) Mesh%nt,  error_L1(1:n_vars)
       WRITE(666,*) "number of elements, errors in l1 norm"
       CLOSE(666)

       FILENAME = 'EndPointsErrorL2_XXXX.dat'
       WRITE(FILENAME(18:21),FMT="(I4.4)") Mesh%nt   
       OPEN(667,FILE=TRIM(FILENAME))
       WRITE(667,*) Mesh%nt,  error_L2(1:n_vars)
       WRITE(667,*) "number of elements, errors in l2 norm"
       CLOSE(667)

       FILENAME = 'EndPointsErrorLinf_XXXX.dat'
       WRITE(FILENAME(20:23),FMT="(I4.4)") Mesh%nt   
       OPEN(668,FILE=TRIM(FILENAME))
       WRITE(668,*) Mesh%nt,  error_linf(1:n_vars)
       WRITE(668,*) "number of elements, errors in linf norm"
       CLOSE(668)


    END IF
 END SUBROUTINE printerrorsEndPoints


  !*--------------------------------------------------------------------------------------
  !*MyPrintFinalSolution prints the numerical solution in a .dat file
  !*--------------------------------------------------------------------------------------
  SUBROUTINE MyPrintFinalSolution(Coeff, Mesh, DATA,kt)
    IMPLICIT NONE
    TYPE(PVar), DIMENSION(:), INTENT(in):: Coeff !*Coefficients of the solution
    TYPE(maillage), INTENT(in)  :: Mesh 
    TYPE(donnees),  INTENT(inout)  :: DATA
    INTEGER, INTENT(IN) :: kt

    TYPE(PVar), DIMENSION(:), ALLOCATABLE :: Val !*Values of the solution

    REAL(DP), DIMENSION(:), ALLOCATABLE :: x !*Abscissae

    !*Loops variables used in the integration of the difference 
    INTEGER :: indt=0 !*Index for the loop on the triangles
    INTEGER :: indi=0 !*Index for the first loop on the DoFs i

    !*To deal with the numeration issues
    INTEGER, DIMENSION(Mesh%nDoFs) :: GlobalNumeration !*GobalNumeration(indi) is the global index of the indi DoF with increasing order of abscissa   
    INTEGER :: orderpolynomials 
    INTEGER :: indsupp !*Loop on the DoFs

    
    TYPE (element) :: e !*For the quick reference to the generic Mesh%e(indt) in the loops

    CHARACTER(LEN=255)   :: FILENAME
    
    CHARACTER(LEN=255)   :: DoFnumber

    TYPE(PVar) :: ValICPrint !*Values of the initial condition printed either analytical or numerical
    CHARACTER(LEN=75) :: perturbationprint, ktprint

    IF( (DATA%test==400 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==401 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==410 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==411 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==420 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==421 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==430 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==431 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==440 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3 .OR. DATA%perturbation==4) ) &
        & .OR. (DATA%test==441 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==415 .AND. DATA%perturbation==0 ) &
        & .OR. (DATA%test==416 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==425 .AND. DATA%perturbation==0 ) &
        & .OR. (DATA%test==426 .AND. (DATA%perturbation==0 &
        &.OR. DATA%perturbation==1 .OR. DATA%perturbation==2 .OR. DATA%perturbation==3) ) &
        & .OR. (DATA%test==435 .AND. DATA%perturbation==0 ) &
        & .OR. (DATA%test==120 .AND. DATA%perturbation==0 ) & 
        & .OR. (DATA%test==121 .AND. DATA%perturbation==0 ) &
        & .OR. (DATA%test==122 .AND. DATA%perturbation==0 ) ) THEN

       !*Allocate the vector of the values
       ALLOCATE(Val(Mesh%nDoFs))
       !*Safe initialization
       Val=0._DP
       !*Passing from the coefficients to the values
       Val=global_Control_to_Cons(Coeff,Mesh) 

       !*Allocate the vector of the abscissae
       ALLOCATE(x(Mesh%nDoFs))
       !*Safe initialization
       x=0._DP
       !*Acquisition
       DO indt=1,Mesh%nt
          e=Mesh%e(indt)
          DO indi=1,e%nsommets
             x(e%nu(indi))=e%coor(indi)
          END DO
       END DO

       !*!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*DO indi=1,Mesh%Ns
       !*   PRINT*, indi, x(indi)
       !*END DO
       !*!!!!!!!!!!!!!!!!!!!!!!!!!

       !*We need now to print x and Val but the order is still not good, we would like to order with increasing abscissa

       !*---------------------------------------------------------------------------------
       !*---------------------------------------------------------------------------------        
       !*---------------------------------------------------------------------------------
       !*Reordering issues
       !*1) Reordering structure
       !*GlobalNumeration, structure needed for the reordering
       !*GobalNumeration(indi) is the global index of the indi DoF with increasing order of abscissa    

       GlobalNumeration=0

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
          orderpolynomials=1 
       CASE(2,3,12) !*B2, P2 (eventually in Gauss-Lobatto points)
          orderpolynomials=2
       CASE(4,5,13) !*P3, B3 (P3 eventually in Gauss-Lobatto points)
          orderpolynomials=3
       CASE(6,7,14) !*B4, P4 (PGL4 eventually in Gauss-Lobatto points)
          orderpolynomials=4
       CASE default
          PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
          STOP
       END SELECT

       !*Vertices
       DO indt=1,Mesh%Nt+1
          GlobalNumeration((indt-1)*orderpolynomials+1)=indt
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*SAFETY CHECK
          !*PRINT*, "Vertex", indt 
          !*PRINT*, "Index in the increasing numeration", indt+(indt-1)*orderpolynomials
          !*PRINT*, "Index in the Global numeration", GlobalNumeration((indt-1)*orderpolynomials+1)
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO



       !*Internal nodes
       indsupp=0
       SELECT CASE(DATA%itype)
       CASE(1,11) !*P1 (eventually in Gauss-Lobatto points)
          !*No internal nodes
       CASE(2,3,12) !*B2, P2 (eventually in Gauss-Lobatto points)
          indsupp=0
          DO indt=1,Mesh%Nt
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+1)=Mesh%Nt+1+indsupp
          END DO
       CASE(4,5,13) !*P3, B3 (P3 eventually in Gauss-Lobatto points)
          indsupp=0
          DO indt=1,Mesh%Nt
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+1)=Mesh%Nt+1+indsupp
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+2)=Mesh%Nt+1+indsupp
          END DO
       CASE(6,7,14) !*P3, B3 (P3 eventually in Gauss-Lobatto points)
          indsupp=0
          DO indt=1,Mesh%Nt
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+1)=Mesh%Nt+1+indsupp
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+2)=Mesh%Nt+1+indsupp
             indsupp=indsupp+1
             GlobalNumeration((indt-1)*orderpolynomials+1+3)=Mesh%Nt+1+indsupp
          END DO



       CASE default
          PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
          STOP
       END SELECT

       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*DO indi=1,Mesh%nDoFs
       !*   PRINT*, indi, GlobalNumeration(indi)
       !*END DO
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*---------------------------------------------------------------------------------
       !*---------------------------------------------------------------------------------
       !*---------------------------------------------------------------------------------

       !*Now we can print


       !*0) Open the file
       SELECT CASE(DATA%test)
       CASE(400) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             !*WRITE(ktprint,*) kt
	          !*FILENAME = TRIM(ADJUSTL('LakeAtRestSwashesDiscreteEquilibrium'//TRIM(ADJUSTL(ktprint))//'.dat'))
             FILENAME = TRIM(ADJUSTL('LakeAtRestSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             !*WRITE(FILENAME,"(A5,I1)") "LakeAtRestSwashesPerturbation", DATA%perturbation
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("LakeAtRestSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT
       CASE(401) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
	          FILENAME = 'LakeAtRestSwashesDiscreteEquilibriumBis.dat'
          CASE(1,2,3)
             !*WRITE(FILENAME,"(A5,I1)") "LakeAtRestSwashesPerturbation", DATA%perturbation
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("LakeAtRestSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT
       CASE(120) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SupercriticalSmoothSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT
       CASE(410) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SupercriticalNotSmoothSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("SupercriticalNotSmoothSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT
       CASE(411) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SupercriticalNotSmoothSwashesDiscreteEquilibriumBis.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("SupercriticalNotSmoothSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
		    CASE DEFAULT
				 PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
				 STOP          
		    END SELECT
       CASE(415) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SupercriticalNotSmoothFrictionSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT
       CASE(416) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SupercriticalNotSmoothFrictionSwashesDiscreteEquilibriumBis.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("SupercriticalNotSmoothFrictionSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
		    CASE DEFAULT
				 PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
				 STOP          
		    END SELECT

       CASE(121) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SubcriticalSmoothSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT
       CASE(420) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SubcriticalNotSmoothSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("SubcriticalNotSmoothSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT

       CASE(421) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SubcriticalNotSmoothSwashesDiscreteEquilibriumBis.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("SubcriticalNotSmoothSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
		    CASE DEFAULT
				 PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
				 STOP          
		    END SELECT

       CASE(425) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SubcriticalNotSmoothFrictionSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT

       CASE(426) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SubcriticalNotSmoothFrictionSwashesDiscreteEquilibriumBis.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("SubcriticalNotSmoothFrictionSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
		    CASE DEFAULT
				 PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
				 STOP          
		    END SELECT

       CASE(122) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('TranscriticalSmoothSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT
       CASE(430) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('TranscriticalNotSmoothSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("TranscriticalNotSmoothSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT

       CASE(431) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('TranscriticalNotSmoothSwashesDiscreteEquilibriumBis.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("TranscriticalNotSmoothSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
		    CASE DEFAULT
				 PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
				 STOP          
		    END SELECT

       CASE(435) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('TranscriticalNotSmoothFrictionSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT




       CASE(440) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SupercriticalFrictionShortChannelSwashesDiscreteEquilibrium.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3,4)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("SupercriticalFrictionShortChannelSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
		    CASE DEFAULT
				 PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
				 STOP          
		    END SELECT

       CASE(441) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
             FILENAME = TRIM(ADJUSTL('SupercriticalFrictionShortChannelSwashesDiscreteEquilibriumBis.dat'))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
          CASE(1,2,3)
             WRITE(perturbationprint,*) DATA%perturbation
             WRITE(FILENAME,*)TRIM(ADJUSTL(TRIM("SupercriticalFrictionShortChannelSwashesPerturbation")//TRIM(ADJUSTL(perturbationprint))//".dat"))
             FILENAME=ADJUSTL(FILENAME)
             FILENAME=TRIM(FILENAME)
		    CASE DEFAULT
				 PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
				 STOP          
		    END SELECT

       CASE default
          PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
          STOP
       END SELECT

       OPEN(700,FILE=TRIM(FILENAME))       


       !*Depending on the case we may want to write different information
       SELECT CASE(DATA%test)
       CASE(400) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
				 WRITE(700,*) "indi"," ","x"," ","h"," ","hu", " ", "h_init", " ", "hu_init"
				 DO indi=1,Mesh%Ns
                ValICPrint=0._DP
                ValICPrint=IC(x(GlobalNumeration(indi)),DATA)
				    WRITE(700,*) indi, x(GlobalNumeration(indi)), &
                & Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2), &
                & ValICPrint%u(1), ValICPrint%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !*SAFETY CHECK
                !*PRINT*, indi, x(GlobalNumeration(indi)), Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				 END DO
          CASE(1,2,3)
				 WRITE(700,*) "indi"," ","x"," ","h"," ","hu", " ", "h_init", " ", "hu_init"
				 DO indi=1,Mesh%Ns
                ValICPrint=0._DP
                ValICPrint=IC(x(GlobalNumeration(indi)),DATA)
				    WRITE(700,*) indi, x(GlobalNumeration(indi)), &
                & Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2), &
                & ValICPrint%u(1), ValICPrint%u(2)
				 END DO
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP
          END SELECT
       CASE(401,410,411,415,416,420,421,425,426,430,431,435,120,121,122) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
				 WRITE(700,*) "indi"," ","x"," ","h"," ","hu", " ", "h_init", " ", "hu_init"
				 DO indi=1,Mesh%Ns
                ValICPrint=0._DP
                ValICPrint=ValICUtils(GlobalNumeration(indi))
				    WRITE(700,*) indi, x(GlobalNumeration(indi)), &
                & Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2), &
                & ValICPrint%u(1), ValICPrint%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !*SAFETY CHECK
                !*PRINT*, indi, x(GlobalNumeration(indi)), Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				 END DO
          CASE(1,2,3)
				 WRITE(700,*) "indi"," ","x"," ","h"," ","hu", " ", "h_init", " ", "hu_init"
				 DO indi=1,Mesh%Ns
                ValICPrint=0._DP
                ValICPrint=ValICUtils(GlobalNumeration(indi))
				    WRITE(700,*) indi, x(GlobalNumeration(indi)), &
                & Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2), &
                & ValICPrint%u(1), ValICPrint%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !*SAFETY CHECK
                !*PRINT*, indi, x(GlobalNumeration(indi)), Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				 END DO
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT
       CASE(440,441) 
          SELECT CASE(DATA%perturbation)
          CASE(0)
				 WRITE(700,*) "indi"," ","x"," ","h"," ","hu", " ", "b", " ", "h_init", " ", "hu_init"
				 DO indi=1,Mesh%Ns
                ValICPrint=0._DP
                ValICPrint=ValICUtils(GlobalNumeration(indi))
				    WRITE(700,*) indi, x(GlobalNumeration(indi)), &
                & Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2), &
                & BathymetryAtDoFsUtils(GlobalNumeration(indi)), &
                & ValICPrint%u(1), ValICPrint%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !*SAFETY CHECK
                !*PRINT*, indi, x(GlobalNumeration(indi)), Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				 END DO
          CASE(1,2,3,4)
				 WRITE(700,*) "indi"," ","x"," ","h"," ","hu", " ", "b", " ", "h_init", " ", "hu_init"
				 DO indi=1,Mesh%Ns
                ValICPrint=0._DP
                ValICPrint=ValICUtils(GlobalNumeration(indi))
				    WRITE(700,*) indi, x(GlobalNumeration(indi)), &
                & Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2), &
                & BathymetryAtDoFsUtils(GlobalNumeration(indi)), &
                & ValICPrint%u(1), ValICPrint%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !*SAFETY CHECK
                !*PRINT*, indi, x(GlobalNumeration(indi)), Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2)
                !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				 END DO
          CASE DEFAULT
		       PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
		       STOP          
          END SELECT






       CASE default
          PRINT*, "Wrong test number for MyPrintFinalSolution, test = ", DATA%test
          STOP
       END SELECT




       CLOSE(700)

       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*OPEN(700,FILE=TRIM(FILENAME))       
       !*READ(700,*)
       !*READ(700,*) indi, x(GlobalNumeration(indi)), Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2)
       !*PRINT*, indi, x(GlobalNumeration(indi)), Val(GlobalNumeration(indi))%u(1), Val(GlobalNumeration(indi))%u(2)
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !*Deallocate the vector of the values
       DEALLOCATE(Val,x)
    END IF
  
  END SUBROUTINE MyPrintFinalSolution

END MODULE postprocessing
