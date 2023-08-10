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



!*---------------------------------------------------------------------------------
!*Arete is the class of the boundaries of the elements which are faces in 3d, segments in 2d, points in 1d, 
!*In this case (1d) it is almost useless
!*We will have a lot of useless fields which just make sense for higher dimension and will not be used/initialized here !*THEY WILL BE MARKED BY USELESS AND IT WILL BE SPECIFIED WHETHER THEY CAN BE SAFELY CANCELED OR NOT !*(NB: If they are used I write GOOD)
!*---------------------------------------------------------------------------------
!*SPOILER: THREE FIELDS ONLY ARE IMPORTANT, THE OTHERS CAN BE DELETED
!*I)bord !*GOOD !*It tells whether the arete is on the boundary of the domain or not
!*NB: 
!*Mesh%edge(1)%bord=T
!*Mesh%edge(Mesh%nsegmt)%bord=Mesh%edge(Mesh%nt+1)%bord=T
!*The others are F
!*II) jt1=-1, jt2 =-1 ! les deux elements de par et d'autre de l'arete. !*GOOD
!*The two (in 1d) elements sharing the arete 
!*NB: The aretes at the boundary of the domain will have just one element !*IN THIS CASE IT IS SET jt1=jt2=that element !*Clearly this happens for the first and the last one
!*Mesh%edge(1)%jt1=Mesh%edge(1)%jt2=1
!*Mesh%edge(nsegmt)%jt1=Mesh%edge(nsegmt)%jt2=Mesh%nt
!*For the rest we have
!*Mesh%edge(indi)%jt1=indi-1
!*Mesh%edge(indi)%jt2=indi
!*---------------------------------------------------------------------------------
MODULE arete_class
  USE PRECISION
  IMPLICIT NONE
  TYPE, PUBLIC:: arete
     !*---SUCCEEDED IN DELETING---
     !*INTEGER:: nsommets, itype, nvertex! nbre de dofs, type element: 1-> P1 !*USELESS ---SUCCEEDED IN DELETING---
     ! 2-> B2
     !3->P2
     ! nombre de sommet dans cet element arete
     !*nsommets !*-10, clearly not used, by logic it should be the number of DoFs in the "edge" which would be 1 in the 1d case
     !*itype -1, clearly not used, by logic it should be the itype 
     !*-1, clearly not used, by logic it should be the number of vertices in the "edge" which would be 1 in the 1d case

     LOGICAL:: bord !*GOOD !*It tells whether the arete is on the boundary of the domain or not
     !*NB: 
     !*Mesh%edge(1)%bord=T
     !*Mesh%edge(Mesh%nsegmt)%bord=Mesh%edge(Mesh%nt+1)%bord=T
     !*The others are F
     INTEGER:: jt1=-1, jt2 =-1 ! les deux elements de par et d'autre de l'arete. !*GOOD
     !*The two elements sharing the arete 
     !*NB: The aretes at the boundary of the domain will have just one element !*IN THIS CASE IT IS SET jt1=jt2=that element !*Clearly this happens for the first and the last one
     !*Mesh%edge(1)%jt1=Mesh%edge(1)%jt2=1
     !*Mesh%edge(nsegmt)%jt1=Mesh%edge(nsegmt)%jt2=Mesh%nt
     !*For the rest we have
     !*Mesh%edge(indi)%jt1=indi-1
     !*Mesh%edge(indi)%jt2=indi

     !*---SUCCEEDED IN DELETING---
     !*INTEGER, DIMENSION(:,:), POINTER :: nu =>Null() ! nu( indice des elements, indice des voisins): on cherche a connaitre le numero local des !*USELESS ---SUCCEEDED IN DELETING--- !*0, clearly not used, by logic it should be the vector of the global indices of the DoFs in the "edge" which would be just the global index of the node corresponding to that "edge" in the 1d case
     ! points communs (puisque le maillage est confome)  aux deux elements sur cette face dans chacun des elements jt1 et jt2
     
     !*---SUCCEEDED IN DELETING---
     !*REAL(dp), DIMENSION(:,:), POINTER :: coor=>Null() ! il s'agit des coordonnees physique des dofs (communs) sur la face (commune) !*USELESS ---SUCCEEDED IN DELETING--- !*NOT EVEN INITIALIZED, when you try to print it you may get a segmentation fault
     !*By logic it should be the coordinates of the DoFs contained by the edge (coordinates of a single point in 1d)


     !*volume and jump_flag USELESS ---SUCCEEDED IN DELETING---
     !*NB: volume and jump_flag were declared without ._DP 
     !*I guess they are useless but in any case...
     !*Anyway
     !*volume and jump_flag are never changed, they stay 0 and 1
     !*I guess they are not used
     !*volume should be the measure of the arete (one point in 1d, so trivially 0)
     !*jump_flag??????????
     !*---SUCCEEDED---!*REAL(dp) :: volume=0._DP  !
     !*---SUCCEEDED---!*REAL(dp) :: jump_flag=1.0_DP !(between 0 and 1 which gives the weight of the edge term
     ! for this one check if there is a discontinuity around and take the maximum value)

     !*USELESS ---SUCCEEDED IN DELETING---
     !*REAL(dp),  DIMENSION(:), POINTER :: n   =>Null()   ! normales exterieures !*USELESS
     !*Actually in 2d it makes sense to define the normal to an arete because it is a segment, in 1d we have a point so it makes no sense !*SEGMENTATION FAULT, NOT EVEN DEFINED

     !*USELESS ---SUCCEEDED IN DELETING---
     !*INTEGER  :: log    ! logique !*Always 0


     !*USELESS ---SUCCEEDED IN DELETING---
     !*The quadrature doesn't make sense for an arete in 1d (point)
!!!!   quadrature de surface !
     !*REAL(dp),   DIMENSION(:,:),POINTER :: quad =>Null()  ! point de quadrature 
     !*REAL(dp),     DIMENSION(:),POINTER :: weight=>Null() ! poids 
     !*INTEGER                        :: nquad  ! nbre de points de quadrature 
!!! quadrature bord (dimension -1)
     !*REAL(dp),   DIMENSION(:,:),POINTER :: quad_1 =>Null()  ! point de quadrature 
     !*REAL(dp),     DIMENSION(:),POINTER :: weight_1=>Null() ! poids 
     !*INTEGER                        :: nquad_1  ! nbre de points de quadrature 
!!!
   !*CONTAINS
     !*PROCEDURE, PUBLIC:: aire=>aire_arete !*USELESS !*---SUCCEEDED IN DELETING---
     !*PROCEDURE, PUBLIC:: quadrature=>quadrature_arete
     !*PROCEDURE, PUBLIC:: normale=>normale_arete !*USELESS !*---SUCCEEDED IN DELETING---
     !FINAL:: clean_arete

  END TYPE arete


CONTAINS

  !*aire_arete and normale_arete are clearly bidimensional oriented
  !*They can be removed. Moreover they involve coor which is not even initialized (if you try to print it you get a segmentation fault)
  
  !*!*-------------------
  !*!*In 2d it would be the length of the boundary element !*USELESS ---SUCCEEDED IN DELETING---
  !*!*------------------- 
  !*REAL(dp) FUNCTION aire_arete(e)
  !*  CLASS(arete), INTENT(in):: e
  !*  REAL(dp), DIMENSION(2):: a
  !*  a= e%coor(:,2)-e%coor(:,1)
  !*
  !*  aire_arete=SQRT(a(1)**2+ a(2)**2 )
  !*END FUNCTION aire_arete


  !*!*-------------------
  !*!*In 2d it would be the normal to the boundary element ---SUCCEEDED IN DELETING---
  !*!*------------------- 
  !*FUNCTION normale_arete(e) RESULT(n)
  !*  CLASS(arete), INTENT(in)::e
  !*  REAL(dp), DIMENSION(2):: n
  !*  INTEGER:: l, k1, k2
  !*  INTEGER, DIMENSION(3), PARAMETER:: ip1=(/2,3,1/)
  !*
  !*  n(1)=e%coor(2,2)-e%coor(2,1)
  !*  n(2)=e%coor(1,1)-e%coor(1,2)
  !*
  !*END FUNCTION normale_arete



  !*SUBROUTINE quadrature_arete(e)
  !*  CLASS(arete), INTENT(inout):: e
  !*  REAL(dp):: w,zo,xo,s
  !*  INTEGER:: nquad
  !*
  !*  !Gaussia formula, exact for polynomials of degree 5
  !*  PRINT*, "quadrature_arete"
  !*  STOP
  !*
  !*  e%nquad=3
  !*  nquad=e%nquad
  !*  ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad))
  !*
  !*  s=SQRT(0.6d0)
  !*  e%quad(1,1)=0.5d0*(1.0d0 - s)
  !*  e%quad(2,1)=0.5d0*(1.0d0 + s)
  !*  e%weight(1) = 5.0d0/18.0d0
  !*  e%quad(1,2)=0.5d0*(1.0d0 + s)
  !*  e%quad(2,2)=0.5d0*(1.0d0 - s)
  !*  e%weight(2) = 5.0d0/18.0d0
  !*  e%quad(1,3)=0.5d0
  !*  e%quad(2,3)=0.5d0
  !*  e%weight(3)=8.0d0/18.0d0
  !*
  !*
  !*
  !*END SUBROUTINE quadrature_arete

  !*SUBROUTINE clean_arete(e)
    !*TYPE(arete):: e
    !*IF (ASSOCIATED(e%nu)) NULLIFY(e%nu) ---SUCCEEDED IN DELETING---
    !*IF (ASSOCIATED(e%coor)) NULLIFY(e%coor) ---SUCCEEDED IN DELETING---
    !*IF (ASSOCIATED(e%quad)) NULLIFY(e%quad)
    !*IF (ASSOCIATED(e%weight)) NULLIFY(e%weight)
  !*END SUBROUTINE clean_arete
END MODULE arete_class
