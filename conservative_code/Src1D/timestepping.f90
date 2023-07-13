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

!*----------------------------------------------------------------------------------------------
!*This module is used to compute the updating quantity that one needs to add to the solution in a single subtimestep (in all the nodes clearly) to perform the iteration of the DeC (RMK: in the single subtimestep of the single iteration)
!*The subroutines of this module are recalled in fin() in the main.
!*fin() performs the updating of the single subtimestep inside the loop on the subtimesteps (inside the iteration loop of the DeC inside the timestepping loop) 
!*----------------------------------------------------------------------------------------------
!*Denoting by m the index of the subtimestep (and by p the index of the iteration), the updating formula performed in fin() reads
!*u_i^m,(p)=u_i^m,(p-1)-1/|C_i|*{\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))+\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}]} (#)
!*----------------------------------------------------------------------------------------------
!*Actually the real updating is performed in fin()
!*We have a loop over the DoFs where we set
!*var%ua(is,k-1)=Var%up(is,k-1)-Var%un(is)*Mesh%aires(is) !*Updating
!*Here we just compute
!*{\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))+\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}]}
!*NB: In reality since we are in 1d we do not have the boundary residual
!*\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}
!*They are not needed here, there are escamotages to avoid them. So we just compute
!*{\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))]}
!*----------------------------------------------------------------------------------------------
MODULE timestepping

  ! In this module are collected the main_update, corresponding to the update of the Flux of the scheme
  ! and the edge_main_update which accounts for the stabilization of the scheme 


  USE overloading
  USE element_class
  USE variable_def
  USE arete_class
  USE param2d
  USE scheme
  USE Model
  USE PRECISION
  USE utils
  IMPLICIT NONE
CONTAINS

  !*-----------------------------------------------------------------------------------
  !*Here we add to Var%un the updating quantity
  !*\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + 
  !*+ dt*\sum_{l=0}^M theta_l^m*\Sum_{K \in K_i} RES_i^K(u^l,(p-1))
  !*to perform the update of the DeC
  !*
  !*So we put in Var%un just: the mass matrix part & the node-element (not boundary) residuals up to the edge stabilization (which would be in theory part of RES_i^K(u^l,(p-1)) )
  !*
  !*NB: This subroutine is called inside fin() i.e. inside the loop on the sutimesteps inside the loop on the iterations of the DeC
  !*NB: IT DEALS WITH A SINGLE SUBTIMESTEP AND IT IS CALLED INSIDE A LOOP ON THE TRIANGLES
  !*SO WE WILL COMPUTE HERE THE CONTRIBUTION OF THE NODES OF THE SINGLE TRIANGLE NAMELY
  !*\int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + 
  !*+ dt*\sum_{l=0}^M theta_l^m*\Sum_{K \in K_i} RES_i^K(u^l,(p-1))
  !*NB: We are not forgetting the sums. Everything we compute is added to its DoF thanks to
  !*Var%un(e%nu(l))=Var%un(e%nu(l))+res(l)
  !*The fact that we add to sth existing performs the sum
  !*Also we do it for every triangle so we have
  !*\Sum_{x_j \in K} here inside
  !*\Sum_{K \in K_i} in the outside loop
  !*
  !*In other words: As specified we deal singularly with each subtimestep so we will get something of dimension Nvars*NDoFs (no subtimesteps)
  !*
  !*NB: Outside we have a loop over the elements K and this subroutine is recalled in this loop so we will have a loop on the DoFs i \in K here to calculate the component^K_i 
  !*We do not initializate Var%un every time so we have component_i in the end
  !*It is recalled on a loop on triangles so we have just to make a loop over the DoFs
  !*-----------------------------------------------------------------------------------
  SUBROUTINE main_update(k, dt, e,Var, DATA,alpha,beta,gamma,n_theta,theta,jt,mesh)
    !*------------------------------------------------------------------------------------
    !*SUPER IMPORTANT: The spatial residuals must be combined in time though the theta coefficients
    !*The strategy used here is the following:
    !*The spatial residuals are sth like
    !*Res^K_i=\int_K d_x f phi_i + some stabilization depending on the scheme depending linearly on u
    !*Note that Res^K_i are linear with respect to the fluc 
    !*We must compute \sum_{l=0}^M theta_l^m* RES_i^K(u^l,(p-1))
    !*SO WE MUST COMPUTE
    !*1)\int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) mass matrix part
    !*2)RES_i^K(u^l,(p-1))
    !*3)dt*\sum_{l=0}^M theta_l^m*RES_i^K(u^l,(p-1))
    !*(As already said the sum \Sum_{K \in K_i} is embedded in the fact that we are calling this inside a loop on the elments)
    !*The crucial point is that since RES_i^K(u^l,(p-1)) is linear with respect to the flux
    !*we can combine the flux through the theta coefficients and then plug the combination into RES_i^K(u^l,(p-1))
    !*In synthesis: The flux is likely to be nonlinear but since the operators in which it is involved in this updating part (RES) are linear we can compute it in all the subtimesteps, combine these evaluations through the thetas corresponding to the processed subtimestep and apply RES only once to this combination
    !*The same for the schemes which involve also u (which will be combined in time though the proper theta and then plugged into RES)  
    !*------------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(dp), INTENT(in):: dt
    REAL(dp),DIMENSION(4,2:5):: alpha
    REAL(dp),DIMENSION(4,4,3:5):: beta, gamma
    INTEGER, INTENT(IN):: jt
    INTEGER, DIMENSION(5):: n_theta
    REAL(dp),DIMENSION(0:4,1:4,1:5):: theta
    TYPE(maillage),INTENT(IN):: mesh
    INTEGER, INTENT(in):: k
    TYPE(element), INTENT(in):: e
    TYPE(variables), INTENT(inout):: Var
    TYPE(donnees), INTENT(in):: DATA
    TYPE(Pvar),DIMENSION(n_dim,e%nsommets):: flux, flux_c
    TYPE(Pvar),DIMENSION(e%nsommets,0:DATA%iordret-1):: u, up, u_p, up_p !*They are PVar in all the nodes (first index), in all the subtimesteps
    TYPE(Pvar),DIMENSION(e%nsommets):: res, difference, uu,uu_c, source, source_c
    INTEGER:: l, lp
    !*NB: We are interested in the updating quantity for the subtimestep k because we are looping
    !*from k=1 to DATA%iordret-1

    DO l=1,e%nsommets
       up(l,:) =Var%up(e%nu(l),:) !*He's copying the PVar of all the nodes (first index) in all the subtimesteps (second index 0,1,...,M=order-1) in up
    ENDDO
    u=up !*Copying the coefficients up in u 
    !*NB: up are coefficients


    !*Here we are computing the coefficients for the mass matrix part
    !*\int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) 
    !*l->j
    !*m->k-1
    DO l=1,e%nsommets
       difference(l)=up(l,k)-up(l,0) 
    ENDDO

    !*Here we are passing from the coefficients to the values of the solution
    DO l=0,DATA%iordret-1
       u_P(:,l) =Control_to_Cons(u(:,l),e)
    ENDDO
    up_P=u_P !*Copying the values u_P in up_P
    
    !*RMK:
    !*-> up are coefficients
    !*-> up_P=u_P are values

    !*Initialization of vectors of PVar
    !*The next variables are referred to the single subtimestep
    !*flux_c(e%nsommets,n_dim) combination thorugh the thetas of the point evaluations of the VALUES of the flux in the nodes of the element coming from all the subtimesteps
    !*uu(e%nsommets) combination thorugh the thetas of the VALUES of the solution in the nodes of the element coming from all the subtimesteps
    !*source_c(e%nsommets) combination thorugh the thetas of the point evaluations of the source in the nodes of the element coming from all the subtimesteps
    flux_c=0._dp; uu=0._dp; source_c=0._dp
    !*BASICALLY DUE TO THE LINEARITY OF RES W.R.T. THE TRIPLE (flux,u,source) we first combine in time through the thetas and then apply res only once to the combination rather then applying it to every subtimestep and combining later

    DO l=1, e%nsommets !*Loop on the DoFs
       DO lp=0, n_theta(DATA%iordret)-1 !*Loop on the subtimesteps that must be combined over the thetas

          !*RMK:
          !*theta coefficients for L^2 in the DeC
          !*theta(p,l,k)=int_0^alpha(l,k) phi_p(s)ds
          !*p->0,1,...,k-1 basis functions
          !*l->1,2,...,k-1 subtimesteps
          !*k order

          !*The third index is referred to the order -> DATA%iordret, fixed
          !*The second index is the one referred to the subtimestep so we must set it to be k-1 (in this case), fixed 
          !*The first index is referred to the basis function in time so it must go with the second index of the coefficients (the first index of the coefficients is referred to the DoFs), IT IS THE CRUCIAL INDEX IN THE COMBINATION

          flux_c(:,l)= flux_c(:,l)+ theta(lp,k,DATA%iordret)* up_P(l,lp)%flux((/e%coor(l)/)) 
          !*Second index of flux referred to the dimensions (1 in this case)
          !*First index of up_P referred to the DoFs, second index to the subtimesteps
          !*NB: flux computed from the values and not from the coefficients

          uu(l)      = uu(l)      + theta(lp,k,DATA%iordret)* up  (l,lp) !*SAME FOR uu !*up are coefficients
          source_c(l)= source_c(l)+ &
          & theta(lp,k,DATA%iordret)*sourcetermfunction(e%coor(l),up_P(l,lp),e%nu(l),DATA%test)
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !*PREVIOUS VERSION WITH REMI'S SOURCE
          !*source_c(l)= source_c(l)+ theta(lp,k,DATA%iordret)* up_P(l,lp)%source((/e%coor(l)/),DATA%test) !*SAME FOR source !*up_P are values
          !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
       ENDDO

       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK 
       !*to make sure that I'm giving as input the global DoF and its coordinate to source
       !*PRINT*, e%nu(l), e%coor(l)
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !*NOW FINALLY WE HAVE THE COMBINATION IN TIME OF f,u and s AND WE CAN APPLY RES ONLY ONCE


#if (1==0)
       !*This is not needed anymore, I guess it belongs to the first version of the operator L^1
       !*It is not even compiled because of #if (1==0)


       ! the pure Jacobi-like does not need this, only the 
       ! the Gauss-Seidel like
       ! version 1
       DO lp=1,k
          flux_c(:,l)= flux_c(:,l)-&
               &GAMMA(lp, DATA%iordret-1,DATA%iordret)*( up_P(l,lp)%flux( (/e%coor(l)/) )&
               & -u_p(l,lp)%flux( (/e%coor(l)/) )   )
          source_c(l)= source_c(l)-&
               &GAMMA(lp, DATA%iordret-1,DATA%iordret)*(sourcetermfunction(e%coor(l),up_P(l,lp),e%nu(l),DATA%test)&
               & -sourcetermfunction(e%coor(l),u_p(l,lp),e%nu(l),DATA%test)   )
               !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !*PREVIOUS VERSION WITH REMI'S SOURCE
               !*source_c(l)= source_c(l)-&
               !*&GAMMA(lp, DATA%iordret-1,DATA%iordret)*( up_P(l,lp)%source( (/e%coor(l)/) ,DATA%test)&
               !*& -u_p(l,lp)%source( (/e%coor(l)/),DATA%test )   )
               !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          uu(l)=uu(l) -GAMMA(lp,DATA%iordret-1,DATA%iordret)*( up(l,lp)-u(l,lp))
       ENDDO
#endif 
    ENDDO

    !*We pass now from the values to the coefficients of the "COMBINED" flux and source
    DO l=1, n_dim
       flux(l,:)=Cons_to_control(flux_c(l,:),e)
    ENDDO
    source = Cons_to_control(source_c,e)

    !*Finally we apply schema to compute the contribution of the element to the DoFs (outside there is a loop on the elements and we are now focusing on the single element) 
    !*NB we can compute RES only once because we already combined the argument through the thetas and it is linear with respect to the argument 
    !*NB: We are entering in this function
    !*flux COEFFICIENTS combined through the thetas
    !*source COEFFICIENTS combined through the thetas
    !*uu COEFFICIENTS, combined through the thetas (because computed from up)
    !*difference variation of the COEFFICIENTS in time
    res=schema( DATA%ischema, e, uu, difference, flux, source, dt, jt,mesh)

    !*As anticipated everything is put into Var%un
    DO l=1, e%nsommets
       Var%un(e%nu(l))=Var%un(e%nu(l))+res(l)
    ENDDO


  END SUBROUTINE main_update


  !*--------------------------------------------------------------------------------------
  !*edge_main_update computes the contribute to the updating of the jump stabilization (CIP/Burman/gradient across the edges)
  !*NB: DIFFERENTLY FROM main_update IT IS NOT CALLED INSIDE A LOOP ON THE ELEMENTS/EDGES SO A GLOBAL SPATIAL LOOP IS NECESSARY HERE INSIDE
  !*Anyway like main_update it is recalled inside a the loop over the subtimesteps in the loop over the iterations in the timestepping loop
  !*--------------------------------------------------------------------------------------
  !*In practice we should make a loop over the DoFs i, look ONCE on each edge f of \tilde{f}_i which is the collection of the edges of the elements of K_i (elements containing i) and compute
  !*\int_f [grad u][grad phi_i]
  !*Add it to Var%un so to add in the end
  !*\sum_{f\in \tilde{f}_i} \int_f [grad u][grad phi_i]
  !*INSTEAD
  !*In practice we make a loop over the edges f (points here), we consider the DoFs "involved" in the edge, the DoFs whose gradient has a jump over f which are basically the DoFs contained in the two elements sharing the edge and we compute
  !*\int_f [grad u][grad phi_i]
  !*Add it to Var%un so to have in the end
  !*\sum_{f\in \tilde{f}_i} \int_f [grad u][grad phi_i]
  !*BASICALLY WE SWITCH THE LOOP ON THE DoFs AND THE LOOP ON THE EDGES
  !*NB: Obviously this is done just for internal edges
  !*--------------------------------------------------------------------------------------
  SUBROUTINE edge_main_update(k,DATA,Mesh,Var,dt,alpha,beta,gamma,n_theta,theta)
    IMPLICIT NONE
    REAL(dp),DIMENSION(4,2:5):: alpha
    REAL(dp),DIMENSION(4,4,3:5):: beta, gamma

    INTEGER, DIMENSION(5):: n_theta
    REAL(dp),DIMENSION(0:4,1:4,1:5):: theta
    INTEGER, INTENT(IN):: k

    TYPE(Maillage), INTENT(in):: Mesh
    TYPE(donnees), INTENT(in):: DATA
    TYPE(variables), INTENT(inout):: Var
    REAL(dp), INTENT(in):: dt
    TYPE(element):: e, e1, e2
    TYPE(arete):: ed
    TYPE(Pvar), DIMENSION(:), ALLOCATABLE :: u1, u2
    TYPE(Pvar),DIMENSION(:,:), ALLOCATABLE:: u, up,ua1,ua2,up1,up2, up_p, u_p
    TYPE(Pvar), DIMENSION(:,:), ALLOCATABLE:: resJ
    INTEGER:: jt1, jt2, iseg, l, lp
    ! question : nothing on Jump for Version 1 ?


    DO iseg=1, Mesh%nsegmt !*Loop on the boundaried between the elements (points here)

       ed=Mesh%edge(iseg) !*Fast reference to the processed edge

       IF (ed%bord)  CYCLE ! we are on the boundary !*The interior penalty can be just made on internal edges, when we are on the bord skip the iteration of the loop
       !*bord is a logical field of arete, it tells whether the arete is on the boundary of the domain or not
       !*NB: 
       !*Mesh%edge(1)%bord=T
       !*Mesh%edge(Mesh%nsegmt)%bord=Mesh%edge(Mesh%nt+1)%bord=T
       !*The others are F
 
       !*jt1 and jt2 are the two elements sharing the arete 
       !*NB: The aretes at the boundary of the domain will have just one element !*IN THIS CASE IT IS SET jt1=jt2=that element !*Clearly this happens for the first and the last one
       !*Mesh%edge(1)%jt1=Mesh%edge(1)%jt2=1
       !*Mesh%edge(nsegmt)%jt1=Mesh%edge(nsegmt)%jt2=Mesh%nt
       !*For the rest we have
       !*Mesh%edge(indi)%jt1=indi-1
       !*Mesh%edge(indi)%jt2=indi
       !*SINCE WE ARE SKIPPING THE BOUNDARY EDGES WE ARE IN THIS LAST CASE
       !*Mesh%edge(indi)%jt1=indi-1
       !*Mesh%edge(indi)%jt2=indi

       jt1 = ed%jt1 !*Quick reference to the index of the first element containing the processed edge
       e1=Mesh%e(jt1) !*Quick reference to the first element containing the processed edge

       jt2 = ed%jt2 !*Quick reference to the index of the second element containing the processed edge
       e2=Mesh%e(jt2) !*Quick reference to the second element containing the processed edge

      !*type_flux is an integer field of the type element which is set to be equal to the integer of the scheme chosen

      !*SUPER IMPORTANT: even if in fin() edge_main_update is called only for DATA%ischema==5 .OR. DATA%ischema==4 it is called also in mood so the following IFs make sense
      !*Anyway note how in the "standard" case where we do not have mood we enter in these IFs

      IF (e1%type_flux==2&!4&
         & .OR.e1%type_flux==1&!5 &
         & .OR.e1%type_flux==4&!5 & !*<--------------IT ENTERS
         & .OR.e1%type_flux==5&!5 & !*<--------------IT ENTERS
         & .OR.e1%type_flux==6 &
         & .OR.e1%type_flux==7 &
         & .OR.e1%type_flux==9) THEN !(with jumps--> Burman'stuff on cell e1 )
!!$     

         IF  (e2%type_flux==2&!4&
            & .OR. e2%type_flux==1&!5 &
            & .OR. e2%type_flux==4&!5 & !*<--------------IT ENTERS
            & .OR. e2%type_flux==5&!5 & !*<--------------IT ENTERS
            & .OR.e2%type_flux==6 &
            & .OR.e2%type_flux==7 &
            & .OR.e2%type_flux==9) THEN !(with jumps--> Burman'stuff on cell e2 )

               !*He performs the same trick as in main_update.
               !*We should compute the contributions of the edge stabilization to the residuals in all the subtimesteps and then combine them through the theta coefficients corresponding to the treated subtimestep 
               !*Since the edge stabilization is linear in u he FIRST combines the u of the different subtimesteps in time through the theta coefficients and THEN he computes the edge stabilization residual just once on the combination
               ALLOCATE(u1(e1%nsommets), u2(e2%nsommets), resJ(e1%nsommets,2) ) 
               ALLOCATE(ua1(e1%nsommets,0:DATA%iordret-1),&
               & ua2(e2%nsommets,0:DATA%iordret-1),&
               & up1(e1%nsommets,0:DATA%iordret-1),&
               & up2(e2%nsommets,0:DATA%iordret-1))
               lp=n_theta(DATA%iordret)-1

               !*Copying of the solution in the two elements sharing the processed edges in all the subtimesteps for fast referencing
               !*NB: ua will not be used, we work with Var%up because this is the phylosophy of the code and then we will update Var%ua (outside this sub obviously) so also here inside ua will not be used 
               !*In the ua and up variables we have two indices
               !*First index -> DoFs of the element
               !*Second index -> subtimesteps 0,1,...,M=order-1
               DO l=0,DATA%iordret-1
                  !*FIRST ELEMENT
                  ua1(:,l)=Var%ua(e1%nu(:),l)
                  up1(:,l)=Var%up(e1%nu(:),l)
                  !*SECOND ELEMENT
                  ua2(:,l)=Var%ua(e2%nu(:),l)
                  up2(:,l)=Var%up(e2%nu(:),l)
               ENDDO

               !*Initialization of the variable u1 which will contain the combination in time through the thetas of the coefficients in the DoFs of the first element sharing the processed edge
               u1=0._dp

               DO l=1,e1%nsommets !*LOOP ON THE DoFs
                  !*TWO ALTERNATIVES, ACTUALLY JUST THE FIRST ONE IS CORRECT, THE OTHER ONE IS IN PRINCIPLE WRONG
#if (1==1)
                  !*1) COMBINATION through the theta coefficients of the coefficients of the DoF l associated to the subtimesteps
                  !interpolate in time
                  DO lp=0, n_theta(DATA%iordret)-1 !*Loop over the subtimesteps !*n_theta(k)=k is the number of basis functions (and so subtimesteps) for each order, basically the order because the basis functions are M+1=accuracy=k so 0,1,...,M
                     u1(l)=u1(l)+theta(lp,k,DATA%iordret)* up1(l,lp) !*Combination through the thetas
                     !*theta has three indices
                     !*theta(p,l,k)=int_0^alpha(l,k) phi_p(s)ds
                     !*p->0,1,...,k-1 basis functions
                     !*l->1,2,...,k-1 subtimesteps
                     !*k order

                     !*The third index is fixed to the order of the scheme DATA%iordret
                     !*The second one is referred to the specific subtimestep we are dealing with i.e. k-1
                     !*The first index moves across the basis functions 
                  ENDDO
#else
                  !*2) NO COMBINATION, for the edge stabilization we simply take, for the subtimestep we are dealing with, the coefficients associated to that sutimestep
                  ! no interpolation in time
                  lp=k !*Subtimestep
                  u1(l)= up1(l,lp) 
#endif
               ENDDO

!*\int_f [grad u][grad phi_i]
  !*Add it to Var%un so to add in the end
  !*\
               !*SAME STORY FOR THE SECOND ELEMENT
               u2=0._dp
               DO l=1,e2%nsommets !*Loop over the DoFs
#if (1==1)
                  !interpolate in time
                  DO lp=0, n_theta(DATA%iordret)-1 !*Loop over the subtimesteps
                     u2(l)=u2(l)+theta(lp,k,DATA%iordret)* up2(l,lp) !*Combination
                  ENDDO
#else
                  ! no interpolation in time
                  lp=k
                  u2(l)= up2(l,lp)
#endif   
               ENDDO



!*ALERT
!*THIS IS NOT EVEN COMPILED, I GUESS IT BELONGS TO THE PREVIOUS VERSION OF THE DeC, WHERE THE OPERATOR L^1 WAS MORE COMPLICATED.
!*ANYWAY IT IS NOT USED
!*I SMELL THAT IT COULD ALSO BE WRONG -> Before the updating ua and up are equal so 
!*up1(l,lp)-ua1(l,lp)
!*I DO NOT DELETE 
!*DAVIDE TOLD ME THAT IT IS IMPORTANT FOR SO I LEAVE IT (reference where it is used: R.Abgrall and D.Torlo - Some preliminary results on a high order asymptotic preserving explicit kinetic scheme)
!*ANYWAY TO BE USED YOU MUST DO STH ON ua OR up OTHERWISE THEY ARE THE SAME
#if (1==0)
              !Gauss-Seidel like
              DO l=1,e1%nsommets
                 DO lp=1,k
                    u1(l)=u1(l)-GAMMA(lp,DATA%iordret-1,DATA%iordret)* (up1(l,lp)-ua1(l,lp))
                 ENDDO
              ENDDO

              DO l=1,e2%nsommets
                 DO lp=1,k
                    u2(l)=u2(l)-GAMMA(lp,DATA%iordret-1,DATA%iordret)* (up2(l,lp)-ua2(l,lp))
                 ENDDO
              ENDDO
#endif


              !*Now that we have in u1 and u2 the combinations through the thetas of the different subtimesteps we can work on them and get the contributes of the edge stabilization to the node residuals through jump which is subroutine

              ! jump of grad u
              CALL jump( DATA%ijump, ed, e1, e2, u1, u2, resJ, DATA%alpha_jump,DATA%alpha_jump2)

              !*Jump takes in input 
              !*-ed, the processed edge (a point in 1d)
              !*-e1 and e2, the two elements sharing that edge 
              !*-u1 and u2, the combinations of the coefficients in the two elements already combined through the theta coefficients
              !*-DATA%alpha_jump and DATA%alpha_jump2, the parameters of the jump stabilization
              !*-resJ which is actually the output, the contribution to the nodal residual of the edge stabilization
             
              !*We add the edge contribution to Var%un of the nodes after multiplying it by dt

              !*For the first element
              DO l=1, e1%nsommets
                 Var%un(e1%nu(l))=Var%un(e1%nu(l))+resJ(l,1) * dt
              ENDDO

              !*For the second element
              DO l=1, e2%nsommets
                 Var%un(e2%nu(l))=Var%un(e2%nu(l))+resJ(l,2) * dt
              ENDDO
 
              !*Clean
              DEALLOCATE(resJ,u1,u2,ua1,ua2,up1,up2)
!!$
           END IF
        END IF
 ENDDO
END SUBROUTINE edge_main_update


!*---------------------------------------------------------------------------
!*test is part of the stabilization technique called "mood"
!*I'm not debugging/checking it but I just report the following ALERT
!*->It's basically all in SP, the real are declared as REAL and not REAL(DP)
!*---------------------------------------------------------------------------
SUBROUTINE test(k_iter,Debug,Var, mesh, DATA, flux_mood) !*ALERT IT WAS MOSTLY IN SP, I'M CHANGING IT INTO DP
  IMPLICIT NONE
  ! tunable parameters:
  REAL(DP), PARAMETER:: cour_max=10000._DP !*Added DP
  REAL(DP) :: eps1, eps2 !*Added DP
  REAL(DP), PARAMETER:: coeff=1._DP !0. !*Added DP
  ! variables to be checked
  INTEGER, PARAMETER:: n_list=1 !3
  INTEGER, DIMENSION(n_list), PARAMETER:: list=(/1/)!,2,3/)
  !
  ! emergency schemes
  INTEGER, PARAMETER:: theflux=1,theflux2=0
  !
  !
  INTEGER, INTENT(in):: k_iter
  TYPE(maillage), INTENT(inout):: mesh
  TYPE(variables), INTENT(in)::  debug, Var
  TYPE(donnees), INTENT(in):: DATA
  INTEGER, INTENT(IN):: flux_mood

  TYPE(arete):: ed
  TYPE(pvar), DIMENSION(:),ALLOCATABLE:: deb
  !REAL(8),DIMENSION(:),ALLOCATABLE:: coefficient!with Remi's version of u2
  TYPE(element):: e, eL, eR
  TYPE(pvar),DIMENSION(:),ALLOCATABLE:: coefficient !!with our version of u2
  TYPE(pvar), DIMENSION(:),ALLOCATABLE:: coefficient1, coefficient2

  TYPE(pvar), DIMENSION(:),ALLOCATABLE:: v_min, v_max,u2_min,u2_max, u2_minloc, u2_maxloc
  REAL(DP),  DIMENSION(:),ALLOCATABLE:: w_min, w_max,vd_min, vd_max !*Added DP
  TYPE(Pvar),DIMENSION(:), ALLOCATABLE::phys, phys_d,phys_control, phys_d_control
  TYPE(Pvar),DIMENSION(:),ALLOCATABLE::temp
  REAL(DP),DIMENSION(:),ALLOCATABLE:: temp_min, temp_max !*Added DP
  INTEGER:: jt, k,  l, lk,diag, is, jseg, lkk,r
  INTEGER:: jt1, jt2
  REAL(DP) :: val_min,val_max,eps,xg,yg,val_min_n, val_max_n !*Added DP
  REAL(DP) :: eps_i, absK !*Added DP
  REAL(DP) :: ratio, alphaR, alphaL, smooth_sensor !*Added DP
  TYPE(pvar),DIMENSION(:),ALLOCATABLE::u1_der, u1_derR, u1_derL
  REAL(DP),DIMENSION(:),ALLOCATABLE::u1_derL_mean, u1_der_mean,u1_derR_mean,coefficient_mean !*Added DP
  REAL(DP),DIMENSION(:),ALLOCATABLE:: vL_min,vL_max,vR_min,vR_max,vL,vR !*Added DP
  INTEGER,PARAMETER :: plateau = 1
  INTEGER,PARAMETER :: NAN_criteria = 2
  INTEGER,PARAMETER :: PAD_criteria = 3
  INTEGER,PARAMETER ::  DMP_B1=4
  INTEGER,PARAMETER :: DMP_nou2=5
  INTEGER,PARAMETER ::DMP_u2_2=6
  INTEGER,PARAMETER ::DMP_u2_1=7

  ! integer, dimension(2,e%nsommets-1):: nuloc
  Mesh%e(:)%diag2= Mesh%e(:)%diag

  Mesh%e(:)%diag=0

  ALLOCATE(v_min(Mesh%nt), v_max(Mesh%nt),w_min(Mesh%nt), w_max(Mesh%nt))
  ALLOCATE(vd_min(Mesh%nt), vd_max(Mesh%nt))
  ALLOCATE(phys(Mesh%ndofs),phys_d(Mesh%ndofs))
  ALLOCATE(phys_control(Mesh%ndofs),phys_d_control(Mesh%ndofs),u2_max(Mesh%nt),u2_min(Mesh%nt))

  ALLOCATE(vL(Mesh%nt),vR(Mesh%nt),vL_min(Mesh%nt),vL_max(Mesh%nt),vR_min(Mesh%nt),vR_max(Mesh%nt))
  ALLOCate(u1_der_mean(Mesh%nt),u1_derL_mean(Mesh%nt),u1_derR_mean(Mesh%nt),coefficient_mean(Mesh%nt))
  !compute physical quantities from interpolated quantities.
  ! is this really usefull for Bezier thanks to the TV property?
  ! If not done : more strict condition since Bezier>0==> positive function
  ! converse not true
  DO jt=1, Mesh%nt
     e=Mesh%e(jt)
     eps1=0._DP!10.**(-3) !*Added DP
     eps2=0._DP!10.**(-4)!%volume !*Added DP
     ALLOCATE(temp(e%nsommets))

     !---------------------------
     !---------------------------
     ! Define the vector of primitives at the initial subtimestep n,0
     temp=Control_to_cons(Var%ua(e%nu,0),e)! before ua. up: we compare with the solution before the update
     DO l=1,e%nsommets
        phys(e%nu(l))=convert_cons2prim(temp(l)) !vector of primitive variables in physical values for n,0
     ENDDO

     !...................................
     ! warning: for convergence study, better to work with phys variables'
     ! because of u2 criterion-> use the right way to compute second derivatives
     ! comment here
     temp=Cons_to_Control(phys(e%nu),e)
     phys_control(e%nu)=Var%ua(e%nu,0)!temp !vector of primitive variables in control values for n,0
     !end comment
     !...................................

     !---------------------------
     !---------------------------
     ! Define the vector of primitives at the current subtimestep n,m
     temp=Control_to_cons(debug%ua(e%nu,k_iter),e)
     DO l=1,e%nsommets
        phys_d(e%nu(l))=convert_cons2prim(temp(l)) !vector of primitive variables in physical values, for n,m
     ENDDO

     !...................................
     ! warning: for convergence study, better to work with phys variables'
     ! because of u2 criterion-> use the right way to compute second derivatives
     ! comment here
     temp=Cons_to_Control(phys_d(e%nu),e)
     phys_d_control(e%nu)=debug%ua(e%nu,k_iter)!temp !vector of primitive variables in control values for n,m
     !end comment
     !...................................
     DEALLOCATE(temp)
  ENDDO


  ListVar: DO lkk=1,n_list 
     lk=list(lkk)


     ! In the following we compute the neighobour cells - but this approach is only good for 1D
     !----------------------------------------------------------------------------------------------------------------
     DO jt = 1,Mesh%nt

        e=Mesh%e(jt)

        !-------------------------------------------------
        !-------------------------------------------------
        ! Definition of the minimum and maximum values of the neighbour elements at time n,0

        !////////////////////////////////////////////////////////////
        !Boundary conditions: Attention you may have to change them according to the test case
#if(1==1)
        !  outflow conditions
        IF (jt > 1) THEN
           eL=Mesh%e(jt-1)
        ELSE
           eL=e
        END IF
        IF (jt < Mesh%nt) THEN
           eR=Mesh%e(jt+1)
        ELSE
           eR=e
        END IF

#else
        ! periodic
        IF (jt > 1) THEN
           eL=Mesh%e(jt-1)
        ELSE
           eL=Mesh%e(Mesh%nt)
        END IF
        IF (jt < Mesh%nt) THEN
           eR=Mesh%e(jt+1)
        ELSE
           eR=Mesh%e(1)
        END IF
#endif
        !////////////////////////////////////////////////////////////

        w_min(jt) = MIN(MINVAL(phys(e%nu)%u(lk)),MINVAL(phys(eL%nu)%u(lk)),MINVAL(phys(eR%nu)%u(lk)))
        w_max(jt) = MAX(MAXVAL(phys(e%nu)%u(lk)),MAXVAL(phys(eL%nu)%u(lk)),MAXVAL(phys(eR%nu)%u(lk)))


        !-------------------------------------------
        !------------------------------------------
        ! Definition of the minimum and maximum values of the neighbour elements' second derivative at time n,0
        ALLOCATE(coefficient(e%nsommets))
        !ALLOCATE(coefficient(e%nsommets-1))
        DO l=1,e%nsommets                 
           coefficient(l)=e%eval_der2(phys_d_control(e%nu), e%x(:,l))
           !coefficient=e%der_sec(phys_control(eL%nu)%u(lk))    
        END DO

        ALLOCATE(coefficient1(eL%nsommets))
        !ALLOCATE(coefficient1(eL%nsommets-1))
        DO l=1,eL%nsommets
           !coefficient1=eL%der_sec(phys_control(eL%nu)%u(lk))           
           coefficient1(l)=eL%eval_der2(phys_d_control(eL%nu), eL%x(:,l))
        END DO

        !ALLOCATE(coefficient2(eR%nsommets-1))
        !coefficient2=eR%der_sec(phys_control(eR%nu)%u(lk))       
        ALLOCATE(coefficient2(eR%nsommets))
        DO l=1,eR%nsommets                 
           coefficient2(l)=eR%eval_der2(phys_d_control(eR%nu), eR%x(:,l))
        END DO

        u2_min(jt)=MIN(MINVAL(coefficient%u(lk)),MINVAL(coefficient1%u(lk)),MINVAL(coefficient2%u(lk)))
        u2_max(jt)=MAX(MAXVAL(coefficient%u(lk)),MAXVAL(coefficient1%u(lk)),MAXVAL(coefficient2%u(lk)))


        ALLOCATE(u1_derL(eL%nsommets))
        DO l=1,eL%nsommets             
           u1_derL(l)=eL%eval_der(phys_d_control(eL%nu),eL%x(:,l))
        END DO

        ALLOCATE(u1_der(e%nsommets))
        DO l=1,e%nsommets             
           u1_der(l)=e%eval_der(phys_d_control(e%nu),e%x(:,l))
        END DO

        ALLOCATE(u1_derR(eR%nsommets))
        DO l=1,eR%nsommets             
           u1_derR(l)=eR%eval_der(phys_d_control(eR%nu),eR%x(:,l))
        END DO

        u1_derL_mean(jt)=sum(u1_derL%u(lk))/REAL(eL%nsommets, DP) !*Added DP
        u1_der_mean(jt)=sum(u1_der%u(lk))/REAL(e%nsommets, DP) !*Added DP
        u1_derR_mean(jt)=sum(u1_derR%u(lk))/REAL(eR%nsommets, DP) !*Added DP
        coefficient_mean(jt)=sum(coefficient%u(lk))/REAL(e%nsommets, DP) !*Added DP

        vL(jt)=u1_der_mean(jt)-e%volume*0.5_DP*coefficient_mean(jt) !*Added DP
        vR(jt)=u1_der_mean(jt)+e%volume*0.5_DP*coefficient_mean(jt) !*Added DP

        vL_min(jt)=min(u1_derL_mean(jt),u1_der_mean(jt))
        vL_max(jt)=max(u1_derL_mean(jt),u1_der_mean(jt))
        vR_min(jt)=min( u1_der_mean(jt),u1_derR_mean(jt))
        vR_max(jt)=max( u1_der_mean(jt),u1_derR_mean(jt))

        DEALLOCATE(coefficient,coefficient1, coefficient2)
        DEALLOCATE(u1_der,u1_derR,u1_derL)

        !=========================================================
        !=========================================================       

     END DO



     !"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
     !"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
     !""""""""""""""" MOOD DETECTION CRITERIA """"""""""""""""""""""""""""""""""""""""""""""""""""
     !"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


     Loop:      DO jt=1,Mesh%nt
        e=Mesh%e(jt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! NaN detection
        !.....................................
        ! at n+1
        IF (.NOT. ALLOCATED(deb)) ALLOCATE(deb(e%nsommets))

        deb=debug%ua(e%nu,k_iter)
        DO l=1, e%nsommets
           DO k=1,n_vars
              IF (deb(l)%u(k).NE.deb(l)%u(k) ) THEN
                 Mesh%e(jt)%diag=NAN_criteria
                 CYCLE loop
              ENDIF
           ENDDO
        ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! PAD test: for fluids, assumes: 1-> density, n_vars-> pressure
        !...................................................................................................................................................
        IF (MINVAL(phys_d(e%nu)%u(1)).LE.0. .OR. MINVAL(phys_d(e%nu)%u(n_vars)).LE.0.) THEN
           Mesh%e(jt)%diag=PAD_criteria

           CYCLE loop
        ENDIF

        val_min=MINVAL(phys_d(e%nu)%u(lk))  ! at n+1
        val_max=MAXVAL(phys_d(e%nu)%u(lk))
        val_min_n=MINVAL(phys(e%nu)%u(lk))  ! at n
        val_max_n=MAXVAL(phys(e%nu)%u(lk))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! plateau detection
        !...............................................
        IF (ABS(w_max(jt)-w_min(jt)) .LT. e%volume**3) THEN
           Mesh%e(jt)%diag= plateau
           CYCLE Loop
        ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Extrema detection DMP+u2
        !.....................................................................
        eps=coeff *MAX(ABS(w_max(jt)-w_min(jt))*eps1,eps2)   ! at n

        IF (val_min.LT.w_min(jt)-eps.OR.val_max.GT.w_max(jt)+eps) THEN
!!$

           IF (e%nsommets .eq. 2) then
              Mesh%e(jt)%diag=DMP_B1
              CYCLE Loop
           ELSE

              IF (vL(jt) .gt. u1_der_mean(jt)) THEN
                 ratio=(vL_max(jt)-u1_der_mean(jt))/(vL(jt)-u1_der_mean(jt))
                 alphaL=min(1._DP, ratio) !*Added DP
              ELSE IF  (vL(jt) .eq. u1_der_mean(jt)) THEN

                 alphaL=1._DP !*Added DP
              ELSE
                 ratio=(vL_min(jt)-u1_der_mean(jt))/(vL(jt)-u1_der_mean(jt))
                 alphaL=min(1._DP,ratio) !*Added DP
              END IF

              IF (vR(jt) .gt. u1_der_mean(jt)) THEN
                 ratio=(vR_max(jt)-u1_der_mean(jt))/(vR(jt)-u1_der_mean(jt))
                 alphaR=min(1._DP, ratio) !*Added DP
              ELSE IF  (vR(jt) .eq. u1_der_mean(jt)) THEN

                 alphaR=1._DP !*Added DP
              ELSE
                 ratio=(vR_min(jt)-u1_der_mean(jt))/(vR(jt)-u1_der_mean(jt))
                 alphaR=min(1._DP,ratio) !*Added DP
              END IF

              smooth_sensor=min(alphaR,alphaL)

              IF (abs(smooth_sensor) .ge. 1._DP-Tiny(1.0_DP)) THEN
                 Mesh%e(jt)%diag=DMP_nou2
            
                 CYCLE Loop
              ELSE
                 Mesh%e(jt)%diag=DMP_u2_2
                 CYCLE Loop
              END IF

           ENDIF


!!!!!!!!!!!!!!

        ELSE
           Mesh%e(jt)%diag=0
           cycle Loop
        END IF
        ! END DO Loop_loc
        !DEALLOCATE(coefficient,coefficient1,coefficient2,u2_minloc,u2_maxloc)
        DEALLOCATE(deb)
     ENDDO Loop

  END DO ListVar
  DEALLOCATE(phys,phys_d, phys_control, phys_d_control)
  DEALLOCATE(v_min,v_max,w_min,w_max,u2_min,u2_max,vd_min,vd_max)
  DEALLOCATE(VL, vR, vL_min, vL_max, vr_min, vr_max)
  DEALLOCATE(u1_derL_mean,u1_der_mean,u1_derR_mean,coefficient_mean)

  !-----------------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------------
  !---------------- Update of the Indicators that allow the scheme switch ---------------------------

  DO jt=1, Mesh%nt

     SELECT CASE(Mesh%e(jt)%diag)
     CASE(0,Plateau,DMP_nou2)
        ! In this case we don't touch the indicator
     CASE(DMP_B1,DMP_u2_1,DMP_u2_2)
        Mesh%e(jt)%type_flux=flux_mood
        ! In this case we switch to a more diffusive scheme
     CASE(PAD_criteria,NAN_criteria)
        Mesh%e(jt)%type_flux=theflux2
        ! In this case we take a first order monotone scheme
     CASE default
        PRINT*, "in test, bad behavior, jt=", jt
        STOP
     END SELECT
     !
  ENDDO
END SUBROUTINE test
END MODULE timestepping








