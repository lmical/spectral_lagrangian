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
!*scheme contains the procedures called in timestepping which involve the computations of the updating part for the DeC
!*We recall here the updating formula of the DeC
!*u_i^m,(p)=u_i^m,(p-1)-1/|C_i|*{\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))+\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}]} (#)
!*
!*The boundary residuals i.e.
!*\Sum_{\Gamma \in \Gamma_i RES_i^\Gamma(u^l,(p-1))}
!*are almost always not present in 1d, in any case they are not computed here 
!*
!*The updating part computed here is
!*\Sum_{K \in K_i} \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[\Sum_{K \in K_i} RES_i^K(u^l,(p-1))]
!*
!*-> schema is called inside a loop on the elements and computes
!*\Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) + dt*\sum_{l=0}^M theta_l^m*[ RES_i^K(u^l,(p-1))]
!*
!*The sum \Sum_{K \in K_i} is implicit in the fact that we add this contribution to each node for each element in the loop
!*
!*NB: In the second piece 
!*dt*\sum_{l=0}^M theta_l^m*[ RES_i^K(u^l,(p-1))]
!*we will not have any combination through the theta coefficients, we pass the coefficients of the flux (and of the source) and the values of u already combined and we rely on the linearity of the operators in which the flux and u are involved 
!*
!*So in schema we compute 
!*\Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0)
!*dt*[ RES_i^K(UC)] where UC means that everything in input here is already combined
!*
!*schema recalls (depending on the parameter e%type_flux) the correct scheme i.e. the specific residuals that we want to use
!*(4)galerkin_new
!*(5)psi_galerkin
!*(-1)lxf
!*NB: The residuals are computed up to the edge stabilization based on the jump of the gradient, which as I explained in edge_main_update in timestepping, it is treated in a "particular" way
!*
!*->jump computes the edge stabilization based on the jump of the gradient
!*
!*
!*----------------------------------------------------------------------------------------
MODULE scheme !1D
  ! In this module all scheme for the flux are collected.
  ! Structure of this module:
  ! - schema -> allows to choose between different scheme for the flux approximation,
  !                 which are coded in this same module and are the following
  !                 a)  : galerkin as such, it is not quadrature free for non linear problems
  !                     : galerkin_new:  this one is quadrature free
  !                 b) galerkin + psi limiting, for galerkin : see above
  !                 c) local Lax-Friedrichs (=Rusanov). use Galerkin as central part, so see a)
  ! - stabilization: Burman's jump
  ! - limiting: generic function to activate the limiting,
  !              which can be done either via conservative variables, or characteristic ones
  USE overloading
  USE element_class
  USE variable_def
  USE arete_class
  USE param2d
  USE Model
  USE PRECISION
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: schema, jump

CONTAINS

  !*--------------------------------------------------------------------------------------
  !*schema computes for each node i of the processed element
  !*\Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0)
  !*dt*[ RES_i^K(UC)] where UC means that everything in input here is already combined
  !*schema recalls (depending on the parameter e%type_flux) the correct scheme i.e. the specific residuals that we want to use
  !*(4)galerkin_new
  !*(5)psi_galerkin
  !*(-1)lxf
  !*NB: The residuals are computed up to the edge stabilization based on the jump of the gradient, which as I explained in edge_main_update in timestepping, it is treated in a "particular" way
  !*--------------------------------------------------------------------------------------
  FUNCTION schema(ischema, e, u, du, flux, source, dt, jt,mesh) RESULT(res)
    ! input:
    ! e: element
    ! du: temporal increment !*OUR (u_j^m,(p-1)-u_j^m,0) already computed
    ! flux: values of the flux at the dofs !*NOT VALUES BUT COEFFICIENTS, already combined through the thetas
    !*u: VALUES of the solution at the DoFs, already combined through the thetas
    ! dt: time step
    ! output:
    ! residual (supg or psi) !*NEIN, probably the function has been recycled from old codes where one had to compute only the residuals now we have also \Sum_{x_j \in K} \int_K phi_i phi_jdx (u_j^m,(p-1)-u_j^m,0) i.e. the mass matrix multiplied by the increment
    !*FURTHER, for what concerns the spatial part we have the following residuals 
    !*(4)galerkin_new
    !*(5)psi_galerkin
    !*(-1)lxf
    !*
    !*RMK
    !*NB: We are entering in this function
    !*flux COEFFICIENTS combined through the thetas
    !*source COEFFICIENTS combined through the thetas
    !*u COEFFICIENTS, combined through the thetas
    !*du variation of the COEFFICIENTS in time
    INTEGER, INTENT(IN):: jt
    TYPE(maillage),INTENT(IN):: mesh
    INTEGER, INTENT(in):: ischema
    TYPE(element),                     INTENT(in):: e
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: u
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: du, source
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux
    REAL(dp)   ,                        INTENT(in):: dt

    TYPE(PVar), DIMENSION(e%nsommets)            :: res

    SELECT CASE (e%type_flux) !*Depending on the scheme we want he calls the proper function
    CASE(4)
       res=galerkin_new(e, u, du, flux, source, dt,jt,mesh)
    CASE(5)
       res=psi_galerkin(e, u, du, flux, source, dt,jt,mesh)
    CASE(-1)
       res=lxf(e, u, du, flux, source, dt,jt,mesh)
    CASE default
       PRINT*, "scheme not defined"
       STOP
    END SELECT
  END FUNCTION schema

  !----------------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------------
  !------------------------ SCHEME FOR FLUX ---------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------------
  !*RMK
  !*flux COEFFICIENTS combined through the thetas
  !*source COEFFICIENTS combined through the thetas
  !*u COEFFICIENTS, combined through the thetas
  !*du variation of the COEFFICIENTS in time





  !*--------------------------------------------------------------------------------------
  !*galerkin computes the galerkin residuals through the divergence theorem (integration by parts in 1d)
  !*\int_e varphi_i *du - dt* \int_e nabla varphi_i . flux + [f*varphi_i]
  !*NB: It is not used anymore, now we use galerkin_new where the integration is performed integrating directly phi div f without passing the derivative on the test function 
  !*--------------------------------------------------------------------------------------
  FUNCTION galerkin(e, u, du, flux, source, dt,jt,mesh) RESULT(res)
    ! input:
    ! e: element
    ! du: temporal increment
    ! flux: values of the flux at the dofs !*NEIN, COEFFICIENTS
    ! dt: time step
    ! output:
    ! \int_e varphi_i *du - dt* \int_e nabla varphi_i . flux !*WRONG
    !*
    !*CORRECT
    !*\int_e varphi_i *du + dt* \int_e div(flux) * varphi_i
    !*OR, applying the divergence theorem (SO, WHAT WE DO HERE)
    !*\int_e varphi_i *du - dt* \int_e nabla varphi_i . flux + [f*varphi_i] <- This is actually the version that we have here
    INTEGER, INTENT(IN):: jt
    TYPE(maillage),INTENT(IN):: mesh
    TYPE(element),                     INTENT(in):: e
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: u
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: du, source
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux
    REAL(dp)  ,                        INTENT(in):: dt

    TYPE(PVar), DIMENSION(e%nsommets)            :: res

    REAL(dp), DIMENSION(e%nsommets):: base
    REAL(dp), DIMENSION(n_dim,e%nsommets):: grad
    REAL(dp), DIMENSION(n_dim):: vit
    REAL(dp), PARAMETER, DIMENSION(n_dim)::xx=0._dp
    TYPE(Pvar),DIMENSION(n_dim):: div_flux

    REAL(dp), DIMENSION(n_vars,n_vars,n_dim):: Jac
    REAL(dp), DIMENSION(n_vars,n_vars):: nnmat
    TYPE(PVar) :: du_loc, flux_loc, source_loc

    REAL(dp)   :: hh, eps
    INTEGER   :: i, iq, l, k
    TYPE(Pvar),DIMENSION(n_dim) :: fluxloc

    TYPE(PVar) :: divflux
    REAL(dp),  DIMENSION(2,2):: xp


    xp(1,1)=0.0_dp; xp(2,1)=1.0_dp; xp(1,2)=1.0_dp; xp(2,2)=0.0_dp

    !*Initialization
    DO l=1, e%nsommets
       res(l)=0._dp
    ENDDO

    !*Loop on the DoFs, i index of the node for which we compute the residual
    DO i=1, e%nsommets
       !*Loop on the quadrature points because we need to perform an integration
       DO iq=1, e%nquad
          !*Compute the basis functions and their gradients in the processed quadrature point 
          DO l=1, e%nsommets
             base(l  )=e%base    (l, e%quad(:,iq))
             grad(:,l)=e%gradient(l, e%quad(:,iq))
          ENDDO


          du_loc  = e%eval_func(du      ,e%quad(:,iq)) !*THIS IS CORRECT because here we have the coefficients
          source_loc  = e%eval_func(source      ,e%quad(:,iq)) !*THIS IS CORRECT because here we have the coefficients



          !*CORRECT
          ! quadrature free !correct
          fluxloc(1)=e%eval_func(flux(1,:),e%quad(:,iq)) !*THIS IS CORRECT: from the coefficients f the flux to the value in the quadrature point


          !*QUADRATURE SUM
          res(i) = res(i) + ( base(i)*du_loc-dt*grad(1,i)*fluxloc(1) -dt*source_loc*base(i) )*e%weight(iq)

          !*base(i)*du_loc CORRECT
          !*-dt*grad(1,i)*fluxloc(1) CORRECT
          !*-dt*source_loc*base(i) CORRECT

       ENDDO ! iq
       res(i) = res(i)* e%volume !*To achieve the integral on the correct domain
    ENDDO ! i

    !*We need to add [f*varphi_i]=DX-SX
    !*RMK:
    !*xp(:,:) are barycentric coordinates
    !*second index referred to the vertices
    !*first index barycentric coordinates of the vertex
    !*Left vertex xp(:,1)=(0,1)
    !*xp(1,1)=0.0_dp; xp(2,1)=1.0_dp; 
    !*Right vertex xp(:,2)=(1,0)
    !*xp(1,2)=1.0_dp; xp(2,2)=0.0_dp
    !*RMK:
    !*x->barycentric coordinates
    !*Barycentric coordinates of the left vertex (0,1)
    !*Barycentric coordinates of the right vertex (1,0)
    !*So the second barycentric coordinate is associated to the left/first vertex, the first one is associated to the right/second vertex 
    !*RMK: To visualize it in a logical way, in the reference interval
    !*x(1)=x even if associated to the second vertex
    !*x(2)=1-x
    eps=-dt ! 
    DO k=1,2
       !*Let's focus on the first cycle
       DO l=1, e%nsommets
          base(l)=e%base(l,xp(:,k)) !*All the basis functions in the left vertex (SX) 
       ENDDO!*NOT USED !*u_diff(i)=(u(i)-vbar(i)) 
       fluxloc=e%eval_func(flux(1,:),xp(:,k)) !*Flux in the left vertex
       res(:) = res(:)+eps*fluxloc(1)*base !*Correct
       eps = -eps !*Next time (DX vertex) it will be added with +
    ENDDO


  END FUNCTION galerkin

  !*--------------------------------------------------------------------------------------
  !*galerkin_new computes the galerkin residuals integrating directly phi div f (and not through the divergence theorem (integration by parts in 1d))
  !*\int_e varphi_i *du - dt* \int_e varphi_i div(flux) 
  !*--------------------------------------------------------------------------------------
  FUNCTION galerkin_new(e, u, du, flux, source, dt,jt,mesh) RESULT(res)
    ! input:
    ! e: element
    ! du: temporal increment !*Coefficients
    ! flux: values of the flux at the dofs !*Coefficients
    ! dt: time step
    ! output:
    ! \int_e varphi_i *du - dt* \int_e nabla varphi_i . flux !*WRONG
    !*
    !*CORRECT
    !*\int_e varphi_i *du + dt* \int_e div(flux) * varphi_i
    
    INTEGER, INTENT(IN):: jt
    TYPE(maillage),INTENT(IN):: mesh
    TYPE(element),                     INTENT(in):: e
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: u
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: du, source
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux
    REAL(dp)  ,                        INTENT(in):: dt

    TYPE(PVar), DIMENSION(e%nsommets)            :: res



    INTEGER   ::  l

    !*RMK
    !*e%mass(i,j)-> local mass matrix int_K phi_i phi_j
    !*e%coeff(i,j) -> another matrix int_K phi_i  d_x phi_j !*NB:DERIVATIVE ON THE SECOND INDEX

    DO l=1, SIZE(flux(1,1)%u) !*Loop on the components
       res(:)%u(l)=MATMUL(e%mass(:,:),du%u(l)-dt*source%u(l)) & !*The index of the columns runs
            &+dt*MATMUL(e%coeff   (:,:), flux(1,:)%u(l) )!& !*The index of the columns runs

            !*coeff_b NEVER USED, IT WAS ALREADY COMMENTED
       !       &+dt*MATMUL(e%coeff_b (:,:), flux(1,:)%u(l) )
    ENDDO
   
  


  END FUNCTION galerkin_new


  !*--------------------------------------------------------------------------------------
  !*Lax Friedrichs = Galerkin + average stabilization
  !*--------------------------------------------------------------------------------------
  FUNCTION lxf(e, u, du, flux, source, dt,jt,mesh) RESULT(phi_lxf)
    ! input:
    ! e: element
    ! du: temporal increment
    ! flux: values of the flux at the dofs
    ! dt: time step
    ! output:
    ! \int_e varphi_i *du - dt* \int_e nabla varphi_i . flux

    INTEGER, INTENT(IN):: jt
    TYPE(maillage),INTENT(IN):: mesh
    TYPE(element),                     INTENT(in):: e
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: u
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: du, source
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux
    REAL(dp)   ,                        INTENT(in):: dt

    TYPE(PVar), DIMENSION(e%nsommets)            :: res, phi_lxf,u_cons
    TYPE(PVar):: u_consL, u_consR
    REAL(dp),DIMENSION(n_dim):: xx=0._dp
    TYPE(PVar), DIMENSION(e%nsommets):: v,u_diff,u_diff2
    TYPE(PVar) ::  ubar
    TYPE(PVar), DIMENSION(e%nsommets)::vbar
    REAL(dp)::   alpha
    INTEGER:: l,i,iq, k

    res= galerkin_new(e, u, du, flux, source, dt,jt,mesh)

    u_cons=Control_to_cons(u,e) !*From coefficient to values

    DO l = 1,n_vars
       ubar%u(l)= SUM(u%u(l))/e%nsommets !*Average of the coefficients at the DoFs
    END DO

    DO l=1,e%nsommets
       vbar(l)=ubar !*Copying ubar in every DoF in the variable vbar
    END DO
    !*DO i=1,e%nsommets
       !*NOT USED !*u_diff(i)=(u(i)-vbar(i)) 
       alpha=ubar%spectral_radius(xx,e%n) 
    !*END DO
    alpha=ubar%spectral_radius(xx,e%n) 
    phi_lxF(:) = res(:) + dt*alpha*(u(:)-vbar(:)) 

  END FUNCTION lxf




  !*--------------------------------------------------------------------------------------
  !*Psi Galerkin = Lax Fridriechs + additional stabilization (limit)
  !*--------------------------------------------------------------------------------------
  FUNCTION psi_galerkin(e, u, du, flux, source, dt,jt, Mesh) RESULT(res)
    ! input:
    ! e: element
    ! du: temporal increment
    ! flux: values of the flux at the dofs
    ! dt: time step
    ! output:
    ! \int_e varphi_i *du - dt* \int_e nabla varphi_i . flux
    !
    INTEGER, INTENT(IN):: jt
    TYPE(maillage),INTENT(IN):: mesh
    TYPE(element),                     INTENT(in):: e
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: u
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: du, source
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux
    REAL(dp)   ,                        INTENT(in):: dt
    TYPE(PVar), DIMENSION(e%nsommets)            :: res, phi_lxF,w
    INTEGER::  k, l

    phi_lxf= lxf(e,u,du,flux, source,dt,jt,mesh) 
    CALL limit(e,u,res,phi_lxf)

  END FUNCTION psi_galerkin


  !----------------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------------
  !------------------------ STABILIZATION ------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------------
!!!!!!!!!! Jump
 !*---------------------------------------------------------------------------------------
 !*Jump computes the stabilization part of the residuals
 !*This is treated in a "particular way", we take in input 
 !*ed, an edge
 !*e1,e2, the two elements sharing the edge
 !*uu1, uu2, the combination through the thetas of the coefficients of the solution in the two elements from all the subtimesteps
 !*alpha stabilization on the first derivative
 !*alpha2 stabilization on the second derivative 
 !*
 !*We call jump in a loop over the (internal) edges
 !*We consider the DoFs "involved" in the edge ed, the DoFs whose gradient has a jump over f=ed which are basically the DoFs contained in the two elements sharing the edge f=ed and we compute
 !*\int_f [grad u][grad phi_i] (later added to Var%un after multiplying by dt)
 !*Basically in edge_main_update we switched the loop on the DoFs with the loop on the edges
 !*---------------------------------------------------------------------------------------
 SUBROUTINE jump(ed ,e1, e2, uu1, uu2, resJ, alpha, alpha2)
    ! e1 is the element before the edge, e2 the one after
    ! e1-> uu1, e2->uu2

    REAL(dp),PARAMETER,DIMENSION(6)::TYPE=(/1._DP,2._DP,2._DP,3._DP,3._DP,4._DP/)
    REAL(dp), INTENT(in) :: alpha, alpha2
    TYPE(arete), INTENT(in):: ed
    TYPE(element),                     INTENT(in):: e1
    TYPE(element),                     INTENT(in):: e2
    TYPE(PVar), DIMENSION(e1%nsommets), INTENT(in):: uu1, uu2
    TYPE(Pvar), DIMENSION(e1%nsommets,2):: resJ
    TYPE(Pvar):: u, ubar1, ubar2
    TYPE(Pvar):: divflux1, divflux2, jumpflux,jumpfluxder2, divdivflux1, divdivflux2
    INTEGER:: is, l
    real(dp),dimension(n_dim):: xx=0._DP
    REAL(dp):: theta, theta2, h,  umax, umin, l2, spectral_radius1, spectral_radius2
    REAL(dp), DIMENSION(n_vars,n_vars,n_dim)::jaco
    REAL(dp), DIMENSION(e1%nsommets):: base1
    REAL(dp),DIMENSION(n_dim,e1%nsommets):: grad1, grad1_der
    REAL(dp), DIMENSION(e2%nsommets):: base2
    REAL(dp),DIMENSION(n_dim,e2%nsommets):: grad2, grad2_der
    REAL(dp),DIMENSION(2), PARAMETER:: x2=(/0._DP,1._DP/),x1=(/1._DP,0._DP/)
    !*Barycentric coordinates of the left vertex (0,1) x2
    !*Barycentric coordinates of the right vertex (1,0) x1

    
    DO l = 1,n_vars
       ubar1%u(l)= SUM(uu1%u(l))/REAL(SIZE(uu1),DP) !*Average in element 1
    END DO
    DO l = 1,n_vars
       ubar2%u(l)= SUM(uu2%u(l))/REAL(SIZE(uu2),DP) !*Average in element 2
    END DO

    spectral_radius1= ubar1%spectral_radius(xx,e1%n) !*Spectral radius in ubar1
    spectral_radius2= ubar2%spectral_radius(xx,e2%n) !*Spectral radius in ubar2

    !*RMK:
    !*ed has two elem(1,0) x1ents which share it
    !*Mesh%edge(indi)%jt1=indi-1 -> e1 is the left element
    !*Mesh%edge(indi)%jt2=indi -> e2 is the right element
    !*EX:
    !*            e1            ed              e2
    !* ------------------------ . -----------------------------
    !*
    !*So... 
    !*For e1 ed is the right edge whith barycentric coordinates (1,0) x1
    !*For e2 dd is the left edge whith barycentric coordinates (0,1) x2
    DO l=1, e1%nsommets !*Loop on the basis functions
       !*DERIVATIVES OF THE BASIS FUNCTIONS
       grad1(:,l) = e1%gradient(l,x1) !*First derivative in the right vertex of e1 OK 
       grad2(:,l) = e2%gradient(l,x2) !*First derivative in the left vertex of e2 OK
       grad1_der(:,l) = e1%gradient2(l,x1) !*Second derivative in the right vertex of e1 OK
       grad2_der(:,l) = e2%gradient2(l,x2) !*Second derivative in the left vertex of e2 OK
    ENDDO

    ! jump of grad u
    !*DERIVATIVES OF THE SOLUTION
    !*eval_der_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the derivative of the solution in y 
    !*eval_der2_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the SECOND derivative of the solution in y 
    divflux1 = e1%eval_der(uu1(:),x1) !*First derivative in the right vertex of e1
    divflux2 = e2%eval_der(uu2(:),x2) !*First derivative in the left vertex of e2
    divdivflux1 = e1%eval_der2(uu1(:),x1) !*Second derivative in the right vertex of e1
    divdivflux2 = e2%eval_der2(uu2(:),x2) !*Second derivative in the left vertex of e2

    !! h=1./((TYPE(e1%itype)*TYPE(e2%itype))*0.5) !* l2!

    h=1._DP/(0.5_DP*(SUM(ABS(grad1(1,:))) + SUM( ABS(grad2(1,:))) ))

    !*Jump of the solution (el2-el1)
    !*NB: The following are vectors 
    jumpflux=(divflux2-divflux1)   !*Jump of the first derivative over the edge
    jumpfluxder2=(divdivflux2-divdivflux1) !*Jump of the second derivative over the edge

    !*"Consistent" stabilization parameters
    theta=alpha*h**2 
    theta2=alpha2*h**4 

    DO l=1, e1%nsommets !*Loop on the DoFs
       ResJ(l,1) = spectral_radius1* &
                   &(jumpflux *(-1._DP)* grad1(1,l) * theta + &
                   & jumpfluxder2*(-grad1_der(1,l))*theta2)
                   !*jumpflux is el2-el1 so we take grad1 and grad1_der with sign -
                   !*We have to be coherent and take Outside - Inside in the two jumps
                   !*1 is inside in both cases  
       ResJ(l,2) = spectral_radius2* &
                   &(jumpflux * grad2(1,l) * theta + &
                   & jumpfluxder2*grad2_der(1,l)*theta2)
                   !*We are coherent also in this case el2-el1 so derivatives taken with sign +
    ENDDO

    !*NB: Most of the treated DoFs receive just one contribution, the one related to the element they belong to, because they are 0 outside that element and that contribution is not counted (the one coming from outside which is null). Just the DoF on the boundary receives the contribution of both the two elements (it actually has a nontrivial jump, it is not 0 on one of the two sides), anyway the two contributions are in ResJ(,1) and ResJ(,2), they will be put together outside this function, when these stabilization residuals will be added to Var%un (after the multiplication by dt)

    !*NB: The computation is very "close" to their definition as residuals
    !*If we consider their contribution to the nodal/cell residual we have
    !*contribution_i^K=int_f [grad u] grad phi_i|_K 
  END SUBROUTINE jump




  !----------------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------------
  !------------------------ LIMITING ---------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------------

  SUBROUTINE limit (e,u,res,phi)
    !using Characteristic variables
    IMPLICIT NONE
    TYPE(element),               INTENT(in):: e
    TYPE(Pvar), DIMENSION(:), INTENT(in):: u
    TYPE(Pvar), DIMENSION(:), INTENT(out):: res
    TYPE(Pvar), DIMENSION(:), INTENT(in)   :: phi
    LOGICAL:: out

    CALL limit_char(e,u,res,phi,out)
    !if NaN detected, then do limitation by variables
    IF (out) THEN
       CALL limit_var(u,res,phi)
    ENDIF
  END SUBROUTINE limit




  SUBROUTINE limit_char(e,u,res,phi,out)
    !using Characteristic variables
    IMPLICIT NONE
    TYPE(Pvar), DIMENSION(:), INTENT(in):: u
    TYPE(element),               INTENT(in):: e
    TYPE(Pvar), DIMENSION(:), INTENT(out):: res
    TYPE(Pvar), DIMENSION(:), INTENT(in)   :: phi
    LOGICAL, INTENT(out):: out
    TYPE(Pvar), DIMENSION(SIZE(phi)) :: phi_ch, res_ch
    REAL(dp)::  l, phi_tot, den,  truc
    REAL(dp), DIMENSION(SIZE(res)):: x
    INTEGER:: i, k
    REAL(dp), DIMENSION(n_Vars,n_Vars) :: EigR, EigL
    TYPE(Pvar) :: ubar, vbar
    TYPE(Pvar), DIMENSION(e%nsommets) :: v
    REAL(dp), DIMENSION(1) ::  v_nn=(/1.0_dp/)
    REAL(dp), PARAMETER:: eps=1.e-9_dp
    out=.FALSE.
    res=0._dp

    v= Control_to_Cons(u,e)
    DO k = 1,n_vars
       ubar%u(k)= SUM(v%u(k))/e%nsommets
    END DO

    EigR = ubar%rvectors(v_nn)
    IF (SUM(ABS(EigR)).NE.SUM(ABS(Eigr)) ) THEN !*Essentially if we get NaN because for Fortran NaN is not equal to NaN. By definition, NAN is not equal to anything, even itself.
       out=.TRUE.
       RETURN
    ENDIF
    EigL = ubar%lvectors(v_nn)


    DO i = 1,SIZE(res)
       phi_ch(i)%u = MATMUL(EigL,phi(i)%u)
    END DO
    !print*, 'phi_ch',maxval(abs(phi_ch(:)%u(1)))

    DO k = 1,n_vars

       phi_tot=SUM(phi_ch(:)%u(k))
       !                 print*, k,((phi_ch(:)%u(1))),phi_tot
       IF (ABS(phi_tot)>TINY(1._dp)) THEN

          x=MAX(phi_ch(:)%u(k)/phi_tot,0._dp)
          den=SUM(x)
          truc=SUM(ABS(phi_ch(:)%u(k))) 

          l=ABS(phi_tot)/truc

          res(:)%u(k)=(1._dp-l)*(Phi_tot*x/den)+l*phi_ch(:)%u(k)
       ELSE
          res(:)%u(k)=0.0_dp 
       ENDIF
    ENDDO
    !read*
    DO i = 1,SIZE(res)
       phi_ch(i)%u=MATMUL(EigR,res(i)%u)
    END DO
    !    print*, phi_ch(:)%u(1)
    res=phi_ch
    !print*
  END SUBROUTINE limit_char

  SUBROUTINE limit_var(u,res,phi)
    !using conserved variables
    IMPLICIT NONE
    TYPE(Pvar), DIMENSION(:), INTENT(in):: u
    TYPE(Pvar), DIMENSION(:), INTENT(out):: res
    TYPE(Pvar), DIMENSION(:), INTENT(in)   :: phi
    REAL(dp):: l1, l2, l,phi_tot, den, umax, umin, truc
    REAL(dp),DIMENSION(SIZE(res)):: x
    INTEGER:: i, k

    DO k = 1,n_vars
       L=0._dp
       res%u(k) = 0.0_dp
       phi_tot=SUM(phi(:)%u(k))
       IF (ABS(phi_tot)>TINY(1.0_dp)) THEN
          x=MAX(phi(:)%u(k)/phi_tot,TINY(1.0_dp))
          den=SUM(x)+TINY(1.0_dp)
          truc=SUM(ABS(phi(:)%u(k)))
          l=ABS(phi_tot)/(SUM(ABS(phi(:)%u(k))))
          res(:)%u(k)= (1.0_dp-l)*phi_tot * x/den+l*phi(:)%u(k)
       ENDIF
    END DO ! k

  END SUBROUTINE limit_var

END MODULE scheme
