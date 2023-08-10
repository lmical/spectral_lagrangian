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
  FUNCTION schema(ischema, e, u, du, flux, source, dt, jt,mesh,flux_velocity, Heta_matrix, source_extra) RESULT(res)
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
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: du, source, source_extra
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux
    REAL(dp)   ,                        INTENT(in):: dt

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*For wb
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux_velocity !*Velocity part of the flux
    REAL(dp), DIMENSION(e%nsommets,e%nsommets):: Heta_matrix !*Matrix made by H_j (H_k+b_k) to integrate H*grad(H)+H*grad(b) 
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    TYPE(PVar), DIMENSION(e%nsommets)            :: res

    SELECT CASE (e%type_flux) !*Depending on the scheme we want he calls the proper function
    CASE(4)
       res=galerkin_new(e, u, du, flux, source, dt,jt,mesh)
    CASE(5)
       res=psi_galerkin(e, u, du, flux, source, dt,jt,mesh)
    CASE(-1)
       res=lxf(e, u, du, flux, source, dt,jt,mesh)
    CASE(14)
       res=galerkin_new_wb(e, u, du, flux_velocity, Heta_matrix, dt, jt, mesh, source_extra)
    CASE(15)
       res=psi_galerkin_wb(e, u, du, flux_velocity, Heta_matrix, dt, jt, mesh, source_extra)
    CASE(-11)
       res=lxf_wb(e, u, du, flux_velocity, Heta_matrix, dt, jt, mesh, source_extra)
    CASE(24)
       res=galerkin_new(e, u, du, flux, source, dt,jt,mesh)
    CASE(25)
       res=psi_galerkin(e, u, du, flux, source, dt,jt,mesh)
    CASE(-21)
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
  !*--------------------------------------------------------------------------------------
  !*Jump computes the stabilization part of the residuals
  !*We call jump in a loop over the (internal) edges
  !*--------------------------------------------------------------------------------------
 SUBROUTINE jump(DATA, ed ,e1, e2, uu1, uu2, resJ, alpha, alpha2,ueta1,ueta2,jump3coeff1, jump3coeff2,&
            & u1k,u2k,A0,b1,b2,alphaconsistency,A1,g1,g2)
    ! e1 is the element before the edge, e2 the one after
    ! e1-> uu1, e2->uu2

    TYPE(donnees), INTENT(in) :: DATA
    REAL(dp), INTENT(in) :: alpha, alpha2
    TYPE(arete), INTENT(in):: ed
    TYPE(element),                     INTENT(in):: e1
    TYPE(element),                     INTENT(in):: e2
    TYPE(PVar), DIMENSION(e1%nsommets), INTENT(in):: uu1, ueta1, jump3coeff1, u1k, g1
    TYPE(PVar), DIMENSION(e2%nsommets), INTENT(in):: uu2, ueta2, jump3coeff2, u2k, g2
    REAL(DP), DIMENSION(2,2,1), INTENT(in):: A0,A1
    REAL(DP), DIMENSION(e1%nsommets), INTENT(in):: b1
    REAL(DP), DIMENSION(e2%nsommets), INTENT(in):: b2
    REAL(dp), INTENT(in) :: alphaconsistency

    TYPE(Pvar), DIMENSION(e1%nsommets,2):: resJ

    REAL(DP), DIMENSION(2,2) :: M0,M1 !*Proper matrix in which to copy A0

    REAL(dp)                        :: spectralradiusatinterface
    real(dp),DIMENSION(2), PARAMETER:: xx=(/1._DP,0._DP/)

    M0=A0(:,:,1)
    M1=A1(:,:,1)

    spectralradiusatinterface=uu1(2)%spectral_radius(xx,e1%n) !*Coefficient BUT also value because boundary of the element ;)

    SELECT CASE (DATA%ijump) !*Depending on the scheme we want he calls the proper function
    !*1,2,3,4,5,6  h defined by Remi
    !*11,12,13,14,15,16  h defined by Mario h=dx
    CASE(1,11) !*Jump of the gradient of the conserved variables
       CALL jump_burman( ed, e1, e2, uu1, uu2, resJ, alpha,alpha2,DATA,spectralradiusatinterface)
    CASE(2,12) !*Jump of the gradient of the conserved variables BUT with eta instead of h
       CALL jump_burman( ed, e1, e2, ueta1, ueta2, resJ, alpha,alpha2,DATA,spectralradiusatinterface) !*ueta instead of uu
    CASE(3,13) !*jump of DUDW*w where u are the conserved variables and w are the entropy variables
       CALL jump_burman( ed, e1, e2, jump3coeff1, jump3coeff2, resJ, alpha,alpha2,DATA,spectralradiusatinterface) !*jumpcoeff instead of uu
    CASE(4,14) !*jump [[ J grad(phi_i) ]] A0 [[ J grad(u) + S ]]
       !*RMK: We directly consider the subtimestep k, no combination through the thetas
       CALL jump_4( ed, e1, e2, u1k, u2k, resJ, alpha,alpha2,M0,b1,b2,alphaconsistency,DATA)      !*<----M0 (|J|^{-1})
    CASE(5,15) !*jump [[ J grad(phi_i) ]] A1 [[ J grad(u) + S ]]
       !*RMK: We directly consider the subtimestep k, no combination through the thetas
       CALL jump_4( ed, e1, e2, u1k, u2k, resJ, alpha,alpha2,M1,b1,b2,alphaconsistency,DATA)      !*<----M1 (rho^{-1})
    CASE(6,16) !*jump of the gradient of the global flux, [[ J grad(phi_i) ]] A1 [[ grad(g) ]]
       CALL jump_4_global_flux( ed, e1, e2, u1k, u2k, resJ, alpha,alpha2,M0,g1,g2,DATA)           !*<----M0 (|J|^{-1})
    CASE(7,17) !*jump of the gradient of the global flux, [[ J grad(phi_i) ]] A1 [[ grad(g) ]]
       CALL jump_4_global_flux( ed, e1, e2, u1k, u2k, resJ, alpha,alpha2,M1,g1,g2,DATA)           !*<----M1 (rho^{-1})
    CASE default
       PRINT*, "jump not defined"
       STOP
    END SELECT



 END SUBROUTINE jump


 !*---------------------------------------------------------------------------------------
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
 SUBROUTINE jump_burman(ed ,e1, e2, uu1, uu2, resJ, alpha, alpha2,DATA,spectralradiusatinterface)
    ! e1 is the element before the edge, e2 the one after
    ! e1-> uu1, e2->uu2

    REAL(dp),PARAMETER,DIMENSION(6)::TYPE=(/1._DP,2._DP,2._DP,3._DP,3._DP,4._DP/)
    REAL(dp), INTENT(in) :: alpha, alpha2
    TYPE(arete), INTENT(in):: ed
    TYPE(element),                     INTENT(in):: e1
    TYPE(element),                     INTENT(in):: e2
    TYPE(PVar), DIMENSION(e1%nsommets), INTENT(in):: uu1
    TYPE(PVar), DIMENSION(e2%nsommets), INTENT(in):: uu2
    TYPE(donnees), INTENT(in) :: DATA
    REAL(dp), INTENT(in) :: spectralradiusatinterface
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

#if(1==1)
    spectral_radius1=spectralradiusatinterface
    spectral_radius2=spectralradiusatinterface
#else
    DO l = 1,n_vars
       ubar1%u(l)= SUM(uu1%u(l))/REAL(SIZE(uu1),DP) !*Average in element 1
    END DO
    DO l = 1,n_vars
       ubar2%u(l)= SUM(uu2%u(l))/REAL(SIZE(uu2),DP) !*Average in element 2
    END DO

    spectral_radius1= ubar1%spectral_radius(xx,e1%n) !*Spectral radius in ubar1
    spectral_radius2= ubar2%spectral_radius(xx,e2%n) !*Spectral radius in ubar2
#endif


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


    SELECT CASE (DATA%ijump) !*Definition of h
    !*1,2,3 h defined by Remi
    !*11,12,13 h defined by Mario h=dx
    CASE(1,2,3) 
       !! h=1./((TYPE(e1%itype)*TYPE(e2%itype))*0.5) !* l2!
       h=1._DP/(0.5_DP*(SUM(ABS(grad1(1,:))) + SUM( ABS(grad2(1,:))) ))
    CASE(11,12,13)  
       h=0.5_DP*(e1%volume+e2%volume) !*dx if the mesh is uniform
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*PRINT*, e1%volume,e2%volume
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!
    CASE default
       PRINT*, "Problems in h in jump_burman"
       STOP
    END SELECT
    

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

                   !*-------------------------------
                   !*SAFETY CHECK
                   !*PRINT*, spectral_radius1* &
                   !*&(jumpflux *(-1._DP)* grad1(1,l) * theta + &
                   !*& jumpfluxder2*(-grad1_der(1,l))*theta2)
                   !*PRINT*, spectral_radius2* &
                   !*&(jumpflux * grad2(1,l) * theta + &
                   !*& jumpfluxder2*grad2_der(1,l)*theta2)
                   !*-------------------------------
                                 

    ENDDO

    !*NB: Most of the treated DoFs receive just one contribution, the one related to the element they belong to, because they are 0 outside that element and that contribution is not counted (the one coming from outside which is null). Just the DoF on the boundary receives the contribution of both the two elements (it actually has a nontrivial jump, it is not 0 on one of the two sides), anyway the two contributions are in ResJ(,1) and ResJ(,2), they will be put together outside this function, when these stabilization residuals will be added to Var%un (after the multiplication by dt)

    !*NB: The computation is very "close" to their definition as residuals
    !*If we consider their contribution to the nodal/cell residual we have
    !*contribution_i^K=int_f [grad u] grad phi_i|_K 
  END SUBROUTINE jump_burman



 !*---------------------------------------------------------------------------------------
 !* [[ J grad(phi_i) ]] A1 [[ J grad(u)-S ]]
 !*---------------------------------------------------------------------------------------
 SUBROUTINE jump_4( ed, e1, e2, u1k, u2k, resJ, alpha,alpha2, M,b1,b2,alphaconsistency,DATA)
    ! e1 is the element before the edge, e2 the one after
    ! e1-> u1k, e2->u2k

    REAL(dp),PARAMETER,DIMENSION(6)::TYPE=(/1._DP,2._DP,2._DP,3._DP,3._DP,4._DP/)
    REAL(dp), INTENT(in) :: alpha, alpha2
    TYPE(arete), INTENT(in):: ed
    TYPE(element),                     INTENT(in):: e1
    TYPE(element),                     INTENT(in):: e2
    TYPE(PVar), DIMENSION(e1%nsommets), INTENT(in):: u1k
    TYPE(PVar), DIMENSION(e2%nsommets), INTENT(in):: u2k
    TYPE(donnees), INTENT(in) :: DATA
    TYPE(Pvar), DIMENSION(e1%nsommets,2):: resJ
    REAL(DP), DIMENSION(2,2), INTENT(in):: M
    REAL(DP), DIMENSION(e1%nsommets), INTENT(in):: b1
    REAL(DP), DIMENSION(e2%nsommets), INTENT(in):: b2
    REAL(dp), INTENT(in) :: alphaconsistency
    TYPE(Pvar):: uinterface
    TYPE(PVar), DIMENSION(e1%nsommets) :: uval1
    TYPE(PVar), DIMENSION(e2%nsommets) :: uval2
    REAL(dp), DIMENSION(n_vars,n_vars,n_dim) :: J
    TYPE(Pvar) :: S1, S2
    REAL(DP), DIMENSION(n_dim):: gradb1, gradb2
    REAL(DP), DIMENSION(n_vars) :: jumpMJgraduS
    INTEGER:: indi
    real(dp) :: x !*For a safety check


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

    !*--------------------------
    !*For a SAFETY CHECK
    TYPE(PVar), DIMENSION(e1%nsommets) :: b1coeff
    TYPE(PVar), DIMENSION(e2%nsommets) :: b2coeff
    TYPE(PVar) :: gradsupp1, gradsupp2
    !*--------------------------


    !* - BASIS FUNCTIONS AT THE INTERFACE
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

    !* - GRAD(U) AT THE INTERFACE (AND ALSO THE SECOND DERIVATIVE EVEN IF I DO NOT KNOW WHAT TO DO WITH IT)

#if(1==1)
    !*eval_der_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the derivative of the solution in y 
    !*eval_der2_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the SECOND derivative of the solution in y 
    divflux1 = e1%eval_der(u1k(:),x1) !*First derivative in the right vertex of e1
    divflux2 = e2%eval_der(u2k(:),x2) !*First derivative in the left vertex of e2
#else
    divflux1=0._dp
    divflux2=0._dp
    DO l=1,n_vars
       divflux1%u(l)=sum(u1k%u(l)*grad1(1,:))
       divflux2%u(l)=sum(u2k%u(l)*grad2(1,:))
    END DO
#endif
    divdivflux1 = e1%eval_der2(u1k(:),x1) !*Second derivative in the right vertex of e1
    divdivflux2 = e2%eval_der2(u2k(:),x2) !*Second derivative in the left vertex of e2


    SELECT CASE (DATA%ijump) !*Definition of h
    !*4,5 h defined by Remi
    !*14,15 h defined by Mario h=dx
    CASE(4,5) 
       !! h=1./((TYPE(e1%itype)*TYPE(e2%itype))*0.5) !* l2!
       h=1._DP/(0.5_DP*(SUM(ABS(grad1(1,:))) + SUM( ABS(grad2(1,:))) ))
    CASE(14,15)  
       h=0.5_DP*(e1%volume+e2%volume) !*dx if the mesh is uniform
    CASE default
       PRINT*, "Problems in h in jump_4"
       STOP
    END SELECT



    !* - U AT THE INTERFACE SO TO GET J AND S
    uinterface=0._DP
    uval1=0._DP
    uval2=0._DP
    uval1(:)=Control_to_Cons(u1k,e1)
    uval2(:)=Control_to_Cons(u2k,e2)
    uinterface=uval2(1)
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFTEY CHECK
    !*PRINT*, "Right", uval2(1)%u
    !*PRINT*
    !*IF( (ABS(uinterface%u(1)-uval2(1)%u(1)).GE. 10._DP**(-15))  .OR. &
    !*   & (ABS(uinterface%u(2)-uval2(1)%u(2)).GE. 10._DP**(-15))  .OR. &
    !*   & (ABS(uinterface%u(1)-uval1(2)%u(1)).GE. 10._DP**(-15))  .OR. &
    !*   & (ABS(uinterface%u(2)-uval1(2)%u(2)).GE. 10._DP**(-15))  ) THEN
    !*   PRINT*, "PROBLEM"   
    !*   PRINT*, uinterface%u(1)-uval2(1)%u(1)
    !*   PRINT*, uinterface%u(2)-uval2(1)%u(2)
    !*   PRINT*, uinterface%u(1)-uval1(2)%u(1) 
    !*   PRINT*, uinterface%u(2)-uval1(2)%u(2)
    !*END IF
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !*Jacobian
    J=0._DP
    J=uinterface%Jacobian( (/e2%coor(1)/) )
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*IF( (ABS(J(1,1,1)-0._DP).GE. 10._DP**(-15)) .OR. &
    !*  & (ABS(J(1,2,1)-1._DP).GE. 10._DP**(-15)) .OR. &
    !*  & (ABS(J(2,1,1)- (-(uinterface%u(2)/uinterface%u(1))**2+grav*uinterface%u(1)) ).GE. 10._DP**(-15)) .OR. &
    !*  & (ABS(J(2,2,1)- 2._DP*uinterface%u(2)/uinterface%u(1) ).GE. 10._DP**(-15))) THEN
    !*   PRINT*, "Problem in Jacobian"
    !*   PRINT*, J(1,1,1)
    !*   PRINT*, J(1,2,1)
    !*   PRINT*, J(2,1,1), -(uinterface%u(2)/uinterface%u(1))**2+grav*uinterface%u(1) 
    !*   PRINT*, J(2,2,1), 2._DP*uinterface%u(2)/uinterface%u(1)
    !*   STOP
    !*END IF
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !*Source
    S1=0._DP
    S2=0._DP
    S1%u(1)=0._DP
    S2%u(1)=0._DP
    gradb1=0._DP
    gradb2=0._DP

    gradb1=SUM( b1(:)*grad1(1,:) ) !*For each component the derivative is the sum of the derivatives of the basis functions in that point multiplied by the respective coefficients 
    gradb2=SUM( b2(:)*grad2(1,:) ) !*For each component the derivative is the sum of the derivatives of the basis functions in that point multiplied by the respective coefficients 

    S1%u(2:2)=-grav*uinterface%u(1)*gradb1 &
    & - grav*uinterface%u(1)*n_manning**2*uinterface%u(2)*ABS(uinterface%u(2))/( uinterface%u(1)**(10._DP/3._DP) )
    S2%u(2:2)=-grav*uinterface%u(1)*gradb2 &
    & - grav*uinterface%u(1)*n_manning**2*uinterface%u(2)*ABS(uinterface%u(2))/( uinterface%u(1)**(10._DP/3._DP) )

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*x=e2%coor(1)
    !*IF((ABS(S1%u(2)-(-grav*uinterface%u(1)* ( &
    !*        & -50._dp*(x-0.4_dp)*2._dp*0.1_dp*exp(-50._dp*(x-0.4_dp)**2) &
    !*        &-300._dp*(x-0.55_dp)*2._dp*0.04_dp*exp(-300._dp*(x-0.55_dp)**2) &
    !*        & ) ) ) .GE. 10._DP**(-7)) &
    !*        & .OR. (ABS(S2%u(2)-(-grav*uinterface%u(1)* ( &
    !*        & -50._dp*(x-0.4_dp)*2._dp*0.1_dp*exp(-50._dp*(x-0.4_dp)**2) &
    !*        &-300._dp*(x-0.55_dp)*2._dp*0.04_dp*exp(-300._dp*(x-0.55_dp)**2) &
    !*        & ) )) .GE. 10._DP**(-7))    ) THEN
    !*           PRINT*, S1%u(2)
    !*           PRINT*, S1%u(2)-(-grav*uinterface%u(1)* ( &
    !*           & -50._dp*(x-0.4_dp)*2._dp*0.1_dp*exp(-50._dp*(x-0.4_dp)**2) &
    !*           &-300._dp*(x-0.55_dp)*2._dp*0.04_dp*exp(-300._dp*(x-0.55_dp)**2) &
    !*           & ) )
    !*           PRINT*, S2%u(2)
    !*           PRINT*, S2%u(2)-(-grav*uinterface%u(1)* ( &
    !*           & -50._dp*(x-0.4_dp)*2._dp*0.1_dp*exp(-50._dp*(x-0.4_dp)**2) &
    !*           &-300._dp*(x-0.55_dp)*2._dp*0.04_dp*exp(-300._dp*(x-0.55_dp)**2) &
    !*           & ) )
    !*END IF
    !*!*SAFETY CHECK, lake at rest
    !*PRINT*, grav*uinterface%u(1)*divflux1%u(1) +grav*uinterface%u(1)*gradb1
    !*PRINT*, grav*uinterface%u(1)*divflux2%u(1) +grav*uinterface%u(1)*gradb2
    !*PRINT*, divflux1%u(1)+gradb1
    !*PRINT*, divflux2%u(1)+gradb2
    !*!*SAFETY CHECK
    !*DO indi=1,e1%nsommets
    !*   b1coeff(indi)=b1(indi)
    !*END DO
    !*DO indi=1,e2%nsommets
    !*   b2coeff(indi)=b2(indi)       
    !*END DO
    !*gradsupp1=e1%eval_der(b1coeff(:),x1)
    !*gradsupp2=e2%eval_der(b2coeff(:),x2)
    !*IF ( (ABS(gradsupp1%u(1)-gradb1(1)) > 10._DP**(-13)) .OR. &
    !*& (ABS(gradsupp2%u(1)-gradb2(1)) > 10._DP**(-13)) ) THEN 
    !*   PRINT*, "1", gradsupp1%u(1)-gradb1
    !*   PRINT*, "2", gradsupp2%u(1)-gradb2
    !*   !*STOP
    !*END IF
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !*"Consistent" stabilization parameters
    theta=alpha*h**2 
    theta2=alpha2*h**4 

    !*Jump of the solution (el2-el1)
    !*NB: The following are vectors 



    !*Jump of J grad(u) - S !*NB: el2-el1 !*ALERT - SOURCE
    jumpMJgraduS=0._DP !*(n_vars)
    jumpMJgraduS=MATMUL(M(:,:), MATMUL(J(:,:,1),divflux2%u-divflux1%u) - (S2%u-S1%u) )


    !*NO SPECTRAL RADIUS
    DO indi=1,e1%nsommets
       ResJ(indi,1)%u = alphaconsistency*theta*(-1._DP)* MATMUL(grad1(1,indi)*(J(:,:,1)), jumpMJgraduS)
    END DO
    DO indi=1,e2%nsommets
       ResJ(indi,2)%u = alphaconsistency*theta*(1._DP)* MATMUL(grad2(1,indi)*(J(:,:,1)), jumpMJgraduS)
    END DO

        

    !*NB: Most of the treated DoFs receive just one contribution, the one related to the element they belong to, because they are 0 outside that element and that contribution is not counted (the one coming from outside which is null). Just the DoF on the boundary receives the contribution of both the two elements (it actually has a nontrivial jump, it is not 0 on one of the two sides), anyway the two contributions are in ResJ(,1) and ResJ(,2), they will be put together outside this function, when these stabilization residuals will be added to Var%un (after the multiplication by dt)

    !*NB: The computation is very "close" to their definition as residuals
    !*If we consider their contribution to the nodal/cell residual we have
    !*contribution_i^K=int_f [grad u] grad phi_i|_K 
  END SUBROUTINE jump_4


 !*---------------------------------------------------------------------------------------
 !*Jump [[ J grad(phi_i) ]] A1 [[ grad(g) ]]
 !*---------------------------------------------------------------------------------------
 SUBROUTINE jump_4_global_flux( ed, e1, e2, u1k, u2k, resJ, alpha,alpha2,M,g1,g2,DATA )
    ! e1 is the element before the edge, e2 the one after
    ! e1-> u1k, e2->u2k

    REAL(dp),PARAMETER,DIMENSION(6)::TYPE=(/1._DP,2._DP,2._DP,3._DP,3._DP,4._DP/)
    REAL(dp), INTENT(in) :: alpha, alpha2
    TYPE(arete), INTENT(in):: ed
    TYPE(element),                     INTENT(in):: e1
    TYPE(element),                     INTENT(in):: e2
    TYPE(PVar), DIMENSION(e1%nsommets), INTENT(in):: u1k
    TYPE(PVar), DIMENSION(e2%nsommets), INTENT(in):: u2k
    TYPE(donnees), INTENT(in) :: DATA
    TYPE(Pvar), DIMENSION(e1%nsommets,2):: resJ
    REAL(DP), DIMENSION(2,2), INTENT(in):: M
    TYPE(PVar), DIMENSION(e1%nsommets), INTENT(in):: g1
    TYPE(PVar), DIMENSION(e2%nsommets), INTENT(in):: g2

    TYPE(Pvar):: uinterface
    TYPE(PVar), DIMENSION(e1%nsommets) :: uval1
    TYPE(PVar), DIMENSION(e2%nsommets) :: uval2
    REAL(dp), DIMENSION(n_vars,n_vars,n_dim) :: J
    TYPE(Pvar) :: S1, S2
    REAL(DP), DIMENSION(n_dim):: gradb1, gradb2
    REAL(DP), DIMENSION(n_vars) :: jumpMgradg
    INTEGER:: indi
    real(dp) :: x !*For a safety check


    TYPE(Pvar):: u
    TYPE(Pvar):: divflux1, divflux2, jumpflux,jumpfluxder2, divdivflux1, divdivflux2
    INTEGER:: is, l
    real(dp),dimension(n_dim):: xx=0._DP
    REAL(dp):: theta, theta2, h,  umax, umin, l2
    REAL(dp), DIMENSION(n_vars,n_vars,n_dim)::jaco
    REAL(dp), DIMENSION(e1%nsommets):: base1
    REAL(dp),DIMENSION(n_dim,e1%nsommets):: grad1, grad1_der
    REAL(dp), DIMENSION(e2%nsommets):: base2
    REAL(dp),DIMENSION(n_dim,e2%nsommets):: grad2, grad2_der
    REAL(dp),DIMENSION(2), PARAMETER:: x2=(/0._DP,1._DP/),x1=(/1._DP,0._DP/)
    !*Barycentric coordinates of the left vertex (0,1) x2
    !*Barycentric coordinates of the right vertex (1,0) x1

    !*--------------------------
    !*For a SAFETY CHECK
    TYPE(PVar), DIMENSION(e1%nsommets) :: b1coeff
    TYPE(PVar), DIMENSION(e2%nsommets) :: b2coeff
    TYPE(PVar) :: gradsupp1, gradsupp2
    !*--------------------------



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

    ! jump of grad g
    !*DERIVATIVES OF THE SOLUTION
    !*eval_der_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the derivative of the solution in y 
    !*eval_der2_element(u,y) takes in input the coefficients of the solution (or a general function) in the DoFs and the barycentric coordinates y and evaluates the SECOND derivative of the solution in y 
    divflux1 = e1%eval_der(g1(:),x1) !*First derivative in the right vertex of e1
    divflux2 = e2%eval_der(g2(:),x2) !*First derivative in the left vertex of e2
    divdivflux1 = e1%eval_der2(g1(:),x1) !*Second derivative in the right vertex of e1
    divdivflux2 = e2%eval_der2(g2(:),x2) !*Second derivative in the left vertex of e2


    SELECT CASE (DATA%ijump) !*Definition of h
    !*6 h defined by Remi
    !*16 h defined by Mario h=dx
    CASE(6,7) 
       !! h=1./((TYPE(e1%itype)*TYPE(e2%itype))*0.5) !* l2!
       h=1._DP/(0.5_DP*(SUM(ABS(grad1(1,:))) + SUM( ABS(grad2(1,:))) ))
    CASE(16,17)  
       h=0.5_DP*(e1%volume+e2%volume) !*dx if the mesh is uniform
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!
       !*SAFETY CHECK
       !*PRINT*, e1%volume,e2%volume
       !*!!!!!!!!!!!!!!!!!!!!!!!!!!
    CASE default
       PRINT*, "Problems in h in jump_burman"
       STOP
    END SELECT
    

    IF (  DATA%iGlobalFlux == 0 ) THEN
        PRINT*, "Stop in jump_burman, in scheme_swwb. Jump 6 or 16 available only with global flux."
        STOP
    END IF

    !*Jump of the solution (el2-el1)
    !*NB: The following are vectors 
    jumpflux=(divflux2-divflux1)   !*Jump of the first derivative over the edge
    jumpfluxder2=(divdivflux2-divdivflux1) !*Jump of the second derivative over the edge, even if not needed now


    !* - U AT THE INTERFACE SO TO GET J
    uinterface=0._DP
    uval1=0._DP
    uval2=0._DP
    uval1(:)=Control_to_Cons(u1k,e1)
    uval2(:)=Control_to_Cons(u2k,e2)
    uinterface=uval2(1)
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFTEY CHECK
    !*PRINT*, "Right", uval2(1)%u
    !*PRINT*
    !*IF( (ABS(uinterface%u(1)-uval2(1)%u(1)).GE. 10._DP**(-15))  .OR. &
    !*   & (ABS(uinterface%u(2)-uval2(1)%u(2)).GE. 10._DP**(-15))  .OR. &
    !*   & (ABS(uinterface%u(1)-uval1(2)%u(1)).GE. 10._DP**(-15))  .OR. &
    !*   & (ABS(uinterface%u(2)-uval1(2)%u(2)).GE. 10._DP**(-15))  ) THEN
    !*   PRINT*, "PROBLEM"   
    !*   PRINT*, uinterface%u(1)-uval2(1)%u(1)
    !*   PRINT*, uinterface%u(2)-uval2(1)%u(2)
    !*   PRINT*, uinterface%u(1)-uval1(2)%u(1) 
    !*   PRINT*, uinterface%u(2)-uval1(2)%u(2)
    !*END IF
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !*Jacobian
    J=0._DP
    J=uinterface%Jacobian( (/e2%coor(1)/) )
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*SAFETY CHECK
    !*IF( (ABS(J(1,1,1)-0._DP).GE. 10._DP**(-15)) .OR. &
    !*  & (ABS(J(1,2,1)-1._DP).GE. 10._DP**(-15)) .OR. &
    !*  & (ABS(J(2,1,1)- (-(uinterface%u(2)/uinterface%u(1))**2+grav*uinterface%u(1)) ).GE. 10._DP**(-15)) .OR. &
    !*  & (ABS(J(2,2,1)- 2._DP*uinterface%u(2)/uinterface%u(1) ).GE. 10._DP**(-15))) THEN
    !*   PRINT*, "Problem in Jacobian"
    !*   PRINT*, J(1,1,1)
    !*   PRINT*, J(1,2,1)
    !*   PRINT*, J(2,1,1), -(uinterface%u(2)/uinterface%u(1))**2+grav*uinterface%u(1) 
    !*   PRINT*, J(2,2,1), 2._DP*uinterface%u(2)/uinterface%u(1)
    !*   STOP
    !*END IF
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !*"Consistent" stabilization parameters
    theta=alpha*h**2 
    theta2=alpha2*h**4 

    jumpMgradg=0._DP
    jumpMgradg=MATMUL(M(:,:), divflux2%u-divflux1%u )


    !*NO SPECTRAL RADIUS
    DO indi=1,e1%nsommets
       ResJ(indi,1)%u = theta*(-1._DP)* MATMUL(grad1(1,indi)*(J(:,:,1)), jumpMgradg)
    END DO
    DO indi=1,e2%nsommets
       ResJ(indi,2)%u = theta*(1._DP)* MATMUL(grad2(1,indi)*(J(:,:,1)),  jumpMgradg)
    END DO

        

    !*NB: Most of the treated DoFs receive just one contribution, the one related to the element they belong to, because they are 0 outside that element and that contribution is not counted (the one coming from outside which is null). Just the DoF on the boundary receives the contribution of both the two elements (it actually has a nontrivial jump, it is not 0 on one of the two sides), anyway the two contributions are in ResJ(,1) and ResJ(,2), they will be put together outside this function, when these stabilization residuals will be added to Var%un (after the multiplication by dt)

    !*NB: The computation is very "close" to their definition as residuals
    !*If we consider their contribution to the nodal/cell residual we have
    !*contribution_i^K=int_f [grad u] grad phi_i|_K 
  END SUBROUTINE jump_4_global_flux




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



  !*--------------------------------------------------------------------------------------
  !*galerkin_new_wb operates like galerkin_new for what concerns the velocity part of the flux 
  !*BUT the hydrostatic part and the source are dealt differently so to guarantee lake at rest 

  !*galerkin_new computes the galerkin residuals integrating directly phi div f (and not through the divergence theorem (integration by parts in 1d))
  !*\int_e varphi_i *du - dt* \int_e varphi_i div(flux) 
  !*--------------------------------------------------------------------------------------
  FUNCTION galerkin_new_wb(e, u, du,  flux_velocity, Heta_matrix, dt,jt,mesh,source_extra) RESULT(res)
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
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: du
    REAL(dp)  ,                        INTENT(in):: dt

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*For wb
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux_velocity !*Velocity part of the flux
    REAL(dp), DIMENSION(e%nsommets,e%nsommets):: Heta_matrix !*Matrix made by H_j (H_k+b_k) to integrate H*grad(H)+H*grad(b) 
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*extra contribution
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: source_extra
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    TYPE(PVar), DIMENSION(e%nsommets)            :: res

    INTEGER   ::  l

    INTEGER :: indi !*For the loop on the DoFs
    INTEGER :: indj !*For the loop on the DoFs
    INTEGER :: indk !*For the loop on the DoFs
    INTEGER :: indc !*For the loop on the components
    INTEGER :: indq !*For the quadrature points

    REAL(DP), DIMENSION(e%nsommets) :: base !*All the basis functions in the quadrature point
    REAL(DP), DIMENSION(n_dim,e%nsommets) :: gradbase !*All the gradients of the basis functions in the quadrature point

    TYPE(PVar), DIMENSION(e%nsommets) :: res_supp

    REAL(DP), DIMENSION(e%nsommets,e%nsommets,e%nsommets,n_dim) :: phiphigradphi
    TYPE(PVar), DIMENSION(e%nsommets) :: res_supp_bis

    TYPE(PVar), DIMENSION(e%nsommets) :: res_supp_tris



    !*RMK
    !*e%mass(i,j)-> local mass matrix int_K phi_i phi_j
    !*e%coeff(i,j) -> another matrix int_K phi_i  d_x phi_j !*NB:DERIVATIVE ON THE SECOND INDEX

    DO l=1, SIZE(flux_velocity(1,1)%u) !*Loop on the components
       res(:)%u(l)=MATMUL(e%mass(:,:),du%u(l)-dt*source_extra%u(l)) & !*The index of the columns runs
            &+dt*MATMUL(e%coeff   (:,:), flux_velocity(1,:)%u(l) )!& !*The index of the columns runs
    ENDDO

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*NB: With respect to galerkin_new we removed the source
    !*DO l=1, SIZE(flux(1,1)%u) !*Loop on the components
    !*   res(:)%u(l)=MATMUL(e%mass(:,:),du%u(l)-dt*source%u(l)) & !*The index of the columns runs
    !*        &+dt*MATMUL(e%coeff   (:,:), flux(1,:)%u(l) )!& !*The index of the columns runs
    !*ENDDO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   
    !*Now we need to add the hydrostatic and source part

    !*I CODED THREE EQUIVALENT ALTERANTIVES
    !*In the first one I do not precompute anything
    !*In the second one I precompute the matrix gradgradphi but here
    !*In the third one I use precomputing structures <- It's the most efficient so I will leave this uncommented and comment the other two
    

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!*1) FIRST ALTERNATIVE: nothing precomputed 
    !*!*I will fill res_supp
    !*!*Safe intialization
    !*base=0._DP
    !*gradbase=0._DP
    !*res_supp=0._DP
    !*DO indq=1,e%nquad !*Loop on the quadrature points
    !*   !*Safe intialization
    !*   base=0._DP
    !*   gradbase=0._DP
    !*   !*Compute the basis functions and the gradients in each DoF
    !*   DO indi=1,e%nsommets 
    !*      base(indi)=e%base(indi,e%quad(:,indq))
    !*      gradbase(:,indi)=e%gradient(indi, e%quad(:,indq))
    !*   END DO
    !*   !*Actual integration
    !*   DO indi=1,e%nsommets
    !*      DO indj=1,e%nsommets
    !*         DO indk=1,e%nsommets
    !*             res_supp(indi)%u(2:2)=res_supp(indi)%u(2:2)+&
    !*             & dt*grav*Heta_matrix(indj,indk)*base(indi)*base(indj)*gradbase(:,indk)*e%weight(indq)
    !*         END DO
    !*      END DO
    !*   END DO
    !*END DO
    !*!*Scaling
    !*res_supp=res_supp*e%volume !*To achieve the integral on the correct domain
    !*DO indi=1,e%nsommets
    !*   res(indi)%u=res(indi)%u+res_supp(indi)%u
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!*2) SECOND ALTERNATIVE: precomputation of the tensor phiphigrad here 
    !*!*I will fill res_supp_bis
    !*!*Safe initialization
    !*phiphigradphi=0._DP
    !*base=0._DP
    !*gradbase=0._DP
    !*DO indq=1,e%nquad
    !*   !*Safe intialization
    !*   base=0._DP
    !*   gradbase=0._DP
    !*   !*Compute the basis functions and the gradients in each DoF
    !*   DO indi=1,e%nsommets 
    !*      base(indi)=e%base(indi,e%quad(:,indq))
    !*      gradbase(:,indi)=e%gradient(indi, e%quad(:,indq))
    !*   END DO
    !*   !*Actual integration
    !*   DO indi=1,e%nsommets
    !*      DO indj=1,e%nsommets
    !*         DO indk=1,e%nsommets
    !*            phiphigradphi(indi,indj,indk,:)=phiphigradphi(indi,indj,indk,:)+&
    !*            & base(indi)*base(indj)*gradbase(:,indk)*e%weight(indq)   
    !*         END DO
    !*      END DO  
    !*   END DO
    !*END DO
    !*res_supp_bis=0._DP
    !*DO indi=1,e%nsommets
    !*   DO indj=1,e%nsommets
    !*      DO indk=1,e%nsommets
    !*         res_supp_bis(indi)%u(2:2)=res_supp_bis(indi)%u(2:2)+&
    !*         & dt*grav*Heta_matrix(indj,indk)*phiphigradphi(indi,indj,indk,:)
    !*      END DO
    !*   END DO  
    !*END DO
    !*res_supp_bis=res_supp_bis*e%volume !*To achieve the integral on the correct domain
    !*!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!*SAFETY SAFETY CHECK, checking the difference between res_supp and res_supp_bis
    !*DO indi=1,e%nsommets
    !*   !*PRINT*, indi, res_supp_bis(indi)%u-res_supp(indi)%u
    !*   !*PRINT*
    !*   IF( (ABS(res_supp_bis(indi)%u(1)-res_supp(indi)%u(1)).GE. 10._DP**(-14)) .OR.&
    !*   & (ABS(res_supp_bis(indi)%u(2)-res_supp(indi)%u(2)).GE. 10._DP**(-14)) ) THEN
    !*      PRINT*, "Azz"
    !*      !*STOP
    !*   END IF
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*DO indi=1,e%nsommets
    !*   res(indi)=res(indi)+res_supp_bis(indi)
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!*3) THIRD ALTERNATIVE: I use the tensor precomputed tensor e%phiphigrad 
    !*!*I will fill res_supp_tris
    res_supp_tris=0._DP
    DO indi=1,e%nsommets
       DO indj=1,e%nsommets
          DO indk=1,e%nsommets
             res_supp_tris(indi)%u(2:2)=res_supp_tris(indi)%u(2:2)+&
             & dt*grav*Heta_matrix(indj,indk)*e%phiphigradphi(indi,indj,indk,:)
          END DO
       END DO  
    END DO
    DO indi=1,e%nsommets
       res(indi)=res(indi)+res_supp_tris(indi)
    END DO

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!*SAFETY CHECK, checking the difference between res_supp, res_supp_bis and res_supp_tris
    !*DO indi=1,e%nsommets
    !*   !*PRINT*, indi, res_supp_bis(indi)%u-res_supp(indi)%u
    !*   !*PRINT*
    !*   IF( (ABS(res_supp_bis(indi)%u(1)-res_supp(indi)%u(1)).GE. 10._DP**(-17)) .OR.&
    !*   & (ABS(res_supp_bis(indi)%u(2)-res_supp(indi)%u(2)).GE. 10._DP**(-17)) .OR.&
    !*   & (ABS(res_supp_tris(indi)%u(1)-res_supp(indi)%u(1)).GE. 10._DP**(-17)) .OR.&
    !*   & (ABS(res_supp_tris(indi)%u(2)-res_supp(indi)%u(2)).GE. 10._DP**(-17))) THEN
    !*      PRINT*, "Azz"
    !*      STOP
    !*   END IF
    !*END DO
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
 
  END FUNCTION galerkin_new_wb


  !*--------------------------------------------------------------------------------------
  !*Lax Friedrichs = Galerkin + average stabilization
  !*--------------------------------------------------------------------------------------
  FUNCTION lxf_wb(e, u, du,  flux_velocity, Heta_matrix, dt,jt,mesh,source_extra) RESULT(phi_lxf)
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
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: du
    REAL(dp)   ,                        INTENT(in):: dt

    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*For wb
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux_velocity !*Velocity part of the flux
    REAL(dp), DIMENSION(e%nsommets,e%nsommets):: Heta_matrix !*Matrix made by H_j (H_k+b_k) to integrate H*grad(H)+H*grad(b) 
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*extra contribution
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: source_extra
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    TYPE(PVar), DIMENSION(e%nsommets)            :: res, phi_lxf,u_cons
    TYPE(PVar):: u_consL, u_consR
    REAL(dp),DIMENSION(n_dim):: xx=0._dp
    TYPE(PVar), DIMENSION(e%nsommets):: v,u_diff,u_diff2
    TYPE(PVar) ::  ubar
    TYPE(PVar), DIMENSION(e%nsommets)::vbar
    REAL(dp)::   alpha
    INTEGER:: l,i,iq, k

    res= galerkin_new_wb(e, u, du,  flux_velocity, Heta_matrix, dt,jt,mesh,source_extra)

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

  END FUNCTION lxf_wb



  !*--------------------------------------------------------------------------------------
  !*Psi Galerkin = Lax Fridriechs + additional stabilization (limit)
  !*--------------------------------------------------------------------------------------
  FUNCTION psi_galerkin_wb(e, u, du,  flux_velocity, Heta_matrix, dt,jt,mesh,source_extra) RESULT(res)
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
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: du
    TYPE(PVar), DIMENSION(n_dim,e%nsommets), INTENT(in):: flux_velocity
    REAL(dp), DIMENSION(e%nsommets,e%nsommets):: Heta_matrix !*Matrix made by H_j (H_k+b_k) to integrate H*grad(H)+H*grad(b) 
    REAL(dp)   ,                        INTENT(in):: dt
    TYPE(PVar), DIMENSION(e%nsommets)            :: res, phi_lxF,w
    INTEGER::  k, l
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*extra contribution
    TYPE(PVar), DIMENSION(e%nsommets), INTENT(in):: source_extra
    !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    phi_lxf= lxf_wb(e, u, du,  flux_velocity, Heta_matrix, dt,jt,mesh,source_extra) 
    CALL limit(e,u,res,phi_lxf)

  END FUNCTION psi_galerkin_wb



END MODULE scheme
