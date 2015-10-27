!  Subroutines for basic 3D  rate-independent
!  elasto-plastic deformation with nonlinear kinematic
!  hardening. 'Armstrong-Frederick' strain
!  hardening is employed.


!==========================SUBROUTINE el_kinematic_viscoplasticity_3Dbasic ==============================
subroutine el_kinematic_viscoplasticity_3Dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint,m_count, n_count,i,j,k,l,a,jj,iter,maxit

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bbar(6,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp(6,length_dof_array)          ! temp matrix when constructing B-bar matrix
    real (prec)  ::  el_vol                            ! element volume
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  ddisp(3,length_coord_array/3)     ! Re-shaped displacement array disp(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12,Y0,c_prop,gamma_prop             ! Material properties
    real (prec)  ::  Iden(3,3),S_tr(3,3),alpha_t(3,3),S_t(3,3),press_t,sig_t(3,3),Sbar_tr
    real (prec)  ::  sdev0(6),sdev(6),press0,press,alpha0(6),alpha(6)
    real (prec)  ::  S_tau(3,3),press_tau,alpha_tau(3,3),sig_tau(3,3)
    real (prec)  ::  deps(3,3),deps0(3,3),tr_deps
    real (prec)  ::  F_tr,F_fun,dFde,error,delta_mage,tol,delta_mage_new,ita,lamda,ditade,dlamdade,beta
    real (prec)  ::  temp1,temp2,temp3,temp4,temp5,temp6,temp7,Ctang(3,3,3,3)
    real (prec)  ::  third,examine,Sbar_tau
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    !     element_properties(3)         Initial yield stress
    !     element_properties(4)         Hardening rate c
    !     element_properties(5)         Nonlinear parameter gamma
    ! Utilities
      third = 1.d0/3.d0
    ! Pass material properties in

       E          = element_properties(1)
       xnu        = element_properties(2)
       Y0         = element_properties(3)
       c_prop     = element_properties(4)
       gamma_prop = element_properties(5)

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))
    ddisp = reshape(dof_increment,(/3,length_dof_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)


    element_residual = 0.d0
    element_stiffness = 0.d0

    el_vol = 0.d0
    dNbardx = 0.d0

! B-bar method
    !     -- Loop over the integration points for B bar method
    do kint = 1, n_points
       call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

    !
        do m_count = 1,n_nodes
           do n_count = 1,3
              dNbardx(m_count,n_count) = dNbardx(m_count,n_count) + &
                                         dNdx(m_count,n_count)*w(kint)*determinant
           end do
        end do
    !   Get the average element volume
        el_vol = el_vol + w(kint)*determinant
    end do

    !   Get the final form of vol_avg_shape function derivatives
       do m_count = 1,n_nodes
         do n_count = 1,3
             dNbardx(m_count,n_count) = dNbardx(m_count,n_count)/el_vol
         end do
       end do
 ! ******************************************************************************************
         !     --  Loop over integration points
    Iden = 0.d0
    do i=1,3
       Iden(i,i)=1.0d0
    end do

    do kint = 1, n_points
        jj = 0 ! jj is used for tracking the state variables
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        B = 0.d0
        temp = 0.d0
        Bbar = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        do m_count = 1,n_nodes
        temp(1:3:1,3*m_count-2) = dNbardx(m_count,1) - dNdx(m_count,1)
        temp(1:3:1,3*m_count-1) = dNbardx(m_count,2) - dNdx(m_count,2)
        temp(1:3:1,3*m_count) = dNbardx(m_count,3) - dNdx(m_count,3)
        end do

        Bbar = B + 1.d0/3.d0 *temp
     !      obtian the state variables at step N
        sdev0  = initial_state_variables(1+jj:6+jj)
        press0 = initial_state_variables(7+jj)
        alpha0 = initial_state_variables(8+jj:13+jj)
     !      Construct them into matrix form for further calculation
           press_t = press0

        do i=1,3
           S_t(i,i)= sdev0(i)
           alpha_t(i,i) = alpha0(i)
        end do
           S_t(1,2) = sdev0(4)
           S_t(2,1) = S_t(1,2)
           S_t(1,3) = sdev0(5)
           S_t(3,1) = S_t(1,3)
           S_t(2,3) = sdev0(6)
           S_t(3,2) = S_t(2,3)

          alpha_t(1,2) = alpha0(4)
          alpha_t(2,1) = alpha_t(1,2)
          alpha_t(1,3) = alpha0(5)
          alpha_t(3,1) = alpha_t(1,3)
          alpha_t(2,3) = alpha0(6)
          alpha_t(3,2) = alpha_t(2,3)


  !**********************************************************************************
    !  Meterial properties
  !**********************************************************************************

    deps = 0.d0
    do i=1,3
       do j=1,3
          do a = 1, n_nodes
             deps(i,j) = deps(i,j) + 0.5d0 *(dNdx(a,j)*ddisp(i,a)+dNdx(a,i)*ddisp(j,a))
          end do
       end do
    end do

   ! Step 1. Compute the deviatoric strain increment
      tr_deps= deps(1,1)+deps(2,2)+deps(3,3)
      deps0 = deps - third*tr_deps*Iden
   ! Step 2. Calculate the 'elastic predictor' S_tr and effectice stress Sbar_tr at the end of the increment
      S_tr = S_t + E/(1+xnu)*deps0
      Sbar_tr = 0.d0
      do i=1,3
        do j=1,3
          Sbar_tr = Sbar_tr + 1.5d0*(S_tr(i,j)-alpha_t(i,j))*(S_tr(i,j)-alpha_t(i,j))
        end do
      end do
      Sbar_tr = dsqrt(Sbar_tr)

   ! Step 3. Check whether this is a elastic or plastic step by yield criterion
      F_tr = Sbar_tr - Y0
    if ( F_tr .LE. 0.d0) then! elastic step
     ! write(*,*) 'elastic in'
      ! return elastic tangents and trial stress as its true stress at the end of the step.
      ! Update the internal variables
       alpha_tau = alpha_t
       S_tau = S_tr
       press_tau = press_t - tr_deps* E/(1.d0-2.d0*xnu)*third

      ! Update the material mechanical tangent modulus
       D = 0.d0
       d44 = 0.5D0*E/(1+xnu)
       d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
       d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
       D(1:3,1:3) = d12
       D(1,1) = d11
       D(2,2) = d11
       D(3,3) = d11
       D(4,4) = d44
       D(5,5) = d44
       D(6,6) = d44
     ! Update stress
       sig_tau = S_tau - press_tau*Iden

    else ! plastic step
    ! write(*,*) 'plastic in '
     ! Step.4 Solve the plastic strain increment
      delta_mage = 1.0d-15
      error = Y0
      tol = 1.d-8 * Y0
      iter = 1
      maxit = 50
     do while(error .GT. tol)

          temp1 = (1.5d0*E/(1+xnu)+c_prop/(1.d0 + gamma_prop*delta_mage))

          temp2 = 0.d0
          do i=1,3
             do j=1,3
              temp2 = temp2+1.5d0* (S_tr(i,j)-alpha_t(i,j)/(1.d0 + gamma_prop*delta_mage)) &
                            *(S_tr(i,j)-alpha_t(i,j)/(1.d0 + gamma_prop*delta_mage))
              end do
          end do
          temp2 = dsqrt(temp2)

          temp3 = 0.d0
          do i=1,3
             do j=1,3
              temp3 = temp3+1.5d0*(S_tr(i,j)-alpha_t(i,j)/(1.d0 + gamma_prop*delta_mage)) &
                            *(alpha_t(i,j)*gamma_prop/(1.d0 + gamma_prop*delta_mage)**2.d0)
             end do
          end do

          F_fun = Y0 +temp1*delta_mage-temp2
          dFde = temp1 - c_prop*delta_mage*gamma_prop/(1.d0 + gamma_prop*delta_mage)**2.d0 &
                 - temp3/temp2
          delta_mage_new = delta_mage - F_fun/dFde
          if (delta_mage_new .LT. 0.d0) then
            delta_mage = delta_mage/10.d0
          else
            delta_mage = delta_mage_new
          end if
          error = abs(F_fun)
          iter = iter + 1

          if (iter .EQ. maxit) then
            write(*,*) 'Warning in solving delta_mage: Max iterations exceeded.'
            stop
          end if

      end do

     ! Step 5. Update the deviatoric stress,current stress and pressure
       temp5 = Y0 + (gamma_prop*Y0+c_prop)*delta_mage
       ita = 1.5d0*E/(1.d0+xnu)*delta_mage/temp5
       temp6 = 1.5d0*E/Y0/(1.d0*xnu)-ita*c_prop/Y0
       lamda = 1.d0/(1.d0+delta_mage*temp6)
       S_tau = lamda*(S_tr+ita*alpha_t)
       sig_tau = S_tau - press_t*Iden + E/3.d0/(1.d0-2.d0*xnu)*tr_deps*Iden
       press_tau = -third*(sig_tau(1,1)+sig_tau(2,2)+sig_tau(3,3))

     ! Step 6. Update the nonlinear kinematic hardening relation: alpha_tau
       temp4 = (1.d0+(gamma_prop+c_prop/Y0)*delta_mage)
       alpha_tau=(alpha_t+delta_mage*c_prop/Y0*S_tau)/temp4
     ! Return the material mechanical tangents Cijkl or Matrix D

      ditade = 1.5d0*E/(1.d0+xnu)*(1.d0/temp5 - (temp5 - Y0)/temp5**2.d0)
      dlamdade = - (lamda**2.d0) *(temp6-delta_mage*c_prop/Y0*ditade)
      temp7 = dlamdade*ita+lamda*ditade
      beta = 1.d0/(Y0+(1.5d0*E/(1.d0+xnu)+c_prop/(1.d0+gamma_prop*delta_mage))*delta_mage)

      do i=1,3
        do j=1,3
          do k =1,3
            do l= 1,3
               Ctang(i,j,k,l) = lamda*E/(1.d0+xnu)*((Iden(i,k)*Iden(j,l)+Iden(j,k)*Iden(i,l))&
               /2.d0- third*Iden(i,j)*Iden(k,l))+third*E/(1.d0-2.d0*xnu)*Iden(i,j)*Iden(k,l) &
               + 1.5d0*beta*E/(1.d0+xnu)/dFde*(dlamdade*S_tr(k,l)*S_tr(i,j)-temp7&
               *alpha_t(k,l)*alpha_t(i,j)/(1.d0+gamma_prop*delta_mage)+ &
               temp7*S_tr(k,l)*alpha_t(i,j)-dlamdade*alpha_t(k,l)* &
               S_tr(i,j)/(1.d0+gamma_prop*delta_mage))
            end do
          end do
        end do
      end do

    ! Construct the D matrix using Cijkl above,not assuming C is symmetric here.
      D = 0.d0
      do i=1,3
        do j=1,3
           D(i,j) = Ctang(i,i,j,j)
        end do
      end do
        do i =1,3
           D(i,4) = Ctang(i,i,1,2)
           D(4,i) = Ctang(1,2,i,i)
        end do
           D(4,4) = Ctang(1,2,1,2)
        do i=1,3
           D(i,5) = Ctang(i,i,1,3)
           D(5,i) = Ctang(1,3,i,i)
        end do
           D(5,5) = Ctang(1,3,1,3)
           D(4,5) = Ctang(1,2,1,3)
           D(5,4) = Ctang(1,3,1,2)
        do i=1,3
           D(i,6) = Ctang(i,i,2,3)
           D(6,i) = Ctang(2,3,i,i)
        end do
           D(6,6) = Ctang(2,3,2,3)
           D(4,6) = Ctang(1,2,2,3)
           D(6,4) = Ctang(2,3,1,2)
           D(5,6) = Ctang(1,3,2,3)
           D(6,5) = Ctang(2,3,1,3)
      write(*,*) 'plastic step'
           write(*,*) ''
      !!!! examine yeild
      Sbar_tau = 0.d0
      do i=1,3
        do j=1,3
          Sbar_tau = Sbar_tau + 1.5d0*(S_tau(i,j)-alpha_tau(i,j))*(S_tau(i,j)-alpha_tau(i,j))
        end do
      end do
      Sbar_tau = dsqrt(Sbar_tau)

      examine = Sbar_tau - Y0
           write(*,*) 'examine'
           write(*,*)   examine
           write(*,*)   ''
    end if
!    Update state variables
     sdev = 0.d0
     alpha = 0.d0
     press = 0.d0
!
        do i=1,3
          sdev(i) = S_tau(i,i)
          alpha(i) = alpha_tau(i,i)
        end do

        sdev(4) = S_tau(1,2)
        sdev(5) = S_tau(1,3)
        sdev(6) = S_tau(2,3)

        alpha(4) = alpha_tau(1,2)
        alpha(5) = alpha_tau(1,3)
        alpha(6) = alpha_tau(2,3)

        press = press_tau

        updated_state_variables(1+jj:6+jj) = sdev
        updated_state_variables(7+jj) = press_tau
        updated_state_variables(8+jj:13+jj) = alpha

        jj = jj + 13

         stress(1) = sig_tau(1,1)
         stress(2) = sig_tau(2,2)
         stress(3) = sig_tau(3,3)
         stress(4) = sig_tau(1,2)
         stress(5) = sig_tau(1,3)
         stress(6) = sig_tau(2,3)
         write(*,*) 's11'
          write(*,*) stress(1)
          write(*,*) ' '

!
       if ( element_identifier == 3001 ) then

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

       else if (element_identifier == 4001) then

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(Bbar),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(Bbar(1:6,1:3*n_nodes)),matmul(D,Bbar(1:6,1:3*n_nodes)))*w(kint)*determinant

       end if

    end do

    return
end subroutine el_kinematic_viscoplasticity_3Dbasic

!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_3d_kinematic_plasticity(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step

    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables

    ! Local Variables
    logical      :: strcmp

    integer      :: n_points,kint,k,m_count, n_count,jj

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bbar(6,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp(6,length_dof_array)          ! temp matrix when constructing B-bar matrix
    real (prec)  ::  el_vol                            ! element volume
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0

    el_vol = 0.d0
    dNbardx = 0.d0
     !     -- Loop over the integration points for B bar method
    do kint = 1, n_points
       call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
    !
        do m_count = 1,n_nodes
           do n_count = 1,3
              dNbardx(m_count,n_count) = dNbardx(m_count,n_count) + &
                                         dNdx(m_count,n_count)*w(kint)*determinant
           end do
        end do
    !   Get the average element volume
        el_vol = el_vol + w(kint)*determinant
    end do

    !   Get the final form of vol_avg_shape function derivatives
       do m_count = 1,n_nodes
         do n_count = 1,3
             dNbardx(m_count,n_count) = dNbardx(m_count,n_count)/el_vol
         end do
       end do
    !     --  Loop over integration points
    do kint = 1, n_points
       jj = 0 ! jj is used for tracking the state variables
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        temp = 0.d0
        Bbar = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
        do m_count=1,n_nodes
        temp(1:3:1,3*m_count-2) = dNbardx(m_count,1) - dNdx(m_count,1)
        temp(1:3:1,3*m_count-1) = dNbardx(m_count,2) - dNdx(m_count,2)
        temp(1:3:1,3*m_count) = dNbardx(m_count,3) - dNdx(m_count,3)
        end do
        Bbar = B + 1.d0/3.d0 *temp

       sdev = updated_state_variables(1+jj:6+jj)
       p= updated_state_variables(7+jj)

        if ( element_identifier == 3001 ) then
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        strain = strain + dstrain
        stress(1:3) = sdev(1:3)-p
        stress(4:6) = sdev(4:6)
        else if ( element_identifier == 4001 ) then
        strain = matmul(Bbar,dof_total)
        dstrain = matmul(Bbar,dof_increment)
        strain = strain + dstrain
        stress(1:3) = sdev(1:3)-p
        stress(4:6) = sdev(4:6)
        end if

        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        jj = jj + 13
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(4)/0.5d0*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(5)/0.5d0*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(6)/0.5d0*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            end if
        end do

    end do

    return
end subroutine fieldvars_3d_kinematic_plasticity

!subroutine nl_kinematic_material(properties,n_props,dstrain,initial_state,updated_state,n_states,stress,D)
!    use Types
!    use ParamIO
!    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
!
!    implicit none
!
!    integer, intent( in ) :: n_props
!    integer, intent( in ) :: n_states
!
!    real (prec), intent( in )         :: properties(n_props)                                                   ! Element number
!
!
!    real (prec), intent( out ) :: stress(6)
!
!
!end subroutine nl_kinematic_material
