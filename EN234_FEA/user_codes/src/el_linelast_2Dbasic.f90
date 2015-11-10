!     Subroutines for basic 2D linear elastic elements



!==========================SUBROUTINE el_linelast_2Dplanestrain ==============================
subroutine el_linelast_2d_planestrain(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
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
    integer      :: n_points,kint,m_count,n_count

    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(3,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bbar(3,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp(3,length_dof_array)          ! temp matrix when constructing B-bar matrix
    real (prec)  ::  el_vol                            ! element volume
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/2,length_coord_array/2/))

    if (n_nodes == 3) n_points = 3
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 8) n_points = 9
    if (n_nodes == 4) n_points = 4

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1.D0+xnu)
    d11 = (1.D0-xnu)*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
    d12 = xnu*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
    D(1,2) = d12
    D(2,1) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d44
    el_vol = 0.d0
    dNbardx = 0.d0
   !     -- Loop over the integration points for B bar method
    do kint = 1, n_points
       call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
    !
        do m_count = 1,n_nodes
           do n_count = 1,2
              dNbardx(m_count,n_count) = dNbardx(m_count,n_count) + &
                                         dNdx(m_count,n_count)*w(kint)*determinant
           end do
        end do
    !   Get the average element volume
        el_vol = el_vol + w(kint)*determinant
    end do
    !   Get the final form of vol_avg_shape function derivatives
       do m_count = 1,n_nodes
         do n_count = 1,2
             dNbardx(m_count,n_count) = dNbardx(m_count,n_count)/el_vol
         end do
       end do

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        B = 0.d0
        temp = 0.d0
        Bbar = 0.d0
        B(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(3,2:2*n_nodes:2) = dNdx(1:n_nodes,1)

        do m_count = 1,n_nodes
        temp(1:2:1,2*m_count-1) = dNbardx(m_count,1) - dNdx(m_count,1)
        temp(1:2:1,2*m_count) = dNbardx(m_count,2) - dNdx(m_count,2)
        end do

        Bbar = B + 1.d0/2.d0 *temp

        if ( element_identifier == 101 ) then
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        stress = matmul(D,strain+dstrain)

        element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
            + matmul(transpose(B(1:3,1:2*n_nodes)),matmul(D,B(1:3,1:2*n_nodes)))*w(kint)*determinant

        else if ( element_identifier == 201 ) then
        strain = matmul(Bbar,dof_total)
        dstrain = matmul(Bbar,dof_increment)
        stress = matmul(D,strain+dstrain)

        element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(Bbar),stress)*w(kint)*determinant

        element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
            + matmul(transpose(Bbar(1:3,1:2*n_nodes)),matmul(D,Bbar(1:3,1:2*n_nodes)))*w(kint)*determinant
        end if
    end do
  
    return
end subroutine el_linelast_2d_planestrain

!==========================SUBROUTINE el_linelast_2Dplanestress ==============================
subroutine el_linelast_2d_planestress(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
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
    integer      :: n_points,kint,m_count,n_count

    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(3,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bbar(3,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp(3,length_dof_array)          ! temp matrix when constructing B-bar matrix
    real (prec)  ::  el_vol                            ! element volume
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.

    x = reshape(element_coords,(/2,length_coord_array/2/))

    if (n_nodes == 3) n_points = 3
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 8) n_points = 9
    if (n_nodes == 4) n_points = 4

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0

    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1.D0+xnu)
    d11 = 1.D0*E/( (1.D0-xnu)*(1.D0+xnu) )
    d12 = xnu*E/( (1.D0-xnu)*(1.D0+xnu) )
    D(1,2) = d12
    D(2,1) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d44

    el_vol = 0.d0
    dNbardx = 0.d0
   !     -- Loop over the integration points for B bar method
    do kint = 1, n_points
       call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
    !
        do m_count = 1,n_nodes
           do n_count = 1,2
              dNbardx(m_count,n_count) = dNbardx(m_count,n_count) + &
                                         dNdx(m_count,n_count)*w(kint)*determinant
           end do
        end do
    !   Get the average element volume
        el_vol = el_vol + w(kint)*determinant
    end do
    !   Get the final form of vol_avg_shape function derivatives
       do m_count = 1,n_nodes
         do n_count = 1,2
             dNbardx(m_count,n_count) = dNbardx(m_count,n_count)/el_vol
         end do
       end do

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        B = 0.d0
        temp = 0.d0
        Bbar = 0.d0
        B(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(3,2:2*n_nodes:2) = dNdx(1:n_nodes,1)

        do m_count = 1,n_nodes
        temp(1:2:1,2*m_count-1) = dNbardx(m_count,1) - dNdx(m_count,1)
        temp(1:2:1,2*m_count) = dNbardx(m_count,2) - dNdx(m_count,2)
        end do

        Bbar = B + 1.d0/2.d0 *temp
        if ( element_identifier == 102 ) then
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        stress = matmul(D,strain+dstrain)

        element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
            + matmul(transpose(B(1:3,1:2*n_nodes)),matmul(D,B(1:3,1:2*n_nodes)))*w(kint)*determinant

        else if ( element_identifier == 202 ) then
        strain = matmul(Bbar,dof_total)
        dstrain = matmul(Bbar,dof_increment)
        stress = matmul(D,strain+dstrain)

        element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(Bbar),stress)*w(kint)*determinant

        element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
            + matmul(transpose(Bbar(1:3,1:2*n_nodes)),matmul(D,Bbar(1:3,1:2*n_nodes)))*w(kint)*determinant
        end if
    end do

    return
end subroutine el_linelast_2d_planestress

!==========================SUBROUTINE el_linelast_2Dplanestrain ==============================
subroutine el_cahn_hilliard_2dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    use Globals, only: TIME,DTIME  ! For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
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
    integer      :: n_points,kint,m_count,n_count

    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bbar(3,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp_1(6)                         ! temp matrix when constructing q matrix
    real (prec)  ::  temp_2(6)                         ! temp matrix when constructing q matrix
    real (prec)  ::  temp_3(6)                         ! temp matrix when constructing q matrix
    real (prec)  ::  q_mat(6)                          !  q matrix
    real (prec)  ::  el_vol                            ! element volume
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  Diff,kappa,theta             ! Material properties
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Diffusive coefficient
    !     element_properties(2)         Kappa
    !     element_properties(3)         theta

    fail = .false.

    x = reshape(element_coords,(/2,length_coord_array/2/))

    if (n_nodes == 3) n_points = 3
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 8) n_points = 9
    if (n_nodes == 4) n_points = 4

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0


    Diff  = element_properties(1)
    kappa = element_properties(2)
    theta = element_properties(3)
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        B = 0.d0
        temp_1 = 0.d0
        temp_2 = 0.d0
        B(1,1:2*n_nodes-1:2) = N(1:n_nodes)
        B(2,2:2*n_nodes:2) = N(1:n_nodes)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(4,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(5,2:2*n_nodes:2) = dNdx(1:n_nodes,1)
        B(6,2:2*n_nodes:2) = dNdx(1:n_nodes,2)

        temp_1 = matmul(B,(dof_total+dof_increment))
        temp_2 = matmul(B,dof_increment)
        temp_3 = matmul(B,(dof_total+theta*dof_increment))
       ! construct q vector
        q_mat =0.d0
        q_mat(1)= temp_1(1) - temp_1(2)*(temp_1(2)**2.d0 -1.d0)
        q_mat(2)= temp_2(2)/ DTIME
        q_mat(3)= -kappa*temp_1(5)
        q_mat(4)= -kappa*temp_1(6)
        q_mat(5)= Diff*temp_3(3)
        q_mat(6)= Diff*temp_3(4)

       ! construct D matrix
        D = 0.d0
        D(1,1)= 1.d0
        D(1,2)= -(3*temp_1(2)**2.d0 -1.d0)
        D(2,2)= 1.d0 / DTIME
        D(3,5)= -kappa
        D(4,6)= D(3,5)
        D(5,3)= theta*Diff
        D(6,4)= D(5,3)


        element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(B),q_mat)*w(kint)*determinant

        element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
            + matmul(transpose(B(1:6,1:2*n_nodes)),matmul(D,B(1:6,1:2*n_nodes)))*w(kint)*determinant

    end do

    return
end subroutine el_cahn_hilliard_2dbasic

!==========================SUBROUTINE fieldvars_linelast_2d_planestrain ==============================
subroutine fieldvars_linelast_2d_planestrain(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
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

    integer      :: n_points,kint,k,m_count,n_count

    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  stress_total(6)                   ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  strain_total(6)                   ! Stress vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(3,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bbar(3,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp(3,length_dof_array)          ! temp matrix when constructing B-bar matrix
    real (prec)  ::  el_vol                            ! element volume
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/2,length_coord_array/2/))

    if (n_nodes == 3) n_points = 3
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 8) n_points = 9
    if (n_nodes == 4) n_points = 4

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1.D0+xnu)
    d11 = (1.D0-xnu)*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
    d12 = xnu*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
    D(1,2) = d12
    D(2,1) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d44

    el_vol = 0.d0
    dNbardx = 0.d0
     !     -- Loop over the integration points for B bar method
    do kint = 1, n_points
       call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
    !
        do m_count = 1,n_nodes
           do n_count = 1,2
              dNbardx(m_count,n_count) = dNbardx(m_count,n_count) + &
                                         dNdx(m_count,n_count)*w(kint)*determinant
           end do
        end do
    !   Get the average element volume
        el_vol = el_vol + w(kint)*determinant
    end do
    !   Get the final form of vol_avg_shape function derivatives
       do m_count = 1,n_nodes
         do n_count = 1,2
             dNbardx(m_count,n_count) = dNbardx(m_count,n_count)/el_vol
         end do
       end do

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        B = 0.d0
        temp = 0.d0
        Bbar = 0.d0
        B(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(3,2:2*n_nodes:2) = dNdx(1:n_nodes,1)
        do m_count = 1,n_nodes
        temp(1:2:1,2*m_count-1) = dNbardx(m_count,1) - dNdx(m_count,1)
        temp(1:2:1,2*m_count) = dNbardx(m_count,2) - dNdx(m_count,2)
        end do

        Bbar = B + 1.d0/2.d0 *temp
        if ( element_identifier == 101 ) then

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        strain = strain + dstrain
        stress = matmul(D,strain)
        else if ( element_identifier == 201 ) then
        strain = matmul(Bbar,dof_total)
        dstrain = matmul(Bbar,dof_increment)
        strain = strain + dstrain
        stress = matmul(D,strain)
        end if
        stress_total(1:2) = stress(1:2)
        stress_total(4) = stress(3)
        stress_total(3) = E*xnu*(strain(1)+strain(2))/(1.D0-2.D0*xnu)/(1.D0+xnu)
        stress_total(5:6) = 0.d0
        strain_total(1:2) = strain(1:2)
        strain_total(4) = (1+xnu)/E*stress_total(4)
        strain_total(3) = 0.d0
        strain_total(5:6) = 0.d0
        p = sum(stress_total(1:3))/3.d0
        sdev = stress_total
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(1)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(2)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(3)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(4)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(5)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(6)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(1)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(2)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(3)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(4)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(5)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(6)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do

    end do
  
    return
end subroutine fieldvars_linelast_2d_planestrain

!==========================SUBROUTINE fieldvars_linelast_2d_planestress ==============================
subroutine fieldvars_linelast_2d_planestress(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only: dNbardx => vol_avg_shape_function_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
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
  
    integer      :: n_points,kint,k,m_count,n_count

    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  stress_total(6)                   ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  strain_total(6)                   ! Stress vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(3,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  Bbar(3,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp(3,length_dof_array)          ! temp matrix when constructing B-bar matrix
    real (prec)  ::  el_vol                            ! element volume
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/2,length_coord_array/2/))

    if (n_nodes == 3) n_points = 3
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 8) n_points = 9
    if (n_nodes == 4) n_points = 4

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0

    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1.D0+xnu)
    d11 = 1.D0*E/( (1.D0-xnu)*(1.D0+xnu) )
    d12 = xnu*E/( (1.D0-xnu)*(1.D0+xnu) )
    D(1,2) = d12
    D(2,1) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d44

    el_vol = 0.d0
    dNbardx = 0.d0
    !     -- Loop over the integration points for B bar method
    do kint = 1, n_points
       call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
    !
        do m_count = 1,n_nodes
           do n_count = 1,2
              dNbardx(m_count,n_count) = dNbardx(m_count,n_count) + &
                                         dNdx(m_count,n_count)*w(kint)*determinant
           end do
        end do
    !   Get the average element volume
        el_vol = el_vol + w(kint)*determinant
    end do
    !   Get the final form of vol_avg_shape function derivatives
       do m_count = 1,n_nodes
         do n_count = 1,2
             dNbardx(m_count,n_count) = dNbardx(m_count,n_count)/el_vol
         end do
       end do

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        B = 0.d0
        temp = 0.d0
        Bbar = 0.d0
        B(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(3,2:2*n_nodes:2) = dNdx(1:n_nodes,1)
        do m_count = 1,n_nodes
        temp(1:2:1,2*m_count-1) = dNbardx(m_count,1) - dNdx(m_count,1)
        temp(1:2:1,2*m_count) = dNbardx(m_count,2) - dNdx(m_count,2)
        end do

        Bbar = B + 1.d0/2.d0 *temp

       if ( element_identifier == 102 ) then

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        strain = strain + dstrain
        stress = matmul(D,strain)
        else if ( element_identifier == 202 ) then
        strain = matmul(Bbar,dof_total)
        dstrain = matmul(Bbar,dof_increment)
        strain = strain + dstrain
        stress = matmul(D,strain)
        end if
        stress_total(1:2) = stress(1:2)
        stress_total(4) = stress(3)
        stress_total(3) = 0.d0
        stress_total(5:6) = 0.d0
        strain_total(1:2) = strain(1:2)
        strain_total(4) = strain(3)
        strain_total(3) = -xnu/E*(stress(1)+stress(2))
        strain_total(5:6) = 0.d0
        p = sum(stress_total(1:3))/3.d0
        sdev = stress_total
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(1)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(2)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(3)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(4)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(5)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress_total(6)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(1)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(2)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(3)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(4)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(5)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'E23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain_total(6)*N(1:n_nodes) &
                *determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_linelast_2d_planestress

