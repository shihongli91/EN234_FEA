!     Subroutines for basic 3D hyperslastic elements



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_hyperslasticity_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    integer      :: n_points,kint,m_count, n_count,i,j,k,l,a

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)                            !
    real (prec)  ::  temp_bar(6,length_dof_array)      ! construct for Bbar
    real (prec)  ::  el_vol                            ! element volume
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  disp(3,length_coord_array/3)     ! Re-shaped displacement array disp(i,a) is ith coord of ath node
    real (prec)  ::  Gshear, Kbulk             ! Material properties
    real (prec)  ::  T_tau(3,3)                ! Cauchy stress tensor
    real (prec)  :: dNbardy(n_nodes,3),temp_p1(3*n_nodes,3*n_nodes),temp_p2(3*n_nodes,3*n_nodes)
    real (prec)  :: F_tau(3,3),F_tau_inv(3,3),detF_tau,dNdy(n_nodes,3),J_bar,dNdyvec(3*n_nodes)
    real (prec)  :: dNdy_aibk(3*n_nodes,3*n_nodes),dNdy_akbi(3*n_nodes,3*n_nodes),Pvec(3*n_nodes),Pmat(3*n_nodes,3*n_nodes)
    real (prec)  :: Fbar_tau(3,3),Fbar_tau_inv(3,3),detFbar_tau,B_CG(3,3), sig(3,3),tr_sig,Idenv(6),B_CGv(6)
    real (prec)  :: B_CG_inv_v(6),temp_D(6,6),IcrossI(6,6),IcrossBinv(6,6),BcrossBinv(6,6),B_bar(6,3*n_nodes)
    real (prec)  :: B_star(9,3*n_nodes),Gmat(6,9),Sigma(3*n_nodes,3*n_nodes),Svec(3*n_nodes),Smat(3*n_nodes,3*n_nodes)
    real (prec)  :: Q_matrix(3*n_nodes,3*n_nodes), P_matrix(3*n_nodes,3*n_nodes),dNbardyvec(3*n_nodes)
    real (prec)  :: temp_B_star(9,3*n_nodes),S(3,length_coord_array/3),B_CG_inv(3,3),Bqq,detB_CG,Iden(3,3)
    real (prec)  :: one,two,three,five,third,zero
    ! Utilities
    zero = 0.d0
    one = 1.d0
    two = 2.d0
    three = 3.d0
    third = 1.d0/3.d0
    five = 5.d0
    Iden =0.d0
    do i=1,3
     Iden(i,i)=1.d0
    end do
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Gshear
    !     element_properties(2)         Bbulk

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))
    disp = reshape(dof_total+dof_increment,(/3,length_dof_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    Gshear = element_properties(1)
    Kbulk  = element_properties(2)

    el_vol = 0.d0
    dNbardy = 0.d0
    J_bar = 0.d0
    temp_p1 = 0.d0
    !     -- Loop over the integration points for B bar method
    do kint = 1, n_points
       call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
    ! Compute the deformation gradient at current step F_tau
    F_tau = Iden
    do i=1,3
       do j= 1,3
          do k = 1,n_nodes
             F_tau(i,j) = F_tau(i,j) + dNdx(k,j)*disp(i,k)
           end do
       end do
    end do

    ! Compute the inverse of the current deformation gradient and J
       call  invert_small(F_tau,F_tau_inv,detF_tau)

    ! Convert the shape function derivatives into deformed configurations dNdy
       dNdy(1:n_nodes,1:3)= matmul(dNdx(1:n_nodes,1:3),F_tau_inv)

    ! Adding contributions to dNbardy and J_bar
        do m_count = 1,n_nodes
           do n_count = 1,3
              dNbardy(m_count,n_count) = dNbardy(m_count,n_count) + &
                                         detF_tau*dNdy(m_count,n_count)*w(kint)*determinant
           end do
        end do
        J_bar = J_bar + detF_tau*w(kint)*determinant
    !   Get the average element volume
        el_vol = el_vol + w(kint)*determinant

    !   Get temp_p1 for the P matrix
         dNdyvec(1:3*n_nodes) = reshape(transpose(dNdy(1:n_nodes,1:3)),(/3*n_nodes/))
         dNdy_aibk = spread(dNdyvec,dim=2,ncopies=3*n_nodes)*spread(dNdyvec,dim=1,ncopies=3*n_nodes)

        do i=1,n_nodes
           Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
           Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
        enddo
           dNdy_akbi = Pmat*transpose(Pmat)

         temp_p1 =temp_p1+ detF_tau*(dNdy_aibk-dNdy_akbi)*w(kint)*determinant

    end do

    !   Get the final form of dNbardy and detF0_tau
        J_bar = J_bar /el_vol
        dNbardy = dNbardy/el_vol/J_bar
        temp_p1 = temp_p1/el_vol/J_bar

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

    ! Compute the deformation gradient at current step F_tau
        F_tau = matmul(disp(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        F_tau = F_tau + Iden

    ! Compute the inverse of the current deformation gradient and J
       call  invert_small(F_tau,F_tau_inv,detF_tau)

    ! Convert the shape function derivatives into deformed configurations dNdy
       dNdy(1:n_nodes,1:3)= matmul(dNdx(1:n_nodes,1:3),F_tau_inv)

    ! Compute F-bar matrix
       Fbar_tau = 0.d0
       Fbar_tau = F_tau*(J_bar/detF_tau)**third

    ! Compute F_bar inverse
       call  invert_small(Fbar_tau,Fbar_tau_inv,detFbar_tau)
    ! detFbar_tau is actually J_bar


    ! Compute the Kirchhoff stress sig and material stiffness tangents D

    ! Left Cauchy-Green tensor
         B_CG = zero
          do i=1,3
             do j=1,3
                do k=1,3
                    B_CG(i,j) = B_CG(i,j) + Fbar_tau(i,k)*Fbar_tau(j,k)
                 enddo
             enddo
          enddo

        Bqq = B_CG(1,1)+B_CG(2,2)+B_CG(3,3)
    ! Cauchy stress T_tau
        T_tau = 0.d0
         do i=1,3
           do j=1,3
              T_tau(i,j) = Gshear/detFbar_tau**(five/three)*(B_CG(i,j) - third*Bqq*Iden(i,j)) + &
                           Kbulk * (detFbar_tau - one)*Iden(i,j)
           end do
         end do
    ! Kirchhoff stress sig
         Sig = detFbar_tau * T_tau
    !  vector form of kirchhoff stress tensor
       stress(1) = Sig(1,1)
       stress(2) = Sig(2,2)
       stress(3) = Sig(3,3)
       stress(4) = Sig(1,2)
       stress(5) = Sig(1,3)
       stress(6) = Sig(2,3)

    ! trace of kirchhoff stress tensor
       tr_sig = stress(1)+stress(2)+stress(3)

    ! Stiffness tangents D
       call  invert_small(B_CG,B_CG_inv,detB_CG)

       Idenv = [Iden(1,1),Iden(2,2),Iden(3,3),Iden(1,2),Iden(1,3),Iden(2,3)]
       B_CGv = [B_CG(1,1),B_CG(2,2),B_CG(3,3),B_CG(1,2),B_CG(1,3),B_CG(2,3)]
       B_CG_inv_v = [B_CG_inv(1,1),B_CG_inv(2,2),B_CG_inv(3,3),B_CG_inv(1,2),B_CG_inv(1,3),B_CG_inv(2,3)]
       IcrossI = spread(Idenv,dim=2,ncopies=6)*spread(Idenv,dim=1,ncopies=6)
       IcrossBinv =  spread(Idenv,dim=2,ncopies=6)*spread(B_CG_inv_v,dim=1,ncopies=6)
       BcrossBinv =  spread(B_CGv,dim=2,ncopies=6)*spread(B_CG_inv_v,dim=1,ncopies=6)

       temp_D = 0.d0
       do i=1,3
         temp_D(i,i)=1.d0
       enddo
       do i=4,6
         temp_D(i,i)=0.5d0
       enddo

       D = Gshear/detFbar_tau**(two/three)*temp_D +Gshear/three/detFbar_tau**(two/three)*(Bqq/three*IcrossBinv - &
           IcrossI-BcrossBinv)+ Kbulk*detFbar_tau*(detFbar_tau-0.5d0)*IcrossBinv
    ! Assemble B_bar
        B=0.d0
        B_bar = 0.d0
        temp_bar = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
        do m_count = 1,n_nodes
        temp_bar(1:3:1,3*m_count-2) = dNbardy(m_count,1) - dNdy(m_count,1)
        temp_bar(1:3:1,3*m_count-1) = dNbardy(m_count,2) - dNdy(m_count,2)
        temp_bar(1:3:1,3*m_count) = dNbardy(m_count,3) - dNdy(m_count,3)
        end do
        B_bar = B + third *temp_bar

    ! Assemble B_Star
        B_star = 0.d0
        temp_B_star = 0.d0
        B_star(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B_star(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B_star(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B_star(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B_star(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B_star(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B_star(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B_star(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B_star(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        do m_count = 1,n_nodes
           temp_B_star(1:3:1,3*m_count-2) = dNbardy(m_count,1) - dNdy(m_count,1)
           temp_B_star(1:3:1,3*m_count-1) = dNbardy(m_count,2) - dNdy(m_count,2)
           temp_B_star(1:3:1,3*m_count)   = dNbardy(m_count,3) - dNdy(m_count,3)
        end do
        B_star = B_star + third * temp_B_star

    ! Assemble Gmat
      Gmat = 0.d0
      Gmat(1,1) = 2.d0*B_CG(1,1)
      Gmat(1,4) = 2.d0*B_CG(1,2)
      Gmat(1,6) = 2.d0*B_CG(1,3)
      Gmat(2,2) = 2.d0*B_CG(2,2)
      Gmat(2,5) = 2.d0*B_CG(1,2)
      Gmat(2,8) = 2.d0*B_CG(2,3)
      Gmat(3,3) = 2.d0*B_CG(3,3)
      Gmat(3,7) = 2.d0*B_CG(1,3)
      Gmat(3,9) = 2.d0*B_CG(2,3)
      Gmat(4,1) = 2.d0*B_CG(1,2)
      Gmat(4,2) = 2.d0*B_CG(1,2)
      Gmat(4,4) = 2.d0*B_CG(2,2)
      Gmat(4,5) = 2.d0*B_CG(1,1)
      Gmat(4,6) = 2.d0*B_CG(2,3)
      Gmat(4,8) = 2.d0*B_CG(1,3)
      Gmat(5,1) = 2.d0*B_CG(1,3)
      Gmat(5,3) = 2.d0*B_CG(1,3)
      Gmat(5,4) = 2.d0*B_CG(2,3)
      Gmat(5,6) = 2.d0*B_CG(3,3)
      Gmat(5,7) = 2.d0*B_CG(1,1)
      Gmat(5,9) = 2.d0*B_CG(1,2)
      Gmat(6,2) = 2.d0*B_CG(2,3)
      Gmat(6,3) = 2.d0*B_CG(2,3)
      Gmat(6,5) = 2.d0*B_CG(1,3)
      Gmat(6,7) = 2.d0*B_CG(1,2)
      Gmat(6,8) = 2.d0*B_CG(3,3)
      Gmat(6,9) = 2.d0*B_CG(2,2)

    ! Compute Sum
    Sigma = 0.d0
    S =reshape(matmul(transpose(B),stress),(/3,length_dof_array/3/))
    do i=1,n_nodes
       Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
       Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
       Svec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
       Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
    end do

     Sigma = Pmat*transpose(Smat)

    ! Compute Qmat
    Q_matrix= 0.d0
    do i=1,n_nodes
           Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
           Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
    enddo
     Q_matrix = Pmat*transpose(Pmat)
    ! Compute Pmat
     P_matrix = 0.d0
     dNbardyvec(1:3*n_nodes) = reshape(transpose(dNbardy(1:n_nodes,1:3)),(/3*n_nodes/))
     temp_p2 = spread(dNbardyvec,dim=2,ncopies=3*n_nodes)*spread(dNbardyvec,dim=1,ncopies=3*n_nodes)
     P_matrix = temp_p1 - temp_p2
    ! Compute element residual

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B_bar),stress) &
                                        *w(kint)*determinant

    ! Compute element stiffness

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B_bar(1:6,1:3*n_nodes)),matmul(D,matmul(Gmat(1:6,1:9),B_star(1:9,1:3*n_nodes)))) &
            *w(kint)*determinant +(-Sigma+tr_sig*third*(P_matrix+Q_matrix))*w(kint)*determinant


    end do

    return
end subroutine el_hyperslasticity_3dbasic



!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_hyperslasticity_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
  
    integer      :: n_points,kint,m_count, n_count,i,j,k,l,a

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  temp_bar(6,length_dof_array)      ! construct for Bbar
    real (prec)  ::  Bbar(6,length_dof_array)          ! strain = Bbar*(dof_total+dof_increment), B-bar method
    real (prec)  ::  temp(6,length_dof_array)          ! temp matrix when constructing B-bar matrix
    real (prec)  ::  el_vol                            ! element volume
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  disp(3,length_coord_array/3)      ! Re-shaped displacement array disp(i,a) is ith coord of ath node
    real (prec)  ::  Gshear, Kbulk                    ! Material properties
    real (prec)  ::  T_tau(3,3),E_tau(3,3)                ! Cauchy stress tensor, lagrangian strain tensor
    real (prec)  :: dNbardy(n_nodes,3),temp_p1(3*n_nodes,3*n_nodes),temp_p2(3*n_nodes,3*n_nodes)
    real (prec)  :: F_tau(3,3),F_tau_inv(3,3),detF_tau,dNdy(n_nodes,3),J_bar,dNdyvec(3*n_nodes)
    real (prec)  :: dNdy_aibk(3*n_nodes,3*n_nodes),dNdy_akbi(3*n_nodes,3*n_nodes),Pvec(3*n_nodes),Pmat(3*n_nodes,3*n_nodes)
    real (prec)  :: Fbar_tau(3,3),Fbar_tau_inv(3,3),detFbar_tau,B_CG(3,3), sig(3,3),tr_sig,Idenv(6),B_CGv(6)
    real (prec)  :: B_CG_inv_v(6),temp_D(6,6),IcrossI(6,6),IcrossBinv(6,6),BcrossBinv(6,6),B_bar(6,3*n_nodes)
    real (prec)  :: B_star(9,3*n_nodes),Gmat(6,9),Sigma(3*n_nodes,3*n_nodes),Svec(3*n_nodes),Smat(3*n_nodes,3*n_nodes)
    real (prec)  :: Q_matrix(3*n_nodes,3*n_nodes), P_matrix(3*n_nodes,3*n_nodes),dNbardyvec(3*n_nodes),p,smises
    real (prec)  :: temp_B_star(9,3*n_nodes),S(3,length_coord_array/3),B_CG_inv(3,3),Bqq,detB_CG,Iden(3,3)
    real (prec)  :: one,two,three,five,third,zero
     ! Utilities
    zero = 0.d0
    one = 1.d0
    two = 2.d0
    three = 3.d0
    third = 1.d0/3.d0
    five = 5.d0
    Iden =0.d0
    do i=1,3
     Iden(i,i)=1.d0
    end do
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Gshear
    !     element_properties(2)         Bbulk


    x = reshape(element_coords,(/3,length_coord_array/3/))
    disp = reshape(dof_total+dof_increment,(/3,length_dof_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0


    Gshear = element_properties(1)
    Kbulk  = element_properties(2)

    el_vol = 0.d0
    dNbardy = 0.d0
    J_bar = 0.d0
    temp_p1 = 0.d0
    !     -- Loop over the integration points for B bar method
    do kint = 1, n_points
       call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
    ! Compute the deformation gradient at current step F_tau
        F_tau = matmul(disp(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
        F_tau = F_tau + Iden

    ! Compute the inverse of the current deformation gradient and J
       call  invert_small(F_tau,F_tau_inv,detF_tau)

    ! Convert the shape function derivatives into deformed configurations dNdy
       dNdy(1:n_nodes,1:3)= matmul(dNdx(1:n_nodes,1:3),F_tau_inv)

    ! Adding contributions to dNbardy and J_bar
        do m_count = 1,n_nodes
           do n_count = 1,3
              dNbardy(m_count,n_count) = dNbardy(m_count,n_count) + &
                                         detF_tau*dNdy(m_count,n_count)*w(kint)*determinant
           end do
        end do
        J_bar = J_bar + detF_tau*w(kint)*determinant
    !   Get the average element volume
        el_vol = el_vol + w(kint)*determinant

    !   Get temp_p1 for the P matrix
         dNdyvec(1:3*n_nodes) = reshape(transpose(dNdy(1:n_nodes,1:3)),(/3*n_nodes/))
         dNdy_aibk = spread(dNdyvec,dim=2,ncopies=3*n_nodes)*spread(dNdyvec,dim=1,ncopies=3*n_nodes)

        do i=1,n_nodes
           Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
           Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
        enddo
           dNdy_akbi = Pmat*transpose(Pmat)

         temp_p1 =temp_p1+ detF_tau*(dNdy_aibk-dNdy_akbi)*w(kint)*determinant

    end do

    !   Get the final form of dNbardy and detF0_tau
        J_bar = J_bar /el_vol
        dNbardy = dNbardy/el_vol/J_bar
        temp_p1 = temp_p1/el_vol/J_bar

     !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

    ! Compute the deformation gradient at current step F_tau
        F_tau = Iden
        do i=1,3
           do j= 1,3
              do k = 1,n_nodes
                 F_tau(i,j) = F_tau(i,j) + dNdx(k,j)*disp(i,k)
               end do
           end do
        end do

    ! Compute the inverse of the current deformation gradient and J
       call  invert_small(F_tau,F_tau_inv,detF_tau)

    ! Convert the shape function derivatives into deformed configurations dNdy
       dNdy(1:n_nodes,1:3)= matmul(dNdx(1:n_nodes,1:3),F_tau_inv)

    ! Compute F-bar matrix
       Fbar_tau = 0.d0
       Fbar_tau = F_tau*(J_bar/detF_tau)**third

    ! Compute F_bar inverse
       call  invert_small(Fbar_tau,Fbar_tau_inv,detFbar_tau)
    ! detFbar_tau is actually J_bar


    ! Compute the Kirchhoff stress sig and material stiffness tangents D

    ! Left Cauchy-Green tensor
         B_CG = zero
          do i=1,3
             do j=1,3
                do k=1,3
                    B_CG(i,j) = B_CG(i,j) + Fbar_tau(i,k)*Fbar_tau(j,k)
                 enddo
             enddo
          enddo

        Bqq = B_CG(1,1)+B_CG(2,2)+B_CG(3,3)
        ! Cauchy stress T_tau
        T_tau = 0.d0
         do i=1,3
           do j=1,3
              T_tau(i,j) = Gshear/detFbar_tau**(five/three)*(B_CG(i,j) - third*Bqq*Iden(i,j)) + &
                           Kbulk * (detFbar_tau - one)*Iden(i,j)
           end do
         end do
        ! Lagrangian strain
         E_tau = 0.5d0*(matmul(transpose(Fbar_tau),Fbar_tau) - Iden)

         strain(1) = E_tau(1,1)
         strain(2) = E_tau(2,2)
         strain(3) = E_tau(3,3)
         strain(4) = E_tau(1,2)
         strain(5) = E_tau(1,3)
         strain(6) = E_tau(2,3)
        ! Cauchy stress
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
         stress(5) = T_tau(1,3)
         stress(6) = T_tau(2,3)

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)

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
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(4)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(5)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'E23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + strain(6)*N(1:n_nodes)*determinant &
                *w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            end if
        end do
 
    end do


!
!        write(*,*) ' smises'
!        write(*,*) nodal_fieldvariables(17,1:n_nodes)
!        write(*,*)  ''
!        write(*,*) 's33'
!        write(*,*) nodal_fieldvariables(12,1:n_nodes)
!        write(*,*)  ''
!
!        stop
    return
end subroutine fieldvars_hyperslasticity_3dbasic


