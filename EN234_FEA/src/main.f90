program en234fea
  use Types
  use ParamIO
  use Globals
  use Controlparameters

  implicit none

  integer :: num_kk
!  Demo codes - basic 3D linear elasticity
!
!  infil = './input_files/linear_elastic_3d.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/linear_elastic_3d.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

! homework3
!      infil(1) = './input_files/Linear_elastic_2d_planestrain.in'
!   outfil(1) = './output_files/Linear_elastic_2d_planestrain.out'
!      infil(2) = './input_files/Linear_elastic_2d_planestress.in'
!   outfil(2) = './output_files/Linear_elastic_2d_planestress.out'
!      infil(3) = './input_files/Holeplate_2d_quad4.in'
!   outfil(3) = './output_files/Holeplate_2d_quad4.out'
!      infil(4) = './input_files/Holeplate_2d_quad8.in'
!   outfil(4) = './output_files/Holeplate_2d_quad8.out'
!      infil(5) = './input_files/Holeplate_2d_tri3.in'
!   outfil(5) = './output_files/Holeplate_2d_tri3.out'
!      infil(6) = './input_files/Holeplate_2d_tri6.in'
!   outfil(6) = './output_files/Holeplate_2d_tri6.out'
!      infil(7) = './input_files/Holeplate_2d_quad4_pstress.in'
!   outfil(7) = './output_files/Holeplate_2d_quad4_pstress.out'
!      infil(8) = './input_files/Holeplate_2d_quad8_pstress.in'
!   outfil(8) = './output_files/Holeplate_2d_quad8_pstress.out'
!      infil(9) = './input_files/Holeplate_2d_tri3_pstress.in'
!   outfil(9) = './output_files/Holeplate_2d_tri3_pstress.out'
!      infil(10) = './input_files/Holeplate_2d_tri6_pstress.in'
!   outfil(10) = './output_files/Holeplate_2d_tri6_pstress.out'

! do num_kk = 1,10
!    open (unit = IOR, file = infil(num_kk), status = 'old', ERR=500)
!    open (UNIT = IOW, FILE = outfil(num_kk), STATUS = 'unknown', ERR=500)


!  None of the files below will work until you write the codes that will use them!
!
!

!
!  Homework 4, crack tip elements and the J integral
!  infil = './input_files/crack_tri6.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/crack_tri6.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
!
!  Homework 5, small-strain B bar element - test with same files as in HW3, but
!  try approaching incompressible limit by making Poisson's ratio close to 0.5
!      infil(1) = './input_files/Linear_elastic_3d.in'
!  outfil(1) = './output_files/Linear_elastic_3d.out'
!      infil(2) = './input_files/Linear_elastic_3d_Bbar.in'
!   outfil(2) = './output_files/Linear_elastic_3d_Bbar.out'
!      infil(3) = './input_files/Holeplate_3d.in'
!   outfil(3) = './output_files/Holeplate_3d.out'
!      infil(4) = './input_files/Holeplate_3d_Bbar.in'
!   outfil(4) = './output_files/Holeplate_3d_Bbar.out'
!      infil(5) = './input_files/Linear_elastic_2d_planestrain.in'
!  outfil(5) = './output_files/Linear_elastic_2d_planestrain.out'
!      infil(6) = './input_files/Linear_elastic_2d_planestrain_Bbar.in'
!   outfil(6) = './output_files/Linear_elastic_2d_planestrain_Bbar.out'
!      infil(7) = './input_files/Linear_elastic_2d_planestress.in'
!   outfil(7) = './output_files/Linear_elastic_2d_planestress.out'
!      infil(8) = './input_files/Linear_elastic_2d_planestress_Bbar.in'
!   outfil(8) = './output_files/Linear_elastic_2d_planestress_Bbar.out'
!      infil(9) = './input_files/Holeplate_2d_tri6_pstress.in'
!   outfil(9) = './output_files/Holeplate_2d_tri6_pstress.out'
!      infil(10) = './input_files/Holeplate_2d_tri6_pstress_Bbar.in'
!   outfil(10) = './output_files/Holeplate_2d_tri6_pstress_Bbar.out'
!      infil(11) = './input_files/Holeplate_2d_tri6.in'
!   outfil(11) = './output_files/Holeplate_2d_tri6.out'
!      infil(12) = './input_files/Holeplate_2d_tri6_Bbar.in'
!   outfil(12) = './output_files/Holeplate_2d_tri6_Bbar.out'

  ! write(*,*) infil
!  do num_kk =1,12
!   write(*,*) 'LOOP'
!    open (unit = IOR, file = infil(num_kk), status = 'old', ERR=500)
!    open (UNIT = IOW, FILE = outfil(num_kk), STATUS = 'unknown', ERR=500)
!!
! Homework 6: small-strain Armstrong-Frederick kinematic hardening model
  infil = './input_files/cyclic_plastic_3d.in'
  open (unit = IOR, file = infil, status = 'old', ERR=500)
  outfil = './Output_files/cyclic_plastic_3d.out'
  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)

! Homework 7, stretch a hyperelastic bar, check stiffness.
!  infil = './input_files/Hyperelastic_bar_stretch.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/hyperelastic_bar_stretch.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
!!
!!  Homework 7, stretch and rotate a hyperelastic bar
!  infil = './input_files/Hyperelastic_stretch_rotate.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/hyperelastic_stretch_rotate.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
!
! Homework 7, stretch a hyperelastic plate with a central hole
!  infil = './input_files/Holeplate_hyperelastic.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/Holeplate_hyperelastic.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!!  Homework 8, solve the 2D Cahn-Hilliard equation
!  infil = './input_files/cahn_hilliard_2d_fine.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/cahn_hilliard_2d_fine.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)


!!  Homework 9, Dynamic fracture with explicit dynamics, finite strain Gurson model.
!  infil = './input_files/notch_fracture_dynamic.in'
!  open (unit = IOR, file = infil, status = 'old', ERR=500)
!  outfil = './Output_files/notch_fracture_dynamic.out'
!  open (UNIT = IOW, FILE = outfil, STATUS = 'unknown', ERR=500)
  write(*,*) 'in'
  call read_input_file
  write(*,*) 'out'
  if (printinitialmesh) call print_initial_mesh
  if (checkstiffness) call check_stiffness(checkstiffness_elementno)

  if (staticstep) then
      call compute_static_step
      if (checkstiffness) call check_stiffness(checkstiffness_elementno)
  endif

  if (explicitdynamicstep) then
   call explicit_dynamic_step
  endif


! end do
   write(6,*) ' Program completed successfully '
   stop
   500 write(6,*) ' Error opening input or output file '




end program en234fea
