C**********************************************************************
C                      USER MATERIAL SUBROUTINE UMAT FOR
C                      -----------------------------------
C                          ISOTROPIC VISCOPLASTICITY
C       -----------------------------------------------------------
C
C            THIS UMAT IS NOT FOR USE 3D VISCOPLASTICITY.
C            ENERGY POWER IS NOT AVAILABEL IN THIS SUBROUTINE
C
C Shihong Li. Nov. 2015
C**********************************************************************
C
C       CONTENTS OF PROPS VECTOR (*USER MATERIAL, CONSTANTS = 8) :
C       -----------------------------------------------------------
C       PROPS(I)
C       --------
C  1   E        ---  young's modulus
C  2   xnu      ---  poisson ratio
C  3   Y        ---  yield limit
C  4   e0       ---  material props
C  5   n        ---  hardening modulus
C  6   edot0    ---  material props
C  7   m        ---  strain rate sensitivity
C
C       STATE VARIABLES (*DEPVAR 2) :
C       ------------------------------
C      STATEV(1)   = Ee          -- EQUIV. SHEAR PLASTIC STRAIN  
C      STATEV(2)   = dEe         -- EQUIV. SHEAR PLASTIC STRAIN RATE

C*********************************************************************
     
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!     WARNING - the aba_param.inc file declares
!        Implicit real*8(a-h,o-z)
!     This means that, by default, any variables with
!     first letter between a-h or o-z are double precision.
!     The rest are integers.
!     Note that this also means that if you type a variable
!     name incorrectly, the compiler won't catch your typo.
!
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


!
!      DDSDDE(NTENS,NTENS)
!         Jacobian matrix of the constitutive model.
!         DDSDDE(I,J) defines the change in the Ith stress component
!         at the end of the time increment caused by an infinitesimal
!         perturbation of the Jth component of the strain increment array.
!         Unless you invoke the unsymmetric equation solution capability
!         for the user-defined material, ABAQUS/Standard will use only
!         the symmetric part of DDSDDE. The symmetric part of the matrix
!         is calculated by taking one half the sum of the matrix and its transpose.


!      STRESS(NTENS)
!         This array is passed in as the stress tensor at the beginning
!         of the increment and must be updated in this routine to be the
!         stress tensor at the end of the increment. If you specified
!         initial stresses (“Initial conditions,” Section 19.2.1), this
!         array will contain the initial stresses at the start of the
!         analysis. The size of this array depends on the value of NTENS
!         as defined below. In finite-strain problems the stress tensor
!         has already been rotated to account for rigid body motion in
!         the increment before UMAT is called, so that only the corotational
!         part of the stress integration should be done in UMAT. The
!         measure of stress used is “true” (Cauchy) stress.
!
!      STATEV(NSTATV)
!         An array containing the solution-dependent state variables.
!         These are passed in as the values at the beginning of the
!         increment unless they are updated in user subroutines USDFLD
!        (“USDFLD,” Section 25.2.39) or UEXPAN (“UEXPAN,” Section 25.2.20),
!        in which case the updated values are passed in. In all cases
!         STATEV must be returned as the values at the end of the increment.
!         The size of the array is defined as described in 
!        “Allocating space” in “User subroutines: overview,” Section 25.1.1.
!
!         In finite-strain problems any vector-valued or tensor-valued
!         state variables must be rotated to account for rigid body
!         motion of the material, in addition to any update in the
!         values associated with constitutive behavior. The rotation
!         increment matrix, DROT, is provided for this purpose.
!
!      SSE, SPD, SCD
!         Specific elastic strain energy, plastic dissipation, and
!         “creep” dissipation, respectively. These are passed in as
!         the values at the start of the increment and should be
!         updated to the corresponding specific energy values at
!         the end of the increment. They have no effect on the solution,
!         except that they are used for energy output.
!
!     Only in a fully coupled thermal-stress analysis
!      RPL
!         Volumetric heat generation per unit time at the end of the increment
!         caused by mechanical working of the material.
!
!     DDSDDT(NTENS)
!          Variation of the stress increments with respect to the temperature.
!
!     DRPLDE(NTENS)
!           Variation of RPL with respect to the strain increments.
!
!     DRPLDT
!           Variation of RPL with respect to the temperature.
!
!     Variables that can be updated
!
!     PNEWDT
!        Ratio of suggested new time increment to the time increment being
!        used (DTIME, see discussion later in this section). This variable
!        allows you to provide input to the automatic time incrementation
!        algorithms in ABAQUS/Standard (if automatic time incrementation is chosen).
!        For a quasi-static procedure the automatic time stepping that ABAQUS/Standard
!        uses, which is based on techniques for integrating standard creep laws
!        (see “Quasi-static analysis,” Section 6.2.5), cannot be controlled from within
!        the UMAT subroutine.
!        PNEWDT is set to a large value before each call to UMAT.
!        If PNEWDT is redefined to be less than 1.0, ABAQUS/Standard must abandon the
!        time increment and attempt it again with a smaller time increment. The
!        suggested new time increment provided to the automatic time integration
!        algorithms is PNEWDT × DTIME, where the PNEWDT used is the minimum value
!        for all calls to user subroutines that allow redefinition of PNEWDT for this
!        iteration.
!        If PNEWDT is given a value that is greater than 1.0 for all calls to user
!        subroutines for this iteration and the increment converges in this iteration,
!        ABAQUS/Standard may increase the time increment. The suggested new time increment
!        provided to the automatic time integration algorithms is PNEWDT × DTIME, where
!        the PNEWDT used is the minimum value for all calls to user subroutines for
!        this iteration.
!        If automatic time incrementation is not selected in the analysis procedure,
!        values of PNEWDT that are greater than 1.0 will be ignored and values of
!        PNEWDT that are less than 1.0 will cause the job to terminate.
!
!    Variables passed in for information
!
!     STRAN(NTENS)
!         An array containing the total strains at the beginning of the increment.
!         If thermal expansion is included in the same material definition, the
!         strains passed into UMAT are the mechanical strains only (that is, the
!         thermal strains computed based upon the thermal expansion coefficient have
!         been subtracted from the total strains). These strains are available for output
!         as the “elastic” strains.
!
!         In finite-strain problems the strain components have been rotated to account for
!         rigid body motion in the increment before UMAT is called and are approximations
!         to logarithmic strain.

!     DSTRAN(NTENS)
!         Array of strain increments. If thermal expansion is included in the same
!         material definition, these are the mechanical strain increments (the total
!         strain increments minus the thermal strain increments).
!
!     TIME(1)
!         Value of step time at the beginning of the current increment.
!
!     TIME(2)
!          Value of total time at the beginning of the current increment.
!
!     DTIME
!        Time increment.
!
!     TEMP
!         Temperature at the start of the increment.
!
!     DTEMP
!         Increment of temperature.
!
!     PREDEF
!        Array of interpolated values of predefined field variables at this point
!        at the start of the increment, based on the values read in at the nodes.
!
!      DPRED
!        Array of increments of predefined field variables.
!
!      CMNAME
!        User-defined material name, left justified. Some internal material models are given names starting with the “ABQ_” character string. To avoid conflict, you should not use “ABQ_” as the leading string for CMNAME.
!
!      NDI
!        Number of direct stress components at this point.
!
!      NSHR
!        Number of engineering shear stress components at this point.
!
!      NTENS
!        Size of the stress or strain component array (NDI + NSHR).
!
!      NSTATV
!         Number of solution-dependent state variables that are associated with
!         this material type (defined as described in “Allocating space” in “User
!         subroutines: overview,” Section 25.1.1).
!
!      PROPS(NPROPS)
!         User-specified array of material constants associated with this user material.
!
!      NPROPS
!         User-defined number of material constants associated with this user material.
!
!      COORDS
!         An array containing the coordinates of this point. These are the current
!         coordinates if geometric nonlinearity is accounted for during the step
!         (see “Procedures: overview,” Section 6.1.1); otherwise, the array contains
!         the original coordinates of the point.
!
!     DROT(3,3)
!          Rotation increment matrix. This matrix represents the increment of rigid
!          body rotation of the basis system in which the components of stress
!          (STRESS) and strain (STRAN) are stored. It is provided so that vector- or
!          tensor-valued state variables can be rotated appropriately in this subroutine:
!          stress and strain components are already rotated by this amount before UMAT
!          is called. This matrix is passed in as a unit matrix for small-displacement
!          analysis and for large-displacement analysis if the basis system for the
!          material point rotates with the material (as in a shell element or when a
!          local orientation is used).
!
!      CELENT
!          Characteristic element length, which is a typical length of a line across
!          an element for a first-order element; it is half of the same typical length
!          for a second-order element. For beams and trusses it is a characteristic length
!          along the element axis. For membranes and shells it is a characteristic length
!          in the reference surface. For axisymmetric elements it is a characteristic length
!          in the  plane only. For cohesive elements it is equal to the constitutive
!          thickness.
!
!      DFGRD0(3,3)
!          Array containing the deformation gradient at the beginning of the increment.
!          See the discussion regarding the availability of the deformation gradient for
!          various element types.
!
!     DFGRD1(3,3)
!            Array containing the deformation gradient at the end of the increment.
!           The components of this array are set to zero if nonlinear geometric effects
!           are not included in the step definition associated with this increment. See
!           the discussion regarding the availability of the deformation gradient for
!           various element types.
!
!      NOEL
!           Element number.
!
!      NPT
!           Integration point number.
!
!      LAYER
!          Layer number (for composite shells and layered solids).
!
!      KSPT
!          Section point number within the current layer.
!
!      KSTEP
!         Step number.
!
!     KINC
!         Increment number.

!      user coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
!      and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT
!
!     Local variables
       real*8  stress_old(ntens),stress_new(ndi,ndi),
     +         Ee_old, dEe_old, Ee_new,dEe_new 
C----------------------------------------------------------------------
C                       INITIALIZATION
C----------------------------------------------------------------------
        stress_old = stress
        call zerom(stress_new)
C----------------------------------------------------------------------
C                       START COMPUTATION
C----------------------------------------------------------------------

!
!       At the start of  an abaqus calculation the state variuables are passed into 
!       Umat with zero values. At this point, the time total_time and step_time both 
!       have a value equal to zero and dtime is equal to 1.d0
!
      
      if ((time(1) .eq. 0.d0) .and. (kstep .eq. 1))  then
         statev(1) = zero
         statev(2) = zero
      end if 
!      Store the value of the state variables at the beginning of the time step
         Ee_old  = statev(1)
         dEe_old = statev(2)      
 
!      Time integration
      call  viscoplastic(
!      input terms
     +                   props,nprops,ndi,ntens,dstran,
     +                   dtime,Ee_old,dEe_old,stress_old,
!      output terms 
     +                   stress, Ee_new,dEe_new,ddsdde)
      

!
!       Update the state variables 
!
        statev(1)  =  Ee_new
        statev(2)  =  dEe_new
!        
      !     write(*,*)  statev(1)
       return 
      end


!**************************************************************************
      subroutine viscoplastic(
!      input terms
     +                   props,nprops,ndi,ntens,dstran,
     +                   dtime,Ee_old,dEe_old,stress_old,
!      output terms 
     +                   stress, Ee_new,dEe_new,ddsdde)
!      time integration procedure       
!***************************************************************************
      implicit none
      
      integer          nprops,ndi,ntens, i,j,k,l
      double precision props(nprops),dstran(ntens),dtime,
     +                 Ee_old, dEe_old,ddsdde(ntens,ntens), 
     +                 stress_old(ntens),Ee_new,dEe_new,stress_new(3,3),
     +                 check, stress(ntens),eps_v,deps0(ntens),sig_v,
     +                 sig_0(ntens),sig_star_vec(ntens),dEe_mid,
     +                 sig_star(ndi,ndi),sig_e, dEe,err,tol,temp_1,
     +                 F_fun, dFde,gamma,Ctang(3,3,3,3)
      double precision E,xnu,Y,e0,edot0,n,m,iter,maxit
      double precision zero,one,half,two,third,Iden(3,3)
      parameter       (zero=0.d0, one=1.d0, half=0.5d0, two=2.d0,
     +                 third= 1.d0/3.d0)
     
        E   = props(1)
        xnu = props(2)
        Y   = props(3)
        e0  = props(4)
        n   = props(5)
        edot0 = props(6)
        m   = props(7)


       call onem(Iden)
!
!     when umat is called for the first iteration in an increment, abaqus does
!     not update the strain, and it passes in a zero strain increment. In this 
!     case do not perform tho ingeration. Compute the elastic stiffness and return.
!
        check = 0.d0
        do i=1,ntens
           check = check + dabs(dstran(i))
        end do
        
        if (check .lt. 1.d-15) then
            dEe_new     = zero
	      Ee_new      = Ee_old
    !   tangent stiffness for elastic material
	      ddsdde(1,1) = 1.d0-xnu
	      ddsdde(1,2) = xnu
	      ddsdde(1,3) = xnu
	      ddsdde(2,1) = xnu
	      ddsdde(2,2) = 1.d0-xnu
	      ddsdde(2,3) = xnu
	      ddsdde(3,1) = xnu
	      ddsdde(3,2) = xnu
	      ddsdde(3,3) = 1.d0-xnu
	      ddsdde(4,4) = 0.5d0*(1.d0-2.d0*xnu)
	      ddsdde(5,5) = ddsdde(4,4)
	      ddsdde(6,6) = ddsdde(4,4)
	      ddsdde = ddsdde*E/( (1.d0+xnu)*(1.d0-2.d0*xnu) )
    !
    !     NOTE: ABAQUS uses engineering shear strains,
    !     i.e. stran(ndi+1) = 2*e_12, etc...
          do i = 1,ntens
	      do j = 1,ntens
	         stress(i) = stress(i) + ddsdde(i,j)*dstran(j)
	      end do
	    end do

          return
        end if

!-----------------------------------------------------------------------
!      For next steps 
!-----------------------------------------------------------------------

!        Construct the deviatoric part of strain and stress vector.
         eps_v = dstran(1)+dstran(2)+dstran(3)
      !
         do i=1,3
           deps0(i) = dstran(i)-third*eps_v
         end do
           
         do i=4,6
          deps0(i) = dstran(i)*half  
         end do
          
         sig_v = stress_old(1)+stress_old(2)+stress_old(3)
         sig_0 = stress_old
    !
         do i=1,3
           sig_0(i) = sig_0(i) - third*sig_v
         end do

    !     Compute the elastic predictors 
         
         sig_star_vec = sig_0 + E/(one+xnu)*deps0
         ! transform it into matrix form for later use
     
        do i=1,3     
          sig_star(i,i) = sig_star_vec(i)
        end do
         !
         sig_star(1,2) = sig_star_vec(4)
         sig_star(2,1) = sig_star_vec(4)
         sig_star(1,3) = sig_star_vec(5)
         sig_star(3,1) = sig_star_vec(5)
         sig_star(2,3) = sig_star_vec(6)
         sig_star(3,2) = sig_star_vec(6)
         
         ! Calculate the Mises stress
      
         sig_e = dsqrt(dot_product(sig_star_vec(1:3),sig_star_vec(1:3))+
     +           two*dot_product(sig_star_vec(4:6),sig_star_vec(4:6)))
     +           *dsqrt(1.5d0)
          
         ! Try to sole the plastic strain increment(Newton-Rapshon iterations)
             
          if (sig_e*edot0 .eq. 0.d0) then
              dEe = 0.d0
          else
            !write(*,*) 'coming in'
            dEe = 1.d-15  ! changing back to 1e-15 when this is not working
            err = Y
            tol = 1.d-6*Y
            iter = 1
            maxit = 50
            do while (err .gt. tol) 
                temp_1 = (one+(Ee_old+dEe)/e0)**(one/n)*
     +                    (dEe/(dtime*edot0))**(one/m)
                F_fun = sig_e/Y - 1.5d0*dEe*E/(Y*(one+xnu))- temp_1
                dFde = -1.5d0*E/Y/(one+xnu) - temp_1*(one/n/
     +                 (Ee_old+dEe+e0) + one/m/dEe)
                dEe_mid = dEe - F_fun/dFde
                err = dsqrt(F_fun*F_fun)
                iter = iter + 1
                
             if (dEe_mid .lt. 0.d0) then
                dEe = dEe/10.d0
             else 
                dEe = dEe_mid
             end if 

          
               if (iter .EQ. maxit) then
                write(*,*) 'Warning in solving dEe: 
     +                      Max iterations exceeded.'
                stop
               end if
                
            end do
      
          end if
             dEe_new = dEe

          if (dEe_new .LT. zero) then
            write(*,*)  'plastic strain increment cannot be negative'
            stop
          end if

    !       Update the total plastic strain 
              Ee_new = Ee_old + dEe_new
            !  write(*,*)  Ee_old       
    !       Update the stress tensor
            stress_new = (one - 1.5d0*E*dEe_new/(one+xnu)/sig_e)*
     +                   sig_star +(sig_v + E*eps_v/(one -two*xnu))
     +                   *Iden*third

!           Update the stress
!
           do i=1,ndi
             stress(i) = stress_new(i,i)
           end do
             stress(ndi+1) = stress_new(1,2)
             stress(ndi+2) = stress_new(1,3)
             stress(ndi+3) = stress_new(2,3)
              
   
    !       Calculathe gamma
            gamma = 1.5d0*E/(one+xnu)/sig_e + (one - 1.5d0*E*dEe_new
     +              /(one+xnu)/sig_e)*(one/n/(e0+Ee_old+dEe_new)
     +              +one/m/dEe_new)
     
    !       Calculate the tangent stiffness
            do i=1,3
              do j=1,3
                do k=1,3
                  do l=1,3
                     Ctang(i,j,k,l) = E/(one+xnu)*(one - 1.5d0*E*dEe_new
     +               /(one+xnu)/sig_e)*(half*(Iden(i,k)*Iden(j,l)+
     +               Iden(j,k)*Iden(i,l))-third*Iden(i,j)*Iden(k,l))+
     +               E/(one+xnu)*9.d0/4.d0*(dEe_new-one/gamma)/(one+xnu)
     +               /sig_e*sig_star(i,j)*sig_star(k,l)/sig_e/sig_e + 
     +               third*E/(one-two*xnu)*Iden(i,j)*Iden(k,l)
                  enddo
                enddo
              enddo
            enddo
            
    !         Construct the ddsdde using Cijkl above,not assuming C is symmetric here.
          ddsdde = 0.d0
          do i=1,3
            do j=1,3
               ddsdde(i,j) = Ctang(i,i,j,j)
            end do
          end do
            do i =1,3
               ddsdde(i,4) = Ctang(i,i,1,2)
               ddsdde(4,i) = Ctang(1,2,i,i)
            end do
               ddsdde(4,4) = Ctang(1,2,1,2)
            do i=1,3
               ddsdde(i,5) = Ctang(i,i,1,3)
               ddsdde(5,i) = Ctang(1,3,i,i)
            end do
               ddsdde(5,5) = Ctang(1,3,1,3)
               ddsdde(4,5) = Ctang(1,2,1,3)
               ddsdde(5,4) = Ctang(1,3,1,2)
            do i=1,3
               ddsdde(i,6) = Ctang(i,i,2,3)
               ddsdde(6,i) = Ctang(2,3,i,i)
            end do
               ddsdde(6,6) = Ctang(2,3,2,3)
               ddsdde(4,6) = Ctang(1,2,2,3)
               ddsdde(6,4) = Ctang(2,3,1,2)
               ddsdde(5,6) = Ctang(1,3,2,3)
               ddsdde(6,5) = Ctang(2,3,1,3)
       return
       end
!----------------------------------------------------------------------
!     Utility subroutines
!----------------------------------------------------------------------
      subroutine zerom(A)
!    construct a 3 by 3 zero matrix
      implicit none
      real*8 A(3,3)
      integer i,j
      
      do i=1,3
        do j=1,3
          A(i,j) = 0.d0
        end do
      end do
      
      return 
      end

!----------------------------------------------------------------------
      subroutine onem(A)
!    constuct a 3 by 3 idenitiy 
      implicit none
      real*8 A(3,3)
      integer i,j
      
      do i=1,3
        do j=1,3
          if (i .eq. j) then
            A(i,j) = 1.d0
          else 
            A(i,j) = 0.d0
          end if
        end do
      end do
      
      return 
      end
