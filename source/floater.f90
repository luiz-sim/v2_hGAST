!----------------------------------------------------------------------
 Subroutine floater_simulator ( run_code, iter )
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   integer, intent(in)  :: run_code, iter
   integer              :: N, nf

!
! run_code = 0: initialization
! run_code = 1: normal run
! run_code = 2: end timestep
! run_code = 3: write results
!

   if     (run_code == 1) then

      call update_velocity_buffer ( iter )
      call calc_wave_force        ( iter )               !diffraction same for all iterations
      call calc_retardation       ( iter )               !retardation history same for all iterations
      call Morison_floater
      call calc_hydro_matrices

   elseif (run_code == 3) then

      call WRITEOUT_floater

   elseif (run_code == 0) then

      write(10,*) 'initialisation of floater simulation'

         N  = 6 *nfloater_hyd

      do nf = 1, nfloater_hyd
         Allocate ( floater(nf)% AC_retardation   (6, N) )
         Allocate ( floater(nf)% AC_total         (6, N) )
         Allocate ( floater(nf)% AM_added_infinity(6, N) )
         Allocate ( floater(nf)% AM_total         (6, N) )
      enddo

      call floater_init
      call diffraction_init
      call retardation_init !!_5  !qq jim
      if  (DragMod  == 1) &
      call Morison_floater_init

   elseif (run_code == 2) then

      Deallocate ( floater )

      if (( Wavemod >= 1 ).and.( Wavemod <= 5 )) then
         Deallocate ( spectral_freq , spectral_ampl    )
         Deallocate ( spectral_phase, spectral_wavenum )
      endif

   else

      write(10,*) 'Error in floater Simulator'
      stop

   endif


 END Subroutine floater_simulator
!----------------------------------------------------------------------
 Subroutine calc_hydro_matrices
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: vel_sqr(6), scalef
   integer :: i, nf, j, N


      N  = nfloater_hyd*6

   do nf = 1, nfloater_hyd
      j  = (nf-1)*6

      floater(nf)% AM_total (1:6,1  :N  ) =  floater(nf)% AM_added_infinity     (1:6,1  :N  )
      floater(nf)% AM_total (1:6,j+1:j+6) =  floater(nf)% AM_total              (1:6,j+1:j+6) &
                                           + floater(nf)% AM_floater_structural (1:6,1  :6  )
      floater(nf)% AC_total (1:6,1  :N  ) =  floater(nf)% AC_retardation        (1:6,1  :N  )
      floater(nf)% AK_total (1:6,1  :6  ) =  floater(nf)% AK_external           (1:6,1  :6  ) &
                                           + floater(nf)% AK_hydrostatic        (1:6,1  :6  ) &
                                           + floater(nf)% AK_floater_structural (1:6,1  :6  )
      floater(nf)% Q_total  (    1  :6  ) =  floater(nf)% AQ_exciting           (    1  :6  ) &
                                           + floater(nf)% AQ_Newman             (    1  :6  ) &
                                           - floater(nf)% AQ_retardation        (    1  :6  ) &
!!                                         + floater(nf)% AQ_morison            (    1  :6  ) &
                                           + floater(nf)% AQ_external           (    1  :6  )

!------ addition of total buoyancy, floater's weight and Moorings_weight
      floater(nf)% Q_total  (          3) =  floater(nf)% Q_total               (3          )         &
                                           + floater(nf)% Buoyancy_total                              &
                                           - floater(nf)% AM_floater_structural (3  ,3      )*gravity &
                                           - floater(nf)% Moorings_weight

!------ addition of gravity moments in case of floater's Xcg offsets wrt centre line [i.e. flooded damage case]
      floater(nf)% Q_total  (          4) =  floater(nf)% Q_total               (4          )         &
                                           - floater(nf)% AM_floater_structural (3  ,4      )*gravity
      floater(nf)% Q_total  (          5) =  floater(nf)% Q_total               (5          )         &
                                           - floater(nf)% AM_floater_structural (3  ,5      )*gravity

!------ Add the External damping
      if     (ITYPE_Damp_ext==1) then

                            scalef =   1.d0
!qq      if (time_h<100.d0) scalef = 100.d0

         floater(nf)% AC_total     (1:6,j+1:j+6) =  floater(nf)% AC_total       (1:6,j+1:j+6) &
                                                  + floater(nf)% AC_external    (1:6,  1:  6) * scalef
         floater(nf)% AQ_damp_quad (1:6        ) = 0.d0;
      elseif (ITYPE_Damp_ext==2) then
         do i = 1, 6
            vel_sqr(i) = floater(nf)% new_velocity(i) * dabs(floater(nf)% new_velocity(i))
         enddo

         floater(nf)% AQ_damp_quad (1:6) = -matmul ( floater(nf)%AC_external (1:6,1:6), vel_sqr(1:6) )
         floater(nf)% Q_total      (1:6) =  floater(nf)% Q_total      (1:6) &
                                          + floater(nf)% AQ_damp_quad (1:6)
      endif
   enddo !nf


 END Subroutine calc_hydro_matrices
!----------------------------------------------------------------------
 Subroutine calc_wave_force ( iter )
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   real(8) :: Lp(6), Lm(6)
   integer :: ifr, iter, nf


   if ( iter > 1 ) return

!------ First order Exciting force
   do    nf = 1, nfloater_hyd
      floater(nf)%AQ_exciting(1:6) = 0.d0;
      floater(nf)%AQ_Newman  (1:6) = 0.d0;

      do ifr= 1, no_spectral_comp  !add all wave components

         floater(nf)%AQ_exciting(    1:6) =                                                   &
         floater(nf)%AQ_exciting(    1:6)                                                     &
 +       floater(nf)%wave_force (ifr,1:6) * spectral_ampl(ifr)                                &
 * dcos( floater(nf)%wave_phase (ifr,1:6) - spectral_freq(ifr)*time_h + spectral_phase(ifr) )
      enddo !ifr
         floater(nf)%AQ_exciting(    1:6) = floater(nf)%AQ_exciting (1:6) * Start_Fact
   enddo !nf


   if (NewmanMod == 0) return
!------ Newman's approximation for 2nd order difference-frequency terms
   do    nf      = 1, nfloater_hyd
         Lp(1:6) = 0.d0;
         Lm(1:6) = 0.d0;
      do ifr     = 1, no_spectral_comp  !add all wave components
         Lp(1:6) = Lp(1:6) + dsqrt ( 2.d0 * floater(nf)%wave_force_Dft_pos (ifr,1:6) )* spectral_ampl (ifr)  &
                           * dcos  (      - spectral_freq(ifr)*time_h                 + spectral_phase(ifr) )

         Lm(1:6) = Lm(1:6) + dsqrt (-2.d0 * floater(nf)%wave_force_Dft_neg (ifr,1:6) )* spectral_ampl (ifr)  &
                           * dcos  (      - spectral_freq(ifr)*time_h                 + spectral_phase(ifr) )
      enddo !ifr
         floater(nf)%AQ_Newman (1:6) = ( Lp(1:6)**2 - Lm(1:6)**2 ) * Start_Fact
   enddo !nf


 END Subroutine calc_wave_force
!----------------------------------------------------------------------
 Subroutine calc_retardation ( iter )
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   integer :: k, iter, nf, nfj,j


   if ( iter > 1 ) return


   do nf = 1, nfloater_hyd

       floater(nf)% AQ_retardation (1:6        ) = 0.d0;

    do nfj= 1, nfloater_hyd

       j=(nfj-1)*6
       floater(nf)% AC_retardation (1:6,j+1:j+6) = 0.5d0 * delta_t * floater(nf)% retardation_function (1:6,j+1:j+6,1)

       floater(nf)% AQ_retardation (1:6        ) =         floater(nf)% AQ_retardation       (1:6)        +  0.5d0 * delta_t * &
                                                   matmul( floater(nf)% retardation_function (1:6,j+1:j+6, retardation_length   ), &
                                                           floater(nfj)% velocity_buffer (1:6, 1)                                    )

       do k = 2, retardation_length-1
       floater(nf)% AQ_retardation(1:6         ) =         floater(nf)% AQ_retardation      (1:6)         +          delta_t * &
                                                   matmul( floater(nf)% retardation_function(1:6,j+1:j+6, retardation_length-k+1), &
                                                           floater(nfj)% velocity_buffer (1:6, k)                                    )
       enddo !k
    enddo !nfj
   enddo !nf


 END Subroutine calc_retardation
!----------------------------------------------------------------------
 Subroutine update_velocity_buffer ( iter )
!----------------------------------------------------------------------

 Use Hydro

   implicit none

   integer :: j, iter, nf


   if ( iter > 1 ) return

   do nf   = 1, nfloater_hyd
      do j = 1, retardation_length-2
         floater(nf)% velocity_buffer(1:6, j                   ) = floater(nf)% velocity_buffer (1:6, j+1)
      enddo
         floater(nf)% velocity_buffer(1:6, retardation_length-1) = floater(nf)% new_velocity    (1:6     )
   enddo !nf


 END Subroutine update_velocity_buffer
!----------------------------------------------------------------------
 Subroutine floater_init
!----------------------------------------------------------------------

 Use Hydro
 Use Paths
 Use Cbeam ! UflInit, IQfl_active_tr

   implicit none

   integer :: ityp, i, j, nf, N
   real(8) :: ex(6)


   N = nfloater_hyd * 6

   open (24,file=trim(file_floater)) !'floater.inp'

   do nf = 1, nfloater_hyd

      read(24,*) !title

      if (nf==2) read(24,*) (RC_hyd(j),j=1,3)
 
      read(24,*) !Draft
      read(24,*)  floater(nf)%Draft

      read(24,*) !Added mass at infinity
      read(24,*) ityp !1:diagonal, 2:symmetrix, 3:full
         call read_matrix (floater(nf)%AM_added_infinity    , 6, N, ityp, 24)

      read(24,*) !Buoyancy
      read(24,*)  floater(nf)%Buoyancy_total

      read(24,*) !Hydrostatic stiffness
      read(24,*) ityp !1:diagonal, 2:symmetrix, 3:full
         call read_matrix (floater(nf)%AK_hydrostatic       , 6, 6, ityp, 24)

      read(24,*) !structural mass of floater
      read(24,*) ityp !1:diagonal, 2:symmetrix, 3:full
         call read_matrix (floater(nf)%AM_floater_structural, 6, 6, ityp, 24)

      read(24,*) !structural stiffness of floater
      read(24,*) ityp !1:diagonal, 2:symmetrix, 3:full
         call read_matrix (floater(nf)%AK_floater_structural, 6, 6, ityp, 24)

      read(24,*) !External Damping type
      read(24,*) ITYPE_Damp_ext     !1: linear external damping, 2: quadratic damping added to rhs

      read(24,*) !External Damping
      read(24,*) ityp !1:diagonal, 2:symmetrix, 3:full
         call read_matrix (floater(nf)%AC_external          , 6, 6, ityp, 24)

      read(24,*) !External stiffness
      read(24,*) ityp !1:diagonal, 2:symmetrix, 3:full
         call read_matrix (floater(nf)%AK_external          , 6, 6, ityp, 24)

      read(24,*) !Moorings weight
      read(24,*)  floater(nf)%Moorings_weight


                  floater(nf)%AQ_external = 0.d0;
!jim
! for linearization
!!    read(24,*) !** External forces

!!    do i=1,6
!!       read(24,*) floater(nf)%AQ_external(i)
!!    enddo


!----- Output information
      if (nf==2) write(10,*  )'R connection'
      if (nf==2) write(10,100) (RC_hyd(j),j=1,3)
      write (10,*)  
      write (10,*) 'Draft'
      write (10,100)   floater(nf)% Draft
      write (10,*)  
      write (10,*) 'AM_added_infinity'
      write (28,'(a)')'**AM_added_infinity'

                 ex(1:3)=0.d0; ex(4:6)=1.d0;
      do i = 1, 6
      write (10,100)  (floater(nf)% AM_added_infinity     (i,j),j=1,N)
      write (28,100)  (floater(nf)% AM_added_infinity     (i,j) &
                       *scal_hyd**(3.0d0+ex(i)+ex(j))           &
                       *rho_sweet/rho_salt                     ,j=1,N)
      enddo
      write (10,*)  
      write (10,*) 'Buoyancy'
      write (28,'(a)') '**Buoyancy'
      write (10,100)   floater(nf)% Buoyancy_total
      write (28,100)   floater(nf)% Buoyancy_total              &
                       *scal_hyd**3.0d0*rho_sweet/rho_salt
      write (10,*)  
      write (10,*) 'AK_hydrostatic'
      write (28,'(a)') '**AK_hydrostatic'
      do i = 1, 6
      write (10,100)  (floater(nf)% AK_hydrostatic        (i,j),j=1,6)
      write (28,100)   floater(nf)% AK_hydrostatic        (i,i) &
                       *scal_hyd**(2.0d0+ex(i)+ex(i))           &
                       *rho_sweet/rho_salt
      enddo
      write (10,*)  
      write (10,*) 'AM_floater_structural'
      write (28,'(a)') '**AM_floater_structural'
      do i = 1, 6
      write (10,100)  (floater(nf)% AM_floater_structural (i,j),j=1,6)
      do j = i, 6
      write (28,100)   floater(nf)% AM_floater_structural (i,j) &
                       *scal_hyd**(3.0d0+ex(i)+ex(j))           &
                       *rho_sweet/rho_salt
      enddo!j
      enddo!i
      write (10,*)  
      write (10,*) 'AK_floater_structural'
      write (28,'(a)') '**AK_floater_structural'
      do i = 1, 6
      write (10,100)  (floater(nf)% AK_floater_structural (i,j),j=1,6)
      write (28,100)   floater(nf)% AK_floater_structural (i,i) &
                       *scal_hyd**(2.0d0+ex(i)+ex(i))           &
                       *rho_sweet/rho_salt
      enddo
      write (10,*)  
      write (10,*) 'ITYPE_Damp_ext'
      write (10,*)  ITYPE_Damp_ext
      write (10,*)  
      write (10,*) 'AC_external'
      write (28,'(a)') '**AC_external'
      do i = 1, 6
      write (10,100)  (floater(nf)% AC_external           (i,j),j=1,6)
      write (28,100)   floater(nf)% AC_external           (i,i) &
                       *scal_hyd**(2.5d0+ex(i)+ex(i))           &
                       *rho_sweet/rho_salt
      enddo
      write (10,*)  
      write (10,*) 'AK_external'
      write (28,'(a)') '**AK_external'
      do i = 1, 6
      write (10,100)  (floater(nf)% AK_external           (i,j),j=1,6)
      do j = i, 6
      write (28,100)   floater(nf)% AK_external           (i,j) &
                       *scal_hyd**(2.0d0+ex(i)+ex(j))           &
                       *rho_sweet/rho_salt
      enddo!j
      enddo!i
      write (10,*)  
      write (10,*) 'Moorings weight-Pretension'
      write (28,'(a)') '**Moorings weight-Pretension'
      write (10,100)   floater(nf)% Moorings_weight
      write (28,100)   floater(nf)% Moorings_weight             &
                       *scal_hyd**3.0d0*rho_sweet/rho_salt
      write (10,*)  
      write (10,*) 'Initial values for the 6 floater dofs'
      do i = 1, 6
      write (10,100)   floater(nf)% UflInit(i)
      enddo
      write (10,*)  
      write (10,*) "Enabled or disabled floater's dofs"
      write (10,*)     floater(nf)% IQfl_active_el
!     write (10,*) 'External Forces'
!     write (10,*)  AQ_external
   enddo !nf


   close(24)

 100 format(100en16.6e3)


 END Subroutine floater_init
!----------------------------------------------------------------------
 Subroutine read_matrix (Matrix, NI, NJ, ITYPE, nfil)
!----------------------------------------------------------------------

   implicit none

   integer :: NI, NJ
   real(8) :: Matrix (NI,NJ)
   integer :: ITYPE, nfil, i, j


   Matrix(:,:) = 0.d0;

   goto (1,2,3), ITYPE

!----- diagonal
 1   if (NI/=NJ) then
        write(*,*)'error in diagonal read_matrix',NI,NJ
        stop
     endif
     do i=1,NI
         read(nfil,*)  Matrix (i,i)
     enddo
   return

!----- symmetric
 2   if (NI/=NJ) then
        write(*,*)'error in symmetric read_matrix',NI,NJ
        stop
     endif
     do i=1,NI
      do j=i,NJ
         read(nfil,*)  Matrix (i,j)
                       Matrix (j,i) = Matrix (i,j)
      enddo
     enddo
   return

!----- full matrix
 3   do i=1,NI
         read(nfil,*) (Matrix (i,j),j=1,NJ)
     enddo
   return


 END Subroutine read_matrix
!----------------------------------------------------------------------
 Subroutine diffraction_init
!----------------------------------------------------------------------

 Use Hydro
 Use Paths

   implicit none

   real(8), allocatable :: inp_frequency   (  :)
   real(8), allocatable :: inp_force_meas  (:,:)
   real(8), allocatable :: inp_force_phas  (:,:)
   real(8), allocatable :: inp_force_dft   (:,:)
   real(8), allocatable :: fy              (  :)
   real(8), allocatable :: d2fy            (  :)
   real(8)              :: f(6),p(6),dr(6)
   real(8)              :: omega, fy_sp
   integer              :: k, no_diff_freq_inp, i,j, ifr, nf, idr
!--for yaw
!! real(8)              :: x_yaw(5),y_yaw(5),ang,y1


                     idr=0
   if (NewmanMod==1) idr=6

   open (24 ,file=trim(file_diffract) ) !'diffract.inp'

      read (24,*) no_diff_freq_inp
                 write(28,'(a)')'**diffract.inp'
                 write(28,*) no_diff_freq_inp

      Allocate ( inp_frequency    (no_diff_freq_inp    ) )
      Allocate ( inp_force_meas   (no_diff_freq_inp,6  ) )
      Allocate ( inp_force_phas   (no_diff_freq_inp,6  ) )
      Allocate ( inp_force_dft    (no_diff_freq_inp,idr) )

      Allocate ( fy               (no_diff_freq_inp    ) )
      Allocate ( d2fy             (no_diff_freq_inp    ) )


   do nf = 1, nfloater_hyd
      Allocate ( floater(nf)%wave_force         (no_spectral_comp,6  ) )
      Allocate ( floater(nf)%wave_phase         (no_spectral_comp,6  ) )
      Allocate ( floater(nf)%wave_force_Dft_pos (no_spectral_comp,idr) )
      Allocate ( floater(nf)%wave_force_Dft_neg (no_spectral_comp,idr) )

      if (nf>1) read(24,*) !leaves a blank line between next and previous floater

      do k = 1, no_diff_freq_inp

         read(24,*) inp_frequency(k),(f(i),i=1,6),(p(i),i=1,6),(dr(i),i=1,idr) ![rad/s], [N or Nm/rho/g/A], [rad], [N or Nm/rho/g/A^2]

         inp_force_meas (k,1:6  ) = f (1:6  )*rho*gravity
         inp_force_phas (k,1:6  ) = p (1:6  )
         inp_force_dft  (k,1:idr) = dr(1:idr)*rho*gravity
                 write(28,'(150e16.6)')     &
                   inp_frequency(k)/scal_hyd**0.5d0   ,&
                  (f (i)           *scal_hyd**2,i=1,3),&
                  (f (i)           *scal_hyd**3,i=4,6),&
                  (p (i)                       ,i=1,6),&
                  (dr(i)                       ,i=1,idr)
      enddo !k

      do i=1,6
!--------- interpolate magnitude of 1st order exciting force
         fy(1:no_diff_freq_inp) = inp_force_meas  (1:no_diff_freq_inp,i)

         call spline1( inp_frequency, fy, d2fy, no_diff_freq_inp )

         do ifr=1, no_spectral_comp
            omega = spectral_freq(ifr)

            if (omega > inp_frequency(no_diff_freq_inp)) then
               floater(nf)%wave_force (ifr,i) = 0.d0
               cycle
            endif

            call splint1( inp_frequency, fy, d2fy, no_diff_freq_inp, omega, fy_sp )
            floater(nf)%wave_force (ifr,i) = fy_sp
         enddo !ifr


!--------- interpolate phase of 1st order exciting force
         fy(1:no_diff_freq_inp) = inp_force_phas  (1:no_diff_freq_inp,i)

         do ifr=1, no_spectral_comp
            omega = spectral_freq(ifr)
           !call lint1  ( xa           , ya, n1, n2              , np              , x    , y     )
            call lint1  ( inp_frequency, fy, 1 , no_diff_freq_inp, no_diff_freq_inp, omega, fy_sp )
            floater(nf)%wave_phase (ifr,i) = fy_sp
         enddo !ifr


!--------- interpolate drift force
         if (NewmanMod==0) cycle
        
         fy(1:no_diff_freq_inp) = inp_force_dft   (1:no_diff_freq_inp,i)

         call spline1( inp_frequency, fy, d2fy, no_diff_freq_inp )

         do ifr=1, no_spectral_comp
            omega = spectral_freq(ifr)

            if (omega > inp_frequency(no_diff_freq_inp)) then
               floater(nf)%wave_force_Dft_pos (ifr,i) = 0.d0
               floater(nf)%wave_force_Dft_neg (ifr,i) = 0.d0
               cycle
            endif

            call splint1( inp_frequency, fy, d2fy, no_diff_freq_inp, omega, fy_sp )

!------------ Split pos vs neg drift forces for Newman's approximation
            if (fy_sp>0.d0) then
               floater(nf)%wave_force_Dft_pos (ifr,i) = fy_sp
               floater(nf)%wave_force_Dft_neg (ifr,i) = 0.d0
            else
               floater(nf)%wave_force_Dft_pos (ifr,i) = 0.d0
               floater(nf)%wave_force_Dft_neg (ifr,i) = fy_sp
            endif
         enddo !ifr
      enddo !i


!------ Change wave force with angle β
!!    do ifr=1,no_spectral_comp
!!       floater(nf)%wave_phase(ifr,2) = floater(nf)%wave_phase(ifr,1)
!!       floater(nf)%wave_phase(ifr,4) = floater(nf)%wave_phase(ifr,5)
!!    enddo


!!    do ifr=1,no_spectral_comp
!!       floater(nf)%wave_force(ifr,2) = floater(nf)%wave_force(ifr,1) * dsin(spectral_angle)
!!       floater(nf)%wave_force(ifr,1) = floater(nf)%wave_force(ifr,1) * dcos(spectral_angle)
!!       floater(nf)%wave_force(ifr,4) =-floater(nf)%wave_force(ifr,5) * dsin(spectral_angle)
!!       floater(nf)%wave_force(ifr,5) = floater(nf)%wave_force(ifr,5) * dcos(spectral_angle)
!!    enddo


!------ Special for FLOTTEK ---
!------ adjast yaw
!!    ang = spectral_angle

!! 21 if (ang<0.d0) then
      
!!      ang=ang+pi2_h
!!      goto 21
      
!!    endif
      
!!----- check
!!    if (ang>pi2_h) then
!!       write(*,*) 'error in yaw init for diffraction'
!!       stop
!!    endif
      
!!    ang = dmod(ang, 90.d0*pi_h/180.d0)
      
!!    write( *,*)'ang 0:90',ang*180.d0/pi_h
!!    write(10,*)'ang 0:90',ang*180.d0/pi_h
      
!!    x_yaw(1) =  0.0d0
!!    x_yaw(2) = 22.5d0*pi_h/180.d0
!!    x_yaw(3) = 45.0d0*pi_h/180.d0
!!    x_yaw(4) = 67.5d0*pi_h/180.d0
!!    x_yaw(5) = 90.0d0*pi_h/180.d0
!!    y_yaw(1) =  0.0d0
!!    y_yaw(3) =  0.0d0
!!    y_yaw(5) =  0.0d0
      
!!    do ifr = 1, no_spectral_comp
!!       y_yaw(2) = floater(nf)%wave_force(ifr,6)
!!       y_yaw(4) =-floater(nf)%wave_force(ifr,6)
      
!!       call lint1 ( x_yaw, y_yaw, 1, 3, 3, ang, y1 )
      
!!       floater(nf)%wave_force(ifr,6) = y1
!!    enddo

!------ write out exciting forces
      if (nf==1) then
       open (204,file='diffract.out')
      
       write(204,101) spectral_angle * 180d0/pi_h
      endif


      do i = 1, no_spectral_comp
         write(204,101)                                               &
          spectral_freq(i)                                           ,&
          (floater(nf)%wave_force        (i,j)/(gravity*rho),j=1,6  ),&
          (floater(nf)%wave_phase        (i,j)              ,j=1,6  ),&
          (floater(nf)%wave_force_Dft_pos(i,j)/(gravity*rho),j=1,idr),&
          (floater(nf)%wave_force_Dft_neg(i,j)/(gravity*rho),j=1,idr)
      enddo !i

      write(204,*)
      write(204,*)
   enddo !nf

   close ( 24)
   close (204)

   Deallocate ( inp_force_meas, inp_force_phas, inp_frequency )
   Deallocate ( fy            , d2fy                          )
   if (NewmanMod==1)&
   Deallocate ( inp_force_dft )

 101 format ( 150f23.9 )


 END Subroutine diffraction_init
!----------------------------------------------------------------------
 Subroutine retardation_init
!----------------------------------------------------------------------

 Use Hydro
 Use Paths

   implicit none

   integer              :: k, num_sigs_inp, i
   real(8), allocatable :: retd_time_inp   (  :), retd_func_inp(:,:)
   real(8), allocatable :: rfy             (  :), d2rfy        (  :)
   real(8), allocatable :: retard_func_tmp (:,:)
   real(8)              :: time_sp, rfy_sp
   integer              :: nf, j, N
   real(8)              :: ex(6)


   N = nfloater_hyd * 6


   open (24,file=trim(file_retard)) !'retard.inp'
   open (25,file='retard.out')


      read(24,*) num_sigs_inp
                 write(28,'(a)')'**retard.inp'
                 write(28,*) num_sigs_inp

                 ex(1:3)=0.d0; ex(4:6)=1.d0;

      Allocate ( retd_time_inp (num_sigs_inp                ) )
      Allocate ( retd_func_inp (num_sigs_inp,36*nfloater_hyd) )


   do nf = 1, nfloater_hyd

      if (nf>1) read(24,*) !leaves a blank line between next and previous floater

      do k=1,num_sigs_inp

         read(24,*) retd_time_inp(k)                       , &  !-1                    --> time
                   (retd_func_inp(k,j),j=1,36*nfloater_hyd)     !-2:37*nfloater_hyd+1  --> retardation functions

                 write(28,'(150e16.6)') retd_time_inp(k) * scal_hyd**0.5d0                               ,&
                  ((retd_func_inp(k,(i-1)*6+j)*scal_hyd**(2.5+ex(i)+ex(j))*rho_sweet/rho_salt,j=1,6),i=1,6)
      enddo

      retardation_length = int((retd_time_inp(num_sigs_inp)-retd_time_inp(1))/delta_t) - 1

      Allocate ( floater(nf)% velocity_buffer      (6   , retardation_length-1) )
      Allocate ( floater(nf)% retardation_function (6, N, retardation_length  ) )


!------ Initialise retardation_function, velocity_buffer
      floater(nf)%retardation_function = 0.d0;
      floater(nf)%velocity_buffer      = 0.d0;

      Allocate ( rfy                 (num_sigs_inp                       ) )
      Allocate ( d2rfy               (num_sigs_inp                       ) )
      Allocate ( retard_func_tmp     (retardation_length,36*nfloater_hyd ) )


      write(* ,*)'retardation_length=',retardation_length
      write(10,*)'retardation_length=',retardation_length


      do i    = 1, 36*nfloater_hyd
         do k = 1, num_sigs_inp
            rfy(k)=retd_func_inp(k,i)
         enddo

         call spline1 ( retd_time_inp, rfy, d2rfy, num_sigs_inp )

         do k = 1, retardation_length
            time_sp = delta_t*(k-1)
            call splint1( retd_time_inp, rfy, d2rfy, num_sigs_inp, time_sp, rfy_sp )
            retard_func_tmp(k,i) = rfy_sp
         enddo
      enddo

      do i = 1, 6
         j = (i-1)*N
         do k=1,retardation_length
            floater(nf)%retardation_function(i,1:N,k) = retard_func_tmp(k, j+1:j+N)
         enddo
      enddo


      do i=1,retardation_length
         write(25,101)                &
          delta_t*(i-1)              ,&
         (retard_func_tmp(i,j),j=1,36*nfloater_hyd)
      enddo

      Deallocate ( rfy, d2rfy, retard_func_tmp )
   enddo !nf

   close(24)
   close(25)

   Deallocate ( retd_time_inp, retd_func_inp   )

101 format ( 1500f23.9 )


 END Subroutine retardation_init
!----------------------------------------------------------------------
 Subroutine retardation_init_5
!----------------------------------------------------------------------

 Use Hydro
 Use Paths

   implicit none

   integer              :: k, num_sigs_inp, i
   real(8), allocatable :: retd_time_inp   (  :), retd_func_inp(:,:)
   real(8), allocatable :: rfy             (  :), d2rfy        (  :)
   real(8), allocatable :: retard_func_tmp (:,:)
   real(8)              :: time_sp, rfy_sp
   integer              :: nf, nfj, j1, j2, j, N
   real(8)              :: z !zero


   N = nfloater_hyd * 6


   open (24,file=trim(file_retard)) !'retard.inp'
   open (25,file='retard.out')


      read(24,*) num_sigs_inp

      Allocate ( retd_time_inp (num_sigs_inp               ) )
      Allocate ( retd_func_inp (num_sigs_inp,5*nfloater_hyd) )


   do nf = 1, nfloater_hyd

      if (nf>1) read(24,*) !leaves a blank line between next and previous floater

      do k=1,num_sigs_inp

         read(24,*) retd_time_inp(k)                      , &  !1--> time
                   (retd_func_inp(k,j),j=1,5*nfloater_hyd)
                  ! retd_func_inp(k,1)                    , &  !2--> Fx
                  ! retd_func_inp(k,2)                    , &  !3--> Fz
                  ! retd_func_inp(k,3)                    , &  !4--> M5
                  ! retd_func_inp(k,4)                    , &  !5--> M51
                  ! retd_func_inp(k,5)                         !6--> M6
      enddo

      retardation_length = int((retd_time_inp(num_sigs_inp)-retd_time_inp(1))/delta_t) - 1

      Allocate ( floater(nf)% velocity_buffer      (6   , retardation_length-1) )
      Allocate ( floater(nf)% retardation_function (6, N, retardation_length  ) )


!------ Initialise retardation_function, velocity_buffer
      floater(nf)%retardation_function = 0.d0;
      floater(nf)%velocity_buffer      = 0.d0;

      Allocate ( rfy                 (num_sigs_inp                      ) )
      Allocate ( d2rfy               (num_sigs_inp                      ) )
      Allocate ( retard_func_tmp     (retardation_length,5*nfloater_hyd ) )


      write(* ,*)'retardation_length=',retardation_length
      write(10,*)'retardation_length=',retardation_length


      do i    = 1, 5*nfloater_hyd
         do k = 1, num_sigs_inp
            rfy(k)=retd_func_inp(k,i)
         enddo

         call spline1 ( retd_time_inp, rfy, d2rfy, num_sigs_inp )

         do k = 1, retardation_length
            time_sp = delta_t*(k-1)
            call splint1( retd_time_inp, rfy, d2rfy, num_sigs_inp, time_sp, rfy_sp )
            retard_func_tmp(k,i) = rfy_sp
         enddo
      enddo

      do nfj = 1, nfloater_hyd
          j1 = (nfj-1)*6
          j2 = (nfj-1)*5
       do i  = 1, retardation_length
          floater(nf)%retardation_function(1,j1+1,i) = retard_func_tmp(i,j2+1)
          floater(nf)%retardation_function(1,j1+5,i) = retard_func_tmp(i,j2+4)
          floater(nf)%retardation_function(2,j1+2,i) = retard_func_tmp(i,j2+1)
          floater(nf)%retardation_function(2,j1+4,i) =-retard_func_tmp(i,j2+4)
          floater(nf)%retardation_function(3,j1+3,i) = retard_func_tmp(i,j2+2)
          floater(nf)%retardation_function(4,j1+2,i) =-retard_func_tmp(i,j2+4)
          floater(nf)%retardation_function(4,j1+4,i) = retard_func_tmp(i,j2+3)
          floater(nf)%retardation_function(5,j1+1,i) = retard_func_tmp(i,j2+4)
          floater(nf)%retardation_function(5,j1+5,i) = retard_func_tmp(i,j2+3)
          floater(nf)%retardation_function(6,j1+6,i) = retard_func_tmp(i,j2+5)
       enddo
      enddo

!------- new write ------------------------
          z=0.d0
          write(122,*  ) num_sigs_inp
       do k  = 1, num_sigs_inp
          write(122,102)                                         &
           retd_time_inp(k)                                     ,&
                retd_func_inp(k,1),z,z,z, retd_func_inp(k,4),z  ,&
             z, retd_func_inp(k,1),z,    -retd_func_inp(k,4),z,z,&
           z,z, retd_func_inp(k,2),                        z,z,z,&
             z,-retd_func_inp(k,4),z,     retd_func_inp(k,3),z,z,&
                retd_func_inp(k,4),z,z,z, retd_func_inp(k,3),z  ,&
           z,z,                    z,z,z, retd_func_inp(k,5)
       enddo
       write(*,*) 'WRITE retard OK'
       stop
102 format(150en16.8e3)
!------- end new write --------------------

      do i=1,retardation_length
         write(25,101)                             &
          delta_t*(i-1)                           ,&
         (retard_func_tmp(i,j),j=1,5*nfloater_hyd)
        ! retard_func_tmp(i,1),&
        ! retard_func_tmp(i,2),&
        ! retard_func_tmp(i,3),&
        ! retard_func_tmp(i,4),&
        ! retard_func_tmp(i,5)
      enddo

      Deallocate ( rfy, d2rfy, retard_func_tmp )
   enddo !nf

   close(24)
   close(25)

   Deallocate ( retd_time_inp, retd_func_inp   )

101 format ( 150f23.9 )


 END Subroutine retardation_init_5
!----------------------------------------------------------------------
 Subroutine WRITEOUT_floater
!----------------------------------------------------------------------

 Use Hydro
 Use Cbeam

   implicit none

   real(8) :: QWT(6), Qrad(6), Qm_hyd(6), Qm_str(6), Qmoor(6), Qdamp(6), Qk_hyd(6), Qk_str(6)
   integer :: i, nf, N


   N  = nfloater_hyd * 6
   do nf = 1, nfloater_hyd
    QWT = 0.d0
    if (nf==1) call  CALC_FLOATLOADS ( QWT )

    QWT   (1:6) = -QWT(1:6)
    Qrad  (1:6) = -matmul( floater(nf)% retardation_function  (1:6,1:N,1), UT1_el(NDFBT_el+1:NDFBT_el+N) ) - floater(nf)%AQ_retardation(1:6)
    Qm_hyd(1:6) = -matmul( floater(nf)% AM_added_infinity     (1:6,1:N  ), UT2_el(NDFBT_el+1:NDFBT_el+N) )
    Qm_str(1:6) = -matmul( floater(nf)% AM_floater_structural (1:6,1:6  ), UT2_el(NDFBT_el+1:NDFBT_el+6) )
    Qk_hyd(1:6) = -matmul( floater(nf)% AK_hydrostatic        (1:6,1:6  ), UT_el (NDFBT_el+1:NDFBT_el+6) )
    Qk_str(1:6) = -matmul( floater(nf)% AK_floater_structural (1:6,1:6  ), UT_el (NDFBT_el+1:NDFBT_el+6) )


!qq --- check Q_truss when ITRUSS_el = 2
!------ Moorings Forcing
    if  ( ITRUSS_el == 2 .and. nf == 1)  &
    call truss2float_loads_uncoupl         !calculate Q_truss

    Qmoor (1:6) =-  matmul( floater(nf)% AK_external          (1:6,1:6  ), UT_el (NDFBT_el+1:NDFBT_el+6) )  &
                 +          floater(nf)% Q_truss              (1:6      )

    if  ( ITRUSS_el == 2 .and. nf == 1)  floater(nf)% Q_truss (1:6      ) = 0.d0;

!---- Damping Forcing
    if      (ITYPE_Damp_ext==1) then
       Qdamp (1:6) = -matmul( floater(nf)%AC_external       (1:6,1:6  ), UT1_el(NDFBT_el+1:NDFBT_el+6) )
    elseif  (ITYPE_Damp_ext==2) then
       Qdamp (1:6) =          floater(nf)%AQ_damp_quad      (1:6      )
    else
       Qdamp (1:6) =  0.d0;
    endif

!--- floater's loading
#ifndef ASCII
#define ASCII 1
#endif

#if   ASCII == 0
    if (nf==1)open (15,file='Floater_force.bin'   , access='append',form='UNFORMATTED')
    if (nf> 1)open (15,file='Floater_force_02.bin', access='append',form='UNFORMATTED')
   write(15)                                     &
#elif ASCII == 1
    if (nf==1)open (15,file='Floater_force.dat'   , access='append')
    if (nf> 1)open (15,file='Floater_force_02.dat', access='append')
   write(15,100)                                 &
#endif
     sngl( TIME                         )       ,&  !time                           ! 1
    (sngl( floater(nf)%AQ_exciting  (i) ),i=1,6),&  !Diffraction            loading ! 2- 7
    (sngl( floater(nf)%AQ_Newman    (i) ),i=1,6),&  !Newman's approximation loading ! 8-13
    (sngl( Qrad                     (i) ),i=1,6),&  !Radiation              loading !14-19
    (sngl( floater(nf)%AQ_morison   (i) ),i=1,6),&  !Morison's eq.          loading !20-25
    (sngl( Qmoor                    (i) ),i=1,6),&  !moorings               loading !26-31
    (sngl( QWT                      (i) ),i=1,6),&  !WT's transfered        loading !32-37
!---- next quantities could be found after the simulation with post process (matrices are constant)
    (sngl( Qdamp                    (i) ),i=1,6),&  !external    damping    loading !38-43
    (sngl( Qk_str                   (i) ),i=1,6),&  !structural  stiffness  loading !44-49
    (sngl( Qk_hyd                   (i) ),i=1,6),&  !hydrostatic stiffness  loading !50-55
    (sngl( Qm_str                   (i) ),i=1,6),&  !structural  inertial   loading !56-61
    (sngl( Qm_hyd                   (i) ),i=1,6)    !hydrostatic add mass   loading !62-67

   close(15)

   enddo !nf

 100  format (150f15.5)


 END Subroutine WRITEOUT_floater
!----------------------------------------------------------------------
 Subroutine WRITEOUT_zita
!----------------------------------------------------------------------

 Use Hydro
 Use Cbeam

   implicit none

   real(8) :: Zita, Zita0, PSI0, PSIB0, X(3), X0(3)


!--- Write time, z_elevation, azimuth angle

   if (ICASE_el == 0) return


    X0 (1:3)     = 0.d0;
    call ZELEV_sea ( X0, Zita0 )

   if (ICASE_el >= 3) then
      X  (1:3)   = UT_el(NDFBT_el+1: NDFBT_el+3);
      call ZELEV_sea ( X , Zita  )
   else
      Zita       = Zita0
   endif

   PSI0          = -UT_el(NDFBT_el+NQSW)
   PSIB0         = dmod (PSI0,PI2)*R2D

#if   ASCII == 0
   open (37,file='z_elevation.bin',access='append',form='UNFORMATTED')
   write(37)                                     &
#elif ASCII == 1
   open (37,file='z_elevation.dat',access='append')
   write(37,100)                                 &
#endif
     sngl(TIME), sngl(Zita0), sngl(Zita), sngl(PSIB0)

   close(37)

 100  format (150f15.5)


 END Subroutine WRITEOUT_zita
!--------------------------------------------------------------------------------
!
!  Subrouine : FLOATLOADS_WT26QS ----------------
!
!  Loads Communication between WT and 6 qs [ICASE_el == 3 or 4]
!
!----------------------------------------------------------------------
 Subroutine FLOATLOADS_WT26Qs
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: Aba      (3,3        )
   real(8) :: Aba0     (3,3,NQS    )
   real(8) :: Rba      (3          )
   real(8) :: Rba0     (3  ,NQS    )
   real(8) :: Ra       (3          )
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )
   real(8) :: AFB2B    (6          )
   integer :: nbod1, nbod2, i, j, NDFPE, iq, nq, iqi
   integer :: nbod, nnsub, nbsub, nn_el, nb_el, nelB, nodB


   if ((ICASE_el /= 3).and.(ICASE_el /= 4)) return


   if (NBODBTWT_el == NBLADE_el) then
    nbod1 = 1
    nbod2 = NBLADE_el
   else
    nbod1 = NBODBTWT_el!shaft or tower
    nbod2 = NBODBTWT_el
   endif


   do nbod  = nbod1, nbod2
!------ Body (B) WT or Jacket
      nnsub = 1
      nbsub = 1

      nn_el = body(nbod)%NNODTGLB_el(nnsub)
      nb_el = body(nbod)%NBODTGLB_el(nbsub)
      nelB  = 1
      nodB  = 1

      NDFPE = subbody(nb_el)%NDFPE_el


!------ Body (A) Floater: Global c.o. System
!------ Aba, Aba0, Rba, Rba0
      Ra  (1:3    ) = UT_el(NDFBT_el+1 : NDFBT_el+3)
      Rba (1:3    ) = transf_mat(nn_el)%R_el (1:3    ) - Ra(1:3)  !transf_mat(nn_el)%A_el (1:3,1:3) x (0,0,Htow0)
      Aba (1:3,1:3) = transf_mat(nn_el)%A_el (1:3,1:3)

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)

         Aba0 (1:3,1:3,nq) = transf_mat(nn_el)%A0_el(iq ,1:3,1:3) !(nbsub,nq,i,j)
         Rba0 (1:3,    nq) = 0.d0
      enddo!iq 
         Rba0 (1:3,   4:6) = transf_mat(nn_el)%R0_el(4:6,1:3    ) !(nbsub,nq,i,j)   !! A0 . (0,0,Htow)??


      call LOC_LOADS_el        ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                 AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                 nb_el  , nn_el  , nelB   , nodB    )


      call LOADS_TRANS_el      ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                 AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                 Aba0   , Aba    , Rba0   , Rba  ,&
                                          nn_el  , 1             ,& !IcalcR
                                 NDFPE  , NDFPEM_el                 )


!------ Assembly to global Matrices
      j    = ACCU%NDFTACC_el(nb_el-1)!1:NDFPE !+subbody(nbB_el)%NDFPNBACC_el(nnB-1)+mdof


      do iq = 1, IQfl
         i  = NDFBT_el + iq 

         AM_el ( i, j+1:j+NDFPE ) = AM_el ( i, j+1:j+NDFPE ) + AMB2B ( IMDOF_el(iq), 1:NDFPE )
         AC_el ( i, j+1:j+NDFPE ) = AC_el ( i, j+1:j+NDFPE ) + ACB2B ( IMDOF_el(iq), 1:NDFPE )
         AK_el ( i, j+1:j+NDFPE ) = AK_el ( i, j+1:j+NDFPE ) + AKB2B ( IMDOF_el(iq), 1:NDFPE )
         AQ_el ( i              ) = AQ_el ( i              ) + AFB2B ( IMDOF_el(iq)          )
      enddo


!------ For all q's
      do iqi = 1, IQfl
         i   = NDFBT_el + iqi

         do iq = 1, QSnode(nn_el)%Tot_el(2)
            nq = QSnode(nn_el)%IORD_el(iq)
            j  = NDFBT_el + nq

            AM_el ( i, j ) = AM_el ( i, j ) + AMB2BNQ ( IMDOF_el(iqi), nq )
            AC_el ( i, j ) = AC_el ( i, j ) + ACB2BNQ ( IMDOF_el(iqi), nq )
            AK_el ( i, j ) = AK_el ( i, j ) + AKB2BNQ ( IMDOF_el(iqi), nq )
         enddo!iqj 
      enddo!iqi

   enddo!nbod


 END Subroutine FLOATLOADS_WT26Qs
!--------------------------------------------------------------------------------
!
!  Subrouine : FLOATLOADS  ----------------
!
!  Loads Communication between flexible Jacket and the 6 q's [ICASE_el == 4]
!
!----------------------------------------------------------------------
 Subroutine FLOATLOADS_JA26Qs
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: Aba      (3,3        )
   real(8) :: Aba0     (3,3,NQS    )
   real(8) :: Rba      (3          )
   real(8) :: Rba0     (3  ,NQS    )
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )
   real(8) :: AFB2B    (6          )
   real(8) :: HTA
   integer :: nbod1, nbod2, i, j, NDFPE, iq, nq, iqi
   integer :: nbod, nbc, nnsub, nbsub, nn_el, nb_el, nelB, nodB, ibcon


   if (ICASE_el /= 4) return

   nbod1  = NBODBTWT_el+1
   nbod2  = NBODBT_el


!--- Body(B), part of the flexible floater
   do nbod  = nbod1, nbod2
      nnsub = 1
      nbsub = 1
      nn_el = body(nbod)%NNODTGLB_el(nnsub)
      nb_el = body(nbod)%NBODTGLB_el(nbsub)
      NDFPE = subbody(nb_el)%NDFPE_el


      do nbc   = 1, boundco (nb_el)%nbcpb
         ibcon = boundco (nb_el)%nbcon (nbc)

         if (ibcon /= 0) cycle

         nelB  = boundco (nb_el)%nel (nbc)
         nodB  = boundco (nb_el)%nod (nbc)
         i     = nodB/NNPE_el
         HTA   = subbody (nb_el)%HTA_el(nelB+i)


!--------- Body (A) Floater: Global c.o. System
!--------- Aba, Aba0, Rba, Rba0
         Rba (1:3    ) = matmul( A_float(1:3,1:3), transf_mat(nn_el)%Rja_el(1:3  )      +&
                                                   transf_mat(nn_el)%Aja_el(1:3,2)* HTA   )
         Aba (1:3,1:3) = transf_mat(nn_el)%A_el (1:3,1:3)

         do iq = 1, QSnode(nn_el)%Tot_el(2)
            nq = QSnode(nn_el)%IORD_el(iq)

            Aba0 (1:3,1:3,nq) = transf_mat(nn_el)%A0_el(iq ,1:3,1:3) !(nbsub,nq,i,j)
            Rba0 (1:3,    nq) = 0.d0
         enddo!iq 

!--------- Rba0 only because the floater rotates
         do i  = 1, 3
            nq = i+3
            Rba0 (1:3,    nq) = matmul( A0_float(1:3,1:3,i), transf_mat(nn_el)%Rja_el(1:3  )      +&
                                                             transf_mat(nn_el)%Aja_el(1:3,2)* HTA   )
         enddo


         call LOC_LOADS_el        ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                    AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                    nb_el  , nn_el  , nelB   , nodB    )


         call LOADS_TRANS_el      ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                    AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                    Aba0   , Aba    , Rba0   , Rba  ,&
                                             nn_el  , 1             ,& !IcalcR
                                    NDFPE  , NDFPEM_el                 )


!--------- Assembly to global Matrices
         j    = ACCU%NDFTACC_el(nb_el-1)!1:NDFPE !+subbody(nbB_el)%NDFPNBACC_el(nnB-1)+mdof


         do iq = 1, IQfl
            i  = NDFBT_el + iq 

            AM_el ( i, j+1:j+NDFPE ) = AM_el ( i, j+1:j+NDFPE ) + AMB2B ( IMDOF_el(iq), 1:NDFPE )
            AC_el ( i, j+1:j+NDFPE ) = AC_el ( i, j+1:j+NDFPE ) + ACB2B ( IMDOF_el(iq), 1:NDFPE )
            AK_el ( i, j+1:j+NDFPE ) = AK_el ( i, j+1:j+NDFPE ) + AKB2B ( IMDOF_el(iq), 1:NDFPE )
            AQ_el ( i              ) = AQ_el ( i              ) + AFB2B ( IMDOF_el(iq)          )
         enddo


!--------- For all q's
         do iqi = 1, IQfl
            i   = NDFBT_el + iqi

            do iq = 1, QSnode(nn_el)%Tot_el(2)
               nq = QSnode(nn_el)%IORD_el(iq)
               j  = NDFBT_el + nq

               AM_el ( i, j ) = AM_el ( i, j ) + AMB2BNQ ( IMDOF_el(iqi), nq )
               AC_el ( i, j ) = AC_el ( i, j ) + ACB2BNQ ( IMDOF_el(iqi), nq )
               AK_el ( i, j ) = AK_el ( i, j ) + AKB2BNQ ( IMDOF_el(iqi), nq )
            enddo!iqj 
         enddo!iqi

      enddo!nbc
   enddo!nbod


 END Subroutine FLOATLOADS_JA26Qs
!--------------------------------------------------------------------------------
!
!  Subrouine : CALC_FLOATLOADS  ----------------
!
!  Loads Communication between bodies (last Subbody with 1st of the next body)
!
!----------------------------------------------------------------------
 Subroutine CALC_FLOATLOADS ( REACTF )
!----------------------------------------------------------------------

 Use Cbeam

   implicit none

   real(8) :: Aba      (3,3        )
   real(8) :: Aba0     (3,3,NQS    )
   real(8) :: Rba      (3          )
   real(8) :: Rba0     (3  ,NQS    )
   real(8) :: Ra       (3          )
   real(8) :: AMB2B    (6,NDFPEM_el)
   real(8) :: ACB2B    (6,NDFPEM_el)
   real(8) :: AKB2B    (6,NDFPEM_el)
   real(8) :: AMB2BNQ  (6,NQS      )
   real(8) :: ACB2BNQ  (6,NQS      )
   real(8) :: AKB2BNQ  (6,NQS      )
   real(8) :: AFB2B    (6          )
   real(8) :: REACTF   (6          )
   integer :: nbod1, nbod2, j, NDFPE, iq, nq
   integer :: nbod, nnsub, nbsub, nn_el, nb_el, nel


   REACTF(1:6) = 0.d0;

   if (NBODBTWT_el.eq.NBLADE_el) then

       nbod1 = 1
       nbod2 = NBLADE_el

   else

       nbod1 = NBODBTWT_el!shaft or tower
       nbod2 = NBODBTWT_el

   endif


   do nbod = nbod1, nbod2
 
!------ Body (B)
      nnsub = 1
      nbsub = 1
      nel   = 1

      nn_el = body(nbod)%NNODTGLB_el(nnsub)
      nb_el = body(nbod)%NBODTGLB_el(nbsub)

      NDFPE  = subbody(nb_el)%NDFPE_el


!------ Body (A) Floater: Global c.o. System

!------ Aba, Aba0, Rba, Rba0
      Ra  (1:3    ) = UT_el(NDFBT_el+1 : NDFBT_el+3)
      Rba (1:3    ) = transf_mat(nn_el)%R_el (1:3    ) - Ra(1:3)
      Aba (1:3,1:3) = transf_mat(nn_el)%A_el (1:3,1:3)

      do iq = 1, QSnode(nn_el)%Tot_el(2)
         nq = QSnode(nn_el)%IORD_el(iq)

         Aba0 (1:3,1:3,nq) = transf_mat(nn_el)%A0_el(iq,1:3,1:3) !(nbsub,nq,i,j)
         Rba0 (1:3,    nq) = 0.d0
      enddo!iq 
         Rba0 (1:3,   4:6) = transf_mat(nn_el)%R0_el(4:6,1:3    ) !(nbsub,nq,i,j)   !! A0 . (0,0,Htow)??


      call LOC_LOADS_el        ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                 AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                 nb_el  , nn_el  , nel    , 1       )!nodB


      call LOADS_TRANS_el      ( AMB2B  , ACB2B  , AKB2B  , AFB2B,&
                                 AMB2BNQ, ACB2BNQ, AKB2BNQ       ,&
                                 Aba0   , Aba    , Rba0   , Rba  ,&
                                          nn_el  , 1             ,& !IcalcR
                                 NDFPE  , NDFPEM_el                 )


      j    = ACCU%NDFTACC_el(nb_el-1)

     !REACTF(1:6) = REACTF(1:6) + matmul ( AMB2B(1:6,1:NDFPE), UT02_el(j+1:j+NDFPE) )
     !REACTF(1:6) = REACTF(1:6) + matmul ( ACB2B(1:6,1:NDFPE), UT01_el(j+1:j+NDFPE) )
     !REACTF(1:6) = REACTF(1:6) + matmul ( AKB2B(1:6,1:NDFPE), UT0_el (j+1:j+NDFPE) )
      REACTF(1:6) = REACTF(1:6) -          AFB2B(1:6        )


     !do iq = 1, QSnode(nn_el)%Tot_el(2)

     !   nq = QSnode(nn_el)%IORD_el(iq)
     !   j  = NDFBT_el + nq

     !   REACTF(1:6) = REACTF(1:6) + AMB2BNQ(1:6,nq) * UT02_el(j)
     !   REACTF(1:6) = REACTF(1:6) + ACB2BNQ(1:6,nq) * UT01_el(j)
     !   REACTF(1:6) = REACTF(1:6) + AKB2BNQ(1:6,nq) * UT0_el (j)

     !enddo!iq 


   enddo!nbod


 END Subroutine CALC_FLOATLOADS
!
!
!
!----------------------------------------------------------------------
 Module Floater_Morison
!----------------------------------------------------------------------
!-- Nodes
 integer,              save :: NNT_mo
 real(8), allocatable, save :: XN_mo     (:,:  )   !(3,NNT_mo)

!-- Slender bodies
 integer,              save :: NBODBT_mo
 real(8), allocatable, save :: D_mo      (:,:  )   !(2,NBODBT_mo)
 real(8), allocatable, save :: Cd_mo     (:,:  )   !(2,NBODBT_mo)
 integer, allocatable, save :: Nod_mo    (:,:  )   !(2,NBODBT_mo)
 integer, allocatable, save :: NZ_mo     (:    )   !(  NBODBT_mo)

!-- Plates
 integer,              save :: Nplates_mo
 real(8), allocatable, save :: D_pl_mo   (:    )   !(  Nplates_mo)
 real(8), allocatable, save :: Cd_pl_mo  (:    )   !(  Nplates_mo)
 integer, allocatable, save :: Nod_pl_mo (:,:  )   !(2,Nplates_mo)


 END Module Floater_Morison
!----------------------------------------------------------------------
 Subroutine Morison_floater
!----------------------------------------------------------------------

 Use Hydro
 Use Cbeam
 Use Floater_Morison

   implicit none

   real(8) :: XG (3), UG(3), UINF(3), Ur(3), Ur_n(3)
   real(8) :: FDrag(3), MDrag(3), Ey(3), RG(3), RL(3)
   real(8) :: XG_nod(3,2), Ey_unit(3), Dy(3), Ur_t(3)
   real(8) :: AINF(3), pdyn
   real(8) :: DL, Ur_n_Meas, CD, CDrag, Diam, C, C1, C2
   real(8) :: Z_calc, Length, Length0, lamda, Ac, Ur_t_Meas, ACL(6,6), LLt(3,3), AImLLT(3,3)
   integer :: NZ, nb_mo, npl_mo, iz, m1, m2, N_start, N_end, Ihalfwet, iq,i,j
   integer :: nf


   do nf = 1, nfloater_hyd
      floater(nf)%AQ_morison (1:6    ) = 0.d0;
      floater(nf)%AC_morison (1:6,1:6) = 0.d0;
      floater(nf)%AK_morison (1:6,1:6) = 0.d0;
   enddo

!qq Valid for morison only applied to 1st floater
   nf = 1

   if (DragMod == 0) return

   Z_calc = tide

   do nb_mo = 1, NBODBT_mo

      m1            = Nod_mo(1,nb_mo)
      m2            = Nod_mo(2,nb_mo)

      XG_nod(1:3,1) = R_float(1:3) + matmul(A_float(1:3,1:3),XN_mo(1:3,m1))
      XG_nod(1:3,2) = R_float(1:3) + matmul(A_float(1:3,1:3),XN_mo(1:3,m2))

!------ check which nod is lower
      if (XG_nod(3,2)>=XG_nod(3,1)) then
         N_start=1
         N_end  =2
      else
         N_start=2
         N_end  =1
      endif

      if (XG_nod(3,N_start)>Z_calc) cycle        !both out of water
                                    Ihalfwet = 0 !indicates if the body is half wet half dry
      if (XG_nod(3,N_end  )>Z_calc) Ihalfwet = 1
      
      Ey(1:3)  = XG_nod(1:3,N_end) - XG_nod(1:3,N_start)  !vector tangential
      Length   = dsqrt(dot_product (Ey,Ey))
      Length0  = Length

      if (Ihalfwet==1) then
!--------- X-X1=λ(X2-X1) ---> λ = (X(3)-X1(3)) / (X2(3)-X1(3)) 
!---------               ---> X(1:2) = X1(1:2) + λ (X2(1:2)-X1(1:2))
         lamda             = (Z_calc - XG_nod(3,N_start))/Ey(3)
         XG_nod(1:2,N_end) = XG_nod(1:2,N_start) + lamda*Ey(1:2)
         XG_nod(  3,N_end) = Z_calc

         Ey(1:3)           = XG_nod(1:3,N_end) - XG_nod(1:3,N_start)
         Length            = dsqrt(dot_product (Ey,Ey))
      endif

      NZ           = NZ_mo(nb_mo)
      Ey_unit(1:3) = Ey(1:3)/Length                       !unit tangential vector
      Dy(1:3)      = Ey(1:3)/dble(NZ)                     !NZ because calculations are at the middle of each segment
      DL           = Length /dble(NZ)


      do iz   = 1, NZ
         C    = dble(2*iz-1)/2.d0
         C2   = C * DL / Length0
         C1   = 1.d0 - C2

         XG (1:3) = XG_nod   (1:3,N_start)    + Dy     (1:3) * C
         RG (1:3) = XG       (1:3)            - R_float(1:3)
         RL (1:3) = matmul ( AT_float(1:3,1:3),  RG    (1:3) )
         UG (1:3) = matmul ( DA_float(1:3,1:3),  RL    (1:3) ) + DR_float(1:3)

         call UINFLOW_sea ( XG, UINF, AINF, pdyn, 2 )

         Ur        (1:3) = UINF(1:3) - UG(1:3)
         Ur_n      (1:3) = Ur(1:3) - dot_product(Ur(1:3),Ey_unit(1:3)) * Ey_unit(1:3)
         Ur_n_Meas       = dsqrt ( dot_product ( Ur_n, Ur_n ) )

         CDrag           = C1*Cd_mo(1,nb_mo) + C2*Cd_mo(2,nb_mo)
         Diam            = C1* D_mo(1,nb_mo) + C2* D_mo(2,nb_mo)
         Ac              = Diam * DL
         CD              = CDrag * rho/2d0 * Ac * Ur_n_Meas 
         FDrag     (1:3) = CD * Ur_n(1:3)

!--------- buoyancy Volume Method check
!!!      FDrag     (  3) = Fdrag(3) + DL * Diam**2/4d0 * PI * rhoGrav

         call EXTEPR ( RG, FDrag, MDrag )

         floater(nf)%AQ_morison(1:3) = floater(nf)%AQ_morison(1:3) + FDrag(1:3)
         floater(nf)%AQ_morison(4:6) = floater(nf)%AQ_morison(4:6) + MDrag(1:3)

!-- new: AC_morison
         do j=1,3
          do i=1,3
             LLt(i,j)  = Ey_unit(i)*Ey_unit(j)
          enddo
         enddo

         ACL = 0.d0

         do iq = 1, 3
            ACL(iq ,iq) = 1.d0
            ACL(1:3,iq+3) = matmul (A0_float(iq,1:3,1:3),RL(1:3))
         enddo
            call DIAGO(3, AImLLt)
            AImLLt = AImLLt - LLt;

            ACL(1:3,1:6 ) = matmul ( AImLLt(1:3,1:3), ACL(1:3,1:6) )
            ACL(1:3,1:6 ) = ACL(1:3,1:6) * CD

            ACL(4  ,1:6 ) = ACL(3  ,1:6)*RG(2) - ACL(2  ,1:6)*RG(3)
            ACL(5  ,1:6 ) = ACL(1  ,1:6)*RG(3) - ACL(3  ,1:6)*RG(1)
            ACL(6  ,1:6 ) = ACL(2  ,1:6)*RG(1) - ACL(1  ,1:6)*RG(2)

         floater(nf)%AC_morison(1:6,1:6) = floater(nf)%AC_morison(1:6,1:6) + ACL(1:6,1:6)
!        floater(nf)%AK_morison(1:6,1:6) = 0.d0;

      enddo !iz
   enddo !nb_mo


   do npl_mo = 1, Nplates_mo

      m1            = Nod_pl_mo(1,npl_mo)
      m2            = Nod_pl_mo(2,npl_mo)

      XG_nod(1:3,1) = R_float(1:3) + matmul(A_float(1:3,1:3),XN_mo(1:3,m1))
      XG_nod(1:3,2) = R_float(1:3) + matmul(A_float(1:3,1:3),XN_mo(1:3,m2))

      if (XG_nod(3,1)>Z_calc) cycle        ! plate out of water
      
      Ey(1:3)      = XG_nod(1:3,1) - XG_nod(1:3,2)
      Length       = dsqrt(dot_product (Ey,Ey))
      Ey_unit(1:3) = Ey(1:3)/Length


         XG (1:3) = XG_nod   (1:3,1)
         RG (1:3) = XG       (1:3  ) - R_float(1:3)
         RL (1:3) = matmul ( AT_float(1:3,1:3),  RG     (1:3) )
         UG (1:3) = matmul ( DA_float(1:3,1:3),  RL     (1:3) ) + DR_float(1:3)

         call UINFLOW_sea ( XG, UINF, AINF, pdyn, 2 )

         Ur        (1:3) = UINF(1:3) - UG(1:3)
         Ur_t      (1:3) = dot_product(Ur(1:3),Ey_unit(1:3)) * Ey_unit(1:3)
         Ur_t_Meas       = dsqrt ( dot_product ( Ur_t, Ur_t ) )

         CDrag           = Cd_pl_mo(npl_mo)
         Ac              =  D_pl_mo(npl_mo)**2 * pi_h/4.d0
         CD              = CDrag * rho/2d0 * Ac * Ur_t_Meas 
         FDrag     (1:3) = CD * Ur_t(1:3)

         call EXTEPR ( RG, FDrag, MDrag )

         floater(nf)%AQ_morison(1:3) = floater(nf)%AQ_morison(1:3) + FDrag(1:3)
         floater(nf)%AQ_morison(4:6) = floater(nf)%AQ_morison(4:6) + MDrag(1:3)

!-- new: AC_morison
         do j=1,3
          do i=1,3
             LLt(i,j)  = Ey_unit(i)*Ey_unit(j)
          enddo
         enddo

         ACL = 0.d0

         do iq = 1, 3
            ACL(iq ,iq) = 1.d0
            ACL(1:3,iq+3) = matmul (A0_float(iq,1:3,1:3),RL(1:3))
         enddo
            call DIAGO(3, AImLLt)
            AImLLt = AImLLt - LLt;

            ACL(1:3,1:6 ) = matmul ( AImLLt(1:3,1:3), ACL(1:3,1:6) )
            ACL(1:3,1:6 ) = ACL(1:3,1:6) * CD

            ACL(4  ,1:6 ) = ACL(3  ,1:6)*RG(2) - ACL(2  ,1:6)*RG(3)
            ACL(5  ,1:6 ) = ACL(1  ,1:6)*RG(3) - ACL(3  ,1:6)*RG(1)
            ACL(6  ,1:6 ) = ACL(2  ,1:6)*RG(1) - ACL(1  ,1:6)*RG(2)

         floater(nf)%AC_morison(1:6,1:6) = floater(nf)%AC_morison(1:6,1:6) + ACL(1:6,1:6)

   enddo !npl_mo


 END Subroutine Morison_floater
!----------------------------------------------------------------------
 Subroutine Morison_floater_init
!----------------------------------------------------------------------

 Use Hydro
 Use Floater_Morison
 Use Paths

   implicit none

   integer              :: nn_mo,nb_mo,npl_mo,n,n1,n2,i, NNTwr_mo
   integer, allocatable :: IORDnod_mo(:) !switch input nodes to 1:NNT_na


   open (1,file=trim(file_morison)) !'morison.inp'
   open (2,file='floater_geom.dat')


!-- Read jacket NODES. Index of node and 3D position wrt global c.o. system
!-- Number of nodes needed and number of maximum index (NNTwr_mo)

   read(1,*)!**NNT_mo,NNTwr_mo
   read(1,*) NNT_mo,NNTwr_mo

   Allocate ( XN_mo     (3,NNT_mo  ),& 
              IORDnod_mo(  NNTwr_mo)  )

   read(1,*)!**nod	X[m]	Y[m]	Z[m]
   do nn_mo = 1, NNT_mo
      read(1,*) n, XN_mo(1,nn_mo), XN_mo(2,nn_mo), XN_mo(3,nn_mo)
      IORDnod_mo (n) = nn_mo
   enddo

!-- Read jacket's BODIES given in undeformed condition by setting 1st and last node
!-- Also set the morison variables D, Cd, NZ

   read(1,*)!**NBODBT_mo
   read(1,*) NBODBT_mo

   Allocate ( Nod_mo (2,NBODBT_mo),&
              NZ_mo  (  NBODBT_mo),&
              D_mo   (2,NBODBT_mo),&
              Cd_mo  (2,NBODBT_mo)  )


   read(1,*)!**nbod	nod1	nod2	D1	D2	Cd1	Cd2	NZ

   do nb_mo = 1, NBODBT_mo
      read(1,*) n,n1,n2, D_mo(1,nb_mo),D_mo(2,nb_mo), Cd_mo(1,nb_mo),Cd_mo(2,nb_mo),NZ_mo(nb_mo)
      Nod_mo(1,nb_mo)=IORDnod_mo(n1)
      Nod_mo(2,nb_mo)=IORDnod_mo(n2)

!-- check 
      if (XN_mo(3,IORDnod_mo(n2)) < XN_mo(3,IORDnod_mo(n1))) then
          write(*,*) 'ERROR, Jacket nod2 must be upper than nod1',n,n1,n2
          stop
      endif

!-- write geometry
      n1 = Nod_mo(1,nb_mo)
      n2 = Nod_mo(2,nb_mo)

      write(2,100) (XN_mo(i,n1),i=1,3)
      write(2,100) (XN_mo(i,n2),i=1,3)
      write(2,100)
      write(2,100)

   enddo


!-- Read plates where drag force is calculated
   read(1,*)!**Nplates_mo
   read(1,*) Nplates_mo

   Allocate ( Nod_pl_mo (2,Nplates_mo),&
              D_pl_mo   (  Nplates_mo),&
              Cd_pl_mo  (  Nplates_mo)  )
                               

   read(1,*)!**npl nod1 nod2  D      Cd   !plate is at nod1, the 2nd nod determins U_n

   do npl_mo = 1, Nplates_mo
      read(1,*) n,n1,n2, D_pl_mo(npl_mo), Cd_pl_mo(npl_mo)
      Nod_pl_mo(1,npl_mo)=IORDnod_mo(n1)
      Nod_pl_mo(2,npl_mo)=IORDnod_mo(n2)
   enddo

   close(1)
   close(2)

   Deallocate ( IORDnod_mo )

 100 format ( 3f15.5 )


 END Subroutine Morison_floater_init
!
!
!
!---------------------------------------------------------------
! 1. M,C,K,Q from main routines need update for more floaters           **** OK
! 2. writefloater loads only for nf=1
! 3. morison             >>      >>
! 4. ICASE_el = 4 not working for NFLOATER_el >1                        **** DOESN'T MATTER 
! 5. QS_equate - initia
! 6. truss only to 1st floater atm                                      **** DOESN'T MATTER
