!----------------------------------------------------------------------
 Subroutine MOORINGS_tr (TTime, int_type, IP, RLX)
!----------------------------------------------------------------------
! Overview / notes for future maintenance:
! - MOORINGS_tr drives a full time step by looping over sub-steps (NSubTime_tr),
!   refreshing connection boundary conditions (set_connect_bc), assembling
!   system matrices (MATRIX_tr, MATRIX_REDUCT_tr) and integrating either
!   statically or via Newmark (Time_Integrate_tr). The logic assumes all
!   body properties such as the unstretched element length (ALENG_tr) are
!   fixed during the iteration; changing them mid-step will require
!   re-deriving any quantities cached in the body_tr and connect_tr arrays.
! - Connection loads are accumulated at the end of the loop and mapped back
!   to the floater/elastic model through truss2float_loads_uncoupl. Any
!   updates to geometric quantities (e.g., dynamic "Length0" adjustments)
!   will need to stay consistent with these downstream mappings.

 use truss
 use Cbeam

   implicit none

   real(8), intent(in) :: TTime, RLX
   integer, intent(in) :: int_type, IP
!-- local vars
   integer :: nc,nb_tr,nod,iloc
   integer :: ICONV_tr, it, nt


   UT_tr (1:NDFT_tr) = UTPP_tr (1:NDFT_tr);
   UT1_tr(1:NDFT_tr) = UTPP1_tr(1:NDFT_tr);
   UT2_tr(1:NDFT_tr) = UTPP2_tr(1:NDFT_tr);

   call set_floater_sway(R_float(2))

   do nt = 1, NSubTime_tr

      UTP_tr (1:NDFT_tr) = UT_tr (1:NDFT_tr);
      UTP1_tr(1:NDFT_tr) = UT1_tr(1:NDFT_tr);
      UTP2_tr(1:NDFT_tr) = UT2_tr(1:NDFT_tr);

      TIME_tr  = TTime - DT_tr * dble(NSubTime_tr-nt)
      ICONV_tr = 0

      call update_repositioning(TIME_tr, DT_tr)

      call set_connect_bc(nt)                                            !-- set the connection vars, based on floater's motion


      do it = 1,ITERMAX_tr
         call MATRIX_tr                                                  !- set the local and the global truss matrices
         call MATRIX_REDUCT_tr (IP)                                      !- reduce the matrix for boundary conditions
         call Time_Integrate_tr(it, int_type, ICONV_tr, MAXERR_tr, RLX)  !- time integration, static solution or newmark

         if (ICONV_tr == 1) goto 1                                       !- convergence check
      enddo

      write(*,*)'Not converged'

 1 enddo !nt
 

!--- Calculate connection loads
   do nc    = 1, NCONNECT_tr
      nb_tr = connect_tr(nc)%IBODCONNECT_tr
      nod   = connect_tr(nc)%INODCONNECT_tr
      iloc  = (nod-1)*NDFPN_tr

      connect_tr(nc)%ALOADB (1:3) = body_tr(nb_tr)%AQLOC_tr(iloc+1:iloc+3)
   enddo !nc


 END Subroutine MOORINGS_tr
!----------------------------------------------------------------------
 Subroutine set_connect_bc (nt)
!----------------------------------------------------------------------
!
! based on Newmark-b method
!

 use truss

   implicit none

   integer,intent(in) :: nt
!-- local variables
   integer :: nc, nb_tr, ni, inod, idof
   real(8) :: UTset (NCONNECT_tr,3)
   real(8) :: UT1set(NCONNECT_tr,3)
   real(8) :: UT2set(NCONNECT_tr,3)
   real(8) :: DT


   if (nt < NSubTime_tr) then

         DT = dble(nt)*DT_tr
      do nc = 1, NCONNECT_tr
         UT2set(nc, 1:3) =   connect_tr(nc)%UTP2_CONNECT_tr(1:3) + &
                           ( connect_tr(nc)% UT2_CONNECT_tr(1:3) - &
                             connect_tr(nc)%UTP2_CONNECT_tr(1:3)    ) *dble(nt)/dble(NSubTime_tr)

         UT1set(nc, 1:3) =   connect_tr(nc)%UTP1_CONNECT_tr(1:3)                                  + &
                             connect_tr(nc)%UTP2_CONNECT_tr(1:3) * (1.0d0-GAMMA_tr) * DT          + &
                             UT2set    (nc,                 1:3) *        GAMMA_tr  * DT

         UTset (nc, 1:3) =   connect_tr(nc)% UTP_CONNECT_tr(1:3)                                  + &
                             connect_tr(nc)%UTP1_CONNECT_tr(1:3)                    * DT          + &
                             connect_tr(nc)%UTP2_CONNECT_tr(1:3) * (0.5d0-BITA_tr)  * DT**2       + &
                             UT2set    (nc,                 1:3) *        BITA_tr   * DT**2
      enddo !nc

   elseif (nt==NSubTime_tr) then

      do nc              = 1, NCONNECT_tr
         UT2set(nc, 1:3) = connect_tr(nc)%UT2_CONNECT_tr(1:3)
         UT1set(nc, 1:3) = connect_tr(nc)%UT1_CONNECT_tr(1:3)
         UTset (nc, 1:3) = connect_tr(nc)% UT_CONNECT_tr(1:3)
      enddo !nc

   endif


   do nc    = 1, NCONNECT_tr
      nb_tr = connect_tr(nc)%IBODCONNECT_tr
      ni    = connect_tr(nc)%INODCONNECT_tr
      inod  = body_tr(nb_tr)%INODE_tr (ni)
      idof  = NDFACC_tr (inod-1)

       UT_tr(idof+1:idof+3) =  UTset(nc, 1:3) - XG_tr(1:3, inod)
      UT1_tr(idof+1:idof+3) = UT1set(nc, 1:3)
      UT2_tr(idof+1:idof+3) = UT2set(nc, 1:3)
   enddo !nc


 END Subroutine set_connect_bc
!----------------------------------------------------------------------
 Subroutine truss2float_loads_uncoupl
!----------------------------------------------------------------------
!
! Add the mooring loads contribution from each connection point,
! tranfered to the floater's referance point [ICASE_el == 3]
! or to the corresponding elastic element    [ICASE_el == 4].
!
! Remark: used only for uncoupled mooring system.
!
!----------------------------------------------------------------------

 use truss
 use Cbeam !A_float, Qtruss
 use Hydro !floater

   implicit none

   real(8) :: RG(3), FF(3), MM(3)
   integer :: nc
   real(8) :: SHAPE (NEQPEM_el,NDFPEM_el)
   real(8) :: SHAPET(NDFPEM_el,NEQPEM_el)
   real(8) :: AFLOC0(NDFPEM_el          )
   real(8) :: RLOC(3)
   real(8) :: AT(3,3)
   real(8) :: AQ(6)
   real(8) :: HTA,HTAL,PHPX,PHPZ,ALLOC
   integer :: nbod, nbsub, nel, nb_el, nn_el, NDFPE, NEQPE, nn, i, nf

   do nf = 1, NFLOATER_el
      floater(nf)% Q_truss (1:6) = 0.d0
   enddo

!qq Valid for moorings attached to 1st floater
   nf = 1

   if    (ICASE_el==3) then

      do nc = 1, NCONNECT_tr
!--------- fault case: 1 mooring line is cut
!        if ((TIME > 1017.49d0).and.(nc == 1)) cycle
!!!      if ((TIME >= 250.00d0).and.(nc == 1)) cycle
       
         RG (1:3) = matmul ( A_float  (1:3,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )
         FF (1:3) = connect_tr(nc)%ALOADB (1:3)
       
         call EXTEPR ( RG, FF, MM )
       
         floater(nf)% Q_truss (1:3) = floater(nf)% Q_truss (1:3) + FF (1:3)
         floater(nf)% Q_truss (4:6) = floater(nf)% Q_truss (4:6) + MM (1:3)
      enddo !nc

   elseif (ICASE_el==4) then

      do nc = 1, NCONNECT_tr
!--------- fault case: 1 mooring line is cut
!        if ((TIME > 1017.49d0).and.(nc == 1)) cycle
!!!      if ((TIME >= 250.00d0).and.(nc == 1)) cycle

         nbod    = connect_tr(nc)%nbod_ja
         nbsub   = 1
         nel     = connect_tr(nc)%nel_ja
         nb_el   = body(nbod)%NBODTGLB_el(nbsub)
         nn_el   = body(nbod)%NNODTGLB_el(nbsub)
         RLOC(1) = connect_tr(nc)%RLOC_ja(1)       !from the start of the element
         RLOC(2) = 0.d0
         RLOC(3) = connect_tr(nc)%RLOC_ja(3)       !from the start of the element

         AT(1:3,1:3) = transf_mat(nn_el)%AT_el(1:3,1:3)

         ALLOC   = subbody(nb_el)%ALENG_el(nel)
         HTA     = connect_tr(nc)%RLOC_ja(2)
         HTAL    = HTA/ALLOC
         PHPX    = subbody(nb_el)%PHIX_el (nel)
         PHPZ    = subbody(nb_el)%PHIZ_el (nel)
         NEQPE   = subbody(nb_el)%NEQPE_el
         NDFPE   = subbody(nb_el)%NDFPE_el

         call SSHAPEFUNC15 (HTAL, ALLOC, PHPX, PHPZ, SHAPE, NDFPE, NEQPE)

         SHAPET  = transpose (SHAPE)

         FF(1:3) = matmul ( AT(1:3,1:3), connect_tr(nc)%ALOADB (1:3) )

         call EXTEPR ( RLOC, FF, MM )

!--------- Using shape functions, add to local element matrices
         AQ (IMDOF_el(1:3)) = FF (1:3)
         AQ (IMDOF_el(4:6)) = MM (1:3)
         AFLOC0             = matmul( SHAPET, AQ )
         nn                 =                                   subbody(nb_el)%NODMB       (nel, 1)!nnod)
         i                  = ACCU%NDFTACC_el ( nb_el-1 )   +   subbody(nb_el)%NDFPNBACC_el(nn-1  )
!--------- STORE_LOC_MATRIX_el
         subbody(nb_el)%AFLOC_el(nel, 1:NDFPE) = subbody(nb_el)%AFLOC_el(nel, 1:NDFPE) + AFLOC0(1:NDFPE)
!--------- ASSEMBLY_MATRIX_el
         AQ_el                  (i+1: i+NDFPE) = AQ_el                  (i+1: i+NDFPE) + AFLOC0(1:NDFPE)
      enddo !nc

   endif


 END Subroutine truss2float_loads_uncoupl
!----------------------------------------------------------------------
 Subroutine calc_connection_UT_tr !UT2_CONNECT_MOORINGS_tr
!----------------------------------------------------------------------

 use truss
 use Cbeam ! [DDR,R]_float, [DDA,A]_float

   implicit none

   real(8) :: SHAPE (NEQPEM_el,NDFPEM_el)
   real(8) :: UT    (NDFPEM_el)
   real(8) :: UT1   (NDFPEM_el)
   real(8) :: UT2   (NDFPEM_el)
   real(8) :: UTT0  (NEQPEM_el)
   real(8) :: UTT1  (NEQPEM_el)
   real(8) :: UTT2  (NEQPEM_el)
   real(8) :: RLOC(3)
   real(8) :: SS(3,6)
   real(8) :: A(3,3), DA(3,3), DDA(3,3), R(3), DR(3), DDR(3)
   real(8) :: HTA,HTAL,HTA0,PHPX,PHPZ,ALLOC
   integer :: nbod, nbsub, nel, nb_el, nn_el, NDFPE, NEQPE
   integer :: nc


   if    (ICASE_el==3) then

      do nc = 1, NCONNECT_tr
         connect_tr(nc)% UT_CONNECT_tr(1:3) =   R_float(1:3) + matmul(   A_float(1:3,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )
         connect_tr(nc)%UT1_CONNECT_tr(1:3) =  DR_float(1:3) + matmul(  DA_float(1:3,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )
         connect_tr(nc)%UT2_CONNECT_tr(1:3) = DDR_float(1:3) + matmul( DDA_float(1:3,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )
      enddo !nc

   elseif (ICASE_el==4) then

      do nc = 1, NCONNECT_tr
             nbod      = connect_tr(nc)%nbod_ja
             nbsub     = 1
             nel       = connect_tr(nc)%nel_ja
             nb_el     = body(nbod)%NBODTGLB_el(nbsub)
             nn_el     = body(nbod)%NNODTGLB_el(nbsub)
             RLOC(1:3) = connect_tr(nc)%RLOC_ja(1:3)       !from the start of the element

           A (1:3,1:3) = transf_mat(nn_el)%A_el  (1:3,1:3)
          DA (1:3,1:3) = transf_mat(nn_el)%DA_el  (1:3,1:3)
         DDA (1:3,1:3) = matmul(A(1:3,1:3),transf_mat(nn_el)%ATDDA_el  (1:3,1:3))
           R (1:3    ) = transf_mat(nn_el)%R_el  (1:3    )
          DR (1:3    ) = transf_mat(nn_el)%DR_el (1:3    )
         DDR (1:3    ) = matmul(A(1:3,1:3),transf_mat(nn_el)%ATDDR_el(1:3))

         call Calc_SS ( SS, RLOC )

         ALLOC   = subbody(nb_el)%ALENG_el(nel)
         HTA0    = subbody(nb_el)%HTA_el  (nel)
         RLOC(2) = RLOC(2) + HTA0
         HTA     = connect_tr(nc)%RLOC_ja(2)
         HTAL    = HTA/ALLOC
         PHPX    = subbody(nb_el)%PHIX_el (nel)
         PHPZ    = subbody(nb_el)%PHIZ_el (nel)
         NEQPE   = subbody(nb_el)%NEQPE_el
         NDFPE   = subbody(nb_el)%NDFPE_el

         call SSHAPEFUNC15 (HTAL, ALLOC, PHPX, PHPZ, SHAPE, NDFPE, NEQPE)
         call LOCAL_UT_1  (nb_el, nel, UT, UT1, UT2)

         UTT0  = matmul (SHAPE, UT )
         UTT1  = matmul (SHAPE, UT1)
         UTT2  = matmul (SHAPE, UT2)

         connect_tr(nc)% UT_CONNECT_tr(1:3) =   R(1:3) +        matmul(   A(1:3,1:3), matmul(SS(1:3,1:6), UTT0(IMDOF_el(1:6))) + RLOC(1:3) )
         connect_tr(nc)%UT1_CONNECT_tr(1:3) =  DR(1:3) +        matmul(  DA(1:3,1:3), matmul(SS(1:3,1:6), UTT0(IMDOF_el(1:6))) + RLOC(1:3) ) &
                                                       +        matmul(   A(1:3,1:3), matmul(SS(1:3,1:6), UTT1(IMDOF_el(1:6)))             )
         connect_tr(nc)%UT2_CONNECT_tr(1:3) = DDR(1:3) +        matmul( DDA(1:3,1:3), matmul(SS(1:3,1:6), UTT0(IMDOF_el(1:6))) + RLOC(1:3) ) &
                                                       + 2.d0 * matmul(  DA(1:3,1:3), matmul(SS(1:3,1:6), UTT1(IMDOF_el(1:6)))             ) &
                                                       +        matmul(   A(1:3,1:3), matmul(SS(1:3,1:6), UTT2(IMDOF_el(1:6)))             )
      enddo !nc

   endif


 END Subroutine calc_connection_UT_tr
!----------------------------------------------------------------------
 Subroutine MOORINGS_UPDATE_tr
!----------------------------------------------------------------------

 use truss
 use Cbeam ![DDR,DR,R]_float, [DDA,DA,A]_float]

   implicit none

   integer :: nc


   UTPP_tr (1:NDFT_tr) = UT_tr (1:NDFT_tr);
   UTPP1_tr(1:NDFT_tr) = UT1_tr(1:NDFT_tr);
   UTPP2_tr(1:NDFT_tr) = UT2_tr(1:NDFT_tr);

   do nc = 1, NCONNECT_tr
      connect_tr(nc)%  UTP_CONNECT_tr(1:3) = connect_tr(nc)%  UT_CONNECT_tr(1:3)
      connect_tr(nc)% UTP1_CONNECT_tr(1:3) = connect_tr(nc)% UT1_CONNECT_tr(1:3)
      connect_tr(nc)% UTP2_CONNECT_tr(1:3) = connect_tr(nc)% UT2_CONNECT_tr(1:3)
   enddo !nc


END Subroutine MOORINGS_UPDATE_tr
!----------------------------------------------------------------------
 Subroutine set_floater_sway(y_value)
!----------------------------------------------------------------------

 use truss

   implicit none

   real(8), intent(in) :: y_value


   y_float = y_value

   if (repos_active .and. .not. y_ref_set) then
      y_ref     = y_value
      y_ref_set = .true.
   endif


 END Subroutine set_floater_sway
!----------------------------------------------------------------------
Subroutine init_repositioning(path_truss)
!----------------------------------------------------------------------

 use truss

   implicit none

  integer :: ios, i, last_slash
  character(len=256) :: line, line_trim
  character(len=*), parameter :: fname = 'repositioning.inp'
  character(len=*), intent(in) :: path_truss
  character(len=256) :: candidate


  repos_active   = .false.
  y_ref_set      = .false.
  N_repos_pts    = 0
  N_winch_elem   = 0
  y_ref          = 0.d0
  y_target_curr  = 0.d0
  y_err_int      = 0.d0

   if (allocated(t_repos   )) deallocate(t_repos)
   if (allocated(y_repos   )) deallocate(y_repos)
   if (allocated(ALENG0_tr )) deallocate(ALENG0_tr)
   if (allocated(ALENG_ctrl)) deallocate(ALENG_ctrl)
   if (allocated(ALENG_min )) deallocate(ALENG_min)
   if (allocated(ALENG_max )) deallocate(ALENG_max)
   if (allocated(winch_elem)) deallocate(winch_elem)

  candidate = fname

!--- try alongside the truss input path first, then fall back to cwd
  last_slash = 0

  do i = len_trim(path_truss), 1, -1
     if (path_truss(i:i) == '/') then
        last_slash = i
        exit
     endif
  enddo

  if (last_slash > 0) then
     candidate = trim(path_truss(1:last_slash)) // fname
  endif

  open (unit=33, file=trim(candidate), status='old', iostat=ios)

  if (ios /= 0) then
     open (unit=33, file=trim(fname), status='old', iostat=ios)
  endif

  if (ios /= 0) then
     write(*,*) 'Repositioning: input file not found, controller disabled'
     return
  endif

!------ find the first non-comment line with controller constants
   do
      read (33,'(A)', iostat=ios) line
      if (ios /= 0) then
         close(33)
         return
      endif

      line_trim = adjustl(line)

      if (len_trim(line_trim) == 0) cycle
      if (line_trim(1:1) == '!') cycle

      read (line_trim,*, iostat=ios) N_repos_pts, repos_vwinch, repos_Kp, repos_Ki
      if (ios /= 0 .or. N_repos_pts <= 0) then
         close(33)
         return
      endif
      exit
   enddo

   allocate ( t_repos(1:N_repos_pts), stat=ios )
   if (ios /= 0) then
      close(33)
      return
   endif

   allocate ( y_repos(1:N_repos_pts), stat=ios )
   if (ios /= 0) then
      deallocate(t_repos)
      close(33)
      return
   endif

   do i = 1, N_repos_pts

      do
         read (33,'(A)', iostat=ios) line

         if (ios /= 0) then
            close(33)
            deallocate(t_repos, y_repos)
            return
         endif

         line_trim = adjustl(line)

         if (len_trim(line_trim) == 0) cycle
         if (line_trim(1:1) == '!') cycle

         read (line_trim,*, iostat=ios) t_repos(i), y_repos(i)
         if (ios /= 0) then
            close(33)
            deallocate(t_repos, y_repos)
            return
         endif

         exit
      enddo

   enddo

   repos_active = .true.

   write(*,*) 'Repositioning: controller active with ', N_repos_pts, ' pts; VWINCH=', repos_vwinch, ' Kp=', repos_Kp, ' Ki=', repos_Ki

   close(33)


 END Subroutine init_repositioning
!----------------------------------------------------------------------
 Subroutine init_winch_elements()
!----------------------------------------------------------------------

 use truss

   implicit none

   integer :: e, ios


   if (.not. repos_active) return

   if (allocated(ALENG0_tr )) deallocate(ALENG0_tr)
   if (allocated(ALENG_ctrl)) deallocate(ALENG_ctrl)
   if (allocated(ALENG_min )) deallocate(ALENG_min)
   if (allocated(ALENG_max )) deallocate(ALENG_max)
   if (allocated(winch_elem)) deallocate(winch_elem)

   allocate ( ALENG0_tr (NBODT_tr), stat=ios )
   if (ios /= 0) then
      repos_active = .false.
      return
   endif

   allocate ( ALENG_ctrl(NBODT_tr), stat=ios )
   if (ios /= 0) then
      repos_active = .false.
      deallocate(ALENG0_tr)
      return
   endif

   allocate ( ALENG_min (NBODT_tr), stat=ios )
   if (ios /= 0) then
      repos_active = .false.
      deallocate(ALENG0_tr, ALENG_ctrl)
      return
   endif

   allocate ( ALENG_max (NBODT_tr), stat=ios )
   if (ios /= 0) then
      repos_active = .false.
      deallocate(ALENG0_tr, ALENG_ctrl, ALENG_min)
      return
   endif

   do e = 1, NBODT_tr
      ALENG0_tr (e) = body_tr(e)%ALENG_tr
      ALENG_ctrl(e) = ALENG0_tr(e)
      ALENG_min (e) = 0.7d0  * ALENG0_tr(e)
      ALENG_max (e) = 1.05d0 * ALENG0_tr(e)
   enddo

   N_winch_elem = min(5, NBODT_tr)

   if (N_winch_elem <= 0) then
      repos_active = .false.
      return
   endif

   allocate ( winch_elem(1:N_winch_elem), stat=ios )
   if (ios /= 0) then
      repos_active = .false.
      return
   endif

   do e = 1, N_winch_elem
      winch_elem(e) = e
   enddo

   write(*,*) 'Repositioning: controlling', N_winch_elem, 'elements starting at index 1'


END Subroutine init_winch_elements
!----------------------------------------------------------------------
 Subroutine update_repositioning(TTIME, DT)
!----------------------------------------------------------------------

 use truss

   implicit none

   real(8), intent(in) :: TTIME
   real(8), intent(in) :: DT

   real(8) :: e, dL_total_dt, dL_step
   real(8) :: w_sum
   integer :: i, idx


   if (.not. repos_active) return
   if (N_winch_elem <= 0) return
   if (.not. y_ref_set) return

!--- interpolate current target
   if (TTIME <= t_repos(1)) then

      y_target_curr = y_repos(1)

   elseif (TTIME >= t_repos(N_repos_pts)) then

      y_target_curr = y_repos(N_repos_pts)

   else

      do i = 1, N_repos_pts-1
         if ( (TTIME >= t_repos(i)) .and. (TTIME <= t_repos(i+1)) ) then
            if (t_repos(i+1) > t_repos(i)) then
               y_target_curr = y_repos(i) + &
                 (y_repos(i+1)-y_repos(i)) * (TTIME - t_repos(i)) / (t_repos(i+1)-t_repos(i))
            else
               y_target_curr = y_repos(i)
            endif
            exit
         endif
      enddo

   endif

!--- control law
   e         = y_target_curr - (y_float - y_ref)
   y_err_int = y_err_int + e * DT

   dL_total_dt = repos_Kp * e + repos_Ki * y_err_int

   if (dL_total_dt >  repos_vwinch) dL_total_dt =  repos_vwinch
   if (dL_total_dt < -repos_vwinch) dL_total_dt = -repos_vwinch

   dL_step = dL_total_dt * DT

   if (dL_step == 0.d0) return

!--- count available segments
   w_sum = 0.d0

   if (dL_step < 0.d0) then
      do i = 1, N_winch_elem
         idx = winch_elem(i)
         if (ALENG_ctrl(idx) > ALENG_min(idx)) w_sum = w_sum + 1.d0
      enddo
   else
      do i = 1, N_winch_elem
         idx = winch_elem(i)
         if (ALENG_ctrl(idx) < ALENG_max(idx)) w_sum = w_sum + 1.d0
      enddo
   endif

   if (w_sum <= 0.d0) return

!--- distribute the change
   do i = 1, N_winch_elem

      idx = winch_elem(i)

      if ( (dL_step < 0.d0) .and. (ALENG_ctrl(idx) <= ALENG_min(idx)) ) cycle
      if ( (dL_step > 0.d0) .and. (ALENG_ctrl(idx) >= ALENG_max(idx)) ) cycle

      ALENG_ctrl(idx) = ALENG_ctrl(idx) + dL_step / w_sum

      if (ALENG_ctrl(idx) < ALENG_min(idx)) ALENG_ctrl(idx) = ALENG_min(idx)
      if (ALENG_ctrl(idx) > ALENG_max(idx)) ALENG_ctrl(idx) = ALENG_max(idx)

      body_tr(idx)%ALENG_tr = ALENG_ctrl(idx)

   enddo


 END Subroutine update_repositioning
!----------------------------------------------------------------------
 Subroutine INIT_tr ( GRAV_in, rho_in, depth_in, DT_in, path )
!----------------------------------------------------------------------

 use truss
 use Cbeam !ICASE_el, body(nbod)%NBODTGLB_el(nbsub), subbody(nb_el)%HTA_el(nel), transf_mat(nn_el)%AT_el(1:3,1:3)
           !          body(nbod)%NNODTGLB_el(nbsub),                             transf_mat(nn_el)%R_el (1:3    )
   implicit none

   real(8)       :: GRAV_in, rho_in, depth_in, DT_in
   integer       :: nb_tr, inod, nc, icon, n1, n2, nbod, nod, k
   real(8)       :: DX(3), ALENGP
   character*256 :: path
!--- for ICASE_el==4
   integer       :: nbsub, nel, nb_el, nn_el
   real(8)       :: Y0(3), AT(3,3), R(3)


   IREDUCT_tr = 1
   pi_tr      = dacos(-1d0)
   GRAV_tr    = GRAV_in
   AROW_tr    = rho_in
   Depth_tr   = depth_in

   open (1,file=trim(path)) !'truss.inp'

!-------------------------
! -Main Parameters
!-------------------------
      read (1,*) !*** Main Parameters ***
      read (1,*) NSubTime_tr
      read (1,*) ITERMAX_tr
      read (1,*) MAXERR_tr
      read (1,*) COEFM_tr
      read (1,*) COEFK_tr
      read (1,*) BITA_tr
      read (1,*) GAMMA_tr
      read (1,*) AMVISCW_tr
      read (1,*) IYN_Morison_tr       !-- IYN_Morison_tr    [0: no hydro loads, 1: Morison's eq.]
      read (1,*) IYN_SeabedInt_tr     !-- IYN_SeabedInt_tr  [O: no Seabed Int , 1: yes          ]
      read (1,*) IWRITE_truss         !-- 1 output elem

   TIME_tr    = 0.d0
   DT_tr      = DT_in / dble(NSubTime_tr)

      write(10,*)
      write(10,*) '*** Truss Init ***'
      write(10,*) 'NSubTime_tr     :',NSubTime_tr
      write(10,*) 'ITERMAX_tr      :',ITERMAX_tr
      write(10,*) 'MAXERR_tr       :',MAXERR_tr
      write(10,*) 'IYN_Morison_tr  :',IYN_Morison_tr
      write(10,*) 'IYN_SeabedInt_tr:',IYN_SeabedInt_tr
      write(10,*)

!-------------------------
! -Nodes co-ordinates
!-------------------------
      read (1,*)
      read (1,*) !*** NODES ***
      read (1,*) NNODE_tr
      read (1,*) !**       Xg             Yg             Zg


!------ Allocate: 1 (NNODE_tr vars)
      Allocate ( NDFACC_tr (0:NNODE_tr) )
      Allocate ( XG_tr     (3,NNODE_tr) )

      do inod = 1,NNODE_tr
         read (1,*) n1,XG_tr(1,inod),XG_tr(2,inod),XG_tr(3,inod)
      enddo

!-------------------------------
! -Elements/Bodies Definition
!-------------------------------
      read (1,*)
      read (1,*) !*** BODIES ***
      read (1,*) NBODT_tr
      read (1,*) !**  nbod nod1 nod2    Length0    Diameter       dens        EA1      EA2      structural damping coeff


!------ Allocate: 2 (NBODT_tr vars)
      Allocate ( body_tr(NBODT_tr) )


      do nb_tr = 1,NBODT_tr
         read (1,*)                    &
           nbod                       ,&   ! Dummy Variable
           body_tr(nb_tr)%INODE_tr(1) ,&   ! Corresponding node 1
           body_tr(nb_tr)%INODE_tr(2) ,&   ! Corresponding node 2
           body_tr(nb_tr)%ALENG_tr    ,&   ! Unstretched Length
           body_tr(nb_tr)%ADIAM_tr    ,&   ! Initial Diameter [m]    [needed for hydrodynamic loading - Morison's eq.]
           body_tr(nb_tr)%DENS_tr     ,&   ! mass density     [kg/L] [needed for inertia      loading]
           body_tr(nb_tr)%WaterFz_tr  ,&   ! weight in water  [N/L]
           body_tr(nb_tr)%EAS1_tr     ,&   ! stiffness for thlipsi [small value]
           body_tr(nb_tr)%EAS2_tr     ,&   ! stiffness for extension
           body_tr(nb_tr)%CDAMP_tr    ,&   ! structural damping constant
           body_tr(nb_tr)%Ca_tr       ,&   ! Added mass coeff ca=cm-1
           body_tr(nb_tr)%Cdnorm_tr   ,&   ! drag coeff for normal     direction
           body_tr(nb_tr)%Cdtang_tr   ,&   ! drag coeff for tangential direction
           body_tr(nb_tr)%Cdfric_tr        !             ****** ATM is no USED!!! ******

!--------- Notes: ALENG_tr is treated as a constant reference length throughout
!          the module. It feeds seabed interaction stiffness (bed_K/bed_C),
!          inertia distribution, and hydrodynamic forcing coefficients. If a
!          future feature needs to update the "Length0"/unstretched length at
!          runtime, all dependent matrices/loads that are preassembled during
!          INIT_tr or cached in body_tr(nb_tr)%AKLOC_tr, ACLOC_tr, AQLOC_tr
!          will have to be rebuilt before the next call to MATRIX_tr.

           body_tr(nb_tr)%AQLOC_tr = 0.d0; ! zero for the 1st writeout at time=0
      enddo


!-------------------------------
! -Concentrated Properties
!-------------------------------
      read (1,*)
      read (1,*) !*** CONCENTRATED PROPERTIES ***
      read (1,*) NCON_tr
      read (1,*) !**   nbod nod  mass    diameter


!------ Allocate: 3 (NCON_tr vars)
      Allocate ( conc_mass_tr(NCON_tr) )
                                    

      do icon = 1, NCON_tr
         read (1,*)                        &
           conc_mass_tr(icon)%IBCON_tr    ,&
           conc_mass_tr(icon)%INCON_tr    ,&
           conc_mass_tr(icon)%AMASSCON_tr ,&
           conc_mass_tr(icon)%DIAMCON_tr
      enddo


!-------------------------------
! -Constraints
!-------------------------------
      read (1,*)
      read (1,*) !*** CONSTRAINTS ***
      read (1,*) NCONSTR_tr
      read (1,*) !**   node        dof1        dof2       dof3

!------ Allocate: 4 (NCONSTR_tr vars)
      Allocate ( constrain_tr(NCONSTR_tr) )


      do nc = 1, NCONSTR_tr
         read (1,*)                              &
           constrain_tr(nc)%INODC_tr            ,&      
          (constrain_tr(nc)%INDXC_tr(k), k=1,3)

         if (nc == 1) cycle
         if (constrain_tr(nc)%INODC_tr < constrain_tr(nc-1)%INODC_tr) then
            write(* ,*)'Constraints must be given with increasing order'
            write(10,*)'Constraints must be given with increasing order'
            stop
         endif
      enddo

!-------------------------------
! -Connections
!-------------------------------
      read (1,*)
      read (1,*) !*** CONNECTIONS ***
      read (1,*) NCONNECT_tr                  !-- read number of connections
      read (1,*) !** nbod        node


!------ Allocate: 5 (NCONNECT_tr vars)
      Allocate ( connect_tr(NCONNECT_tr) )


      do nc = 1, NCONNECT_tr

         connect_tr(nc)%ALOADB(1:3) = 0.d0;

         if     (ICASE_el == 3) then
            read (1,*)                                   &
              connect_tr(nc)%IBODCONNECT_tr             ,&
              connect_tr(nc)%INODCONNECT_tr             ,&
             (connect_tr(nc)%RLOC_CONNECT_tr(k), k=1,3)
         elseif (ICASE_el == 4) then
            read (1,*)                                   &
              connect_tr(nc)%IBODCONNECT_tr             ,&
              connect_tr(nc)%INODCONNECT_tr             ,&
             (connect_tr(nc)%RLOC_CONNECT_tr(k), k=1,3) ,&
              connect_tr(nc)%nbod_ja                    ,&
              connect_tr(nc)%nel_ja
         endif
      enddo !nc

   close (1)


!--- Modify the connection nodes, based on the floater's initial positions
   call set_init_truss_pos


!--- Set RLOC_ja: local distance between the connection and the elastic(beam) body 
!---              wrt the start of the element
   if (ICASE_el == 4) then
      do nc = 1, NCONNECT_tr
!--------- truss bodies
         nbod = connect_tr(nc)%IBODCONNECT_tr
         nod  = connect_tr(nc)%INODCONNECT_tr
         n1   = body_tr(nbod)%INODE_tr (nod)

!--------- beam bodies
         nbod        = connect_tr(nc)%nbod_ja
         nbsub       = 1
         nel         = connect_tr(nc)%nel_ja
         nb_el       = body(nbod)%NBODTGLB_el(nbsub)
         nn_el       = body(nbod)%NNODTGLB_el(nbsub)
         Y0(1:3    ) = 0.d0
         Y0(  2    ) = subbody(nb_el)%HTA_el(nel)
         AT(1:3,1:3) = transf_mat(nn_el)%AT_el(1:3,1:3)
         R (1:3    ) = transf_mat(nn_el)%R_el (1:3    )
!---------
!--------- RG=Ro+A(Yo+hta+Su) <->  Yo = At(RG-Ro)-hta
!---------    RG : the truss node
!---------    Ro, A, AT: jacket matrices
!---------    Su : initialy zero [no deflections]
!---------    hta: {0,hta,0} distance along beam axid of all previous elements of the body
!---------    Yo : {x,y  ,z} offsets wrt to the element [CALCULATED now]
!---------
         connect_tr(nc)%RLOC_ja(1:3) = matmul(  AT(1:3,1:3), XG_tr(1:3,n1) - R(1:3)  ) - Y0(1:3)

         write(10,*)'nc     ', nc
         write(10,*)'RLOC_ja', connect_tr(nc)%RLOC_ja(1:3)
      enddo !nc
   endif


!--- initial values for UT_CONNECT_tr, UT1_CONNECT_tr, UT2_CONNECT_tr
   do nc = 1, NCONNECT_tr
      nbod = connect_tr(nc)%IBODCONNECT_tr
      nod  = connect_tr(nc)%INODCONNECT_tr
      n1   = body_tr(nbod)%INODE_tr (nod)
      connect_tr(nc)%UT_CONNECT_tr (1:3) = XG_tr(1:3,n1)
      connect_tr(nc)%UT1_CONNECT_tr(1:3) = 0.d0
      connect_tr(nc)%UT2_CONNECT_tr(1:3) = 0.d0
   enddo !nc


!--- NDFT_tr variables
   NDFACC_tr(0) = 0
   do inod = 1,NNODE_tr
      NDFACC_tr(inod) = NDFACC_tr(inod-1) + NDFPN_tr
   enddo

   NDFT_tr = NDFACC_tr(NNODE_tr)

   write(*,*) 'NDFT_tr=',NDFT_tr


!--- Allocate: 6 (NDFT_tr vars)
   Allocate ( UT_tr  (NDFT_tr), UT1_tr  (NDFT_tr), UT2_tr  (NDFT_tr) )
   Allocate ( UTP_tr (NDFT_tr), UTP1_tr (NDFT_tr), UTP2_tr (NDFT_tr) )
   Allocate ( UTPP_tr(NDFT_tr), UTPP1_tr(NDFT_tr), UTPP2_tr(NDFT_tr) )

   Allocate ( AM_tr     (NDFT_tr,NDFT_tr) )
   Allocate ( AC_tr     (NDFT_tr,NDFT_tr) )
   Allocate ( AK_tr     (NDFT_tr,NDFT_tr) )
   Allocate ( AQ_tr     (NDFT_tr        ) )
   Allocate ( INDSYSB_tr(NDFT_tr        ) )


!--- Initial zeroing of displacement dofs
   UT_tr    = 0.d0;   UTP_tr   = 0.d0;   UTPP_tr  = 0.d0;
   UT1_tr   = 0.d0;   UTP1_tr  = 0.d0;   UTPP1_tr = 0.d0;
   UT2_tr   = 0.d0;   UTP2_tr  = 0.d0;   UTPP2_tr = 0.d0;


!--- Optional repositioning controller
   call init_repositioning(path)

   if (repos_active) call init_winch_elements()


!------ calculate total length
      ALENGP  = 0.d0
   do nb_tr   = 1, NBODT_tr
      n1      = body_tr(nb_tr)%INODE_tr(1)
      n2      = body_tr(nb_tr)%INODE_tr(2)
      DX(1:3) = XG_tr(1:3,n2)-XG_tr(1:3,n1)
      ALENGP  = ALENGP + dsqrt(DX(1)**2+DX(2)**2+DX(3)**2)
   enddo

   write(*,*)'Ltot',ALENGP


 END Subroutine INIT_tr
!----------------------------------------------------------------------
 Subroutine FINIT_tr
!----------------------------------------------------------------------
 use truss

   implicit none


   Deallocate ( NDFACC_tr, XG_tr )
   Deallocate ( body_tr      )
   Deallocate ( conc_mass_tr )
   Deallocate ( constrain_tr )
   Deallocate ( connect_tr   )
   Deallocate ( UT_tr  , UTP_tr  , UTPP_tr  , AM_tr ,&
                UT1_tr , UTP1_tr , UTPP1_tr , AC_tr ,&
                UT2_tr , UTP2_tr , UTPP2_tr , AK_tr ,&
                                              AQ_tr     )
   Deallocate ( INDSYSB_tr )

   if (allocated(t_repos   )) deallocate(t_repos)
    if (allocated(y_repos   )) deallocate(y_repos)
    if (allocated(ALENG0_tr )) deallocate(ALENG0_tr)
    if (allocated(ALENG_ctrl)) deallocate(ALENG_ctrl)
    if (allocated(ALENG_min )) deallocate(ALENG_min)
    if (allocated(ALENG_max )) deallocate(ALENG_max)
    if (allocated(winch_elem)) deallocate(winch_elem)


 END Subroutine FINIT_tr
!----------------------------------------------------------------------
 Subroutine MATRIX_tr
!----------------------------------------------------------------------
 use truss

   implicit none

   integer :: nb_tr, ni,nj, il,jl, i,j ,iloc,jloc, inod,jnod


!-- Local Matrices
   do nb_tr                   = 1,NBODT_tr
      body_tr(nb_tr)%AMLOC_tr = 0.d0;
      body_tr(nb_tr)%ACLOC_tr = 0.d0;
      body_tr(nb_tr)%AKLOC_tr = 0.d0;
      body_tr(nb_tr)%AQLOC_tr = 0.d0;
   enddo

   call INERTMAT_tr
   call STIFFMAT_tr
   call DAMPSMAT_tr    !- called after M,K for the Rayleigh damping
   call QFORCMAT_tr
   call CONCMAT_tr
   call SEABED_Int_tr  !- take into account the seabed interaction

   AM_tr = 0.d0;
   AC_tr = 0.d0;
   AK_tr = 0.d0;
   AQ_tr = 0.d0;

!-- Assembly local matrices to global 
   do nb_tr = 1, NBODT_tr

      do ni = 1, NNPE_tr
         inod  = body_tr(nb_tr)%INODE_tr (ni)

         do il = 1,NDFPN_tr
            iloc  = (ni-1)*NDFPN_tr    + il
            i     = NDFACC_tr (inod-1) + il

            do nj = 1,NNPE_tr
               jnod  = body_tr(nb_tr)%INODE_tr (nj)

               do jl = 1,NDFPN_tr

                  jloc = (nj-1)*NDFPN_tr    + jl
                  j    = NDFACC_tr (jnod-1) + jl

                  AM_tr(i,j) = AM_tr(i,j) + body_tr(nb_tr)%AMLOC_tr(iloc,jloc)
                  AC_tr(i,j) = AC_tr(i,j) + body_tr(nb_tr)%ACLOC_tr(iloc,jloc)
                  AK_tr(i,j) = AK_tr(i,j) + body_tr(nb_tr)%AKLOC_tr(iloc,jloc)

               enddo !jl
            enddo !nj

                  AQ_tr(i  ) = AQ_tr(i  ) + body_tr(nb_tr)%AQLOC_tr(iloc     )

         enddo !il
      enddo !ni
   enddo !nb_tr


 END Subroutine MATRIX_tr
!----------------------------------------------------------------------
 Subroutine SEABED_Int_tr
!----------------------------------------------------------------------

 use truss

   implicit none

   integer :: nb_tr, i, n, nod
   real(8) :: Zoff, bed_K, bed_C, damp_scale
   real(8) :: ZN(NNPE_tr), UN(NNPE_tr)
   real(8) :: k33,k36,k63,k66,f3,f6, zop1,zop2, L,L1,L2, dd
   real(8) :: c33,c36,c63,c66,       uz1 ,uz2


   if ( IYN_SeabedInt_tr /= 1 ) return

   Zoff          = 0.10d0  !-- distance from sea bed to simulate the seabed contact
   damp_scale    = 1.00d0  !-- scale damping compared to stiffness
   bed_K         = 1.d0/Zoff**2 * body_tr(1)%WaterFz_tr !nb_tr=1
   bed_C         = bed_K * damp_scale
   dd            = Depth_tr - Zoff

   do nb_tr      = 1, NBODT_tr
!------ find nodal displacements
      do nod     = 1, NNPE_tr
         n       = body_tr(nb_tr)%INODE_tr(nod)
         i       = NDFACC_tr (n-1) + 3                        !-- 3rd equation
!--------- Displacement of the node wrt the start of the spring at -depth+Zoff
         ZN(nod) = -(XG_tr(3,n) + UT_tr(i) + Depth_tr - Zoff) !-- positive upwards for Z>-depth+Zoff, so -
         UN(nod) = -UT1_tr(i)
      enddo !nod

      L          = body_tr(nb_tr)%ALENG_tr
      zop1       = -ZN(1) - dd
      zop2       = -ZN(2) - dd
      uz1        = -UN(1)
      uz2        = -UN(2)

      if    (ZN(1)<=0.d0.and.ZN(2)<=0.d0) then
         cycle                                                !no seabed contact
      elseif(ZN(1)>=0.d0.and.ZN(2)>=0.d0) then
         L1      = 0.d0                                       !both nodes of the element are inside the buffer area
         L2      = L
      elseif(ZN(2)<=0.d0                ) then
         L1      = 0.d0                                       !only the 1st node is inside the buffer area
         L2      = -ZN(1)*L/(ZN(2)-ZN(1))
      elseif(ZN(1)<=0.d0                ) then
         L1      = -ZN(1)*L/(ZN(2)-ZN(1))                     !only the 2nd node is inside the buffer area
         L2      = L
      endif
      
!------ K
      k33 =  (bed_K*(4*dd*L*(L1 - L2)*(3*L**2 + L1**2 + L1*L2 + L2**2 - 3*L*(L1 + L2)) + 12*L**3*(L1 - L2)*zop1 + 4*L*(L1**3 - L2**3)*(3*zop1 - 2*zop2)                                                   - &
              3*(L1**4 - L2**4)*(zop1 - zop2) - 6*L**2*(L1 - L2)*(L1 + L2)*(3*zop1 - zop2)))/(6.d0*L**3)                                                                                                  + &
             (bed_C*(12*L**3*(L1 - L2)*uz1 + 4*L*(L1**3 - L2**3)*(3*uz1 - 2*uz2) - 3*(L1**4 - L2**4)*(uz1 - uz2) - 6*L**2*(L1 - L2)*(L1 + L2)*(3*uz1 - uz2)))/(12.*L**3)
      k36 = -(bed_K*(-(L1**2*(6*L**2*(dd + zop1) + 3*L1**2*(zop1 - zop2) - 4*L*L1*(dd + 2*zop1 - zop2))) + L2**2*(6*L**2*(dd + zop1) + 3*L2**2*(zop1 - zop2) - 4*L*L2*(dd + 2*zop1 - zop2))))/(6.d0*L**3) + &
             (bed_C*(6*L**2*(L1 - L2)*(L1 + L2)*uz1 + 3*(L1**4 - L2**4)*(uz1 - uz2) - 4*L*(L1**3 - L2**3)*(2*uz1 - uz2)))/(12.*L**3)
      k63 =   k36
      k66 =  (bed_K*(4*dd*L*(L1**3 - L2**3) + 4*L*(L1**3 - L2**3)*zop1 - 3*(L1**4 - L2**4)*(zop1 - zop2)))/(6.d0*L**3)                                                                                    + &
             (bed_C*(4*L*(L1**3 - L2**3)*uz1 - 3*(L1**4 - L2**4)*(uz1 - uz2)))/(12.d0*L**3)
!------ C
      c33 =  (bed_C*(4*dd*L*(L1 - L2)*(3*L**2 + L1**2 + L1*L2 + L2**2 - 3*L*(L1 + L2)) + 12*L**3*(L1 - L2)*zop1 + 4*L*(L1**3 - L2**3)*(3*zop1 - 2*zop2) -                                                   &
               3*(L1**4 - L2**4)*(zop1 - zop2) - 6*L**2*(L1 - L2)*(L1 + L2)*(3*zop1 - zop2)))/(12.d0*L**3)
      c36 = -(bed_C*(-(L1**2*(6*L**2*(dd + zop1) + 3*L1**2*(zop1 - zop2) - 4*L*L1*(dd + 2*zop1 - zop2))) + L2**2*(6*L**2*(dd + zop1) + 3*L2**2*(zop1 - zop2) - 4*L*L2*(dd + 2*zop1 - zop2))))/(12.d0*L**3)
      c63 =   c36
      c66 =  (bed_C*(4*dd*L*(L1**3 - L2**3) + 4*L*(L1**3 - L2**3)*zop1 - 3*(L1**4 - L2**4)*(zop1 - zop2)))/(12.d0*L**3)
!------ Q
      f3  =  (bed_K*(-(L1*(12*L**3*(dd + zop1)**2 - 6*L**2*L1*(dd + zop1)*(dd + 3*zop1 - 2*zop2) - 3*L1**3*(zop1 - zop2)**2 + 4*L*L1**2*(zop1 - zop2)*(2*dd + 3*zop1 - zop2)))                                                   + &
              L2*(12*L**3*(dd + zop1)**2 - 6*L**2*L2*(dd + zop1)*(dd + 3*zop1 - 2*zop2) - 3*L2**3*(zop1 - zop2)**2 + 4*L*L2**2*(zop1 - zop2)*(2*dd + 3*zop1 - zop2))))/(12.d0*L**3)                                              + &
             (bed_C*(2*dd*L*(L1 - L2)*(-6*L**2*uz1 + 3*L*(L1 + L2)*(2*uz1 - uz2) + 2*(L1**2 + L1*L2 + L2**2)*(-uz1 + uz2)) + 12*L**3*(-L1 + L2)*uz1*zop1 + 3*(L1**4 - L2**4)*(uz1 - uz2)*(zop1 - zop2)                           + &
              6*L**2*(L1 - L2)*(L1 + L2)*(3*uz1*zop1 - uz2*zop1 - uz1*zop2) - 4*L*(L1**3 - L2**3)*(3*uz1*zop1 - 2*uz2*zop1 - 2*uz1*zop2 + uz2*zop2)))/(12.d0*L**3)
      f6  =  (bed_K*(-(L1**2*(6*L**2*(dd + zop1)**2 - 8*L*L1*(dd + zop1)*(zop1 - zop2) + 3*L1**2*(zop1 - zop2)**2)) + L2**2*(6*L**2*(dd + zop1)**2 - 8*L*L2*(dd + zop1)*(zop1 - zop2) + 3*L2**2*(zop1 - zop2)**2)))/(12.d0*L**3) + &
             (bed_C*(2*dd*L*(3*L*(-L1**2 + L2**2)*uz1 + 2*(L1**3 - L2**3)*(uz1 - uz2)) + 6*L**2*(-L1**2 + L2**2)*uz1*zop1 + 3*(L1**4 - L2**4)*(uz1 - uz2)*(-zop1 + zop2) + 4*L*(L1**3 - L2**3)*(2*uz1*zop1 - uz2*zop1 - uz1*zop2)))/(12.d0*L**3)


!------ Assembly to element matrices
      body_tr(nb_tr)%AKLOC_tr(3,3) = body_tr(nb_tr)%AKLOC_tr(3,3) + k33
      body_tr(nb_tr)%AKLOC_tr(3,6) = body_tr(nb_tr)%AKLOC_tr(3,6) + k36
      body_tr(nb_tr)%AKLOC_tr(6,3) = body_tr(nb_tr)%AKLOC_tr(6,3) + k63
      body_tr(nb_tr)%AKLOC_tr(6,6) = body_tr(nb_tr)%AKLOC_tr(6,6) + k66
      body_tr(nb_tr)%ACLOC_tr(3,3) = body_tr(nb_tr)%ACLOC_tr(3,3) + c33
      body_tr(nb_tr)%ACLOC_tr(3,6) = body_tr(nb_tr)%ACLOC_tr(3,6) + c36
      body_tr(nb_tr)%ACLOC_tr(6,3) = body_tr(nb_tr)%ACLOC_tr(6,3) + c63
      body_tr(nb_tr)%ACLOC_tr(6,6) = body_tr(nb_tr)%ACLOC_tr(6,6) + c66
      body_tr(nb_tr)%AQLOC_tr(3  ) = body_tr(nb_tr)%AQLOC_tr(3  ) + f3
      body_tr(nb_tr)%AQLOC_tr(6  ) = body_tr(nb_tr)%AQLOC_tr(6  ) + f6
   enddo !nb_tr


 END Subroutine SEABED_Int_tr
!----------------------------------------------------------------------
 Subroutine SEABED_Int_tr1
!----------------------------------------------------------------------

 use truss

   implicit none

   integer :: nb_tr,i, n,nod, iloc(NNPE_tr)
   real(8) :: Z2, bed_K0, bed_K, bed_C, ZN, damp_scale


   if ( IYN_SeabedInt_tr /= 1 ) return

   Z2         = 0.10d0  !-- distance from sea bed to simulate the seabed contact
   damp_scale = 1.00d0  !-- scale damping compared to stiffness
   iloc(1)    = 3
   iloc(2)    = 6


   do nb_tr = 1, NBODT_tr
      bed_K0 = 1.d0/Z2**2 * body_tr(nb_tr)%ALENG_tr * body_tr(nb_tr)%WaterFz_tr /2.d0

      do nod = 1, NNPE_tr
         n   = body_tr(nb_tr)%INODE_tr(nod)
         i   = NDFACC_tr (n-1) + 3                        !-- 3rd equation

!--- Displacement of the node wrt the start of the spring at -depth+z2
         ZN  = -(XG_tr(3,n) + UT_tr(i) + Depth_tr - Z2)   !-- positive upwards for Z>-depth+Z2, so -

         if (ZN<=0.d0) cycle

            bed_K   = bed_K0
            bed_C   = bed_K * damp_scale * Z2

         body_tr(nb_tr)%ACLOC_tr(iloc(nod),iloc(nod)) = body_tr(nb_tr)%ACLOC_tr(iloc(nod),iloc(nod)) + bed_C
!        body_tr(nb_tr)%AKLOC_tr(iloc(nod),iloc(nod)) = body_tr(nb_tr)%AKLOC_tr(iloc(nod),iloc(nod)) + bed_K
!        body_tr(nb_tr)%AQLOC_tr(iloc(nod)          ) = body_tr(nb_tr)%AQLOC_tr(iloc(nod)          ) + bed_K * ZN    - bed_C * UT1_tr(i)
         body_tr(nb_tr)%AKLOC_tr(iloc(nod),iloc(nod)) = body_tr(nb_tr)%AKLOC_tr(iloc(nod),iloc(nod)) + bed_K * 2d0 * ZN
         body_tr(nb_tr)%AQLOC_tr(iloc(nod)          ) = body_tr(nb_tr)%AQLOC_tr(iloc(nod)          ) + bed_K * ZN**2 - bed_C * UT1_tr(i)
      enddo

   enddo !nb_tr


 END Subroutine SEABED_Int_tr1
!----------------------------------------------------------------------
 Subroutine SEABED_Int_tr0
!----------------------------------------------------------------------

 use truss

   implicit none

   integer :: nb_tr,i, n,nod, iloc(NNPE_tr)
   real(8) :: Z1, Z2, bed_K0, bed_K, bed_C, ZN, rat, damp_scale


   if ( IYN_SeabedInt_tr /= 1 ) return

   Z2         = 0.10d0  !-- distance from sea bed to simulate the seabed contact
   Z1         = Z2/2.d0 !-- distance to linearly vary the stiffness
   damp_scale =  1.00d0 !-- scale damping compared to stiffness
   iloc(1)    = 3
   iloc(2)    = 6


   do nb_tr = 1, NBODT_tr
      bed_K0 = 1.d0/Z2    * body_tr(nb_tr)%ALENG_tr * body_tr(nb_tr)%WaterFz_tr /2.d0
!!    bed_K0 = 1.d0/Z2**2 * body_tr(nb_tr)%ALENG_tr * body_tr(nb_tr)%WaterFz_tr /2.d0

      do nod = 1, NNPE_tr
         n   = body_tr(nb_tr)%INODE_tr(nod)
         i   = NDFACC_tr (n-1) + 3                        !-- 3rd equation

!--- Displacement of the node wrt the start of the spring at -depth+z2
         ZN  = -(XG_tr(3,n) + UT_tr(i) + Depth_tr - Z2)   !-- positive upwards for Z>-depth+Z2, so -

         if (ZN<=0.d0) cycle

         if (ZN>=Z2-Z1) then
            bed_K   = bed_K0
            bed_C   = bed_K * damp_scale
         else
            rat     = ZN/(Z2-Z1)
            bed_K   = bed_K0 * rat
            bed_C   = bed_K * damp_scale
         endif

         body_tr(nb_tr)%ACLOC_tr(iloc(nod),iloc(nod)) = body_tr(nb_tr)%ACLOC_tr(iloc(nod),iloc(nod)) + bed_C
         body_tr(nb_tr)%AKLOC_tr(iloc(nod),iloc(nod)) = body_tr(nb_tr)%AKLOC_tr(iloc(nod),iloc(nod)) + bed_K
         body_tr(nb_tr)%AQLOC_tr(iloc(nod)          ) = body_tr(nb_tr)%AQLOC_tr(iloc(nod)          ) + bed_K * ZN    - bed_C * UT1_tr(i)
   !!    body_tr(nb_tr)%AKLOC_tr(iloc(nod),iloc(nod)) = body_tr(nb_tr)%AKLOC_tr(iloc(nod),iloc(nod)) + bed_K * 2d0 * ZN
   !!    body_tr(nb_tr)%AQLOC_tr(iloc(nod)          ) = body_tr(nb_tr)%AQLOC_tr(iloc(nod)          ) + bed_K * ZN**2 - bed_C * UT1_tr(i)
      enddo

   enddo !nb_tr


 END Subroutine SEABED_Int_tr0
!----------------------------------------------------------------------
 Subroutine BOUNDC_tr(IP)
!----------------------------------------------------------------------

 use truss

   implicit none

   integer :: nc, inod, ndf, il, i, j, IP


   do nc = 1, NCONSTR_tr 

      inod  = constrain_tr(nc)%INODC_tr
      ndf   = NDFACC_tr (inod-1)
         
      do il = 1, NDFPN_tr
          
         if (constrain_tr(nc)%INDXC_tr(il) /= 1) cycle

         i = ndf + il

         do j = 1, NDFT_tr
            AM_tr(i,j) = 0.d0
            AC_tr(i,j) = 0.d0
            AK_tr(i,j) = 0.d0
         enddo

         if (IP==0) then
            AK_tr(i,i) = 1.d0
            AQ_tr(i  ) = -UT_tr (i)
         else
            AM_tr(i,i) = 1.d0
            AQ_tr(i  ) = -UT2_tr(i)
         endif

      enddo !il
   enddo !nc


 END Subroutine BOUNDC_tr
!------------------------------------------------------------------
 Subroutine INERTMAT_tr
!-----------------------------------------------------------------------

 use truss

   implicit none

   real(8) :: D1, D2
   integer :: nb_tr,il,jl,ik,nj,j,nnod,i


   do nb_tr = 1, NBODT_tr
 
      D1 = body_tr(nb_tr)%DENS_tr * body_tr(nb_tr)%ALENG_tr / 3.d0        !-- Analytic integration
      D2 = D1/2.d0

      do i = 1, 3
         body_tr(nb_tr)%AMLOC_tr(  i,  i) = D1
         body_tr(nb_tr)%AMLOC_tr(3+i,3+i) = D1
         body_tr(nb_tr)%AMLOC_tr(  i,3+i) = D2
         body_tr(nb_tr)%AMLOC_tr(3+i,  i) = D2
      enddo


!-- Transfer the previous iteration value to the RHS
      do il = 1, NDFPE_tr

         jl = 0
         do ik = 1, NNPE_tr

            nnod  = body_tr(nb_tr)%INODE_tr(ik)

            do nj = 1, NDFPN_tr
               jl = jl + 1 
               j  = NDFACC_tr(nnod-1) + nj
               body_tr(nb_tr)%AQLOC_tr(il) =   body_tr(nb_tr)%AQLOC_tr(il)             &
                                             - body_tr(nb_tr)%AMLOC_tr(il,jl)*UT2_tr(j)
            enddo !nj
         enddo !ik
      enddo !il

   enddo !nb_tr


 END Subroutine INERTMAT_tr
!------------------------------------------------------------------
 Subroutine DAMPSMAT_tr
!-----------------------------------------------------------------------

 use truss

   implicit none

   real(8) :: AMAT(6,6), UNT3(3,3)
   real(8) :: TAUTAU(3,3), TAUDU (3,3), DU(3), DL(3),TAU(3), RLEN, Cdamp
   real(8) :: TERM1(3,6),TERM2(3,6),FORC(3)
   integer :: nb_tr,n1,n2,i1,i2,il,jl,ik,nj,j,nnod,i

!-- set AMAT, UNT3
   AMAT (:,:) =  0.d0;
   AMAT (1,1) = -1.d0
   AMAT (1,4) =  1.d0
   AMAT (2,2) = -1.d0
   AMAT (2,5) =  1.d0
   AMAT (3,3) = -1.d0
   AMAT (3,6) =  1.d0

   call DIAGO (3,UNT3)

!-- Rayleigh damping
   do nb_tr = 1, NBODT_tr
      Cdamp = body_tr(nb_tr)%CDAMP_tr

      n1    = body_tr(nb_tr)%INODE_tr(1)
      n2    = body_tr(nb_tr)%INODE_tr(2)
      i1    = NDFACC_tr (n1-1)
      i2    = NDFACC_tr (n2-1)

!-- Rayleigh damping added to previous damping
       body_tr(nb_tr)%ACLOC_tr(1:NDFPE_tr,1:NDFPE_tr) =   body_tr(nb_tr)%ACLOC_tr(1:NDFPE_tr,1:NDFPE_tr)            &
                                                        + body_tr(nb_tr)%AKLOC_tr(1:NDFPE_tr,1:NDFPE_tr) * COEFK_tr &
                                                        + body_tr(nb_tr)%AMLOC_tr(1:NDFPE_tr,1:NDFPE_tr) * COEFM_tr


!-- Transfer the previous iteration value to the RHS
      do il = 1, NDFPE_tr
        jl  = 0

        do ik = 1, NNPE_tr

           nnod = body_tr(nb_tr)%INODE_tr (ik)

           do nj = 1, NDFPN_tr
              jl = jl + 1 
              j  = NDFACC_tr(nnod-1) + nj
              body_tr(nb_tr)%AQLOC_tr(il) =  body_tr(nb_tr)%AQLOC_tr(il   ) &
                                           - body_tr(nb_tr)%ACLOC_tr(il,jl)*UT1_tr(j)
           enddo !nj
        enddo !ik
      enddo !il


!-- Katsa damping
!-- Find DL, L, Tau
      DL  (1:3)=  XG_tr(1:3, n2  ) -  XG_tr(1:3, n1  ) &
                 +UT_tr(i2+1:i2+3) -  UT_tr(i1+1:i1+3)
      DU  (1:3)= UT1_tr(i2+1:i2+3) - UT1_tr(i1+1:i1+3)

      RLEN     = dsqrt(dot_product(DL,DL))
      TAU(1:3) = DL(1:3)/RLEN

      do j = 1, 3
       do i = 1, 3
          TAUTAU(i,j) = TAU(i)*TAU(j)
          TAUDU (i,j) = TAU(i)*DU (j)
       enddo
      enddo

      TERM1(1:3,1:6) = Cdamp / RLEN * matmul(dot_product(TAU(1:3),DU(1:3)) * UNT3(1:3,1:3) + TAUDU(1:3,1:3), &
                                             matmul(UNT3(1:3,1:3)-TAUTAU(1:3,1:3), AMAT(1:3,1:6)) )
      TERM2(1:3,1:6) = Cdamp        * matmul( TAUTAU(1:3,1:3)    , AMAT(1:3,1:6) )
      FORC (1:3    ) = Cdamp        * matmul( TAUTAU(1:3,1:3), DU(1:3) )

      body_tr(nb_tr)%AKLOC_tr(1:3,1:6) = body_tr(nb_tr)%AKLOC_tr(1:3,1:6) - TERM1(1:3,1:6)
      body_tr(nb_tr)%ACLOC_tr(1:3,1:6) = body_tr(nb_tr)%ACLOC_tr(1:3,1:6) - TERM2(1:3,1:6)
      body_tr(nb_tr)%AQLOC_tr(1:3    ) = body_tr(nb_tr)%AQLOC_tr(1:3    ) + FORC (1:3    )
      body_tr(nb_tr)%AKLOC_tr(4:6,1:6) = body_tr(nb_tr)%AKLOC_tr(4:6,1:6) + TERM1(1:3,1:6)
      body_tr(nb_tr)%ACLOC_tr(4:6,1:6) = body_tr(nb_tr)%ACLOC_tr(4:6,1:6) + TERM2(1:3,1:6)
      body_tr(nb_tr)%AQLOC_tr(4:6    ) = body_tr(nb_tr)%AQLOC_tr(4:6    ) - FORC (1:3    )
   enddo !nb_tr
     

 END Subroutine DAMPSMAT_tr
!------------------------------------------------------------------
 Subroutine STIFFMAT_tr
!-----------------------------------------------------------------------

 use truss

   implicit none

   real(8) :: AMAT(3,6),Trans_AMAT(6,3), UNT3(3,3), DL(3)    , TAU (3)  , TAUTAU(3,3)
   real(8) ::            TENSION  (6)  , TEM1(3,3), TEM2(3,3), STFM(6,6)
   real(8) :: ALENG, EA, RLEN, STRAIN !!, STRAIN1, STRAIN2, rat
   integer :: nb_tr,n1,n2,i1,i2,i,j


!-- set AMAT, Transpose(AMAT), UNT3
   AMAT (:,:) =  0.d0;
   AMAT (1,1) = -1.d0
   AMAT (1,4) =  1.d0
   AMAT (2,2) = -1.d0
   AMAT (2,5) =  1.d0
   AMAT (3,3) = -1.d0
   AMAT (3,6) =  1.d0
   Trans_AMAT = transpose(AMAT)

   call DIAGO (3,UNT3)

   do nb_tr = 1, NBODT_tr

      ALENG = body_tr(nb_tr)%ALENG_tr

      n1 = body_tr(nb_tr)%INODE_tr (1)
      n2 = body_tr(nb_tr)%INODE_tr (2)
      i1 = NDFACC_tr (n1-1)
      i2 = NDFACC_tr (n2-1)

!-- Find DL, L, Tau
      DL  (1:3)= XG_tr(1:3, n2  ) - XG_tr(1:3, n1  ) &
                +UT_tr(i2+1:i2+3) - UT_tr(i1+1:i1+3)

      RLEN     = dsqrt(dot_product(DL,DL))
      TAU(1:3) = DL(1:3)/RLEN

!-- Find Strain, EA
      STRAIN   = RLEN/ALENG - 1.d0
  
      EA       = body_tr(nb_tr)%EAS2_tr
      if (STRAIN < 0d0) &
      EA       = body_tr(nb_tr)%EAS1_tr

!!    STRAIN1    = 0.000d0
!!    STRAIN2    = 0.001d0

!!    if     (STRAIN > STRAIN2) then
!!       EA  = body_tr(nb_tr)%EAS2_tr
!!    elseif (STRAIN < STRAIN1) then
!!       EA  = body_tr(nb_tr)%EAS1_tr
!!    else
!!       rat = (STRAIN-STRAIN1)/(STRAIN2-STRAIN1)
!!       EA  = body_tr(nb_tr)%EAS1_tr*(1.d0-rat) &
!!            +body_tr(nb_tr)%EAS2_tr*(     rat)
!!    endif

!-- Find Tension
      TENSION(1:6) = EA*STRAIN*matmul(trans_AMAT(1:6,1:3),TAU(1:3))

!-- Find Stiffness matrix
      do j = 1, 3
       do i = 1, 3
          TAUTAU(i,j)=TAU(i)*TAU(j)
       enddo
      enddo

      TEM1(1:3,1:3) = EA/RLEN*STRAIN*(UNT3(1:3,1:3)-TAUTAU(1:3,1:3));
      TEM2(1:3,1:3) = EA/ALENG      *               TAUTAU(1:3,1:3) ;
      STFM(1:6,1:6) = matmul(trans_AMAT(1:6,1:3),matmul(TEM1(1:3,1:3)+TEM2(1:3,1:3), AMAT(1:3,1:6)))

      body_tr(nb_tr)%AKLOC_tr(1:NDFPE_tr,1:NDFPE_tr) =  body_tr(nb_tr)%AKLOC_tr(1:NDFPE_tr,1:NDFPE_tr) &
                                                       +               STFM    (1:NDFPE_tr,1:NDFPE_tr)
      body_tr(nb_tr)%AQLOC_tr(1:NDFPE_tr           ) =  body_tr(nb_tr)%AQLOC_tr(1:NDFPE_tr           ) &
                                                       -               TENSION (1:NDFPE_tr           )
   enddo !nb_tr


 END Subroutine STIFFMAT_tr
!------------------------------------------------------------------
 Subroutine STIFFMAT_tr_OLD
!-----------------------------------------------------------------------

 use truss

   implicit none

   real(8) :: a1,a2
   integer :: nb_tr,n1,n2,i1,i2,il,jl,i
   real(8) :: AMAT(6,6), C(6), C0(6),    &
              X_XP(6)  , X(6), XP(6), EA, TRACT, ALENG
   real(8) :: STRAIN, STRAIN1, STRAIN2, rat, ALENG_N_2, ALENG_0_2


   do nb_tr=1,NBODT_tr

      ALENG = body_tr(nb_tr)%ALENG_tr

      n1 = body_tr(nb_tr)%INODE_tr (1)
      n2 = body_tr(nb_tr)%INODE_tr (2)
      i1 = NDFACC_tr (n1-1)
      i2 = NDFACC_tr (n2-1)

!(3.57 or 3.141) different variables order. Here is (x1,y1,z1,x2,y2,z2) and not (x1,x2,y1,y2,z1,2) as in book
      AMAT = 0.d0;
      do i=1,3
         AMAT (  i,  i) = 0.25d0
         AMAT (3+i,3+i) = 0.25d0
         AMAT (  i,3+i) =-0.25d0
         AMAT (3+i,  i) =-0.25d0
      enddo

!(3.48)
      X    (1:3) = XG_tr(1:3,n1)  
      X    (4:6) = XG_tr(1:3,n2) 
!(3.49)
      XP   (1:3) = UT_tr(i1+1:i1+3)
      XP   (4:6) = UT_tr(i2+1:i2+3)
!(3.47)
      X_XP (1:6) = X (1:6) + XP(1:6)
!(3.56)
      C0        = 4d0 * matmul (AMAT,X   ) !CMAT0
!(3.142)
      C         = 4d0 * matmul (AMAT,X_XP) !CMAT

!!    do i=1,3
!!       C0 (i)   =  X   (i) - X   (i+3)
!!       C  (i)   =  X_XP(i) - X_XP(i+3)
!!       C0 (i+3) = -C0  (i)
!!       C  (i+3) = -C   (i)
!!    enddo


!(3.8) Green's strain e = (ln^2-lo^2)/(2 Lo^2)
   ALENG_N_2  = C(1)**2+C(2)**2+C(3)**2
   ALENG_0_2  = body_tr(nb_tr)%ALENG_tr**2
   STRAIN     = (ALENG_N_2-ALENG_0_2) / (2.d0*ALENG_0_2)
   STRAIN1    = 0.000d0
   STRAIN2    = 0.001d0

   if     (STRAIN < STRAIN1) then
      EA  = body_tr(nb_tr)%EAS1_tr
   elseif (STRAIN > STRAIN2) then
      EA  = body_tr(nb_tr)%EAS2_tr
   else
      rat = (STRAIN-STRAIN1)/(STRAIN2-STRAIN1)
      EA  = body_tr(nb_tr)%EAS1_tr*(1.d0-rat) &
           +body_tr(nb_tr)%EAS2_tr*(     rat)
   endif

!(2.8)
   TRACT = EA*STRAIN

!(3.64) Kt = Kt1 + Kt2 + Ktσ
!(3.69) Ktσ
      a1 = TRACT / ALENG * 4d0
      body_tr(nb_tr)%AKLOC_tr(1:NDFPE_tr,1:NDFPE_tr) = a1 * AMAT(1:NDFPE_tr,1:NDFPE_tr)

!(3.93, 3.94) Kt1' = Kt1 + Kt2a + (Kt2)T + Kt2b
      a2 = EA / ALENG**3
      do il=1,NDFPE_tr
       do jl=1,NDFPE_tr
         body_tr(nb_tr)%AKLOC_tr(il,jl) =  body_tr(nb_tr)%AKLOC_tr(il,jl) + a2*C(il)*C(jl)
       enddo
      enddo

!(3.92)
      do il=1,NDFPE_tr
         body_tr(nb_tr)%AQLOC_tr(il) = body_tr(nb_tr)%AQLOC_tr(il) - TRACT / ALENG * C(il)
      enddo

   enddo !nb_tr


 END Subroutine STIFFMAT_tr_OLD
!------------------------------------------------------------------
!     Subroutine: QFORCMAT_tr
!
!     distributed (mooring elements)
!     1. gravitational loads
!     2. buoyancy
!     3. normal  drag
!     4. tangent drag
!     5. accelaration loads
!
!-----------------------------------------------------------------------
 Subroutine QFORCMAT_tr
!-----------------------------------------------------------------------

 use truss

   implicit none


   real(8) :: ASECT0,AQ,DX(3),ALENGP,RE_norm ,&
              TT_Meas, UREL_n_Meas, &
              CMASS,CDRAG,CDNORM,CDTANG, Ca
   real(8) :: SHAPE      (3,6), UINF  (3), AINF  (3), TT       (3),&
              SHAPET     (6,3), UREL  (3), AREL  (3), TT_prod(3,3),&
              SHAPETAMat (6,3), UREL_t(3), UREL_n(3), XDMAT    (6),&
              AMat0      (6,6), AINF_n(3), AREL_n(3), X0       (3),&
              AF0        (6  ), AC  (3,3), AM  (3,3), F        (3), pdyn
   real(8) :: UT      (NDFPE_tr )         ,&
              UT1     (NDFPE_tr )         ,&
              UT2     (NDFPE_tr )         ,&
!!            UTT     (NEQPE_tr )         ,&
              UTT1    (NEQPE_tr )         ,&
              UTT2    (NEQPE_tr )
!!            SHAPE   (NEQPE_tr,NDFPE_tr) ,&
!!            SHAPET  (NDFPE_tr,NEQPE_tr)
   real(8) :: HTAL, PA01, PA02
   integer :: nb_tr,i,j ,n,nod !!n1,n2,i1,i2,
   integer :: NX_hyd, k
   integer :: i1(NNPE_tr),i2(NNPE_tr)


   NX_hyd = 1
   i1(1)  = 1
   i2(1)  = 3
   i1(2)  = 4
   i2(2)  = 6

   do nb_tr = 1, NBODT_tr

!--- Gravitational loads and bouyancy (Analytic integration)
                           AQ    = -body_tr(nb_tr)%WaterFz_tr  *  body_tr(nb_tr)%ALENG_tr / 2.d0
      body_tr(nb_tr)%AQLOC_tr(3) =  body_tr(nb_tr)%AQLOC_tr(3) + AQ
      body_tr(nb_tr)%AQLOC_tr(6) =  body_tr(nb_tr)%AQLOC_tr(6) + AQ


!--- Morison's equation Hydrodynamic loads
      if ( IYN_Morison_tr == 0 ) cycle


!--- Coordinates of element nodes
      do nod = 1, NNPE_tr
         n   = body_tr(nb_tr)%INODE_tr(nod)
         i   = NDFACC_tr (n-1)
         UT    (i1(nod):i2(nod)) = UT_tr  (i+1:i+3)
         UT1   (i1(nod):i2(nod)) = UT1_tr (i+1:i+3)
         UT2   (i1(nod):i2(nod)) = UT2_tr (i+1:i+3)
         XDMAT (i1(nod):i2(nod)) = XG_tr(1:3,n) + UT(i1(nod):i2(nod))
      enddo


!--- Surface of the element
      ASECT0      = pi_tr*body_tr(nb_tr)%ADIAM_tr**2/4.d0 !*(body_tr(nb_tr)%ADIAMEQ / body_tr(nb_tr)%ADIAM_tr)

!--- Unit vector in direction tangent to the element
      TT (1:3)    = XDMAT (4:6) - XDMAT (1:3)
      TT_Meas     = dsqrt (dot_product(TT,TT))
      TT (1:3)    = TT (1:3)/TT_Meas

      do j=1,3
       do i=1,3
          TT_prod(i,j) = TT(i) * TT(j)
       enddo
      enddo

!--- Current Length
      DX (1:3)    = XDMAT(4:6)-XDMAT(1:3);
      ALENGP      = dsqrt(dot_product(DX,DX)) / dble(NX_hyd)


      do k = 1, NX_hyd

         HTAL     = dble(2*k-1)/2.d0 / dble(NX_hyd) !* DL/ALLOC

         PA01     = 1.d0 - HTAL
         PA02     = HTAL

!--- Coordinates of the middle of current section
         X0 (1:3) = PA01 * XDMAT(1:3) + PA02 * XDMAT(4:6)

!--- Shape function 
         call SHAPEFUNC6 ( HTAL, SHAPE, NDFPE_tr, NEQPE_tr )
         SHAPET   = transpose (SHAPE)
       !!UTT      = matmul (SHAPE, UT )
         UTT1     = matmul (SHAPE, UT1)
         UTT2     = matmul (SHAPE, UT2)


!--- GLOBAL U, A due to wave / current
         call UINFLOW_sea ( X0, UINF, AINF, pdyn, 1 )

!!       UREL  (1:3) = UINF(1:3) - 0.5d0*(UT1_tr(i1+1:i1+3)+UT1_tr(i2+1:i2+3))
         UREL  (1:3) = UINF(1:3) - UTT1(1:3)
         UREL_t(1:3) = dot_product(UREL(1:3),TT(1:3)) * TT(1:3)
         UREL_n(1:3) = UREL(1:3) - UREL_t(1:3)
         UREL_n_Meas = dsqrt ( dot_product ( UREL_n, UREL_n ) )
         RE_norm     = AROW_tr * UREL_n_Meas * body_tr(nb_tr)%ADIAM_tr / AMVISCW_tr


!--- Normal Drag Force 
!!       call CDRAGNORM  (RE_norm, CDNORM, pi_tr)
         call DIAGO (3,AC)

         CDNORM  = body_tr(nb_tr)%Cdnorm_tr
         CDRAG   = CDNORM * AROW_tr/2.d0 * body_tr(nb_tr)%ADIAM_tr * UREL_n_Meas * ALENGP
         AC      = CDRAG  * (AC - TT_prod);
         F (1:3) = CDRAG  * UREL_n(1:3)

         SHAPETAMat                       = matmul( SHAPET    , AC    )
         AMat0                            = matmul( SHAPETAMat, SHAPE )
         AF0                              = matmul( SHAPET    , F     )
         body_tr(nb_tr)%ACLOC_tr(1:6,1:6) = body_tr(nb_tr)%ACLOC_tr(1:6,1:6) + AMat0(1:6,1:6);
         body_tr(nb_tr)%AQLOC_tr(1:6    ) = body_tr(nb_tr)%AQLOC_tr(1:6    ) + AF0  (1:6    );


!--- Tangential Drag Force
!!       call CDRAGTANG  (RE_norm, AMVISCW_tr, CDTANG, pi_tr)
         call DIAGO (3,AC)

         CDTANG  = body_tr(nb_tr)%Cdtang_tr
         CDRAG   = CDTANG * ALENGP
         AC      = CDRAG  * TT_prod;
         F (1:3) = CDRAG  * UREL_t(1:3)

         SHAPETAMat                       = matmul( SHAPET    , AC    )
         AMat0                            = matmul( SHAPETAMat, SHAPE )
         AF0                              = matmul( SHAPET    , F     )
         body_tr(nb_tr)%ACLOC_tr(1:6,1:6) = body_tr(nb_tr)%ACLOC_tr(1:6,1:6) + AMat0(1:6,1:6);
         body_tr(nb_tr)%AQLOC_tr(1:6    ) = body_tr(nb_tr)%AQLOC_tr(1:6    ) + AF0  (1:6    );


!--- Added mass component 1     Ca * ρ * Vo * Arel_n, Ca=Cm-1
         AINF_n(1:3) = AINF(1:3) - dot_product(AINF(1:3),TT(1:3)) * TT(1:3)
!!!      AREL  (1:3) = AINF(1:3) - 0.5d0*(UT2_tr(i1+1:i1+3)+UT2_tr(i2+1:i2+3))
         AREL  (1:3) = AINF(1:3) - UTT2(1:3)
         AREL_n(1:3) = AREL(1:3) - dot_product(AREL(1:3),TT(1:3)) * TT(1:3)

         Ca          = body_tr(nb_tr)%Ca_tr
         CMASS       = Ca * AROW_tr*ASECT0*ALENGP

         call DIAGO (3,AM)
         AM          = CMASS * (AM - TT_prod);
         F (1:3)     = CMASS * AREL_n(1:3)

         SHAPETAMat                       = matmul( SHAPET    , AM    )
         AMat0                            = matmul( SHAPETAMat, SHAPE )
         AF0                              = matmul( SHAPET    , F     )
         body_tr(nb_tr)%AMLOC_tr(1:6,1:6) = body_tr(nb_tr)%AMLOC_tr(1:6,1:6) + AMat0(1:6,1:6);
         body_tr(nb_tr)%AQLOC_tr(1:6    ) = body_tr(nb_tr)%AQLOC_tr(1:6    ) + AF0  (1:6    );


!--- Added mass component 2    (not relative) ρ * Vo * A_n
         CMASS       = AROW_tr*ASECT0*ALENGP
         F (1:3)     = CMASS * AINF_n(1:3)

         AF0                              = matmul( SHAPET    , F     )
         body_tr(nb_tr)%AQLOC_tr(1:6    ) = body_tr(nb_tr)%AQLOC_tr(1:6    ) + AF0  (1:6    );

      enddo !k

   enddo !nb_tr


 END Subroutine QFORCMAT_tr
!------------------------------------------------------------------
!
!  A.Extra Inertia
!
!  B.Extra Forces
!     1. gravitational loads
!     2. buoyancy               
!     3. drag
!     4. accelaration loads     
!
!------------------------------------------------------------------
 Subroutine CONCMAT_tr
!-----------------------------------------------------------------------

 use truss

   implicit none

   real(8) :: AMASS,AMASSW, DMET, X0(3),URELMET, ASURF, AVOLU,RE_norm,CDCON, &
              CDRAG,CMASS1,CMASS2, Ca
   real(8) :: UREL(3), AREL(3), UINF(3), AINF(3), pdyn
   integer :: icon,nbod,inod,nnod,ii,i,n,iloc


!--Inertia
   do icon  = 1, NCON_tr

      nbod  = conc_mass_tr(icon)%IBCON_tr
      inod  = conc_mass_tr(icon)%INCON_tr
      AMASS = conc_mass_tr(icon)%AMASSCON_tr

      nnod  = body_tr(nbod)%INODE_tr(inod)
      ii    = NDFACC_tr (nnod-1)
      iloc  = NDFPN_tr* (inod-1)

      do n  = 1, 3
         i  = iloc+n
         body_tr(nbod)%AMLOC_tr(i,i) = body_tr(nbod)%AMLOC_tr(i,i) + AMASS
         body_tr(nbod)%AQLOC_tr(i  ) = body_tr(nbod)%AQLOC_tr(i  ) - AMASS*UT2_tr(ii+n)
      enddo

   enddo !icon


!--Forces
   do icon      = 1, NCON_tr

      nbod      = conc_mass_tr(icon)%IBCON_tr
      inod      = conc_mass_tr(icon)%INCON_tr
      DMET      = conc_mass_tr(icon)%DIAMCON_tr
      AMASS     = conc_mass_tr(icon)%AMASSCON_tr

      nnod      = body_tr(nbod)%INODE_tr(inod) 
      ii        = NDFACC_tr (nnod-1)
      iloc      = NDFPN_tr* (inod-1)

      ASURF     = 1.d0/4.d0*pi_tr* DMET**2
      AVOLU     = 4.d0/3.d0*pi_tr*(DMET/2.d0)**3
      AMASSW    = AROW_tr*AVOLU


!--- Buoyancy and Gravity Forces
      i         = iloc+3
      body_tr(nbod)%AQLOC_tr(i) = body_tr(nbod)%AQLOC_tr(i) - AMASS  * GRAV_tr  &
                                                            + AMASSW * GRAV_tr

!--- Morison's Equation Forces
      if ( IYN_Morison_tr == 0 ) cycle

      X0(1:3)   = XG_tr(1:3,nnod) + UT_tr(ii+1:ii+3) 

      call UINFLOW_sea ( X0, UINF, AINF, pdyn, 1 )

      UREL(1:3) = UINF(1:3) - UT1_tr(ii+1:ii+3)
      AREL(1:3) = AINF(1:3) - UT2_tr(ii+1:ii+3)
      URELMET   = dsqrt (UREL(1)**2+UREL(2)**2+UREL(3)**2)

!--drag coefficient of a sphere
      RE_norm   = AROW_tr*URELMET*DMET/AMVISCW_tr

      call SPHERE_CDRAG (RE_norm, CDCON)

      open (1,file='re.dat',access='append')
        write (1,*) TIME_tr,RE_norm,CDCON
      close (1)

      CDCON     = 0.5d0
      CDRAG     = 0.5d0*AROW_tr*CDCON*ASURF*URELMET

!--added mass coefficient of a sphere
      Ca        = 0.5d0
      CMASS1    = Ca  * AROW_tr * AVOLU
      CMASS2    =       AROW_tr * AVOLU

      do n      = 1, 3
         i      = iloc+n

         body_tr(nbod)%ACLOC_tr(i,i) = body_tr(nbod)%ACLOC_tr(i,i) + CDRAG
         body_tr(nbod)%AMLOC_tr(i,i) = body_tr(nbod)%AMLOC_tr(i,i) + CMASS1
         body_tr(nbod)%AQLOC_tr(i  ) = body_tr(nbod)%AQLOC_tr(i  ) + CDRAG  * UREL(n) &
                                                                   + CMASS1 * AREL(n) &
                                                                   + CMASS2 * AINF(n)
      enddo
   enddo !icon


 END Subroutine CONCMAT_tr
!--------------------------------------------------------------------------------
!
!  Subroutine : SPHERE CDRAG  ------------------
!
!  Normal drag coefficient               
! 
!--------------------------------------------------------------------------------
 Subroutine SPHERE_CDRAG (REN, Sphere_CD)
!-----------------------------------------------------------------------

   implicit none

   real(8) :: REN, Sphere_CD, REN5, REN263


   if     (REN<=0.1d0) then
      Sphere_CD = 0.d0
   elseif (REN<=2.d0) then
      Sphere_CD = 24d0/REN
   elseif (REN<=1.d6) then
      REN5   = REN/5d0
      REN263 = REN/263000d0

      Sphere_CD = 24d0/REN + 2.6d0*REN5/(1d0+REN5**1.52d0)                &
                           + 0.411d0*REN263**(-7.94d0)/(1+REN263**(-8d0)) &
                           + REN**0.8d0/461000d0
   else
      write(*,*)'Re>10**6'
      stop      
   endif


 END Subroutine SPHERE_CDRAG
!--------------------------------------------------------------------------------
!
!  Subroutine : CDRAGNORM  ------------------
!
!  Normal drag coefficient               
! 
!--------------------------------------------------------------------------------
 Subroutine CDRAGNORM (REN, CDNORM, pi)
!-----------------------------------------------------------------------

   implicit none

   real(8) :: REN, CDNORM, SS, pi


   if     (REN == 0.d0) then
      CDNORM = 0.d0
   elseif (REN.le.1.d0) then
      pi     = dacos(-1d0)
      SS     = -0.077215665d0 + dlog(8.d0/REN)
      CDNORM = 8.d0*pi/(REN*SS)*(1.d0-0.87d0/SS**2)
   elseif (REN.gt.1.d0.and.REN.lt.30.d0) then 
      CDNORM = 1.45d0 + 8.55d0*REN**(-0.9d0)
   else
      CDNORM = 1.1d0 + 4.d0*REN**(-0.5d0)
   endif


 END Subroutine CDRAGNORM
!--------------------------------------------------------------------------------
!
!  Subroutine : CDRAGTANG  ------------------
!
!  Tangential drag coefficient               
! 
!--------------------------------------------------------------------------------
 Subroutine CDRAGTANG (REN, AMVISCW, CDTANG, pi)
!--------------------------------------------------------------------------------

   implicit none

   real(8) :: REN, AMVISCW, CDTANG, pi


   CDTANG = pi*AMVISCW*(0.55*REN**0.5d0 + 0.084d0*REN**(2.d0/3.d0))


 END Subroutine CDRAGTANG
!--------------------------------------------------------------------------------
!
!  Subroutine :WRITELEM_tr         -------------
!
!--------------------------------------------------------------------------------
 Subroutine WRITELEM_tr (NTIME)
!-----------------------------------------------------------------------

 use truss

   implicit none

   real(8)      :: R1(3), R2(3), RR(3), fi, Tens1, Tens2, Tens1hor, Tens2hor
   integer      :: NTIME,nb_tr,n1,n2,i1,i2,nbod,inod,ndf,nc,iloc,nod,j, IWR
   character*1  :: CNUM(5),CNUM2(2)
   character*80 :: outfil, outfil2


!- Write the position of the truss elements every timestep
!- Write the forces of each truss element

                                            IWR=0
   if (IWRITE_truss==1 .or. TIME_tr<=DT_tr) IWR=1

   call INT_2_CHAR ( 5, CNUM, NTIME )

#ifndef ASCII
#define ASCII 1
#endif

   if (IWR==1) then
#if   ASCII == 0
      outfil = 'elem'//CNUM(1)//CNUM(2)//CNUM(3)//CNUM(4)//CNUM(5)//'.bin'
      open (11,file=outfil,form='UNFORMATTED')
#elif ASCII == 1
      outfil = 'elem'//CNUM(1)//CNUM(2)//CNUM(3)//CNUM(4)//CNUM(5)//'.dat'
      open (11,file=outfil)
#endif
   endif

   do nb_tr    = 1, NBODT_tr
      n1       = body_tr(nb_tr)%INODE_tr(1)
      n2       = body_tr(nb_tr)%INODE_tr(2)
      i1       = NDFACC_tr (n1-1)
      i2       = NDFACC_tr (n2-1)
      R1 (1:3) = XG_tr(1:3,n1) + UT_tr(i1+1:i1+3)
      R2 (1:3) = XG_tr(1:3,n2) + UT_tr(i2+1:i2+3)
      RR (1:3) = R2(1:3)-R1(1:3)
      fi       = dacos( dsqrt( (RR(1)**2+RR(2)**2) / (RR(1)**2+RR(2)**2+RR(3)**2) ) ) * (180.d0/pi_tr)
      Tens1    = dsqrt(body_tr(nb_tr)%AQLOC_tr(1)**2 + body_tr(nb_tr)%AQLOC_tr(2)**2 + body_tr(nb_tr)%AQLOC_tr(3)**2)
      Tens2    = dsqrt(body_tr(nb_tr)%AQLOC_tr(4)**2 + body_tr(nb_tr)%AQLOC_tr(5)**2 + body_tr(nb_tr)%AQLOC_tr(6)**2)
      Tens1hor = dsqrt(body_tr(nb_tr)%AQLOC_tr(1)**2 + body_tr(nb_tr)%AQLOC_tr(2)**2)
      Tens2hor = dsqrt(body_tr(nb_tr)%AQLOC_tr(4)**2 + body_tr(nb_tr)%AQLOC_tr(5)**2)

      if (IWR==1) then
#if   ASCII == 0
         write (11    ) (sngl(R1(j)),j=1,3)
         write (11    ) (sngl(R2(j)),j=1,3)
#elif ASCII == 1
         write (11,101) (     R1(j) ,j=1,3)
         write (11,101) (     R2(j) ,j=1,3)
#endif
      endif

      call INT_2_CHAR ( 2, CNUM2, nb_tr )

#if   ASCII == 0
      outfil2 = 'bodytr_'//CNUM2(1)//CNUM2(2)//'.bin'
      open (12,file=outfil2,access='append',form='UNFORMATTED' )
      write(12)                                       &
#elif ASCII == 1
      outfil2 = 'bodytr_'//CNUM2(1)//CNUM2(2)//'.dat'
      open (12,file=outfil2,access='append')
      write (12,100)                                  &
#endif
       sngl(TIME_tr)                                 ,& ! 1
      (sngl(-body_tr(nb_tr)%AQLOC_tr(iloc)),iloc=1,6),& ! 2-7
       sngl(Tens1)                                   ,& ! 8    nod1
       sngl(Tens2)                                   ,& ! 9    nod2
       sngl(fi)                                      ,& !10
       sngl(Tens1hor)                                ,& !11    nod1
       sngl(Tens2hor)                                   !12    nod2
      close(12)
   enddo !nb_tr

   if (IWR==1) &
   close (11)

RETURN


!--- Write position, velocity and acceleration of each Connection
#if   ASCII == 0
   outfil = 'pos_connect'//CNUM(1)//CNUM(2)//CNUM(3)//CNUM(4)//CNUM(5)//'.bin'
   open (12,file=outfil,form='UNFORMATTED' )
#elif ASCII == 1
   outfil = 'pos_connect'//CNUM(1)//CNUM(2)//CNUM(3)//CNUM(4)//CNUM(5)//'.dat'
   open (12,file=outfil)
#endif

   do nc      = 1, NCONNECT_tr
      nbod    = connect_tr(nc)%IBODCONNECT_tr
      nod     = connect_tr(nc)%INODCONNECT_tr
      inod    = body_tr(nbod)%INODE_tr(nod)
      ndf     = NDFACC_tr (inod -1)
      R2(1:3) = XG_tr(1:3,inod) + UT_tr(ndf+1:ndf+3)

      call INT_2_CHAR ( 2, CNUM2, nc )

#if   ASCII == 0
      outfil = 'deform_connect'//CNUM2(1)//CNUM2(2)//'.bin'
      open  (11,file=outfil,access='append',form='UNFORMATTED' )
      write (11)                      &
#elif ASCII == 1
      outfil = 'deform_connect'//CNUM2(1)//CNUM2(2)//'.dat'
      open  (11,file=outfil,access='append')
      write (11,100)                  &
#endif
       sngl(TIME_tr      )           ,& ! 1
       sngl(R2(1)        )           ,& ! 2
       sngl(R2(2)        )           ,& ! 3
       sngl(R2(3)        )           ,& ! 4
       sngl(UT1_tr(ndf+1))           ,& ! 5
       sngl(UT1_tr(ndf+2))           ,& ! 6
       sngl(UT1_tr(ndf+3))           ,& ! 7
       sngl(UT2_tr(ndf+1))           ,& ! 8
       sngl(UT2_tr(ndf+2))           ,& ! 9
       sngl(UT2_tr(ndf+3))           ,& !10
       sngl(connect_tr(nc)%ALOADB(1)),& !11
       sngl(connect_tr(nc)%ALOADB(2)),& !12
       sngl(connect_tr(nc)%ALOADB(3))   !13
      close (11)

#if   ASCII == 0
      write (12)                                &
#elif ASCII == 1
      write (12,100)                            &
#endif
       sngl(connect_tr(nc)% UT_CONNECT_tr(1:3)),&
       sngl(connect_tr(nc)%UT1_CONNECT_tr(1:3)),&
       sngl(connect_tr(nc)%UT2_CONNECT_tr(1:3))
   enddo !nc

   close(12)

 100  format (150f15.5)
 101  format (150f15.8)


 END Subroutine WRITELEM_tr
!----------------------------------------------------------------------
 Subroutine Time_Integrate_tr (it, itype, ICONV, ERR, RLX)
!Subroutine Time_Integrate_tr (AM,AC,AK,AQ,INDSYSB,NDFT,NFDT0, &
!                              UT,UT1,UT2,UTP,UTP1,UTP2,       &
!                              TIME,DT,BITA,GAMMA, RLXs,       &
!                              it, itype, ICONV, ERR)
!----------------------------------------------------------------------

 use truss

   implicit None

   real(8) :: UPRE  (NDFT0_tr)
   real(8) :: UTPRE (NDFT0_tr)
   integer :: Ipivot(NDFT0_tr)
   real(8) :: c1, c2, c3, c4, c5, RESDT,RESDTM, ERR, RLX
   integer :: it, itype, ICONV, i, is, info(2), N, imax


!-- Form the system's matrices, depending on the solution type.
   if     (itype == 2) then
              !-- Newmark-b method
                   c1 = (0.5d0 - BITA_tr )  * DT_tr**2
                   c2 = (1.0d0 - GAMMA_tr)  * DT_tr
                   c3 =  1.d0    / (BITA_tr * DT_tr**2)
                   c4 = GAMMA_tr / (BITA_tr * DT_tr   )
                   c5 = GAMMA_tr            * DT_tr
                   N  = NDFT0_tr

      UPRE  (1:N    ) = UTP_tr (INDSYSB_tr(1:N))  +  UTP1_tr(INDSYSB_tr(1:N))*DT_tr  +  UTP2_tr(INDSYSB_tr(1:N)) * c1
      UTPRE (1:N    ) = UTP1_tr(INDSYSB_tr(1:N))                                     +  UTP2_tr(INDSYSB_tr(1:N)) * c2
      AK_tr (1:N,1:N) = AK_tr(1:N,1:N)       + &
                        AC_tr(1:N,1:N) * c4  + &
                        AM_tr(1:N,1:N) * c3
      AQ_tr (1:N    ) = AQ_tr(1:N) + matmul( AM_tr(1:N,1:N),              c3 * (UPRE(1:N)-UT_tr(INDSYSB_tr(1:N))) + UT2_tr(INDSYSB_tr(1:N)) ) &
                                   - matmul( AC_tr(1:N,1:N), UTPRE(1:N) - c4 * (UPRE(1:N)-UT_tr(INDSYSB_tr(1:N))) - UT1_tr(INDSYSB_tr(1:N)) )
   elseif (itype == 1 ) then
              !-- Static solution
                   c3 = 0.d0
                   c5 = 0.d0
                   N  = NDFT0_tr
      UPRE  (1:N    ) = 0.d0;
      UTPRE (1:N    ) = 0.d0;
   endif


!-- call routines to solve the linear system using LU decomposition
   call DGETRF (     NDFT0_tr, NDFT0_tr, AK_tr, NDFT_tr, Ipivot,                 info(1))
   call DGETRS ('N', NDFT0_tr,   1     , AK_tr, NDFT_tr, Ipivot, AQ_tr, NDFT_tr, info(2))

   if (info(1)>0) write(*,*) 'LU Factorization    ',info(1)
   if (info(2)>0) write(*,*) 'LU Back Substitution',info(2)


!-- calculate the maximum and the mean errors, perform convergence check
   RESDT  = 0.d0
   RESDTM = 0.d0
   imax   = 1

   do i=1, NDFT0_tr
      RESDT  = RESDT + dabs(AQ_tr(i))
      if (RESDTM<dabs(AQ_tr(i))) then
         imax   = i
         RESDTM = dabs(AQ_tr(i))
      endif
   enddo

   RESDT = RESDT/dble(NDFT0_tr)

   if (RESDTM < ERR) ICONV = 1

   write (*,*) '      it = ', it, '  RESDT = ', RESDTM
!   write (* ,1) it, RESDT, RESDTM, imax, UT_tr(INDSYSB_tr(imax)) + AQ_tr(imax)
!   write (10,1) it, RESDT, RESDTM, imax, UT_tr(INDSYSB_tr(imax)) + AQ_tr(imax)

!1  format (2x,'it = ', i2,4x,'  RESDT =',2(e23.16,2x),i5,  4x,e23.16)


!-- Substitute the solution for all active dofs
   do i  = 1, NDFT0_tr
      is = INDSYSB_tr (i)
   
      UT_tr  (is) = UT_tr (is) + AQ_tr(i) * RLX
      UT2_tr (is) = c3 * ( UT_tr(is)-UPRE(i) )
      UT1_tr (is) = UTPRE(i) + c5 * UT2_tr (is)
   enddo


 END Subroutine Time_Integrate_tr
!-----------------------------------------------------------------------
 Subroutine set_init_truss_pos
!-----------------------------------------------------------------------
!
!-- Modify the connection nodes, based on the floater's initial positions
!
!-----------------------------------------------------------------------

 use truss
 use Cbeam !UflInit
 use Hydro

   implicit none

   real(8) :: AQ0(3,3), AQ01(3,3), AQ02(3,3), AQ03(3,3), A(3,3)
   integer :: nc, nbod,nod,n1,i, nf

!qq Valid for moorings attached to 1st floater
   nf = 1

   do nc   = 1, NCONNECT_tr
      nbod = connect_tr(nc)%IBODCONNECT_tr
      nod  = connect_tr(nc)%INODCONNECT_tr
      n1   = body_tr(nbod)%INODE_tr (nod)

      call DIAGO (3,A)

      do i = 1, 3
         call ROT_MATRIX ( i , floater(nf)% UflInit (i+3), AQ0, AQ01, AQ02, AQ03 )
         A = matmul(A,AQ0)
      enddo

      write(*,*)'xcon_pre',XG_tr(1:3,n1)
      XG_tr(1:3,n1) = floater(nf)% UflInit (1:3) + matmul( A(1:3,1:3), connect_tr(nc)%RLOC_CONNECT_tr(1:3) )
      write(*,*)'xcon_aft',XG_tr(1:3,n1)
   enddo !nc


 END Subroutine set_init_truss_pos
!----------------------------------------------------------------------
 Subroutine SHAPEFUNC6 ( HTAL, SHAPE,  NDFPE, NEQPE )
!----------------------------------------------------------------------

   implicit none

   integer :: NEQPE, NDFPE
   real(8) :: SHAPE (NEQPE, NDFPE), HTAL


   SHAPE      = 0.d0;
   SHAPE(1,1) = 1.d0 - HTAL
   SHAPE(2,2) = 1.d0 - HTAL
   SHAPE(3,3) = 1.d0 - HTAL
   SHAPE(1,4) = HTAL
   SHAPE(2,5) = HTAL
   SHAPE(3,6) = HTAL


 END Subroutine SHAPEFUNC6
!----------------------------------------------------------------------
 Subroutine MATRIX_REDUCT_tr (IP)
!----------------------------------------------------------------------

 use truss

   implicit none

   integer :: IP, nc, inod, ndf, il,jl, i, ndf0,dndf,ii,jj,ndf1


   if (IREDUCT_tr == 0) then

!-- Boundary Conditions
      call BOUNDC_tr (IP)              !- set the boundary conditions

      NDFT0_tr = NDFT_tr
      do i = 1, NDFT0_tr
         INDSYSB_tr(i) = i
      enddo

      return
   endif


!---------------------
!-- Keep Lines --
!---------------------
   ndf0 = 1
   ii   = 0

   do nc = 1,NCONSTR_tr 

!-- keep previous non zero dofs
      inod = constrain_tr(nc)%INODC_tr
      ndf  = NDFACC_tr (inod-1)
      dndf=ndf-ndf0+1

      if (dndf>0) then
         AM_tr     (ii+1:ii+dndf, 1:NDFT_tr) = AM_tr(ndf0:ndf, 1:NDFT_tr)
         AC_tr     (ii+1:ii+dndf, 1:NDFT_tr) = AC_tr(ndf0:ndf, 1:NDFT_tr)
         AK_tr     (ii+1:ii+dndf, 1:NDFT_tr) = AK_tr(ndf0:ndf, 1:NDFT_tr)
         AQ_tr     (ii+1:ii+dndf           ) = AQ_tr(ndf0:ndf           )
         do i=1,dndf
            INDSYSB_tr(ii+i) = ndf0+i-1
         enddo
         ii = ii + dndf
      endif

!-- Keep free b.c.
      do il=1,NDFPN_tr
         if (constrain_tr(nc)%INDXC_tr(il) == 1) cycle
         ii   = ii + 1
         ndf1 = ndf + il
         AM_tr     (ii, 1:NDFT_tr) = AM_tr(ndf1, 1:NDFT_tr)
         AC_tr     (ii, 1:NDFT_tr) = AC_tr(ndf1, 1:NDFT_tr)
         AK_tr     (ii, 1:NDFT_tr) = AK_tr(ndf1, 1:NDFT_tr)
         AQ_tr     (ii           ) = AQ_tr(ndf1           )
         INDSYSB_tr(ii           ) =       ndf1
      enddo
         
      ndf0 = ndf+NDFPN_tr+1
   enddo !nc

!-- keep last non zero dofs
      ndf  = NDFT_tr
      dndf=ndf-ndf0+1

      if (dndf>0) then
         AM_tr     (ii+1:ii+dndf, 1:NDFT_tr) = AM_tr(ndf0:ndf, 1:NDFT_tr)
         AC_tr     (ii+1:ii+dndf, 1:NDFT_tr) = AC_tr(ndf0:ndf, 1:NDFT_tr)
         AK_tr     (ii+1:ii+dndf, 1:NDFT_tr) = AK_tr(ndf0:ndf, 1:NDFT_tr)
         AQ_tr     (ii+1:ii+dndf           ) = AQ_tr(ndf0:ndf           )
         do i=1,dndf
            INDSYSB_tr(ii+i) = ndf0+i-1
         enddo
         ii = ii + dndf
      endif

!-- Set NDFT0_tr
   NDFT0_tr = ii


!---------------------
!-- Keep Columns --
!---------------------
   ndf0 = 1
   jj   = 0

   do nc = 1,NCONSTR_tr 

!-- keep previous non zero dofs
      inod = constrain_tr(nc)%INODC_tr
      ndf  = NDFACC_tr (inod-1)
      dndf=ndf-ndf0+1

      if (dndf>0) then
         AM_tr     (1:NDFT0_tr, jj+1:jj+dndf) = AM_tr(1:NDFT0_tr, ndf0:ndf)
         AC_tr     (1:NDFT0_tr, jj+1:jj+dndf) = AC_tr(1:NDFT0_tr, ndf0:ndf)
         AK_tr     (1:NDFT0_tr, jj+1:jj+dndf) = AK_tr(1:NDFT0_tr, ndf0:ndf)
         jj = jj + dndf
      endif

!-- Keep free b.c.
      do jl=1,NDFPN_tr
         if (constrain_tr(nc)%INDXC_tr(jl) == 1) cycle
         jj   = jj + 1
         ndf1 = ndf + jl
         AM_tr     (1:NDFT0_tr, jj) = AM_tr(1:NDFT0_tr, ndf1)
         AC_tr     (1:NDFT0_tr, jj) = AC_tr(1:NDFT0_tr, ndf1)
         AK_tr     (1:NDFT0_tr, jj) = AK_tr(1:NDFT0_tr, ndf1)
      enddo
         
      ndf0 = ndf+NDFPN_tr+1
   enddo !nc

!-- keep last non zero dofs
      ndf  = NDFT_tr
      dndf=ndf-ndf0+1

      if (dndf>0) then
         AM_tr     (1:NDFT0_tr, jj+1:jj+dndf) = AM_tr(1:NDFT0_tr, ndf0:ndf)
         AC_tr     (1:NDFT0_tr, jj+1:jj+dndf) = AC_tr(1:NDFT0_tr, ndf0:ndf)
         AK_tr     (1:NDFT0_tr, jj+1:jj+dndf) = AK_tr(1:NDFT0_tr, ndf0:ndf)
      endif


 END Subroutine MATRIX_REDUCT_tr
