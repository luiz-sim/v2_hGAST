!---------------------------------------------------------------------------------------------------------------------------------------------------!
!   National Technical University of Athens                                              Giannis Serafeim                                           !
!   School of Mechanical Engineering                                                     Dip. Mechanical Engineering - Phd Student                  !
!   Fluids Section                                                                       seraf@fluid.mech.ntua.gr                                   !
!   Laboratory of Aerodynamics                                                           Sunday, 09 00:53 Jun 2019 - Athens, Greece                 !
!---------------------------------------------------------------------------------------------------------------------------------------------------!
program loads_create

   implicit none
   integer(4)::i, j, k, l, m                               ! counters
   integer(4)::step_sta, step_end                          ! steps of start and end
   real(8)::ukn, coe                                       ! undefine value and loads coeffision
   real(8)::fx, fy, fz, mx, my, mz                         ! loads fx, fy, fz, mx, my and mz
   real(8)::fx_max, fy_max, fz_max, mx_max, my_max, mz_max ! maximum loads fx, fy, fz, mx, my and mz
   real(8)::max_loads(1:30,1:6)                            ! matix of maximum loads for each cross-section
   character(len=2)::num2                                  ! name whith 2 characters
   character(len=300)::path_dir                            ! path direct
   character(len=300)::dlc_dir                             ! dlc direct

   ! read path dir
   !--->
      read(*,"(a)") path_dir
      read(*,"(a)") dlc_dir
      read(*,*) coe, step_sta, step_end
   !--->

   ! bin2ascii loads
   !--->
      call system("cd "//trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//" ; /home/seraf/hGAST/bin2ascii.out<bin2ascii.inp")
   !--->

   do k=1,30
      call numbering2(k,num2)

      ! find maximum loads
      !--->
         open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/loads001_"//num2//".dat") ; open(unit=2, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/loads002_"//num2//".dat") ; open(unit=3, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/loads003_"//num2//".dat")
            do l=1,(step_sta-1)
               read(1,*) ; read(2,*) ; read(3,*)
            end do
            fx_max=0.0 ; fy_max=0.0 ; fz_max=0.0 ; mx_max=0.0 ; my_max=0.0 ; mz_max=0.0
            do l=step_sta,step_end
               read(1,*) ukn, ukn, fx, fy, fz, mx, my, mz
               if (fx>fx_max) fx_max=fx ; if (fy>fy_max) fy_max=fy ; if (fz>fz_max) fz_max=fz ; if (mx>mx_max) mx_max=mx ; if (my>my_max) my_max=my ; if (mz>mz_max) mz_max=mz
               read(2,*) ukn, ukn, fx, fy, fz, mx, my, mz
               if (fx>fx_max) fx_max=fx ; if (fy>fy_max) fy_max=fy ; if (fz>fz_max) fz_max=fz ; if (mx>mx_max) mx_max=mx ; if (my>my_max) my_max=my ; if (mz>mz_max) mz_max=mz
               read(3,*) ukn, ukn, fx, fy, fz, mx, my, mz
               if (fx>fx_max) fx_max=fx ; if (fy>fy_max) fy_max=fy ; if (fz>fz_max) fz_max=fz ; if (mx>mx_max) mx_max=mx ; if (my>my_max) my_max=my ; if (mz>mz_max) mz_max=mz
            end do
         close (3) ; close (2) ; close (1)
         max_loads(k,1)=1000.0*coe*fx ; max_loads(k,2)=1000.0*coe*fy ; max_loads(k,3)=1000.0*coe*fz ; max_loads(k,4)=1000.0*coe*mx ; max_loads(k,5)=1000.0*coe*my ; max_loads(k,6)=1000.0*coe*mz
      !--->

      ! find number of loads which check in stresses analysis
      !--->
         open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/loads001_"//num2//".dat") ; open(unit=2, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/loads002_"//num2//".dat") ; open(unit=3, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/loads003_"//num2//".dat")
            do l=1,(step_sta-1)
               read(1,*) ; read(2,*) ; read(3,*)
            end do
            m=0
            do l=step_sta,step_end
               read(1,*) ukn, ukn, fx, fy, fz, mx, my, mz
               if ((fx>0.975*fx_max).or.(fy>0.975*fy_max).or.(fz>0.975*fz_max).or.(mx>0.900*mx_max).or.(my>0.975*my_max).or.(mz>0.975*mz_max)) m=m+1
               read(2,*) ukn, ukn, fx, fy, fz, mx, my, mz
               if ((fx>0.975*fx_max).or.(fy>0.975*fy_max).or.(fz>0.975*fz_max).or.(mx>0.900*mx_max).or.(my>0.975*my_max).or.(mz>0.975*mz_max)) m=m+1
               read(3,*) ukn, ukn, fx, fy, fz, mx, my, mz
               if ((fx>0.975*fx_max).or.(fy>0.975*fy_max).or.(fz>0.975*fz_max).or.(mx>0.900*mx_max).or.(my>0.975*my_max).or.(mz>0.975*mz_max)) m=m+1
            end do
         close (3) ; close (2) ; close (1)
      !--->

      ! write the loads which will be check in stresses analysis
      !--->
         open(unit=4, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/max_loa_"//num2//".dat")
         open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/loads001_"//num2//".dat") ; open(unit=2, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/loads002_"//num2//".dat") ; open(unit=3, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/loads003_"//num2//".dat")
            write(4,*) m
            do l=1,(step_sta-1)
               read(1,*) ; read(2,*) ; read(3,*)
            end do
            do l=step_sta,step_end
               read(1,*) ukn, ukn, fx, fy, fz, mx, my, mz
               if ((fx>0.975*fx_max).or.(fy>0.975*fy_max).or.(fz>0.975*fz_max).or.(mx>0.900*mx_max).or.(my>0.975*my_max).or.(mz>0.975*mz_max)) write(4,"(6(f25.5,3x))") 1000.0*coe*fx, 1000.0*coe*fy, 1000.0*coe*fz, 1000.0*coe*mx, 1000.0*coe*my, 1000.0*coe*mz
               read(2,*) ukn, ukn, fx, fy, fz, mx, my, mz
               if ((fx>0.975*fx_max).or.(fy>0.975*fy_max).or.(fz>0.975*fz_max).or.(mx>0.900*mx_max).or.(my>0.975*my_max).or.(mz>0.975*mz_max)) write(4,"(6(f25.5,3x))") 1000.0*coe*fx, 1000.0*coe*fy, 1000.0*coe*fz, 1000.0*coe*mx, 1000.0*coe*my, 1000.0*coe*mz
               read(3,*) ukn, ukn, fx, fy, fz, mx, my, mz
               if ((fx>0.975*fx_max).or.(fy>0.975*fy_max).or.(fz>0.975*fz_max).or.(mx>0.900*mx_max).or.(my>0.975*my_max).or.(mz>0.975*mz_max)) write(4,"(6(f25.5,3x))") 1000.0*coe*fx, 1000.0*coe*fy, 1000.0*coe*fz, 1000.0*coe*mx, 1000.0*coe*my, 1000.0*coe*mz
            end do
         close (3) ; close (2) ; close (1)
         close (4)
      !---->

   end do

   call system("/home/seraf/hGAST/post_hGAST.dir/pre_smalt.exe<"//trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/path.inp")

   ! write results of maximum loads
   !--->
      open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/max_loads.dat")
         write(1,"(a)") "num   Fx[N]   Fy[N]   Fz[N]   Mx[N/m]   My[N/m]   Mz[N/m]"
         do k=1,30
            write(1,"(i3,3x,6(f25.5,3x))") k, (max_loads(k,m), m=1,6)
         end do
      close (1)
   !--->

end program loads_create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine numbering2(i,num2)

   implicit none
   integer(4)::i          ! counter
   character(len=2)::num2 ! name with 2 characters

   ! select cases
   !--->
      select case (i)
         case (1)  ; num2="01" ; case (2)  ; num2="02" ; case (3)  ; num2="03" ; case (4)  ; num2="04" ; case (5)  ; num2="05" ; case (6)  ; num2="06" ; case (7)  ; num2="07" ; case (8)  ; num2="08" ; case (9)  ; num2="09" ; case (10)    ; num2="10"
         case (11) ; num2="11" ; case (12) ; num2="12" ; case (13) ; num2="13" ; case (14) ; num2="14" ; case (15) ; num2="15" ; case (16) ; num2="16" ; case (17) ; num2="17" ; case (18) ; num2="18" ; case (19) ; num2="19" ; case (20)    ; num2="20"
         case (21) ; num2="21" ; case (22) ; num2="22" ; case (23) ; num2="23" ; case (24) ; num2="24" ; case (25) ; num2="25" ; case (26) ; num2="26" ; case (27) ; num2="27" ; case (28) ; num2="28" ; case (29) ; num2="29" ; case default ; num2="30"
      end select
   !--->

   return
end subroutine numbering2
!---------------------------------------------------------------------------------------------------------------------------------------------------!
!                                                                       -END-                                                                       !
!---------------------------------------------------------------------------------------------------------------------------------------------------!
