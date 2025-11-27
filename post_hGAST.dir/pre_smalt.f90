!---------------------------------------------------------------------------------------------------------------------------------------------------!
!   National Technical University of Athens                                              Giannis Serafeim                                           !
!   School of Mechanical Engineering                                                     Dip. Mechanical Engineering - Phd Student                  !
!   Fluids Section                                                                       seraf@fluid.mech.ntua.gr                                   !
!   Laboratory of Aerodynamics                                                           Sunday, 09 01:47 Jun 2019 - Athens, Greece                 !
!---------------------------------------------------------------------------------------------------------------------------------------------------!
program pre_smalt

   implicit none
   integer(4)::i, j, k, l, m                                                            ! counters
   integer(4)::num_layers, num_steps                                                    ! number of layers and number of steps which check
   integer(4)::kp(1:30,1:16), kp_nose, kp_6(1:30), kp_7(1:30), kp_10(1:30), kp_11(1:30) ! key-points matrix, key-point of nose and impermanent key-points of kp-6, kp-7, kp-10 and kp-11
   real(8)::ukn, th, d, d_min                                                           ! undefined variable and two help values
   real(8)::radius(1:30), x_swep(1:30), z_pre(1:30)                                     ! radius, x-sweep and z-prebent
   real(8)::chord(1:30), twist_or(1:30), twist_re(1:30), ply(1:30)                      ! chords, twists of default blade, twist of retwist and ply angle
   real(8)::caps_move(1:30), coeff(1:30)                                                ! cap's move and coeffision of laminate thikness
   real(8)::triax(1:30,1:9), biax(1:30,1:9), uniax(1:30,1:9), balsa(1:30,1:9)           ! laminate thikness for triax, biax, uniax and balsa
   real(8)::pro(1:19,1:4), air(1:109,1:2,1:30)                                          ! material's properties and airfoils geometry
   real(8)::x, z                                                                        ! x and z axis
   real(8)::gra1, gra2, dth(1:30)                                                       ! gradient-1, gradient-2 and angle between initial and move web-a
   real(8)::s, s_max, s_max1, s_max2                                                    ! legth and maximum legth
   real(8)::fx, fy, fz, mx, my, mz                                                      ! loads in hGAST system: fx, fy, fz, mx, my, mz
   real(8)::f_x_max(1:30), f_y_max(1:30), f_z_max(1:30)                                 ! maximum loads fx, fy and fz
   real(8)::m_x_max(1:30), m_y_max(1:30), m_z_max(1:30)                                 ! maximum loads mx, my and mz
   real(8)::t_max, n_max, p_max, mat(1:108,1:7)                                         ! help values for tsai-wu, normal, perimeter and general
   real(8)::twu_max(1:30), nor_max(1:30), per_max(1:30)                                 ! maximum values for tsai-wu, normal and perimeter
   real(8)::loa_max(1:30,1:6)                                                           ! loads which give the maximum tsai-wu for each cross-section
   character(len=2)::num                                                                ! name of number
   character(len=10)::name_part                                                         ! name of cross-section part
   character(len=300)::path_dir                                                         ! path direct
   character(len=300)::dlc_dir                                                          ! dlc direct

   read(*,"(a)") path_dir
   read(*,"(a)") dlc_dir

   twu_max=0.0 ; nor_max=0.0 ; per_max=0.0
   call system("rm -r "//trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/other.dir ; mkdir "//trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/other.dir")

   ! read blade_construction data
   !--->
      open(unit=1, file=trim(adjustl(path_dir))//"/data.dir/blade_construction.dat")
         read(1,*)
         do i=1,30
            read(1,*) ukn, radius(i), x_swep(i), z_pre(i), chord(i), twist_or(i), twist_re(i), ply(i), caps_move(i), coeff(i)
         end do
      close (1)
   !--->

   ! read properties data
   !--->
      open(unit=1, file=trim(adjustl(path_dir))//"/data.dir/properties.dat")
         read(1,*)
         do i=1,19
            read(1,*) (pro(i,j), j=1,4)
         end do
      close (1)
   !--->

   ! read laminate tickness for each part
   !--->
      do k=1,9
         select case (k)
            case (1) ; name_part="tail_v"  ; case (2) ; name_part="tail_a"   ; case (3)     ; name_part="tail_b"
            case (4) ; name_part="tail_c"  ; case (5) ; name_part="trailing" ; case (6)     ; name_part="caps"
            case (7) ; name_part="leading" ; case (8) ; name_part="nose"     ; case default ; name_part="webs"
         end select
         open(unit=1, file=trim(adjustl(path_dir))//"/data.dir/laminate.dir/"//trim(adjustl(name_part))//".dat")
            read(1,*)
            do i=1,30
               read(1,*) ukn, triax(i,k), biax(i,k), uniax(i,k), balsa(i,k)
            end do
         close (1)
      end do
   !--->

   ! read airfoil's geometry
   !--->
      do i=1,30
         call numbering2(i,num)
         open(unit=1, file=trim(adjustl(path_dir))//"/data.dir/airfoils.dir/airfoil_"//num//".dat")
            read(1,*)
            do j=1,109
               read(1,*) ukn, (air(j,k,i), k=1,2)
            end do
         close (1)
      end do
   !--->

   ! create key-points for each cross-section
   !--->

      ! read keypoints
      !--->
         open(unit=1, file=trim(adjustl(path_dir))//"/data.dir/keypoints.dat")
            read(1,*)
            do i=1,30
               read(1,*) kp(i,1), kp(i,2), kp(i,3), kp(i,4), kp(i,5), kp_6(i), kp_7(i), kp(i,8), kp(i,9), kp_10(i), kp_11(i), kp(i,12), kp(i,13), kp(i,14), kp(i,15), kp(i,16)
            end do
         close (1)
      !--->

      ! move caps and calculate the angle rotation
      !--->
         do i=1,30
            kp(i,6)=kp_6(i)+caps_move(i)   ; if (kp(i,6)<=(kp(i,5)+1)) kp(i,6)=kp(i,5)+1
            kp(i,7)=kp_7(i)+caps_move(i)   ; if (kp(i,7)>=(kp(i,8)-1)) kp(i,7)=kp(i,8)-1
            kp(i,10)=kp_10(i)+caps_move(i) ; if (kp(i,10)<=(kp(i,9)+1)) kp(i,10)=kp(i,9)+1
            kp(i,11)=kp_11(i)+caps_move(i) ; if (kp(i,11)>=(kp(i,12)-1)) kp(i,11)=kp(i,12)-1
            gra1=(air(kp_7(i),1,i)-air(kp_10(i),1,i))/(air(kp_7(i),2,i)-air(kp_10(i),2,i))
            gra2=(air(kp(i,7),1,i)-air(kp(i,10),1,i))/(air(kp(i,7),2,i)-air(kp(i,10),2,i))
            dth(i)=(180.0/(4.0*atan(1.0)))*atan((gra2-gra1)/(1.0+gra2*gra1))
         end do
      !--->
   !--->

   f_x_max=0.0 ; f_y_max=0.0 ; f_z_max=0.0
   m_x_max=0.0 ; m_y_max=0.0 ; m_z_max=0.0
   do i=1,30
      call numbering2(i,num)

      ! write "geometry.dat"
      !--->
         open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/other.dir/geometry.dat")
            write(1,"(a)") "**number of data points y(n).z(n)"
            write(1,"(a)") "109"
            write(1,"(a)") "**points y(i).z(i)"
            write(1,"(a)") "1"
            do j=1,109
               th=(twist_or(i)+twist_re(i)-90.0)*4.0*atan(1.0)/180.0
               x=chord(i)*(air(j,1,i)*cos(th)-air(j,2,i)*sin(th))
               z=chord(i)*(air(j,1,i)*sin(th)+air(j,2,i)*cos(th))
               write(1,"(i3,3x,f15.7,3x,f15.7)") j, x, z
            end do
         close (1)
      !--->

      ! write structury.dat"
      !--->
         open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/other.dir/structury.dat")
            write(1,"(a)") "** input for multiple section [glass-polyester] airfoil with a sigle web"
            write(1,"(a)") "4   10   0"
            write(1,"(a)") "** define materials & prop's..."
            write(1,"(i3,3x,5(f25.5,3x))") 1, pro(10,1), pro(1,1), pro(2,1), pro(4,1), pro(7,1)
            write(1,"(5(f25.5,3x))")          pro(9,1), pro(3,1), pro(6,1), pro(8,1), pro(5,1)
            write(1,"(5(f25.5,3x))")          pro(11,1), pro(12,1), pro(13,1), pro(14,1), pro(15,1)
            write(1,"(4(f25.5,3x))")          pro(16,1), pro(17,1), pro(18,1), pro(19,1)
            write(1,*)
            write(1,"(i3,3x,5(f25.5,3x))") 2, pro(10,2), pro(1,2), pro(2,2), pro(4,2), pro(7,2)
            write(1,"(5(f25.5,3x))")          pro(9,2), pro(3,2), pro(6,2), pro(8,2), pro(5,2)
            write(1,"(5(f25.5,3x))")          pro(11,2), pro(12,2), pro(13,2), pro(14,2), pro(15,2)
            write(1,"(4(f25.5,3x))")          pro(16,2), pro(17,2), pro(18,2), pro(19,2)
            write(1,*)
            write(1,"(i3,3x,5(f25.5,3x))") 3, pro(10,3), pro(1,3), pro(2,3), pro(4,3), pro(7,3)
            write(1,"(5(f25.5,3x))")          pro(9,3), pro(3,3), pro(6,3), pro(8,3), pro(5,3)
            write(1,"(5(f25.5,3x))")          pro(11,3), pro(12,3), pro(13,3), pro(14,3), pro(15,3)
            write(1,"(4(f25.5,3x))")          pro(16,3), pro(17,3), pro(18,3), pro(19,3)
            write(1,*)
            write(1,"(i3,3x,5(f25.5,3x))") 4, pro(10,4), pro(1,4), pro(2,4), pro(4,4), pro(7,4)
            write(1,"(5(f25.5,3x))")          pro(9,4), pro(3,4), pro(6,4), pro(8,4), pro(5,4)
            write(1,"(5(f25.5,3x))")          pro(11,4), pro(12,4), pro(13,4), pro(14,4), pro(15,4)
            write(1,"(4(f25.5,3x))")          pro(16,4), pro(17,4), pro(18,4), pro(19,4)
            write(1,*)
            write(1,"(a)") "** laminate configurations..."

            ! tail-v
            !--->
               num_layers=0 ; if (triax(i,1)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,1)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,1)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,1)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 1, num_layers
               num_layers=1
               if (triax(i,1)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,1), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,1)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,1),  "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,1)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,1), "0.0" ; num_layers=num_layers+1 ; end if
               if (balsa(i,1)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 4 ", coeff(i)*balsa(i,1), "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,1)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,1), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,1)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,1),  "0.0" ; num_layers=num_layers+1 ; end if
               if (triax(i,1)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,1), "0.0" ; num_layers=num_layers+1 ; end if
            !--->

            ! tail-a
            !--->
               num_layers=0 ; if (triax(i,2)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,2)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,2)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,2)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 2, num_layers
               num_layers=1
               if (triax(i,2)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,1), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,2)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,2),  "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,2)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,3), "0.0" ; num_layers=num_layers+1 ; end if
               if (balsa(i,2)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 4 ", coeff(i)*balsa(i,4), "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,2)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,3), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,2)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,2),  "0.0" ; num_layers=num_layers+1 ; end if
               if (triax(i,2)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,1), "0.0" ; num_layers=num_layers+1 ; end if
            !--->

            ! tail-b
            !--->
               num_layers=0 ; if (triax(i,3)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,3)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,3)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,3)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 3, num_layers
               num_layers=1
               if (triax(i,3)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,3), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,3)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,3),  "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,3)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,3), "0.0" ; num_layers=num_layers+1 ; end if
               if (balsa(i,3)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 4 ", coeff(i)*balsa(i,3), "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,3)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,3), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,3)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,3),  "0.0" ; num_layers=num_layers+1 ; end if
               if (triax(i,3)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,3), "0.0" ; num_layers=num_layers+1 ; end if
            !--->

            ! tail-c
            !--->
               num_layers=0 ; if (triax(i,4)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,4)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,4)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,4)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 4, num_layers
               num_layers=1
               if (triax(i,4)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,4), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,4)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,4),  "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,4)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,4), "0.0" ; num_layers=num_layers+1 ; end if
               if (balsa(i,4)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 4 ", coeff(i)*balsa(i,4), "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,4)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,4), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,4)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,4),  "0.0" ; num_layers=num_layers+1 ; end if
               if (triax(i,4)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,4), "0.0" ; num_layers=num_layers+1 ; end if
            !--->

            ! trailing
            !--->
               num_layers=0 ; if (triax(i,5)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,5)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,5)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,5)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 5, num_layers
               num_layers=1
               if (triax(i,5)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,5), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,5)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,5),  "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,5)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,5), "0.0" ; num_layers=num_layers+1 ; end if
               if (balsa(i,5)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 4 ", coeff(i)*balsa(i,5), "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,5)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,5), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,5)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,5),  "0.0" ; num_layers=num_layers+1 ; end if
               if (triax(i,5)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,5), "0.0" ; num_layers=num_layers+1 ; end if
            !--->

            ! caps
            !--->
               num_layers=0 ; if (triax(i,6)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,6)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,6)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,6)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 6, num_layers
               num_layers=1
               if (triax(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 1 ", coeff(i)*triax(i,6), "0.0"   ; num_layers=num_layers+1 ; end if
               if (biax(i,6)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 2 ", coeff(i)*biax(i,6),  "0.0"   ; num_layers=num_layers+1 ; end if
               if (uniax(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,f9.5)") num_layers, " 3 ", coeff(i)*uniax(i,6),  ply(i) ; num_layers=num_layers+1 ; end if
               if (balsa(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 4 ", coeff(i)*balsa(i,6), "0.0"   ; num_layers=num_layers+1 ; end if
               if (uniax(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,f9.5)") num_layers, " 3 ", coeff(i)*uniax(i,6),  ply(i) ; num_layers=num_layers+1 ; end if
               if (biax(i,6)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 2 ", coeff(i)*biax(i,6),  "0.0"   ; num_layers=num_layers+1 ; end if
               if (triax(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 1 ", coeff(i)*triax(i,6), "0.0"   ; num_layers=num_layers+1 ; end if
            !--->

            ! leading
            !--->
               num_layers=0 ; if (triax(i,7)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,7)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,7)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,7)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 7, num_layers
               num_layers=1
               if (triax(i,7)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,7), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,7)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,7),  "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,7)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,7), "0.0" ; num_layers=num_layers+1 ; end if
               if (balsa(i,7)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 4 ", coeff(i)*balsa(i,7), "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,7)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,7), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,7)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,7),  "0.0" ; num_layers=num_layers+1 ; end if
               if (triax(i,7)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,7), "0.0" ; num_layers=num_layers+1 ; end if
            !--->

            ! nose
            !--->
               num_layers=0 ; if (triax(i,8)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,8)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,8)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,8)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 8, num_layers
               num_layers=1
               if (triax(i,8)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,8), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,8)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,8),  "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,8)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,8), "0.0" ; num_layers=num_layers+1 ; end if
               if (balsa(i,8)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 4 ", coeff(i)*balsa(i,8), "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,8)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,8), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,8)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,8),  "0.0" ; num_layers=num_layers+1 ; end if
               if (triax(i,8)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,8), "0.0" ; num_layers=num_layers+1 ; end if
            !--->

            ! caps
            !--->
               num_layers=0 ; if (triax(i,6)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,6)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,6)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,6)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 9, num_layers
               num_layers=1
               if (triax(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 1 ", coeff(i)*triax(i,6),  "0.0"  ; num_layers=num_layers+1 ; end if
               if (biax(i,6)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 2 ", coeff(i)*biax(i,6),   "0.0"  ; num_layers=num_layers+1 ; end if
               if (uniax(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,f9.5)") num_layers, " 3 ", coeff(i)*uniax(i,6), -ply(i) ; num_layers=num_layers+1 ; end if
               if (balsa(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 4 ", coeff(i)*balsa(i,6),   "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,f9.5)") num_layers, " 3 ", coeff(i)*uniax(i,6), -ply(i) ; num_layers=num_layers+1 ; end if
               if (biax(i,6)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 2 ", coeff(i)*biax(i,6),   "0.0"  ; num_layers=num_layers+1 ; end if
               if (triax(i,6)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)")    num_layers, " 1 ", coeff(i)*triax(i,6),  "0.0"  ; num_layers=num_layers+1 ; end if
            !--->

            ! webs
            !--->
               num_layers=0 ; if (triax(i,9)>0.000000000000001) num_layers=num_layers+2 ; if (biax(i,9)>0.000000000000001) num_layers=num_layers+2 ; if (uniax(i,9)>0.000000000000001) num_layers=num_layers+2 ; if (balsa(i,9)>0.000000000000001) num_layers=num_layers+1 ; write(1,"(i3,3x,i3)") 10, num_layers
               num_layers=1
               if (triax(i,9)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,9), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,9)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,9),  "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,9)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,9), "0.0" ; num_layers=num_layers+1 ; end if
               if (balsa(i,9)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 4 ", coeff(i)*balsa(i,9), "0.0" ; num_layers=num_layers+1 ; end if
               if (uniax(i,9)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 3 ", coeff(i)*uniax(i,9), "0.0" ; num_layers=num_layers+1 ; end if
               if (biax(i,9)>0.000000000000001) then  ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 2 ", coeff(i)*biax(i,9),  "0.0" ; num_layers=num_layers+1 ; end if
               if (triax(i,9)>0.000000000000001) then ; write(1,"(i3,3x,a,3x,f25.5,3x,a)") num_layers, " 1 ", coeff(i)*triax(i,9), "0.0" ; num_layers=num_layers+1 ; end if
            !--->

            write(1,*)
            write(1,*) "1"
            write(1,"(a)")               "Section 1 17"
            write(1,"(a)")               "Skin 15"
            write(1,"(a,3x,i3,3x,i3,a)") "  1 ", kp(i,1),  kp(i,2),  " 1 " ! tail-v 
            write(1,"(a,3x,i3,3x,i3,a)") "  2 ", kp(i,2),  kp(i,3),  " 2 " ! tail-a
            write(1,"(a,3x,i3,3x,i3,a)") "  3 ", kp(i,3),  kp(i,4),  " 3 " ! tail-b
            write(1,"(a,3x,i3,3x,i3,a)") "  4 ", kp(i,4),  kp(i,5),  " 4 " ! tail-c
            write(1,"(a,3x,i3,3x,i3,a)") "  5 ", kp(i,5),  kp(i,6),  " 5 " ! trailing
            write(1,"(a,3x,i3,3x,i3,a)") "  6 ", kp(i,6),  kp(i,7),  " 6 " ! caps
            write(1,"(a,3x,i3,3x,i3,a)") "  7 ", kp(i,7),  kp(i,8),  " 7 " ! leading
            write(1,"(a,3x,i3,3x,i3,a)") "  8 ", kp(i,8),  kp(i,9),  " 8 " ! nose
            write(1,"(a,3x,i3,3x,i3,a)") "  9 ", kp(i,9),  kp(i,10), " 7 " ! leaading
            write(1,"(a,3x,i3,3x,i3,a)") " 10 ", kp(i,10), kp(i,11), " 9 " ! caps
            write(1,"(a,3x,i3,3x,i3,a)") " 11 ", kp(i,11), kp(i,12), " 5 " ! trailing
            write(1,"(a,3x,i3,3x,i3,a)") " 12 ", kp(i,12), kp(i,13), " 4 " ! tail-c
            write(1,"(a,3x,i3,3x,i3,a)") " 13 ", kp(i,13), kp(i,14), " 3 " ! tail-b
            write(1,"(a,3x,i3,3x,i3,a)") " 14 ", kp(i,14), kp(i,15), " 2 " ! tail-a
            write(1,"(a,3x,i3,3x,i3,a)") " 15 ", kp(i,15), kp(i,16), " 1 " ! tail-v

            write(1,"(a)")               "Webs 2"
            write(1,"(a,3x,i3,3x,i3,a)") " 16 ", kp(i,7),  kp(i,10), " 10 " ! webs
            write(1,"(a,3x,i3,3x,i3,a)") " 17 ", kp(i,11), kp(i,6),  " 10 " ! webs
            write(1,"(a)")               "Loops 3"
            write(1,"(a)")               " 1   4 "
            write(1,"(a)")               " 7   8   9   -16"
            write(1,"(a)")               " 2   4 "
            write(1,"(a)")               " 6   16   10   17 "
            write(1,"(a)")               " 3   11 "
            write(1,"(a)")               " 1   2   3   4   5   -17   11   12   13   14   15"
         close (1)
      !--->

      open(unit=2, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/max_loa_"//num//".dat")
         read(2,*) num_steps
         t_max=0.0   ; n_max=0.0   ; p_max=0.0
         do k=1,num_steps

            ! read fx, fy, fz, mx, my, mz
            !--->
               read(2,*) fx, fy, fz, mx, my, mz
               if (abs(fx)>abs(f_x_max(i))) f_x_max(i)=abs(fx) ; if (abs(mx)>abs(m_x_max(i))) m_x_max(i)=abs(mx)
               if (abs(fy)>abs(f_y_max(i))) f_y_max(i)=abs(fy) ; if (abs(my)>abs(m_y_max(i))) m_y_max(i)=abs(my)
               if (abs(fz)>abs(f_z_max(i))) f_z_max(i)=abs(fz) ; if (abs(mz)>abs(m_z_max(i))) m_z_max(i)=abs(mz)
            !--->

            ! write "loads.dat"
            !--->
               open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/other.dir/loads.dat")
                  write(1,*) fy
                  write(1,*) fz
                  write(1,*) fx
                  write(1,*) mx
                  write(1,*) mz
                  write(1,*) my
               close (1)
            !--->

            ! run smalt
            !--->
               call system("/home/seraf/hGAST/post_hGAST.dir/smalt.exe <"//trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/path.inp")
            !--->

            ! check stresses results
            !--->
               open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/other.dir/stresses_tswu.dat")
                  do l=1,108
                     read(1,*) (mat(l,m), m=1,7)
                  end do
               close (1)
               t_max=maxval(mat(:,:))

               open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/other.dir/stresses_norm.dat")
                  do l=1,108
                     read(1,*) (mat(l,m), m=1,7)
                  end do
               close (1)
               n_max=maxval(mat(:,:))

               open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/other.dir/stresses_peri.dat")
                  do l=1,108
                     read(1,*) (mat(l,m), m=1,7)
                  end do
               close (1)
               p_max=maxval(mat(:,:))
            !--->

            if (t_max>twu_max(i)) then
               loa_max(i,1)=fx ; loa_max(i,2)=fy ; loa_max(i,3)=fz ; loa_max(i,4)=mx ; loa_max(i,5)=my ; loa_max(i,6)=mz ; twu_max(i)=t_max
            end if
            if (n_max>nor_max(i)) nor_max(i)=n_max
            if (p_max>per_max(i)) per_max(i)=p_max
         end do
      close (2)
   end do

   ! write results
   !--->
      open(unit=1, file=trim(adjustl(path_dir))//"/"//trim(adjustl(dlc_dir))//"/max_twu_nor_per.dat")
         write(1,"(a)") "num   tsai-wu[-]   normal[N/m^2]   perimeter[N/m^2]   !   Fx[N]   Fy[N]   Fz[N]   Mx[N/m]   My[N/m]   Mz[N/m] ---> loads which give maximum tsai-wu"
         do i=1,30
            write(1,"(i3,3x,3(f25.5,3x),a,6(f25.5,3x))") i, twu_max(i), nor_max(i), per_max(i), " ! ", loa_max(i,1)/f_x_max(i), loa_max(i,2)/f_y_max(i), loa_max(i,3)/f_z_max(i), loa_max(i,4)/m_x_max(i), loa_max(i,5)/m_y_max(i), loa_max(i,6)/m_z_max(i)
         end do
      close (1)
   !--->

end program pre_smalt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine numbering2(i,num)

   implicit none
   integer(4)::i         ! counter
   character(len=2)::num ! name with 2 characters

   !select case
   !--->
      select case(i)
         case (1) ; num="01" ; case (6)  ; num="06" ; case (11) ; num="11" ; case (16) ; num="16" ; case (21) ; num="21" ; case (26)    ; num="26"
         case (2) ; num="02" ; case (7)  ; num="07" ; case (12) ; num="12" ; case (17) ; num="17" ; case (22) ; num="22" ; case (27)    ; num="27"
         case (3) ; num="03" ; case (8)  ; num="08" ; case (13) ; num="13" ; case (18) ; num="18" ; case (23) ; num="23" ; case (28)    ; num="28"
         case (4) ; num="04" ; case (9)  ; num="09" ; case (14) ; num="14" ; case (19) ; num="19" ; case (24) ; num="24" ; case (29)    ; num="29"
         case (5) ; num="05" ; case (10) ; num="10" ; case (15) ; num="15" ; case (20) ; num="20" ; case (25) ; num="25" ; case default ; num="30"
      end select
   !--->

   return
end subroutine numbering2
!---------------------------------------------------------------------------------------------------------------------------------------------------!
!                                                                       -END-                                                                       !
!---------------------------------------------------------------------------------------------------------------------------------------------------!
