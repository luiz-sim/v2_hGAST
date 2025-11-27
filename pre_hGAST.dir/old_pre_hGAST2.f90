!---------------------------------------------------------------------------------------------------------------------------------------------------!
!   National Technical University of Athens                                              Giannis Serafeim                                           !
!   School of Mechanical Engineering                                                     Dip. Mechanical Engineering - Phd Student                  !
!   Fluids Section                                                                       seraf@fluid.mech.ntua.gr                                   !
!   Laboratory of Aerodynamics                                                           Tuesday, 04 22:41 Jun 2019 - Athens, Greece                !
!---------------------------------------------------------------------------------------------------------------------------------------------------!
program pre_hGAST

   implicit none
   integer(4)::i, j, k                                                                  ! counters
   integer(4)::num_layers                                                               ! number of layers
   integer(4)::kp(1:30,1:16), kp_nose, kp_6(1:30), kp_7(1:30), kp_10(1:30), kp_11(1:30) ! key-points matrix, key-point of nose and impermanent key-points of kp-6, kp-7, kp-10 and kp-11
   integer(4)::caps_move(1:30)
   real(8)::ukn, th, d, d_min                                                           ! undefined variable and three help values
   real(8)::radius(1:30), x_swep(1:30), z_pre(1:30)                                     ! radius, x-sweep and z-prebent
   real(8)::chord(1:30), twist_or(1:30), twist_re(1:30), ply(1:30)                      ! chords, twists of default blade, twist of retwist and ply angle
   real(8)::caps_mov(1:30)                                                              ! caps move
   real(8)::coeff(1:30)                                                                 ! coeffision of laminate thikness
   real(8)::triax(1:30,1:9), biax(1:30,1:9), uniax(1:30,1:9), balsa(1:30,1:9)           ! laminate thikness for triax, biax, uniax and balsa
   real(8)::pro(1:19,1:4), air(1:109,1:2,1:30)                                          ! material's properties and airfoils geometry
   real(8)::x, z                                                                        ! x and z axis
   real(8)::gra1, gra2, dth(1:30)                                                       ! gradient-1, gradient-2 and angle between initial and move web-a
   character(len=2)::num                                                                ! name of number
   character(len=10)::name_part                                                         ! name of cross-section part
   character(len=300)::path_dir                                                         ! path direct

   read(*,"(a)") path_dir
   call system("rm -r "//trim(adjustl(path_dir))//"/other.dir ; mkdir "//trim(adjustl(path_dir))//"/other.dir")

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

         open(unit=1, file=trim(adjustl(path_dir))//"/data.dir/caps_rotation.txt")
            write(1,"(a)") "num   radius[m]   rotataion[deg]"
            do i=1,30
               write(1,"(i3,3x,2(f15.7,3x))") i, radius(i), dth(i)
            end do
         close (1)
      !--->
   !--->

   do i=1,30
      ! write "geometry.dat"
      !--->
         open(unit=1, file=trim(adjustl(path_dir))//"/other.dir/geometry.dat")
            write(1,"(a)") "**number of data points y(n).z(n)"
            write(1,"(a)") "109"
            write(1,"(a)") "**points y(i).z(i)"
            write(1,"(a)") "1"
            th=(twist_or(i)+twist_re(i)-90.0)*4.0*atan(1.0)/180.0
            do j=109,1,-1
               x=chord(i)*(air(j,1,i)*cos(th)-air(j,2,i)*sin(th))
               z=chord(i)*(air(j,1,i)*sin(th)+air(j,2,i)*cos(th))
               write(1,"(i3,3x,f15.7,3x,f15.7)") (110-j), x, z
            end do
         close (1)
      !--->
!if (i==7) call system("cd "//trim(adjustl(path_dir))//"/other.dir/ ; mv geometry.dat geom7.txt")
      ! write structury.dat"
      !--->
         open(unit=1, file=trim(adjustl(path_dir))//"/other.dir/structury.dat")
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
            write(1,"(a,3x,i3,3x,i3,a)") "  1 ", (110-kp(i,16)), (110-kp(i,15)), " 1 " ! tail-v 
            write(1,"(a,3x,i3,3x,i3,a)") "  2 ", (110-kp(i,15)), (110-kp(i,14)), " 2 " ! tail-a
            write(1,"(a,3x,i3,3x,i3,a)") "  3 ", (110-kp(i,14)), (110-kp(i,13)), " 3 " ! tail-b
            write(1,"(a,3x,i3,3x,i3,a)") "  4 ", (110-kp(i,13)), (110-kp(i,12)), " 4 " ! tail-c
            write(1,"(a,3x,i3,3x,i3,a)") "  5 ", (110-kp(i,12)), (110-kp(i,11)), " 5 " ! trailing
            write(1,"(a,3x,i3,3x,i3,a)") "  6 ", (110-kp(i,11)), (110-kp(i,10)), " 6 " ! caps
            write(1,"(a,3x,i3,3x,i3,a)") "  7 ", (110-kp(i,10)), (110-kp(i,9)),  " 7 " ! leading
            write(1,"(a,3x,i3,3x,i3,a)") "  8 ", (110-kp(i,9)),  (110-kp(i,8)),  " 8 " ! nose
            write(1,"(a,3x,i3,3x,i3,a)") "  9 ", (110-kp(i,8)),  (110-kp(i,7)),  " 7 " ! leaading
            write(1,"(a,3x,i3,3x,i3,a)") " 10 ", (110-kp(i,7)),  (110-kp(i,6)),  " 9 " ! caps
            write(1,"(a,3x,i3,3x,i3,a)") " 11 ", (110-kp(i,6)),  (110-kp(i,5)),  " 5 " ! trailing
            write(1,"(a,3x,i3,3x,i3,a)") " 12 ", (110-kp(i,5)),  (110-kp(i,4)),  " 4 " ! tail-c
            write(1,"(a,3x,i3,3x,i3,a)") " 13 ", (110-kp(i,4)),  (110-kp(i,3)),  " 3 " ! tail-b
            write(1,"(a,3x,i3,3x,i3,a)") " 14 ", (110-kp(i,3)),  (110-kp(i,2)),  " 2 " ! tail-a
            write(1,"(a,3x,i3,3x,i3,a)") " 15 ", (110-kp(i,2)),  (110-kp(i,1)),  " 1 " ! tail-v

            write(1,"(a)")               "Webs 2"
            write(1,"(a,3x,i3,3x,i3,a)") " 16 ", (110-kp(i,10)), (110-kp(i,7)),  " 10 " ! webs
            write(1,"(a,3x,i3,3x,i3,a)") " 17 ", (110-kp(i,6)),  (110-kp(i,11)), " 10 " ! webs
            write(1,"(a)")               "Loops 3"
            write(1,"(a)")               " 1   4 "
            write(1,"(a)")               " 7   8   9   -16"
            write(1,"(a)")               " 2   4 "
            write(1,"(a)")               " 6   16   10   17 "
            write(1,"(a)")               " 3   11 "
            write(1,"(a)")               " 1   2   3   4   5   -17   11   12   13   14   15"
         close (1)
      !--->

      call system("/home/seraf/hGAST/pre_hGAST.dir/smalt.exe <"//trim(adjustl(path_dir))//"/path.inp")
      call system("less "//trim(adjustl(path_dir))//"/other.dir/section_matrix_hGAST.dat >> "//trim(adjustl(path_dir))//"/other.dir/matrix_hGAST.dat")
   end do

   call system("cd /home/seraf/hGAST/pre_hGAST.dir/ ; ./model_create.exe <"//trim(adjustl(path_dir))//"/path.inp")

end program pre_hGAST
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
