#
# --- hGAST -  GenUVP - VPM ---
#
#DIR0* :contains the source
#DIR   :contains the objects and the modules
DIR0_str   = ./source
DIR0_ctl   = ./source/ctrl
EXE0       = hGAST

   GITVER := $(shell git describe --dirty --always --tags)
#             DCONTROLLER_TYPE: 0:NREL , 1:DLL, 2:DTU_v1, 3:DTU_v2, 33:DTU_v2.3
#             OS              : 0:LINUX, 1:WINDOWS
#             ASCII           : 0:*.bin, 1:*.dat
   DIRECT  = -DCONTROLLER_TYPE=3  -DCONTROLLER_DLL=5 -DASCII=0 -DHAVE_OMP -DHAVE_MODAL -DOS=0 -DVER=\"$(GITVER)\" #-DHAVE_MPI -DHAVE_GENUVP
dbg ?=0
cmp ?=0
ifeq ($(cmp), 1)
#--- gfortran ---
   FC      = gfortran
   FLAGS0  = $(DIRECT) -llapack -lblas -fopenmp -cpp -ffree-line-length-none #-Wline-truncation -fwhole-file #-std=f2008 #-Wall -Wline-truncation -fwhole-file 
   mod     = -J
ifeq ($(dbg), 1)
#--- debug mode
#  FLAGS   = $(FLAGS0) -finit-real=nan -finit-integer=n #-ffpe-trap=invalid,zero,overflow,underflow,precision,denormal
   FLAGS   = $(FLAGS0) -fimplicit-none -Wall -Wextra -fcheck=all -fbacktrace #-pedantic -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file #-O2
   FLAGS1  = $(FLAGS)
   DIR_str = $(DIR0_str)/debug_gfo
   EXE     = $(EXE0)_gfo_debug
else
#--- release mode
   FLAGS   = $(FLAGS0) -O3 -march=native -fexternal-blas #-ffast-math -fforce-mem -fforce-add -funroll-all-loops #-ffree-line-length-0
   FLAGS1  = $(FLAGS)
   DIR_str = $(DIR0_str)/release_gfo
   EXE     = $(EXE0)_gfo
endif #dbg

else

#--- ifort ---
   FC      = ifort#mpif90
   FLAGS0  = $(DIRECT) -mkl -fpp #-qopenmp -static
   mod     = -module
ifeq ($(dbg), 1)
#--- debug mode
   FLAGS   = $(FLAGS0) -C -traceback -fpe0 #-g
   FLAGS1  = $(FLAGS)  -warn all #,errors # -warn general -warn interfaces -warn errors -ftrapuv #-warn stderrors
   DIR_str = $(DIR0_str)/debug_ifo
   EXE     = $(EXE0)_debug
else
#--- release mode
   FLAGS   = $(FLAGS0) -O3 -parallel
   FLAGS1  = $(FLAGS) #-warn all
   DIR_str = $(DIR0_str)/release_ifo
   EXE     = $(EXE0)
endif #dbg
endif #cmp
   VPATH   = $(DIR_str)
   OBJ_str = main_el.o raft.o      hydro.o    gelast.o    modal.o     init_el.o controls.o    basica_el.o \
             math_el.o librar_el.o eigen_el.o rotate_el.o rotmat_el.o inflow.o  writeout_el.o


.PHONY        : all clean cleanall complete install prepare remake uninstall

complete      : prepare $(EXE) install

$(EXE)        : $(OBJ_str)
		$(FC) $(patsubst %.o,$(DIR_str)/%.o,$(OBJ_str)) \
		$(mod) $(DIR_str) $(FLAGS) -o $(EXE)


main_el.o     : $(DIR0_str)/main_el.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

gelast.o      : $(DIR0_str)/gelast.f90 $(DIR0_str)/main_el.f90 $(DIR0_str)/truss.f90 $(DIR0_str)/module_truss.f90 $(DIR0_str)/truss_coupled.f90 #$(DIR0_str)/flap.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

modal.o       : $(DIR0_str)/modal.f90 $(DIR0_str)/main_el.f90 $(DIR0_str)/gast2modal_sb.f90 $(DIR0_str)/loads_integr.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

init_el.o     : $(DIR0_str)/init_el.f90 $(DIR0_str)/main_el.f90 $(DIR0_str)/jacket_init.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@
#		$(FC) -c $(FLAGS1) -C -traceback $(mod) $(DIR_str) $< -o $(DIR_str)/$@

basica_el.o   : $(DIR0_str)/basica_el.f90 $(DIR0_str)/main_el.f90 $(DIR0_str)/jacket.f90 $(DIR0_str)/stifcon.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

controls.o    : $(DIR0_str)/controls.f90 $(DIR0_str)/main_el.f90 $(DIR0_ctl)/oc3_control.f90    \
                                                                 $(DIR0_ctl)/controls3*.f90     \
                                                                 $(DIR0_ctl)/10MW_control*      \
                                                                 $(DIR0_ctl)/dtu_control_v2.f90 \
                                                                 $(DIR0_ctl)/dll_control.f90    \
                                                                 $(DIR0_ctl)/dll_interface.f90  \
                                                                 $(DIR0_str)/lib_ctrl.f90       \
                                                                 Makefile
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

math_el.o     : $(DIR0_str)/math_el.f90 $(DIR0_str)/main_el.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

librar_el.o   : $(DIR0_str)/librar_el.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

eigen_el.o    : $(DIR0_str)/eigen_el.f90 $(DIR0_str)/main_el.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

rotate_el.o   : $(DIR0_str)/rotate_el.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

rotmat_el.o   : $(DIR0_str)/rotmat_el.f90 $(DIR0_str)/main_el.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

raft.o        : $(DIR0_str)/raft.f90 $(DIR0_str)/main_el.f90 $(DIR0_str)/foilfs.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

inflow.o      : $(DIR0_str)/inflow.f90 $(DIR0_str)/anemos.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

writeout_el.o : $(DIR0_str)/writeout_el.f90 $(DIR0_str)/main_el.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@

hydro.o       : $(DIR0_str)/hydro.f90 $(DIR0_str)/main_el.f90 $(DIR0_str)/floater.f90 $(DIR0_str)/substruct.f90 $(DIR0_str)/uinflow_sea.f90
		$(FC) -c $(FLAGS1) $(mod) $(DIR_str) $< -o $(DIR_str)/$@


clean         :
		@rm -r -f $(DIR_str)/                             $(EXE)
		@reset
cleanall      :
		@rm -r -f $(DIR0_str)/debug* $(DIR0_str)/release* $(EXE0)*
		@reset
install       :
		@ln -s -f $(shell pwd)/$(EXE) ~/bin
prepare       :
		@mkdir -p $(DIR_str)
		@mkdir -p ~/bin
remake        : clean complete
uninstall     :
		@rm ~/bin/$(EXE0)*
all           :
		make -f Makefile cmp=0 dbg=0
		make -f Makefile cmp=0 dbg=1
#		make -f Makefile cmp=1 dbg=0
#		make -f Makefile cmp=1 dbg=1

