rem mklvars32.bat [C:\Program Files\Intel\Compiler\11.0\061\fortran\mkl\tools\environment] not needed!
rem open fortran build environment command line
rem adjust time
cls
del *.obj *.mod *.pdb *.ilk hGAST.exe

rem                                          DCONTROLLER_TYPE:0:NREL, 1:DLL, 2:DTU_v1, 3:DTU_v2
ifort /exe:hGAST.exe /F100000000 /fast /fpp /DCONTROLLER_TYPE=0 /DCONTROLLER_DLL=2 /DASCII=1 /DOS=1 /DVER=0 /static     ^
   main_el.f90 raft.f90      hydro.f90    gelast.f90                  init_el.f90 controls.f90    basica_el.f90  ^
   math_el.f90 librar_el.f90 eigen_el.f90 rotate_el.f90 rotmat_el.f90 inflow.f90  writeout_el.f90                ^
   /link mkl_intel_c.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib

rem modal.f90  /DHAVE_MODAL
rem /libs:static /MT /Qmkl:sequential
rem /warn:all /traceback /debug:all
