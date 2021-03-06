rem Before executing this script, start a command window, and start the ms vc mode compile like this
rem cmd /k ""C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat"" amd64
rem then execute this script from the root directory of shyft
pause
@if exist D:\Anaconda\64 set PYTHONROOT=D:\Anaconda\64
@if exist C:\Anaconda set PYTHONROOT=C:\Anaconda
@if exist C:\Andaconda\64 set PYTHONROOT=C:\Anaconda\64
@echo Assuming anaconda 64 is installed at this location: %PYTHONROOT%
@echo Generating the swig interface api.py
cd api\python
nmake /F nmakefile.swig_run
cd ..\..\shyft
@echo Compiling sources
cl /I".." /I"../.." /I"../../armadillo/include" /I"../../dlib" /I"../../boost" /I"%PYTHONROOT%\Include" /I"%PYTHONROOT%\lib\site-packages\numpy\core\include" /D "SWIG_PYTHON_INTERPRETER_NO_DEBUG"  /D "ARMA_DONT_PRINT_ERRORS" /D "ARMA_USE_CXX11" /D "BOOSTSERIAL" /D "BOOST_THREAD_USE_DLL" /D "BOOST_LIB_DIAGNOSTIC=1" /D "BOOST_NUMERIC_ODEINT_CXX11" /D  "BOOST_ALL_DYN_LINK=1" /D "_WINDOWS" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_WINDLL" /D "_MBCS" /D "_CRT_SECURE_NO_WARNINGS" /WX- /Zc:forScope /Gd /MD /EHsc /GS- /W1 /Gy /Zc:wchar_t /bigobj /Zi /Gm- /Ox /fp:precise /Fd"vc120.pdb"  /c  ../core/utctime_utilities.cpp  ../api/python/__init___wrap.cpp ../api/python/pt_gs_k_wrap.cpp ../api/python/pt_ss_k_wrap.cpp  ../core/sceua_optimizer.cpp ../core/dream_optimizer.cpp ../api/api.cpp
@echo Linking ___init__.pyd
link /OUT:api\___init__.pyd /Libpath:..\..\shyft-data\blaslapack /libpath:%PYTHONROOT%\libs /MANIFEST /NXCOMPAT /PDB:"_api.pdb" /DYNAMICBASE "blas_win64_MT.lib" "lapack_win64_MT.lib" "kernel32.lib" "user32.lib"  __init___wrap.obj  utctime_utilities.obj sceua_optimizer.obj dream_optimizer.obj api.obj /IMPLIB:"_api.lib" /DEBUG /DLL /MACHINE:X64 /OPT:REF /INCREMENTAL:NO /SUBSYSTEM:CONSOLE /MANIFESTUAC:"level='asInvoker' uiAccess='false'" /ManifestFile:"_api.pyd.intermediate.manifest" /OPT:ICF /ERRORREPORT:PROMPT /NOLOGO /TLBID:1 
@echo Linking _pt_gs_k.pyd
link /OUT:api\_pt_gs_k.pyd /Libpath:..\..\shyft-data\blaslapack /libpath:%PYTHONROOT%\libs /MANIFEST /NXCOMPAT /PDB:"_pt_gs_k.pdb" /DYNAMICBASE "blas_win64_MT.lib" "lapack_win64_MT.lib" "kernel32.lib" "user32.lib"  pt_gs_k_wrap.obj  utctime_utilities.obj sceua_optimizer.obj dream_optimizer.obj api.obj /IMPLIB:"_pt_gs_k.lib" /DEBUG /DLL /MACHINE:X64 /OPT:REF /INCREMENTAL:NO /SUBSYSTEM:CONSOLE /MANIFESTUAC:"level='asInvoker' uiAccess='false'" /ManifestFile:"_pt_gs_k.pyd.intermediate.manifest" /OPT:ICF /ERRORREPORT:PROMPT /NOLOGO /TLBID:1 
@echo Linking _pt_ss_k.pyd
link /OUT:api\_pt_ss_k.pyd /Libpath:..\..\shyft-data\blaslapack /libpath:%PYTHONROOT%\libs /MANIFEST /NXCOMPAT /PDB:"_pt_gs_k.pdb" /DYNAMICBASE "blas_win64_MT.lib" "lapack_win64_MT.lib" "kernel32.lib" "user32.lib"  pt_ss_k_wrap.obj  utctime_utilities.obj sceua_optimizer.obj dream_optimizer.obj api.obj /IMPLIB:"_pt_ss_k.lib" /DEBUG /DLL /MACHINE:X64 /OPT:REF /INCREMENTAL:NO /SUBSYSTEM:CONSOLE /MANIFESTUAC:"level='asInvoker' uiAccess='false'" /ManifestFile:"_pt_ss_k.pyd.intermediate.manifest" /OPT:ICF /ERRORREPORT:PROMPT /NOLOGO /TLBID:1 
cd ..
@echo done