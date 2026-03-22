@echo off

REM This batch file compiles the lpsolve libraries with the Microsoft Visual C/C++ compiler under Windows

set src=../lp_MDO.c ../commonlib.c ../colamd/colamd.c ../lp_rlp.c ../lp_crash.c ../lp_etaPFI.c ../lp_Hash.c ../lp_lib.c ../lp_wlp.c ../lp_rlpt.c ../lp_matrix.c ../lp_mipbb.c ../lp_MPS.c ../lp_presolve.c ../lp_price.c ../lp_pricePSE.c ../lp_report.c ../lp_scale.c ../lp_simplex.c ../lp_SOS.c ../lp_utils.c ../lp_wlpt.c ../yacc_read.c

set c=cl

rc lpsolve.rc
%c% -I.. -I../colamd /LD /MD /O2 /Zp8 /Gz -D_WINDLL -D_USRDLL -DWIN32 -DYY_NEVER_INTERACTIVE -DPARSER_CPLEX -DPARSER_LP %src% lpsolve.res lp_solve.def -olpsolve5.dll
rem /link /LINK50COMPAT

if exist a.obj del a.obj
%c% -I.. -I../colamd /MT /O2 /Zp8 /Gd /c -DWIN32 -DYY_NEVER_INTERACTIVE -DPARSER_CPLEX -DPARSER_LP %src%
lib *.obj /OUT:liblpsolve5.lib

if exist a.obj del a.obj
%c% -I.. -I../colamd /MTd /O2 /Zp8 /Gd /c -DWIN32 -DYY_NEVER_INTERACTIVE -DPARSER_CPLEX -DPARSER_LP %src%
lib *.obj /OUT:liblpsolve5d.lib

if exist *.obj del *.obj
