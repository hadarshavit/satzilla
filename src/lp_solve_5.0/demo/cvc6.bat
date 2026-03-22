@echo off

REM This batch file compiles the demo program with the Microsoft Visual C/C++ compiler under Windows

set src=../lp_MDO.c ../commonlib.c ../colamd/colamd.c ../fortify.c ../lp_rlp.c ../lp_crash.c ../lp_etaPFI.c ../lp_Hash.c ../lp_lib.c ../lp_wlp.c ../lp_rlpt.c ../lp_matrix.c ../lp_mipbb.c ../lp_MPS.c ../lp_presolve.c ../lp_price.c ../lp_pricePSE.c ../lp_report.c ../lp_scale.c ../lp_simplex.c demo.c ../lp_SOS.c ../lp_utils.c ../lp_wlpt.c ../yacc_read.c
set c=cl

%c% -I.. -I../colamd /O2 /Zp8 /Gd -DWIN32 -DYY_NEVER_INTERACTIVE -DPARSER_CPLEX -DPARSER_LP %src% -odemo.exe

if exist *.obj del *.obj
