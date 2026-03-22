@echo off

REM This batch file compiles the demo program with the gnu gcc compiler under DOS/Windows

REM There are two ways to do that: use the lpsolve code directly or use that static library.

rem link lpsolve code with application
set src=../lp_MDO.c ../commonlib.c ../colamd/colamd.c ../fortify.c ../lp_rlp.c ../lp_crash.c ../lp_etaPFI.c ../lp_Hash.c ../lp_lib.c ../lp_wlp.c ../lp_rlpt.c ../lp_matrix.c ../lp_mipbb.c ../lp_MPS.c ../lp_presolve.c ../lp_price.c ../lp_pricePSE.c ../lp_report.c ../lp_scale.c ../lp_simplex.c ../lp_SOS.c ../lp_utils.c ../lp_wlpt.c ../yacc_read.c

rem statically link lpsolve library
rem set src=../lpsolve5/liblpsolve5.a

set c=gcc

%c% -I.. -I../colamd -O3 -DYY_NEVER_INTERACTIVE -DPARSER_CPLEX -DPARSER_LP demo.c %src% -odemo.exe
