@echo off

REM This batch file compiles the lp_solve driver program with the GNU gcc compiler under Windows


set src=../lp_MDO.c ../commonlib.c ../colamd/colamd.c ../lp_rlp.c ../lp_crash.c ../lp_etaPFI.c ../lp_Hash.c ../lp_lib.c ../lp_wlp.c ../lp_rlpt.c ../lp_matrix.c ../lp_mipbb.c ../lp_MPS.c ../lp_presolve.c ../lp_price.c ../lp_pricePSE.c ../lp_report.c ../lp_scale.c ../lp_simplex.c lp_solve.c ../lp_SOS.c ../lp_utils.c ../lp_wlpt.c ../yacc_read.c

set c=gcc

%c% -I.. -I../colamd -O3 -DBFP_CALLMODEL=__stdcall -DYY_NEVER_INTERACTIVE -DPARSER_CPLEX -DPARSER_LP %src% -olp_solve.exe
