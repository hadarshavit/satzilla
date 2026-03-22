@echo off

REM This batch file compiles the lpsolve libraries with the GNU gcc compiler under Windows


set src=../lp_MDO.c ../commonlib.c ../colamd/colamd.c ../lp_rlp.c ../lp_crash.c ../lp_etaPFI.c ../lp_Hash.c ../lp_lib.c ../lp_wlp.c ../lp_rlpt.c ../lp_matrix.c ../lp_mipbb.c ../lp_MPS.c ../lp_presolve.c ../lp_price.c ../lp_pricePSE.c ../lp_report.c ../lp_scale.c ../lp_simplex.c ../lp_SOS.c ../lp_utils.c ../lp_wlpt.c ../yacc_read.c

set c=gcc

rem rc lpsolve.rc
%c% -I.. -I../colamd -O3 -shared -enable-stdcall-fixup -D_USRDLL -DWIN32 -DYY_NEVER_INTERACTIVE -DPARSER_CPLEX -DPARSER_LP %src% lp_solve.def -olpsolve5.dll

if exist a.o del a.o
%c% -I.. -I../colamd -O3 -c -DYY_NEVER_INTERACTIVE -DPARSER_CPLEX -DPARSER_LP %src%
ar rv liblpsolve5.a *.o

if exist *.o del *.o
