#ifndef HEADER_myblas
#define HEADER_myblas

/* ************************************************************************ */
/* BLAS function interface with local and external loadable versions        */
/* ************************************************************************ */

#define BASE         1
#define UseMacroVector
#define LoadableBlasLib


/* ************************************************************************ */
/* Include necessary libraries                                              */
/* ************************************************************************ */
#include "commonlib.h"
#ifdef LoadableBlasLib
  #ifdef WIN32
    #include <windows.h>
  #else
    #include <dlfcn.h>
  #endif
#endif


#ifdef __cplusplus
extern "C" {
#endif


/* ************************************************************************ */
/* BLAS functions                                                           */
/* ************************************************************************ */

#ifndef BLAS_CALLMODEL
  #ifdef WIN32
    #define BLAS_CALLMODEL _cdecl
  #else
    #define BLAS_CALLMODEL
  #endif
#endif

typedef void   (BLAS_CALLMODEL BLAS_dscal_func) (int *n, double *da, double *dx, int *incx);
typedef void   (BLAS_CALLMODEL BLAS_dcopy_func) (int *n, double *dx, int *incx,  double *dy, int *incy);
typedef void   (BLAS_CALLMODEL BLAS_daxpy_func) (int *n, double *da, double *dx, int *incx,  double *dy, int *incy);
typedef void   (BLAS_CALLMODEL BLAS_dswap_func) (int *n, double *dx, int *incx,  double *dy, int *incy);
typedef double (BLAS_CALLMODEL BLAS_ddot_func)  (int *n, double *dx, int *incx,  double *dy, int *incy);
typedef int    (BLAS_CALLMODEL BLAS_idamax_func)(int *n, double *x,  int *is);
typedef void   (BLAS_CALLMODEL BLAS_dload_func) (int *n, double *da, double *dx, int *incx);
typedef double (BLAS_CALLMODEL BLAS_dnormi_func)(int *n, double *x);

#ifndef __WINAPI
# ifdef WIN32
#  define __WINAPI WINAPI
# else
#  define __WINAPI
# endif
#endif

void init_BLAS();
MYBOOL is_nativeBLAS();
MYBOOL load_BLAS(char *libname);
MYBOOL unload_BLAS();

/* ************************************************************************ */
/* User-callable BLAS definitions (C base 1)                                */
/* ************************************************************************ */
void dscal ( int n, double da,  double *dx, int incx );
void dcopy ( int n, double *dx, int incx,   double *dy, int incy );
void daxpy ( int n, double da,  double *dx, int incx,   double *dy, int incy );
void dswap ( int n, double *dx, int incx,   double *dy, int incy );
REAL ddot  ( int n, double *dx, int incx,   double *dy, int incy );
int  idamax( int n, double *x,  int is );
void dload ( int n, double da,  double *dx, int incx );
REAL dnormi( int n, double *x );


/* ************************************************************************ */
/* Locally implemented BLAS functions (C base 0)                            */
/* ************************************************************************ */
void BLAS_CALLMODEL my_dscal ( int *n, double *da, double *dx, int *incx );
void BLAS_CALLMODEL my_dcopy ( int *n, double *dx, int *incx,  double *dy, int *incy );
void BLAS_CALLMODEL my_daxpy ( int *n, double *da, double *dx, int *incx,  double *dy, int *incy );
void BLAS_CALLMODEL my_dswap ( int *n, double *dx, int *incx,  double *dy, int *incy );
REAL BLAS_CALLMODEL my_ddot  ( int *n, double *dx, int *incx,  double *dy, int *incy );
int  BLAS_CALLMODEL my_idamax( int *n, double *x,  int *is );
void BLAS_CALLMODEL my_dload ( int *n, double *da, double *dx, int *incx );
REAL BLAS_CALLMODEL my_dnormi( int *n, double *x );


/* ************************************************************************ */
/* Function pointers for external BLAS library (C base 0)                   */
/* ************************************************************************ */
BLAS_dscal_func  *BLAS_dscal;
BLAS_dcopy_func  *BLAS_dcopy;
BLAS_daxpy_func  *BLAS_daxpy;
BLAS_dswap_func  *BLAS_dswap;
BLAS_ddot_func   *BLAS_ddot;
BLAS_idamax_func *BLAS_idamax;
BLAS_dload_func  *BLAS_dload;
BLAS_dnormi_func *BLAS_dnormi;


/* ************************************************************************ */
/* Subvector and submatrix access routines (Fortran compatibility)          */
/* ************************************************************************ */
#ifdef UseMacroVector
  #define subvec(item) (item - 1)
#else
  int subvec( int item );
#endif

int submat( int nrowb, int row, int col );
int posmat( int nrowb, int row, int col );


/* ************************************************************************ */
/* Randomization functions                                                  */
/* ************************************************************************ */
void randomseed(int *seeds);
void randomdens( int n, REAL *x, REAL r1, REAL r2, REAL densty, int *seeds);
void ddrand( int n, REAL *x, int incx, int *seeds );


#ifdef __cplusplus
}
#endif

#endif
