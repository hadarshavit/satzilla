
#include <stdlib.h>
#include <stdio.h>
/*#include <memory.h>*/
#include <string.h>
#include <math.h>
#include "myblas.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif

/* ************************************************************************ */
/* Initialize BLAS interfacing routines                                     */
/* ************************************************************************ */
MYBOOL mustinitBLAS = TRUE;
#ifdef WIN32
  HINSTANCE hBLAS = NULL;
#else
  void      *hBLAS = NULL;
#endif


/* ************************************************************************ */
/* Define the BLAS interfacing routines                                     */
/* ************************************************************************ */

void init_BLAS()
{
  if(mustinitBLAS) {
    load_BLAS(NULL);
    mustinitBLAS = FALSE;
  }
}

MYBOOL is_nativeBLAS()
{
#ifdef LoadableBlasLib
  return( (MYBOOL) (hBLAS == NULL) );
#else
  return( TRUE );
#endif
}

MYBOOL load_BLAS(char *libname)
{
  MYBOOL result = TRUE;

#ifdef LoadableBlasLib
  if(hBLAS != NULL) {
  #ifdef WIN32
    FreeLibrary(hBLAS);
  #else
    dlclose(hBLAS);
  #endif
    hBLAS = NULL;
  }
#endif

  if(libname == NULL) {
    if(!mustinitBLAS && is_nativeBLAS())
      return( FALSE );
    BLAS_dscal = my_dscal;
    BLAS_dcopy = my_dcopy;
    BLAS_daxpy = my_daxpy;
    BLAS_dswap = my_dswap;
    BLAS_ddot = my_ddot;
    BLAS_idamax = my_idamax;
    BLAS_dload = my_dload;
    BLAS_dnormi = my_dnormi;
    if(mustinitBLAS)
      mustinitBLAS = FALSE;
  }
  else {
#ifdef LoadableBlasLib
  #ifdef WIN32
   /* Get a handle to the Windows DLL module. */
    hBLAS = LoadLibrary(libname);

   /* If the handle is valid, try to get the function addresses. */
    if(hBLAS == NULL) {
      unload_BLAS();
      result = FALSE;
    }
    else {
      BLAS_dscal  = (BLAS_dscal_func *)  GetProcAddress(hBLAS, "dscal");
      BLAS_dcopy  = (BLAS_dcopy_func *)  GetProcAddress(hBLAS, "dcopy");
      BLAS_daxpy  = (BLAS_daxpy_func *)  GetProcAddress(hBLAS, "daxpy");
      BLAS_dswap  = (BLAS_dswap_func *)  GetProcAddress(hBLAS, "dswap");
      BLAS_ddot   = (BLAS_ddot_func *)   GetProcAddress(hBLAS, "ddot");
      BLAS_idamax = (BLAS_idamax_func *) GetProcAddress(hBLAS, "idamax");
#if 0
      BLAS_dload  = (BLAS_dload_func *)  GetProcAddress(hBLAS, "dload");
      BLAS_dnormi = (BLAS_dnormi_func *) GetProcAddress(hBLAS, "dnormi");
#endif
    }
  #else
   /* First standardize UNIX .SO library name format. */
    char blasname[260], *ptr;

    strcpy(blasname, libname);
    if((ptr = strrchr(libname, '/')) == NULL)
      ptr = libname;
    else
      ptr++;
    blasname[(int) (ptr - libname)] = 0;
    if(strncmp(ptr, "lib", 3))
      strcat(blasname, "lib");
    strcat(blasname, ptr);
    if(strcmp(blasname + strlen(blasname) - 3, ".so"))
      strcat(blasname, ".so");

   /* Get a handle to the module. */
    hBLAS = dlopen(blasname, RTLD_LAZY);

   /* If the handle is valid, try to get the function addresses. */
    if(hBLAS == NULL) {
      unload_BLAS();
      result = FALSE;
    }
    else {
      BLAS_dscal  = (BLAS_dscal_func *)  dlsym(hBLAS, "dscal");
      BLAS_dcopy  = (BLAS_dcopy_func *)  dlsym(hBLAS, "dcopy");
      BLAS_daxpy  = (BLAS_daxpy_func *)  dlsym(hBLAS, "daxpy");
      BLAS_dswap  = (BLAS_dswap_func *)  dlsym(hBLAS, "dswap");
      BLAS_ddot   = (BLAS_ddot_func *)   dlsym(hBLAS, "ddot");
      BLAS_idamax = (BLAS_idamax_func *) dlsym(hBLAS, "idamax");
#if 0
      BLAS_dload  = (BLAS_dload_func *)  dlsym(hBLAS, "dload");
      BLAS_dnormi = (BLAS_dnormi_func *) dlsym(hBLAS, "dnormi");
#endif
    }
  #endif
#endif
    /* Do validation */
    if(result &&
       ((BLAS_dscal  == NULL) ||
        (BLAS_dcopy  == NULL) ||
        (BLAS_daxpy  == NULL) ||
        (BLAS_dswap  == NULL) ||
        (BLAS_ddot   == NULL) ||
        (BLAS_idamax == NULL) ||
        (BLAS_dload  == NULL) ||
        (BLAS_dnormi == NULL))
      ) {
      load_BLAS(NULL);
      result = FALSE;
    }
  }
  return( result );
}
MYBOOL unload_BLAS()
{
  return( load_BLAS(NULL) );
}


/* ************************************************************************ */
/* Now define the unoptimized local BLAS functions                          */
/* ************************************************************************ */
void daxpy( int n, double da, double *dx, int incx, double *dy, int incy)
{
  dx++;
  dy++;
  BLAS_daxpy( &n, &da, dx, &incx, dy, &incy);
}
/*void BLAS_CALLMODEL my_daxpy( int n, double da, double *dx, int incx, double *dy, int incy)*/
void BLAS_CALLMODEL my_daxpy( int *_n, double *_da, double *dx, int *_incx, double *dy, int *_incy)
{

/* constant times a vector plus a vector.
   uses unrolled loops for increments equal to one.
   jack dongarra, linpack, 3/11/78.
   modified 12/3/93, array[1] declarations changed to array[*] */

  int      i, ix, iy, m, mp1;
  register double rda;
  double   da = *_da;
  int      n = *_n, incx = *_incx, incy = *_incy;

  if (n <= 0) return;
  if (da == 0.0) return;

  dx--;
  dy--;
  ix = 1;
  iy = 1;
  if (incx < 0)
     ix = (-n+1)*incx + 1;
  if (incy < 0)
     iy = (-n+1)*incy + 1;
  rda = da;

/* CPU intensive loop; option to do pointer arithmetic */
#if defined DOFASTMATH
  {
    double *xptr, *yptr;
    for (i = 1, xptr = dx + ix, yptr = dy + iy;
         i <= n; i++, xptr += incx, yptr += incy)
      (*yptr) += rda*(*xptr);
    return;
  }
#endif

  if (incx==1 && incy==1) goto x20;

/* code for unequal increments or equal increments not equal to 1 */
  for (i = 1; i<=n; i++) {
     dy[iy]+= rda*dx[ix];
     ix+= incx;
     iy+= incy;
  }
  return;

/*  code for both increments equal to 1 */

/*  clean-up loop */
x20:
  m = n % 4;
  if (m == 0) goto x40;
  for (i = 1; i<=m; i++)
     dy[i]+= rda*dx[i];
  if(n < 4) return;
x40:
  mp1 = m + 1;
  for (i = mp1; i<=n; i=i+4) {
    dy[i]+= rda*dx[i];
    dy[i + 1]+= rda*dx[i + 1];
    dy[i + 2]+= rda*dx[i + 2];
    dy[i + 3]+= rda*dx[i + 3];
  }
}


/* ************************************************************************ */
void dcopy( int n, double *dx, int incx, double *dy, int incy)
{
  dx++;
  dy++;
  BLAS_dcopy( &n, dx, &incx, dy, &incy);
}

/*void BLAS_CALLMODEL my_dcopy (int n, double *dx, int incx, double *dy, int incy)*/
void BLAS_CALLMODEL my_dcopy (int *_n, double *dx, int *_incx, double *dy, int *_incy)
{

/* copies a vector, x, to a vector, y.
   uses unrolled loops for increments equal to one.
   jack dongarra, linpack, 3/11/78.
   modified 12/3/93, array[1] declarations changed to array[*] */

  int      i, ix, iy, m, mp1;
  int      n = *_n, incx = *_incx, incy = *_incy;

  if (n<=0) return;

  dx--;
  dy--;
  ix = 1;
  iy = 1;
  if (incx<0)
    ix = (-n+1)*incx + 1;
  if (incy<0)
    iy = (-n+1)*incy + 1;


/* CPU intensive loop; option to do pointer arithmetic */
#if defined DOFASTMATH
  {
    double *xptr, *yptr;
    for (i = 1, xptr = dx + ix, yptr = dy + iy;
         i <= n; i++, xptr += incx, yptr += incy)
      (*yptr) = (*xptr);
    return;
  }
#endif

  if (incx==1 && incy==1) goto x20;

/* code for unequal increments or equal increments not equal to 1 */

  for (i = 1; i<=n; i++) {
    dy[iy] = dx[ix];
    ix+= incx;
    iy+= incy;
  }
  return;

/* code for both increments equal to 1 */

/* version with fast machine copy logic (requires memory.h or string.h) */
x20:
#if defined DOFASTMATH
  memcpy(&dy[1], &dx[1], n*sizeof(double));
  return;
#endif

  m = n % 7;
  if (m == 0) goto x40;
  for (i = 1; i<=m; i++)
     dy[i] = dx[i];
  if (n < 7) return;

x40:
  mp1 = m + 1;
  for (i = mp1; i<=n; i=i+7) {
     dy[i] = dx[i];
     dy[i + 1] = dx[i + 1];
     dy[i + 2] = dx[i + 2];
     dy[i + 3] = dx[i + 3];
     dy[i + 4] = dx[i + 4];
     dy[i + 5] = dx[i + 5];
     dy[i + 6] = dx[i + 6];
  }
}


/* ************************************************************************ */

void dscal (int n, double da, double *dx, int incx)
{
  dx++;
  BLAS_dscal (&n, &da, dx, &incx);
}

/*void BLAS_CALLMODEL my_dscal (int n, double da, double *dx, int incx)*/
void BLAS_CALLMODEL my_dscal (int *_n, double *_da, double *dx, int *_incx)
{

/* Multiply a vector by a constant.

     --Input--
        N  number of elements in input vector(s)
       DA  double precision scale factor
       DX  double precision vector with N elements
     INCX  storage spacing between elements of DX

     --Output--
       DX  double precision result (unchanged if N.LE.0)

     Replace double precision DX by double precision DA*DX.
     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX. */

  int      i, ix, m, mp1;
  register double rda;
  double   da = *_da;
  int      n = *_n, incx = *_incx;

  if (n <= 0)
    return;
  rda = da;

  dx--;

/* Optionally do fast pointer arithmetic */
#if defined DOFASTMATH
  {
    double *xptr;
    for (i = 1, xptr = dx + 1; i <= n; i++, xptr += incx)
      (*xptr) *= rda;
    return;
  }
#endif

  if (incx == 1)
    goto x20;

/* Code for increment not equal to 1 */
  ix = 1;
  if (incx < 0)
    ix = (-n+1)*incx + 1;
  for(i = 1; i <= n; i++, ix += incx)
    dx[ix] *= rda;
  return;

/* Code for increment equal to 1. */
/* Clean-up loop so remaining vector length is a multiple of 5. */
x20:
  m = n % 5;
  if (m == 0) goto x40;
  for( i = 1; i <= m; i++)
    dx[i] *= rda;
  if (n < 5)
    return;
x40:
  mp1 = m + 1;
  for(i = mp1; i <= n; i += 5) {
    dx[i]   *= rda;
    dx[i+1] *= rda;
    dx[i+2] *= rda;
    dx[i+3] *= rda;
    dx[i+4] *= rda;
  }
}


/* ************************************************************************ */

double ddot (int n, double *dx, int incx, double *dy, int incy)
{
  dx++;
  dy++;
  return( BLAS_ddot (&n, dx, &incx, dy, &incy) );
}

/*double BLAS_CALLMODEL my_ddot (int n, double *dx, int incx, double *dy, int incy)*/
double BLAS_CALLMODEL my_ddot (int *_n, double *dx, int *_incx, double *dy, int *_incy)
{

/* forms the dot product of two vectors.
   uses unrolled loops for increments equal to one.
   jack dongarra, linpack, 3/11/78.
   modified 12/3/93, array[1] declarations changed to array[*] */

/*  long double dtemp; */
  register double dtemp;
  int      i, ix, iy, m, mp1;
  int      n = *_n, incx = *_incx, incy = *_incy;

  dtemp = 0.0;
  if (n<=0)
    return( (double) dtemp);

  dx--;
  dy--;
  ix = 1;
  iy = 1;
  if (incx<0)
     ix = (-n+1)*incx + 1;
  if (incy<0)
     iy = (-n+1)*incy + 1;

/* CPU intensive loop; option to do pointer arithmetic */

#if defined DOFASTMATH
  {
    double *xptr, *yptr;
    for (i = 1, xptr = dx + ix, yptr = dy + iy;
         i <= n; i++, xptr += incx, yptr += incy)
      dtemp+= (*yptr)*(*xptr);
    return(dtemp);
  }
#endif

  if (incx==1 && incy==1) goto x20;

/* code for unequal increments or equal increments not equal to 1 */

  for (i = 1; i<=n; i++) {
     dtemp+= dx[ix]*dy[iy];
     ix+= incx;
     iy+= incy;
  }
  return(dtemp);

/* code for both increments equal to 1 */

/* clean-up loop */

x20:
  m = n % 5;
  if (m == 0) goto x40;
  for (i = 1; i<=m; i++)
     dtemp+= dx[i]*dy[i];
  if (n < 5) goto x60;

x40:
  mp1 = m + 1;
  for (i = mp1; i<=n; i=i+5)
     dtemp+= dx[i]*dy[i] + dx[i + 1]*dy[i + 1] +
             dx[i + 2]*dy[i + 2] + dx[i + 3]*dy[i + 3] + dx[i + 4]*dy[i + 4];

x60:
  return(dtemp);
}


/* ************************************************************************ */

void dswap( int n, double *dx, int incx, double *dy, int incy )
{
  dx++;
  dy++;
  BLAS_dswap( &n, dx, &incx, dy, &incy );
}

/*void BLAS_CALLMODEL my_dswap( int n, double *dx, int incx, double *dy, int incy )*/
void BLAS_CALLMODEL my_dswap( int *_n, double *dx, int *_incx, double *dy, int *_incy )
{
  int      i, ix, iy, m, mp1, ns;
  double   dtemp1, dtemp2, dtemp3;
  int      n = *_n, incx = *_incx, incy = *_incy;

  if (n <= 0) return;

  dx--;
  dy--;
  ix = 1;
  iy = 1;
  if (incx < 0)
    ix = (-n+1)*incx + 1;
  if (incy < 0)
    iy = (-n+1)*incy + 1;

/* CPU intensive loop; option to do pointer arithmetic */
#if defined DOFASTMATH
  {
    double *xptr, *yptr;
    for (i = 1, xptr = dx + ix, yptr = dy + iy;
         i <= n; i++, xptr += incx, yptr += incy) {
      dtemp1 = (*xptr);
     (*xptr) = (*yptr);
     (*yptr) = dtemp1;
    }
    return;
  }
#endif

  if (incx == incy) {
    if (incx <= 0) goto x5;
    if (incx == 1) goto x20;
    goto x60;
  }

/* code for unequal or nonpositive increments. */
x5:
  for (i = 1; i<=n; i++) {
     dtemp1 = dx[ix];
     dx[ix] = dy[iy];
     dy[iy] = dtemp1;
     ix+= incx;
     iy+= incy;
  }
  return;

/* code for both increments equal to 1.
   clean-up loop so remaining vector length is a multiple of 3. */
x20:
  m = n % 3;
  if (m == 0) goto x40;
  for (i = 1; i<=m; i++) {
     dtemp1 = dx[i];
     dx[i] = dy[i];
     dy[i] = dtemp1;
  }
  if (n < 3) return;

x40:
  mp1 = m + 1;
  for (i = mp1; i<=n; i=i+3) {
     dtemp1 = dx[i];
     dtemp2 = dx[i+1];
     dtemp3 = dx[i+2];
     dx[i] = dy[i];
     dx[i+1] = dy[i+1];
     dx[i+2] = dy[i+2];
     dy[i] = dtemp1;
     dy[i+1] = dtemp2;
     dy[i+2] = dtemp3;
  }
  return;

/* code for equal, positive, non-unit increments. */
x60:
  ns = n*incx;
  for (i = 1; i<=ns; i=i+incx) {
     dtemp1 = dx[i];
     dx[i] = dy[i];
     dy[i] = dtemp1;
  }

}


/* ************************************************************************ */

void dload (int n, double da, double *dx, int incx)
{
  dx++;
  BLAS_dload (&n, &da, dx, &incx);
}

/*void BLAS_CALLMODEL my_dload (int n, double da, double *dx, int incx)*/
void BLAS_CALLMODEL my_dload (int *_n, double *_da, double *dx, int *_incx)
{
/* copies a scalar, a, to a vector, x.
   uses unrolled loops when incx equals one.

   To change the precision of this program, run the change
   program on dload.f
   Alternatively, to make a single precision version append a
   comment character to the start of all lines between sequential
      precision > double
   and
      end precision > double
   comments and delete the comment character at the start of all
   lines between sequential
      precision > single
   and
      end precision > single
   comments.  To make a double precision version interchange
    the append and delete operations in these instructions. */

  int      i, ix, m, mp1;
  double   da = *_da;
  int      n = *_n, incx = *_incx;

  if (n<=0) return;
  dx--;
  if (incx==1) goto x20;

/* code for incx not equal to 1 */

  ix = 1;
  if (incx<0)
     ix = (-n+1)*incx + 1;
  for (i = 1; i<=n; i++) {
     dx[ix] = da;
     ix+= incx;
  }
  return;

/* code for incx equal to 1 and clean-up loop */

x20:
  m = n % 7;
  if (m == 0) goto x40;
  for (i = 1; i<=m; i++)
     dx[i] = da;
  if (n < 7) return;

x40:
  mp1 = m + 1;
  for (i = mp1; i<=n; i=i+7) {
     dx[i] = da;
     dx[i + 1] = da;
     dx[i + 2] = da;
     dx[i + 3] = da;
     dx[i + 4] = da;
     dx[i + 5] = da;
     dx[i + 6] = da;
  }
}

/* ************************************************************************ */
int idamax( int n, double *x, int is )
{
  x++;
  return ( BLAS_idamax( &n, x, &is ) );
}

/*int BLAS_CALLMODEL my_idamax( int n, double *x, int is )*/
int BLAS_CALLMODEL my_idamax( int *_n, double *x, int *_is )
{
  double xmax, xtest;
  int    i, ii, imax = 0;
  int    n = *_n, is = *_is;

  if((n < 1) || (is <= 0))
    return(imax);
  imax = 1;
  if(n == 1)
    return(imax);

#if defined DOFASTMATH
  xmax = fabs(*x);
  for (i = 2, x += is; i <= n; i++, x += is) {
    xtest = fabs(*x);
    if (xtest > xmax) {
      xmax = xtest;
      imax = i;
    }
  }
#else
  x--;
  ii = 1;
  xmax = fabs(x[ii]);
  for (i = 2, ii+ = is; i <= n; i++, ii+ = is) {
    xtest = fabs(x[ii]);
	  if (xtest > xmax) {
      xmax = xtest;
		  imax = i;
    }
  }
#endif
  return(imax);
}


/* ************************************************************************ */
double dnormi( int n, double *x )
{
  x++;
  return( BLAS_dnormi( &n, x ) );
}

/*double BLAS_CALLMODEL my_dnormi( int n, double *x )*/
double BLAS_CALLMODEL my_dnormi( int *_n, double *x )
{
/* ===============================================================
   dnormi  returns the infinity-norm of the vector x.
   =============================================================== */
   int      j;
   register double hold;
   int      n = *_n;

   x--;
   hold = 0.0;
/*   for(j = 1; j <= n; j++) */
   for(j = n; j > 0; j--)
     hold = MAX( hold, fabs(x[j]) );

   return( hold );
}


/* ************************************************************************ */
/* Subvector and submatrix access routines (Fortran compatibility)          */
/* ************************************************************************ */

#ifndef UseMacroVector
int  subvec( int item)
{
  return( item-1 );
}
#endif

int submat( int nrowb, int row, int col)
{
  return( nrowb*(col-1) + subvec(row) );
}

int posmat( int nrowb, int row, int col)
{
  return( submat(nrowb, row, col)+BASE );
}

/* ************************************************************************ */
/* Randomization functions                                                  */
/* ************************************************************************ */

void randomseed(int seeds[])
/* Simply create some default seed values */
{
  seeds[1] = 123456;
  seeds[2] = 234567;
  seeds[3] = 345678;
}

void randomdens( int n, double *x, double r1, double r2, double densty, int *seeds )
{
/* ------------------------------------------------------------------
   random  generates a vector x[*] of random numbers
   in the range (r1, r2) with (approximate) specified density.
   seeds[*] must be initialized before the first call.
   ------------------------------------------------------------------ */

  int   i;
  double  *y;

  y = (double *) malloc(sizeof(*y) * (n+1));
  ddrand( n, x, 1, seeds );
  ddrand( n, y, 1, seeds );

  for (i = 1; i<=n; i++) {
     if (y[i] < densty)
        x[i] = r1  +  (r2 - r1) * x[i];
     else
        x[i] = 0.0;
  }
  free(y);
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void ddrand( int n, double *x, int incx, int *seeds )
{

/* ------------------------------------------------------------------
   ddrand fills a vector x with uniformly distributed random numbers
   in the interval (0, 1) using a method due to  Wichman and Hill.

   seeds[1..3] should be set to integer values
   between 1 and 30000 before the first entry.

   Integer arithmetic up to 30323 is required.

   Blatantly copied from Wichman and Hill 19-January-1987.
   14-Feb-94. Original version.
   30 Jun 1999. seeds stored in an array.
   30 Jun 1999. This version of ddrand.
   ------------------------------------------------------------------ */

  int    ix, xix;

  if (n < 1) return;

  for (ix = 1; ix<=1+(n-1)*incx; ix=ix+incx) {
     seeds[1]     = 171*(seeds[1] % 177) -  2*(seeds[1]/177);
     seeds[2]     = 172*(seeds[2] % 176) - 35*(seeds[2]/176);
     seeds[3]     = 170*(seeds[3] % 178) - 63*(seeds[3]/178);

     if (seeds[1] < 0) seeds[1] = seeds[1] + 30269;
     if (seeds[2] < 0) seeds[2] = seeds[2] + 30307;
     if (seeds[3] < 0) seeds[3] = seeds[3] + 30323;

	 x[ix]  = ((double) seeds[1])/30269.0 +
             ((double) seeds[2])/30307.0 +
             ((double) seeds[3])/30323.0;
     xix    = (int) x[ix];
	 x[ix]  = fabs(x[ix] - xix);
   }

}

