#ifndef HEADER_lp_GLPKLU
#define HEADER_lp_GLPKLU

/*  Modularized simplex inverse modules - w/interface for lp_solve v5.0+
   ----------------------------------------------------------------------------------
    Author:        Kjell Eikland
    Contact:       kjell.eikland@broadpark.no
    License terms: LGPL.

    Requires:      lp_lib.h   glplib.h      glplib1a.c
                                            glplib2.c
                                            glplib4.c
                              glpluf.h      glpluf.c
                              glpinv.h      glpinv.c

    Release notes:
    v1.0    1 February 2003        First implementation of the GLPK interface.

   ---------------------------------------------------------------------------------- */


/* Include libraries for this inverse system */
#include "lp_types.h"
#include "glpinv.h"


/* GLPKLU matrix defines */
#define MATINDEXBASE             1  /* Inversion engine index start for arrays */

/*#define SafeMemLU*/
/*#define DeltaSize     0*/
/*#define DeltaPivot    0*/
#define LU_START_SIZE        10000  /* Start size of LU; realloc'ed if needed */
#define INVDELTAROWS             1  /* Additional rows inserted at the top    */
#define DEF_MAXPIVOT            42  /* Maximum number of pivots before reinversion */
#define BASE                     1

/* typedef */ struct _INVrec
{
  int       status;                 /* Last operation status code */
  int       dimcount;               /* The actual number of LU rows/columns */
  int       dimalloc;               /* The allocated LU rows/columns size */
  int       user_colcount;          /* The number of user LU columns */
  INV       *LU;
  int       col_enter;              /* The full index of the entering column */
  int       col_leave;              /* The full index of the leaving column */
  int       col_pos;                /* The B column to be changed at the next update using data in value[.]*/
  REAL      *value;
  REAL      *pcol;                  /* Reference to the elimination vector */
  REAL      theta_enter;            /* Value of the entering column theta */

  int       max_Bsize;              /* The largest B matrix of user variables */
  int       max_colcount;           /* The maximum number of user columns in LU */
  int       max_LUsize;             /* The largest NZ-count of LU-files generated */
  int       num_refact;             /* Number of times the basis was refactored */
  double    time_refactstart;       /* Time since start of last refactorization-pivots cyle */
  double    time_refactnext;        /* Time estimated to next refactorization */
  int       num_pivots;             /* Number of pivots since last refactorization */
  int       num_singular;           /* The total number of singular updates */
  REAL      extraD;                 /* The dual objective function offset for the current inverse */
  MYBOOL    is_dirty;               /* Specifies if a column is incompletely processed */
  MYBOOL    force_refact;           /* Force refactorization at the next opportunity */
  MYBOOL    set_Bidentity;          /* Force B to be the identity matrix at the next refactorization */
} /* INVrec */;


#ifdef __cplusplus
namespace GLPKLU
extern "C" {
#endif

/* Put function headers here */
#include "lp_BFP.h"

#ifdef __cplusplus
 }
#endif

#endif
