
/*  Modularized simplex basis factorization module - w/interface for lp_solve v5.0+
   ----------------------------------------------------------------------------------
    Author:        Kjell Eikland
    Contact:       kjell.eikland@broadpark.no
    License terms: LGPL.

    Requires:      lp_lib.h   glplib.h      glplib1.c
                                            glplib2.c
                              glpluf.h      glpluf.c
                              glpinv.h      glpinv.c

    Release notes:
    v1.0    1 February 2003     First implementation of the GLPK interface.
    v1.0.1  23 May 2004         Moved mustrefact() function into the BFP structure.
                                BFPVERSION incremented to 2.

   ---------------------------------------------------------------------------------- */


/* Generic include libraries */
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

/* Include libraries for this factorization system */
#include "lp_lib.h"
#include "commonlib.h"
#include "lp_glpkLU.h"
#include "glpinv.h"

/* Include routines common to factorization implementations */
#include "lp_BFP1.c"
#include "lp_BFP2.c"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif


/* MUST MODIFY */
char * BFP_CALLMODEL bfp_name(void)
{
  return( "GLPKLU v4.4" );
}


/* MUST MODIFY */
void BFP_CALLMODEL bfp_resize(lprec *lp, int newsize)
{
  INVrec *lu;

  lu = lp->invB;

  /* Increment dimensionality since we put the objective row at the top */
  newsize = newsize + INVDELTAROWS;
  lu->dimalloc = newsize;

  /* Allocate index tracker arrays, LU matrices and various work vectors */
  lu->value = (REAL *) realloc(lu->value, (newsize + BASE) * sizeof(*(lu->value)));

  /* Data specific to the factorization engine */
  if(lu->LU != NULL)
    glp_inv_delete(lu->LU);
  if(newsize > 0)
    lu->LU = glp_inv_create(newsize, bfp_pivotmax(lp));
  else
    lu->LU = NULL;
  lu->dimcount = newsize;

}


/* MUST MODIFY */
void BFP_CALLMODEL bfp_free(lprec *lp)
{
  INVrec *lu;

  lu = lp->invB;
  if(lu == NULL)
    return;

  /* General arrays */
  FREE(lu->value);

  /* Data specific to the factorization engine */
  glp_inv_delete(lu->LU);

  FREE(lu);
  lp->invB = NULL;
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_nonzeros(lprec *lp, MYBOOL maximum)
{
  if(maximum == TRUE)
    return(lp->invB->max_LUsize);
  else if(maximum == AUTOMATIC)
    return(lp->invB->max_Bsize);
  else {
    int count;
    count = lp->invB->LU->luf->nnz_v + lp->invB->LU->luf->nnz_v;
    return(count);
  }
}


/* MUST MODIFY (or ignore) */
int BFP_CALLMODEL bfp_memallocated(lprec *lp)
{
  return(0);
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_preparefactorization(lprec *lp)
{
  INVrec *lu;

  lu = lp->invB;

  /* Finish any outstanding business */
  if(lu->is_dirty == AUTOMATIC)
    bfp_finishfactorization(lp);

  /* Reset additional indicators */
  lu->force_refact = FALSE;
  lu->num_refact++;
  lu->user_colcount = 0;
  lu->col_pos = 0;

  /* Signal that we are reinverting */
  lu->is_dirty = AUTOMATIC;

 /* Set time of start of current refactorization cycle */
  lu->time_refactstart = timeNow();
  lu->time_refactnext  = 0;

  return(0);

}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_factorize(lprec *lp, int uservars, int Bsize, MYBOOL *usedpos)
{
  static double piv_tol[1+3] = { 0.00, 0.10, 0.30, 0.70 };
  int    n, i;
  int    singularities = 0;

 /* Update dimensions */
  lp->invB->max_Bsize = my_max(lp->invB->max_Bsize, Bsize+(1+lp->rows-uservars));

 /* Refactorize the entire basis matrix */
  for(n = 1; n <= 3; n++) {
    lp->invB->LU->luf->piv_tol = piv_tol[n];
    i = glp_inv_decomp(lp->invB->LU, lp, lp->get_basiscolumn);
    if(i == 0) break;
  }
  if(i == 0)
    lp->invB->user_colcount = uservars;
  else
    lp->report(lp, CRITICAL, "bfp_factorize: Could not refactorize basis (error code %d)\n", i);

  lp->invB->num_singular += singularities;    /* The total number of singular updates */
  return(singularities);
}


/* MUST MODIFY */
MYBOOL BFP_CALLMODEL bfp_finishupdate(lprec *lp, MYBOOL changesign)
/* Was addetacol() in versions of lp_solve before 4.0.1.8 - KE */
{
  int    i, status;
  INVrec *lu;

  lu = lp->invB;
  if(!lu->is_dirty)
    return( FALSE );

/*  if(changesign) */  /* Not reliable; other operations overwrite internal GLPK data */
  {
    lp->get_lpcolumn(lp, lu->col_enter, lu->value+INVDELTAROWS, NULL, NULL);
    lu->value[0] = 0;
    glp_inv_ftran(lu->LU, lu->value, TRUE);
  }

  i = lu->col_pos+INVDELTAROWS;
  status = glp_inv_update(lu->LU, i);

  /* Update mapping arrays */
  if(lu->col_leave > lu->dimcount-INVDELTAROWS)
    lu->user_colcount--;
  if(lu->col_enter > lu->dimcount-INVDELTAROWS)
    lu->user_colcount++;
  lu->num_pivots++;
  lu->col_pos = 0;
  if(lu->is_dirty != AUTOMATIC)
    lu->is_dirty = FALSE;
  lp->justInverted = FALSE;

  /* Check for update errors and do any corrective action */
  if(status != 0) {
    lu->num_singular++;    /* The total number of singular updates */
    /* Error codes are:

       1 - the adjacent basis matrix is structurally singular, since after
           changing the j-th column of the matrix V by the new column (see
           the algorithm below) the case k1 > k2 occured;
       2 - the factorization is inaccurate, since after transforming the
           k2-th row of the matrix U = P*V*Q, the diagonal element u[k2,k2]
           is zero or close to zero;
       3 - maximal number of updates has been reached;
       4 - overflow of the sparse vector area.

       In all cases, the remedial action is to refactorize/reinvert! */

    lp->invert(lp, INITSOL_USEZERO);
  }
  return( TRUE );

} /* bfp_finishupdate */


void BFP_CALLMODEL bfp_ftran_local(lprec *lp, REAL *pcol, MYBOOL prepareupdate)
{
  int    n;
  INVrec *lu;

  lu = lp->invB;
  if(!lu->LU->valid)
    lp->invert(lp, INITSOL_USEZERO);

  n = lu->dimcount;
#ifdef SafeMemLU
  MEMCPY(lu->value+INVDELTAROWS, pcol, lu->dimcount);
  glp_inv_ftran(lu->LU, lu->value, prepareupdate);
  MEMCPY(pcol, lu->value+INVDELTAROWS, lu->dimcount);
#else
  glp_inv_ftran(lu->LU, (pcol-INVDELTAROWS), prepareupdate);
#endif

}


void BFP_CALLMODEL bfp_ftran_normal(lprec *lp, REAL *pcol, int *nzidx)
{
  REAL   roundzero = lp->epsmachine;

  bfp_ftran_local(lp, pcol, FALSE);

  /* Transfer and round small values to zero */
  if(roundzero != 0) {
    int i, n;
    n = lp->invB->dimcount;
    for(i = 1; i <= n; i++)
      my_roundzero(pcol[i-INVDELTAROWS], roundzero);
  }
}


void BFP_CALLMODEL bfp_ftran_prepare(lprec *lp, REAL *pcol, int *nzidx)
{
  bfp_ftran_local(lp, pcol, TRUE);
}


void BFP_CALLMODEL bfp_btran_normal(lprec *lp, REAL *row, int *nzidx)
{
  int    n;
  INVrec *lu;
  REAL   roundzero = lp->epsmachine;

  lu = lp->invB;
  if(!lu->LU->valid)
    lp->invert(lp, INITSOL_USEZERO);

  n = lu->dimcount;
#ifdef SafeMemLU
  MEMCPY(lu->value+INVDELTAROWS, row, lu->dimcount);
  glp_inv_btran(lu->LU, lu->value);
  MEMCPY(row, lu->value+INVDELTAROWS, lu->dimcount);
#else
  glp_inv_btran(lu->LU, (row-INVDELTAROWS));
#endif

  /* Round small values to zero */
  if(roundzero != 0) {
    int i;
    for(i = 1; i <= n; i++)
      my_roundzero(row[i-INVDELTAROWS], roundzero);
  }
}
