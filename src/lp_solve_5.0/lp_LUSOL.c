
/*  Modularized simplex basis factorization module - w/interface for lp_solve v5.0+
   ----------------------------------------------------------------------------------
    Author:        Kjell Eikland
    Contact:       kjell.eikland@broadpark.no
    License terms: LGPL.

    Requires:      lusol.h, lp_lib.h, myblas.h

    Release notes:
    v2.0.0  1 March 2004        First implementation of the LUSOL v2.0 C translation.
    v2.0.1  1 April 2004        Added singularity recovery and fast/reuse update logic.
    v2.0.2  23 May 2004         Moved mustrefact() function into the BFP structure.
                                BFPVERSION incremented to 2.

   ---------------------------------------------------------------------------------- */

/* Generic include libraries */
#include <malloc.h>
#include <string.h>
#include "lp_lib.h"

/* Include libraries for this factorization system */
#include "myblas.h"
#include "commonlib.h"
#include "lp_LUSOL.h"
#include "lusol.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif


/* Include routines common to factorization engine implementations */
#include "lp_BFP1.c"
#include "lp_BFP2.c"


/* MUST MODIFY */
char * BFP_CALLMODEL bfp_name(void)
{
  return( "LUSOL v2.0" );
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
  allocREAL(lp, &(lu->value), newsize+MATINDEXBASE, AUTOMATIC);

  /* Data specific to the factorization engine */
  if(lu->LUSOL != NULL) {
    if(newsize > 0)
      LUSOL_sizeto(lu->LUSOL, newsize, newsize, 0);
    else {
      LUSOL_free(lu->LUSOL);
      lu->LUSOL = NULL;
    }
  }
  else if(newsize > 0) {
    int  asize;
    REAL bsize;
    
    lu->LUSOL = LUSOL_create(NULL, 0, LUSOL_PIVOT_TPP, bfp_pivotmax(lp)*0);
#ifdef LUSOL_UseBLAS
/*    if(fileSearchPath("PATH", "myBLAS.DLL", NULL) && load_BLAS("myBLAS")) */
    if(is_nativeBLAS() && load_BLAS("myBLAS"))
      lp->report(lp, NORMAL, "Optimized BLAS was successfully loaded for LUSOL.\n");
#endif

    /* Try to minimize memory allocation if we have a large number of unit columns */
    bsize = (REAL) lp->get_nonzeros(lp);
    if(newsize > lp->columns)
      bsize += newsize;
    else  
      bsize = bsize/lp->columns*newsize;
    /* Add a "reasonable" delta to allow for B and associated factorizations 
       that are denser than average; this makes reallocations less frequent.
       Values between 1.2 and 1.5 appear to be reasonable. */
    asize = (int) (bsize*MAX_DELTAFILLIN*1.3333);
    LUSOL_sizeto(lu->LUSOL, newsize, newsize, asize);
  }
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
  LUSOL_free(lu->LUSOL);

  FREE(lu);
  lp->invB = NULL;
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_nonzeros(lprec *lp, MYBOOL maximum)
{
  INVrec *lu;

  lu = lp->invB;
  if(maximum == TRUE)
    return(lu->max_LUsize);
  else if(maximum == AUTOMATIC)
    return(lu->max_Bsize);
  else
    return(lu->LUSOL->luparm[LUSOL_IP_NONZEROS_L0]+lu->LUSOL->luparm[LUSOL_IP_NONZEROS_U0]);
/*    return(lu->LUSOL->luparm[LUSOL_IP_NONZEROS_ROW]); */
}


/* MUST MODIFY (or ignore) */
int BFP_CALLMODEL bfp_memallocated(lprec *lp)
{
  int      mem;
  LUSOLrec *LUSOL = lp->invB->LUSOL;

  mem = sizeof(REAL) * (LUSOL->lena+LUSOL->maxm+LUSOL_RP_LASTITEM);
  mem += sizeof(int) * (2*LUSOL->lena+5*LUSOL->maxm+5*LUSOL->maxn+LUSOL_IP_LASTITEM);
  if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TCP)
    mem += sizeof(REAL) * LUSOL->maxn + 2*sizeof(REAL)*LUSOL->maxn;
  else if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TRP)
    mem += sizeof(REAL) * LUSOL->maxn;
  if(!LUSOL->luparm[LUSOL_IP_KEEPLU])
    mem += sizeof(REAL) * LUSOL->maxn;
  return( mem );
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_preparefactorization(lprec *lp)
{
  INVrec *lu;

  lu = lp->invB;
  
  /* Finish any outstanding business */
  if(lu->is_dirty == AUTOMATIC)
    bfp_finishfactorization(lp);

  /* Clear the existing LU matrices - specific for the factorization engine */
  LUSOL_clear(lu->LUSOL, TRUE);

  /* Reset additional indicators */
  lu->force_refact = FALSE;
  lu->num_refact++;
  lu->user_colcount = 0;
  lu->col_pos = 0;

  /* Signal that we are refactorizing */
  lu->is_dirty = AUTOMATIC;

  /* Set time of start of current refactorization cycle */
  lu->time_refactstart = timeNow();
  lu->time_refactnext  = 0;
  
  return(0);

}


/* LOCAL HELPER ROUTINE - Replace a basis column with corresponding slack */
int bfp_LUSOLsetcolumn(lprec *lp, int posnr, int colnr)
{
  int nz, inform;

  nz = lp->get_lpcolumn(lp, colnr, lp->invB->LUSOL->w + INVDELTAROWS, NULL, NULL);
  inform = LUSOL_replaceColumn(lp->invB->LUSOL, posnr, lp->invB->LUSOL->w);
  return( inform );
}


/* LOCAL HELPER ROUTINE - force the basis to be the identity matrix */
int bfp_LUSOLidentity(lprec *lp, int *rownum)
{
  int    i, nz;
  INVrec *invB = lp->invB;

  /* Reset the factorization engine */
  LUSOL_clear(invB->LUSOL, TRUE);

  /* Add the basis columns */
  lp->invB->set_Bidentity = TRUE;
  for(i = 1; i <= invB->dimcount; i++) {
    nz = lp->get_basiscolumn(lp, i, rownum, invB->value);
    LUSOL_loadColumn(invB->LUSOL, rownum, i, invB->value, nz);
  }
  lp->invB->set_Bidentity = FALSE;

  /* Factorize */
  LU1FAC(invB->LUSOL, &i);

  return( i );
}


/* LOCAL HELPER ROUTINE */
int bfp_LUSOLfactorize(lprec *lp, MYBOOL *usedpos, int *rownum, int *singular)
{
  int    i, j, nz;
  INVrec *invB = lp->invB;

  /* Handle normal, presumed nonsingular case */
  if(singular == NULL) {

  /* Optionally do a symbolic minimum degree ordering; 
     not that slack variables should not be processed */
/*#define UsePreprocessMDO*/
#ifdef UsePreprocessMDO
    int *mdo;
    mdo = bfp_createMDO(lp, usedpos, lp->rows, TRUE);
    if(mdo != NULL) {
      for(i = 1; i <= lp->rows; i++)
        lp->set_basisvar(lp, i, mdo[i]);
      FREE(mdo);
    }
#endif    

    /* Reset the factorization engine */
    LUSOL_clear(invB->LUSOL, TRUE);

    /* Add the basis columns in the original order */
    for(i = 1; i <= invB->dimcount; i++) {
      nz = lp->get_basiscolumn(lp, i, rownum, invB->value);
      LUSOL_loadColumn(invB->LUSOL, rownum, i, invB->value, nz);
      if((i > INVDELTAROWS) && (lp->var_basic[i-INVDELTAROWS] > lp->rows))
        lp->invB->user_colcount++;
    }

    /* Factorize */
    LU1FAC(invB->LUSOL, &i);
  }

  /* Handle case where a column may be singular */
  else {
    LLrec *map;

    /* Reset the factorization engine */
    i = bfp_LUSOLidentity(lp, rownum);

    /* Build map of available columns */
    nz = createLink(lp->rows, &map, NULL);
    for(i = 1; i <= lp->rows; i++) {
      if(lp->var_basic[i] <= lp->rows)
        removeLink(map, i);
    }

    /* Rebuild the basis, column by column, while skipping slack columns */
    j = firstActiveLink(map);
    for(i = 1; i <= lp->rows; i++) {
      if(lp->var_basic[i] <= lp->rows)
        continue;
      nz = bfp_LUSOLsetcolumn(lp, j+INVDELTAROWS, lp->var_basic[i]);
      if(nz == LUSOL_INFORM_LUSUCCESS)
        lp->invB->user_colcount++;
      else {
        nz = bfp_LUSOLsetcolumn(lp, j+INVDELTAROWS, i);
        lp->set_basisvar(lp, i, i);
      }
      j = nextActiveLink(map, j);
    }

    /* Sort the basis list */
    MEMCOPY(rownum, lp->var_basic, lp->rows+1);
    sortByINT(lp->var_basic, rownum, lp->rows, 1, TRUE);

  }

  return( i );
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_factorize(lprec *lp, int uservars, int Bsize, MYBOOL *usedpos)
{
  int  kcol, inform, *rownum = NULL;
  int  singularities = 0;
  LUSOLrec *LUSOL = lp->invB->LUSOL;

 /* Set dimensions and create work array */
  lp->invB->max_Bsize = my_max(lp->invB->max_Bsize, Bsize+(1+lp->rows-uservars));
  kcol = lp->invB->dimcount;
  LUSOL->m = kcol;
  LUSOL->n = kcol;
  allocINT(lp, &rownum, kcol+1, FALSE);

 /* Reload B and factorize */
  inform = bfp_LUSOLfactorize(lp, usedpos, rownum, NULL);

 /* Do some checks */
#ifdef Paranoia
  if(uservars != lp->invB->user_colcount) {
    lp->report(lp, SEVERE, "bfp_factorize: User variable count reconciliation failed\n");
    return( singularities );
  }
#endif

  /* Check result */
  if(inform != LUSOL_INFORM_LUSUCCESS) {
    if(inform == LUSOL_INFORM_LUUNSTABLE) {
      lp->report(lp, NORMAL, "bfp_factorize: Factorization %d failed;\n%s\n",
                     lp->invB->num_refact, LUSOL_informstr(LUSOL, inform));
      inform = bfp_LUSOLfactorize(lp, usedpos, rownum, &singularities);
    }
    else {
      lp->report(lp, DETAILED, "bfp_factorize: %d singularities at refactorization %d, iteration %d.\n",
                               LUSOL->luparm[LUSOL_IP_SINGULARITIES], lp->invB->num_refact,
                               lp->total_iter+lp->current_iter);

      /* Make sure we do not tighten factorization pivot criteria too often, and simply
         accept the substitution of slack columns into the basis */ 
      if((lp->invB->num_singular+1) % TIGHTENAFTER == 0) {

        /* Otherwise first see if we should tighten the pivot threshold parameters,
           or alternatively switch to a more powerful pivoting model */
        if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TPP) {
          if(LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] > 10)
            LUSOL_tightenpivot(LUSOL);
          else {
            LUSOL_setpivotmodel(LUSOL, LUSOL_PIVOT_TRP);
            lp->report(lp, IMPORTANT, "bfp_factorize: Switched to %s to enhance stability!\n",
                                      LUSOL_pivotLabel(LUSOL));
          }
        }
        else if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVOT_TRP) {
          if(LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] > 1.1)
            LUSOL_tightenpivot(LUSOL);
          else
            lp->report(lp, IMPORTANT, "bfp_factorize: Extremely hard numerics; but will not tighten more\n");
        }
      }

      /* Try to restore a non-singular basis by substituting singular columns with slacks */
      while(inform == LUSOL_INFORM_LUSINGULAR) {
        kcol = LUSOL->luparm[LUSOL_IP_SINGULARINDEX];
        kcol = LUSOL_findColumnPosition(LUSOL, kcol);
        lp->set_basisvar(lp, kcol-INVDELTAROWS, kcol-INVDELTAROWS);
#if 0  /* Optimistic version that assumes that update following singularity is Ok      */
        inform = bfp_LUSOLsetcolumn(lp, kcol, kcol);
#else  /* Pessimistic version that assumes that update following singularity is not Ok */
        inform = bfp_LUSOLfactorize(lp, usedpos, rownum, NULL);
#endif
        singularities++;
      }
    }
  }

  /* Clean up before returning */
  FREE(rownum);

  lp->invB->num_singular += singularities;    /* The total number of singular updates */
  return(singularities);
}


/* MUST MODIFY */
MYBOOL BFP_CALLMODEL bfp_finishupdate(lprec *lp, MYBOOL changesign)
/* Was addetacol() in versions of lp_solve before 4.0.1.8 - KE */
{
  int      i, k, kcol;
  REAL     DIAG, VNORM;
  INVrec   *lu;
  LUSOLrec *LUSOL = lp->invB->LUSOL;

  lu = lp->invB;
  if(!lu->is_dirty)
    return( FALSE );
  if(lu->is_dirty != AUTOMATIC)
    lu->is_dirty = FALSE;

  /* Perform the update */
  k = lu->col_pos+INVDELTAROWS;
  lu->num_pivots++;
  if(lu->col_leave > lu->dimcount-INVDELTAROWS)
    lu->user_colcount--;
  if(lu->col_enter > lu->dimcount-INVDELTAROWS)
    lu->user_colcount++;
  kcol = lu->col_pos;
  lu->col_pos = 0;

  /* Do standard update */  
#ifdef LUSOLSafeFastUpdate      /* NB! Defined in lusol.h */
  if(!changesign)
    /* Execute the update using data prepared earlier */
    LU8RPC(LUSOL, LUSOL_UPDATE_OLDNONEMPTY, LUSOL_UPDATE_USEPREPARED,
           k, NULL, NULL, &i, &DIAG, &VNORM);
  else
#endif
  {
    /* Retrieve the data for the entering column (base 0) */
    i = lp->get_lpcolumn(lp, lu->col_enter, lu->value+INVDELTAROWS, NULL, NULL);
    lu->value[0] = 0;
    /* Execute the update */
    LU8RPC(LUSOL, LUSOL_UPDATE_OLDNONEMPTY, LUSOL_UPDATE_NEWNONEMPTY,
           k, lu->value, NULL, &i, &DIAG, &VNORM);
  }
  lp->justInverted = FALSE;
  
  if(i == LUSOL_INFORM_LUSUCCESS) {
    
    /* Check if we should refactorize based on accumulation of fill-in */
    DIAG  = LUSOL->luparm[LUSOL_IP_NONZEROS_L]+LUSOL->luparm[LUSOL_IP_NONZEROS_U];
    VNORM = LUSOL->luparm[LUSOL_IP_NONZEROS_L0]+LUSOL->luparm[LUSOL_IP_NONZEROS_U0];
#if 0    
    /* This is Michael Saunder's fixed parameter */
    VNORM *= MAX_DELTAFILLIN;
#else
    /* This is Kjell Eikland's dynamic error accumulation measure */
    VNORM *= pow(MAX_DELTAFILLIN, pow((0.5*LUSOL->nelem/VNORM), 0.25));
#endif
    lu->force_refact = (MYBOOL) ((DIAG > VNORM) && (lu->num_pivots > 20));

#if 0
    /* Additional KE logic to reduce maximum pivot count based on the density of B */
    if(!lu->force_refact) {
      VNORM = lp->rows+1;
      VNORM = 1.0 - pow((REAL) LUSOL->nelem/VNORM/VNORM, 0.2);
      lu->force_refact = (MYBOOL) (lu->num_pivots > VNORM*bfp_pivotmax(lp));
    }
#endif    
  }
  
  /* Handle errors */
  else {
    k = DETAILED;
    lp->report(lp, k, "bfp_finishupdate: Failed at iteration %d, pivot %d;\n%s\n",
                   lp->total_iter+lp->current_iter, lu->num_pivots, LUSOL_informstr(LUSOL, i));
    if(i == LUSOL_INFORM_ANEEDMEM) {       /* To compress used memory and realloc, if necessary */
      lp->invert(lp, INITSOL_USEZERO);
      if(i != LUSOL_INFORM_LUSUCCESS)
        lp->report(lp, NORMAL, "bfp_finishupdate: Insufficient memory at iteration %d;\n%s\n",
                       lp->total_iter+lp->current_iter, LUSOL_informstr(LUSOL, i));
    }
    else if(i == LUSOL_INFORM_RANKLOSS) {  /* To fix rank loss and clear cumulative errors */
#if 0
      /* This is test code to do pivot in slack BEFORE refactorization (pessimistic approach);
        assumes that LUSOL returns correct information about the source of the singularity */
/*      kcol = LUSOL->luparm[LUSOL_IP_SINGULARINDEX];*/
      kcol = LUSOL_findColumnPosition(LUSOL, kcol+INVDELTAROWS);
      lp->set_basisvar(lp, kcol-INVDELTAROWS, kcol-INVDELTAROWS);
#endif
      lp->invert(lp, INITSOL_USEZERO);
      i = LUSOL->luparm[LUSOL_IP_INFORM];
      if(i != LUSOL_INFORM_LUSUCCESS)
        lp->report(lp, NORMAL, "bfp_finishupdate: Recovery attempt unsuccessful at iteration %d;\n%s\n",
                       lp->total_iter+lp->current_iter, LUSOL_informstr(LUSOL, i));
      else
        lp->report(lp, k, "bfp_finishupdate: Correction or recovery was successful.\n");
    }
  }
  return( (MYBOOL) (i == LUSOL_INFORM_LUSUCCESS) );

} /* bfp_finishupdate */


/* MUST MODIFY */
void BFP_CALLMODEL bfp_ftran_normal(lprec *lp, REAL *pcol, int *nzidx)
{
  int    i;
  INVrec *lu;

  lu = lp->invB;

  /* Do the LUSOL ftran */
  i = LUSOL_ftran(lu->LUSOL, pcol-INVDELTAROWS, nzidx, FALSE);
  if(i != LUSOL_INFORM_LUSUCCESS) {
    lu->status = BFP_STATUS_ERROR;
    lp->report(lp, NORMAL, "bfp_ftran_normal: Failed at iteration %d, pivot %d;\n%s\n",
                   lp->total_iter+lp->current_iter, lu->num_pivots, LUSOL_informstr(lu->LUSOL, i));
  }
}


/* MAY MODIFY */
void BFP_CALLMODEL bfp_ftran_prepare(lprec *lp, REAL *pcol, int *nzidx)
{
  int    i;
  INVrec *lu;

  lu = lp->invB;

  /* Do the LUSOL ftran */
  i = LUSOL_ftran(lu->LUSOL, pcol-INVDELTAROWS, nzidx, TRUE);
  if(i != LUSOL_INFORM_LUSUCCESS) {
    lu->status = BFP_STATUS_ERROR;
    lp->report(lp, NORMAL, "bfp_ftran_prepare: Failed at iteration %d, pivot %d;\n%s\n",
                   lp->total_iter+lp->current_iter, lu->num_pivots, LUSOL_informstr(lu->LUSOL, i));
  }
}


/* MUST MODIFY */
void BFP_CALLMODEL bfp_btran_normal(lprec *lp, REAL *prow, int *nzidx)
{
  int    i;
  INVrec *lu;

  lu = lp->invB;

  /* Do the LUSOL btran */
  i = LUSOL_btran(lu->LUSOL, prow-INVDELTAROWS, nzidx);
  if(i != LUSOL_INFORM_LUSUCCESS) {
    lu->status = BFP_STATUS_ERROR;
    lp->report(lp, NORMAL, "bfp_btran_normal: Failed at iteration %d, pivot %d;\n%s\n",
                   lp->total_iter+lp->current_iter, lu->num_pivots, LUSOL_informstr(lu->LUSOL, i));
  }
}

