
#include "denselib.h"
#include "lu_sparse.h"
/*#include "lumod.h"*/

#ifdef FORTIFY
# include "lp_fortify.h"
#endif


/* ---------------------------------------------------------------------------------- */
/* Section containing inverse management routines                                     */
/* ---------------------------------------------------------------------------------- */

/* Include routines common to inverse implementations */
#include "lp_BFP1.c"
#include "lp_BFP2.c"


/* MUST MODIFY */
char * BFP_CALLMODEL inv_name(void)
{
  return( "LUMOD v1.0" );
}


/* DON'T MODIFY */
MYBOOL BFP_CALLMODEL bfp_init(lprec *lp, int size, int delta)
{
  INVrec *lu;

  lp->invB = (INVrec *) calloc(1, sizeof(*(lp->invB)));
  lu = lp->invB;

  bfp_resize(lp, size);
  bfp_pivotalloc(lp, size+delta);

  lu->max_Bsize = 0;       /* The largest B matrix of user variables */
  lu->max_colcount = 0;    /* The maximum number of user columns in LU */
  lu->max_LUsize = 0;      /* The largest NZ-count of LU-files generated */
  lu->num_refact = 0;

  bfp_preparefactorization(lp);
  lu->extraP = MAXINT;
  lu->extraD = lp->infinite;
  lu->num_refact = 0;

  return(TRUE);
}


/* DON'T MODIFY */
void BFP_CALLMODEL bfp_pivotalloc(lprec *lp, int newsize)
{
  INVrec *lu;

  lu = lp->invB;

  newsize += BASE;
  lu->colindex = (int *) realloc(lu->colindex, newsize * sizeof(*(lu->colindex)));
  lu->theta = (REAL *) realloc(lu->theta, newsize * sizeof(*(lu->theta)));
}


/* MUST MODIFY */
void BFP_CALLMODEL bfp_resize(lprec *lp, int newsize)
{
  INVrec *lu;

  lu = lp->invB;

  /* Increment dimensionality since we put the objective row at the top */
  newsize = newsize + LU_DELTAROWS;
  lu->dimalloc = newsize;

  /* Allocate index tracker arrays, LU matrices and various work vectors */
  lu->value = (REAL *) realloc(lu->value, (newsize + BASE) * sizeof(*(lu->value)));

  /* Data specific to the inversion engine */
  lu->L = createMatrix(newsize, newsize, newsize);
  lu->U = createMatrix(newsize, newsize, newsize);
  lu->x = (double *) realloc(lu->x, (newsize + BASE) * sizeof(*(lu->x)));
  lu->y = (double *) realloc(lu->y, (newsize + BASE) * sizeof(*(lu->y)));
  lu->w = (double *) realloc(lu->w, (newsize + BASE) * sizeof(*(lu->w)));
  lu->dimcount = 0;

}


/* MUST MODIFY */
void BFP_CALLMODEL bfp_free(lprec *lp)
{
  INVrec *lu;

  lu = lp->invB;
  if(lu == NULL)
    return;

  /* General arrays */
  free(lu->colindex);
  free(lu->theta);
  free(lu->value);

  /* Data specific to the inversion engine */
  freeMatrix(lu->L);
  freeMatrix(lu->U);
  free(lu->x);
  free(lu->y);
  free(lu->w);

  free(lu);
  lp->invB = NULL;
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_preparefactorization(lprec *lp)
{
  int i, size;
  INVrec *lu;

  lu = lp->invB;

  /* Finish any outstanding business */
  bfp_invertend(lp);

  /* Clear the existing LU matrices */
  if(lu->num_refact>0)
    for(i = 1; i <= lu->dimcount; i++) {
      clearVector(lu->L->list[i-1], 0, 0);
      clearVector(lu->U->list[i-1], 0, 0);
    }

  /* Reset additional indicators */
  lu->num_refact++;
  lu->user_colcount = 0;

  /* Signal that we are reinverting */
  lu->is_dirty = AUTOMATIC;

  /* Re-initialize LU to identity (I) */
  size = lu->dimalloc;
  MEMCLEAR(lu->x, size + BASE);
  MEMCLEAR(lu->y, size + BASE);
  lu->dimcount = 0;
  lu->colchange = 0;
  for(i = 1; i<=size; i++) {
    lu->x[i-1] = 0;
    lu->y[i-1] = 0;
    lu->x[i]   = 1;
    lu->y[i]   = 1;
    lu->colindex[i] = i - LU_DELTAROWS;  /* Set column index to lp_solve index (0=dummy) */
    LUmod(1, i,                          /* Add row & column; current row size is <i> */
          i, i,                          /* Indexes of the row and colum to be added */
          lu->L, lu->U,                  /* The L and U matrices to be updated */
          lu->x, lu->y,                  /* The row and column vectors to be added */
          lu->w);                        /* Work vector */
    lu->dimcount++;
  }

  return(0);

}


/* DON'T MODIFY */
void BFP_CALLMODEL bfp_finishfactorization(lprec *lp)
{
  INVrec *lu;

  lu = lp->invB;

  lu->max_colcount = my_max(lu->max_colcount, bfp_colcount(lp));
  lu->max_LUsize = my_max(lu->max_LUsize, bfp_nonzeros(lp, FALSE));

  /* Signal that we done reinverting */
  lu->is_dirty = FALSE;
  lp->justInverted = TRUE;
  lp->doInvert = FALSE;

  /* Recompute the RHS ( Ref. lp_solve inverse logic and Chvatal p. 121 ) */
  if(inv_colcount(lp) > 0)
    recompute_solution(lp, AUTOMATIC);
  lp->doRebase = FALSE;

  /* Store information about the current inverse */
  lu->extraP = lp->Extrap;
  lu->extraD = lp->Extrad;
  lu->num_pivots = 0;

}


/* DON'T MODIFY */
int BFP_CALLMODEL bfp_colcount(lprec *lp)
{
  return(lp->invB->user_colcount);
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_nonzeros(lprec *lp, MYBOOL maximum)
{
  if(maximum == TRUE)
    return(lp->invB->max_LUsize);
  else if(maximum == AUTOMATIC)
    return(lp->invB->max_Bsize);
  else
    return(NZcountMatrix(lp->invB->L)+NZcountMatrix(lp->invB->U));
}


/* MUST MODIFY (or ignore) */
int BFP_CALLMODEL bfp_memallocated(lprec *lp)
{
  return(0);
}


/* DON'T MODIFY */
int BFP_CALLMODEL bfp_usercolumn(lprec *lp, int datacolumn)
/* Find user column index of a given data column */
{
  if(datacolumn >= lp->invB->dimcount)
    return( -1 );
  else
    return( lp->invB->colindex[datacolumn+LU_DELTAROWS] );
}


/* DON'T MODIFY */
int BFP_CALLMODEL bfp_datacolumn(lprec *lp, int usercolumn)
/* Get data column containing a particular user column */
{
  int i;
  INVrec *lu;

  lu = lp->invB;

  for(i = 1; i <= lu->dimcount-LU_DELTAROWS; i++)
    if(inv_usercolumn(lp, i) == usercolumn)
      return(i);
  return(-1);
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_preparepivot(lprec *lp, int row_nr, int col_nr, REAL *pcol, MYBOOL *frow)
{
  int  i, rowalt;
  REAL hold, test;

  fsolve(lp, col_nr, pcol, lp->epsvalue, 1.0, TRUE);

  if(frow == NULL) {
    if(fabs(pcol[row_nr]) < lp->epsvalue)
      return( -1 );
    else
      return( row_nr );
  }

  /* Find largest coefficient row available to be pivoted;
     better numerical stability, but requires more CPU time - KE */
  row_nr  = lp->rows + 1;
  rowalt = 0;
  hold  = -lp->infinite;
  for(i = 1; i <= lp->rows; i++)
    /* Only examine rows that currently contain slacks;
       (i.e. the positions that are not already taken by another pivot) */
    if(!frow[i]) {

    test = fabs(pcol[i]);
    if(test > hold) {                  /* Largest absolute value */
      hold = test;
      row_nr = i;
    }
    if(/*!lp->scaling_used && */
       (rowalt == 0) && (test == 1))   /* Unit row */
      rowalt = i;

  }

  if((row_nr > lp->rows) || (hold < lp->epsvalue))
    row_nr = -1;
  else if(rowalt != 0 && (fabs(log10(hold)) <= 1))
    row_nr = rowalt;

  /* Nothing to prepare, simply return the ingoing target column to be replaced */
  return( row_nr );
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_pivotrow(lprec *lp, int datacolumn)
/* Get user pivot-out row for a particular data column */
{
  return( datacolumn );
}


/* DON'T MODIFY */
LREAL BFP_CALLMODEL bfp_putcolumn(lprec *lp, int row_nr, int col_nr, REAL *pcol)
/* Was condensecol() in versions of lp_solve before 4.0.1.8 - KE */
{
  LREAL  pivValue;
  INVrec *lu;

  lu = lp->invB;

  /* Store the incoming pivot value for RHS update purposes */
  lu->colindex[0] = col_nr;  /* The index of the new data column */
  lu->colchange   = row_nr;  /* The basis column to be replaced */
  if(pcol == NULL)
    pivValue = 0;
  else
    pivValue = pcol[row_nr];
  lu->theta[0] = pivValue;

  /* Save reference to the elimination vector */
  lu->pcol = pcol;

  /* Set completion status; but hold if we are reinverting */
  if(lu->is_dirty != AUTOMATIC)
    lu->is_dirty = TRUE;

  return( pivValue );
}


/* SUPPORT FUNCTION FOR FACTORIZATION */
void updateCounts(lprec *lp, MYBOOL *usedpos, int rownr, int colnr)
{
  usedpos[rownr] = AUTOMATIC;
  usedpos[colnr] = AUTOMATIC;
}


/* MUST MODIFY */
int BFP_CALLMODEL bfp_factorize(lprec *lp, int uservars, MYBOOL *usedpos)
{
  MATitem *matentry;
  REAL    *pcol, hold;
  int     *colnum, *rownum;
  int     *mdo = NULL;
  int     k, kk, kkk, i, j, v, rownr, colnr, numit;
  int     singularities = 0;

 /* Check if there is anyting to do */
  lp->invB->max_Bsize = my_max(lp->invB->max_Bsize, Bsize);
  if(uservars == 0)
    return(singularities);

 /* Allocate other necessary working arrays */
  pcol   = (REAL *) calloc(lp->rows + 1, sizeof(*pcol));
  colnum = (int *) calloc(lp->columns + 1, sizeof(*colnum));
  rownum = (int *) calloc(lp->rows + 1, sizeof(*rownum));

 /* Get net row and column entry counts for basic user variables
    (include OF row data, exclude eliminated rows due to slacks) */
  /* Get (net/un-eliminated) row and column entry counts for basic user variables */
  for(j = 1; j <= lp->columns; j++) {

    /* If it is a basic user variable...*/
    if(usedpos[lp->rows+j]) {
      i = lp->matA->col_end[j - 1];
      kk = lp->matA->col_end[j];

      /* Count relevant non-zero values */
      for(matentry = lp->matA->mat + i, numit = 0; i < kk; i++, matentry++, numit++) {
        k = (*matentry).row_nr;

        /* Check for objective row phase 1 adjustments;
           increment count only if necessary following modifyOF1() test */
        if((numit == 0) && !usedpos[0]) {
          if(k == 0)
            hold = (*matentry).value;
          else
            hold = 0;
          if(!modifyOF1(lp, lp->rows+j, &hold, 1.0)) {
            if(k == 0)
              continue;
          }
          else if(k > 0) {
            numit++;
            colnum[j]++;
            rownum[0]++;
          }
        }

       /* Exclude pre-eliminated rows due to basic slacks;
          this is a characteristic of the eta-model that presumes an
          initial basis with all slacks; i.e. an identity matrix */
        if(!usedpos[k]) {
          numit++;
          colnum[j]++;
          rownum[k]++;
        }
      }
    }
  }

 /* Determine the number of remaining pivots and allocate memory
    for the minimum degree ordering column index array */
  k = lp->rows;
  mdo = (int *) calloc(k + 1, sizeof(*mdo));

  /* Fill the mdo[] array with remaining full-pivot basic user variables */
  kk = 0;
  for(j = 1; j <= lp->columns; j++) {
    i = lp->rows + j;
    if(usedpos[i] == TRUE) {
      kk++;
      mdo[kk] = i;
    }
  }
  mdo[0] = kk;
  if(kk == 0)
    goto Cleanup;

 /* Calculate the approximate minimum degree column ordering */
#ifndef UseLegacyOrdering
/*  i = getMDO(lp, usedpos, mdo, NULL, TRUE); */
  i = getMDO(lp, usedpos, mdo, NULL, FALSE);
  if(i != 0) {
    report(lp, CRITICAL, "inv_refactorize: Internal error %d in minimum degree ordering routine", i);
    goto Cleanup;
  }
#endif

 /* Loop over all non-slack user basis columns, finding an appropriate
    unused slack column position to put the user column into */

  for(i = 1; i <= kk; i++) {

  /* Get the entering variable */
    colnr = mdo[i];

  /* Solve for and eliminate the entering column / variable */
    rownr = bfp_preparepivot(lp, 0, colnr, pcol, usedpos);
    if(rownr <= 0) {
     /* This column is singular!  Just skip it, leaving one of the
        slack variables basic in its place... (Source: Geosteiner changes!) */
      if(lp->spx_trace)
        report(lp, DETAILED, "inv_refactorize: Skipped singular column %d\n", colnr);
      singularities++;
      continue;
    }
    else {
      simpleiteration(lp, rownr, colnr, pcol);
    }
   /* Update occupancy states for pivot column/row */
    updateCounts(lp, usedpos, rownr, colnr);

   /* Check timeout and user abort again */
    if(userabort(lp, -1))
      goto Cleanup;
  }


Cleanup:
  if(mdo != NULL)
    free(mdo);
  if(pcol != NULL) {
    free(pcol);
    free(colnum);
    free(rownum);
  }

  return(singularities);
}


/* DON'T MODIFY */
REAL BFP_CALLMODEL bfp_pivotRHS(lprec *lp, LREAL theta, int row_nr, REAL *pcol)
/* Was rhsmincol(), ie. "rhs minus column" in versions of lp_solve before 4.0.1.8 - KE */
{
  int    i;
  LREAL  f;
  REAL   roundzero = lp->epsvalue;
  INVrec *lu;

  lu = lp->invB;

  if(row_nr > lp->rows + 1) {
    if (lp->spx_trace)
      report(lp, IMPORTANT, "inv_pivotRHS(1..%d): Called with invalid row %d; numerical instability!\n",
                            lp->rows, row_nr);
    lp->spx_status = NUMFAILURE;
    return( 0 );
  }

  if(pcol == NULL) {

    if(theta != 0) {
      for(i = 0; i <= lp->rows; i++) {
        f = lp->rhs[i] - theta*lu->pcol[i];
        my_round(f, roundzero);
        lp->rhs[i] = f;
      }
    }
    f = lu->theta[0];

  }
  else {

    for(i = 0; i<= lp->rows; i++) {
      f = pcol[i];
      if(f == 0) continue;
      lp->rhs[i] -= theta*f;
      my_round(lp->rhs[i], roundzero);
    }
    f = 0;
  }

  if(row_nr > 0)
    lp->rhs[row_nr] = theta;

  return( f );

}


/* MUST MODIFY */
void BFP_CALLMODEL bfp_finishupdate(lprec *lp, MYBOOL changesign)
/* Was addetacol() in versions of lp_solve before 4.0.1.8 - KE */
{
  int  i;
  INVrec *lu;

  lu = lp->invB;
  if(!lu->is_dirty) return;

  obtain_column(lp, lu->colindex[0], lu->value+LU_DELTAROWS, NULL);
  lu->value[0] = 0;
  if(changesign) {
    for(i = 1; i <= lu->dimcount; i++)
       lu->value[i] *= -1;
  }

  i = lu->colchange+LU_DELTAROWS;
  LUmod(2, lu->dimcount,
           0, i,
           lu->L, lu->U,
           lu->y, lu->value,
           lu->w);

  /* Update mapping arrays */
  if(lu->colindex[i] > lu->dimcount-LU_DELTAROWS)
    lu->user_colcount--;
  lu->colindex[i] = lu->colindex[0];
  if(lu->colindex[i] > lu->dimcount-LU_DELTAROWS)
    lu->user_colcount++;
  lu->theta[i] = lu->theta[0];
  lu->num_pivots++;
  lu->colchange = 0;
  if(lu->is_dirty != AUTOMATIC)
    lu->is_dirty = FALSE;
  lp->justInverted = FALSE;

} /* bfp_finishupdate */


void BFP_CALLMODEL bfp_ftran_normal(lprec *lp, REAL *pcol, REAL roundzero)
{
  int    i, n;
  INVrec *lu;

  lu = lp->invB;

  n = lu->dimcount;
  Lprod (1, n, lu->L, (pcol-LU_DELTAROWS), lu->w );
  Usolve(1, n, lu->U, lu->w);
  MEMCPY(pcol, (lu->w+LU_DELTAROWS), n);

  /* Transfer and round small values to zero */
  if(roundzero != 0)
  for(i = 1; i <= n; i++)
    my_round(pcol[i-LU_DELTAROWS], roundzero);

}
void BFP_CALLMODEL bfp_ftran_prepare(lprec *lp, REAL *pcol, REAL roundzero)
{
  bfp_ftran_normal(lp, pcol, roundzero);
}

void BFP_CALLMODEL bfp_btran_normal(lprec *lp, REAL *row, REAL roundzero)
{
  int    i, n;
  INVrec *lu;

  lu = lp->invB;

  n = lu->dimcount;
  MEMCPY((lu->w+LU_DELTAROWS), row, n);
  Usolve(2, n, lu->U, lu->w);
  Lprod (2, n, lu->L, lu->w, (row-LU_DELTAROWS));

  /* Round small values to zero */
  if(roundzero != 0)
    for(i = 1; i <= n; i++)
      my_round(row[i-LU_DELTAROWS], roundzero);
}


void BFP_CALLMODEL bfp_btran_double(lprec *lp, REAL *prow, REAL proundzero, REAL *drow, REAL droundzero)
{
  if(prow != NULL)
    btran(lp, prow, proundzero);
  if(drow != NULL)
    btran(lp, drow, droundzero);
}

