
/*  Modularized external language interface module - w/interface for lp_solve v5.0+
    This module is the implemenation of the Ampl subset language named GNU MathProg.
   ----------------------------------------------------------------------------------
    Author:        Kjell Eikland
    Contact:       kjell.eikland@broadpark.no
    License terms: LGPL.

    Template used: lp_XLI2.c
    Requires:      lp_lib.h

    Release notes:
    v1.0.0  15 June 2004        First implementation.
    v1.0.1   9 July 2004        Added print hook routines to print via lp_solve 
                                report routine

   ---------------------------------------------------------------------------------- */

/* Generic include libraries */
#include <malloc.h>
#include <string.h>
#include "lp_lib.h"

/* Include libraries for this language system */
#include <math.h>
#include "glplib.h"
#include "glplpx.h"
#include "glpmpl.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif


/* Include routines common to language interface implementations */
#include "lp_XLI1.c"


/* CAN MODIFY */
char * XLI_CALLMODEL xli_name(void)
{
  return( "XLI_MathProg v1.0" );
}

static int _MathProg_hook_fn(lprec *lp, char *msg)
{
  char *str;

  if ((str = malloc(strlen(msg) + 1 + 1)) != NULL) {
    sprintf(str, "%s\n", msg);
    lp->report(lp, NORMAL, str);
    free(str);
    return(TRUE);
  }
  return(FALSE);
}

static int _MathProg_print_hook_fn(void *lp, char *msg)
{
  return(_MathProg_hook_fn((lprec *)lp, msg));
}

static int _MathProg_fault_hook_fn(void *lp, char *msg)
{
  return(_MathProg_hook_fn((lprec *)lp, msg));
}

/* CAN MODIFY */
MYBOOL XLI_CALLMODEL xli_readmodel(lprec *lp, char *model, char *data, char *options, int verbose)
{
  MPL    *mpl;
  int    ret, m, n, i, ii, j, *ndx, len, kind, type;
  double lb, ub, *val;
  char   *output = NULL;
  MYBOOL status = FALSE;

  /* create and initialize the translator database */
  mpl = mpl_initialize();

  /* set hook for print function */
  lib_set_fault_hook((void *)lp, &_MathProg_fault_hook_fn);
  lib_set_print_hook((void *)lp, &_MathProg_print_hook_fn);

  /* read model section and optional data section */
  ret = mpl_read_model(mpl, model, data != NULL);
  if (ret == 4) goto done;
  insist(ret == 1 || ret == 2);

  /* read data section, if necessary */
  if (data != NULL) {
    insist(ret == 1);
    ret = mpl_read_data(mpl, data);
    if (ret == 4) goto done;
    insist(ret == 2);
  }

  /* generate model */
  ret = mpl_generate(mpl, output);
  if (ret == 4) goto done;
  insist(ret == 3);

  /* set problem name and temporary row storage order */
  lp->set_lp_name(lp, mpl_get_prob_name(mpl));
  lp->set_add_rowmode(lp, TRUE);

  /* build columns (variables) */
  n = mpl_get_num_cols(mpl);
  for (j = 1; j <= n; j++) {

    /* set column name */
    lp->set_col_name(lp, j, mpl_get_col_name(mpl, j));

    /* set column kind */
    kind = mpl_get_col_kind(mpl, j);
    switch (kind) {
      case MPL_NUM: break;
      case MPL_INT:
      case MPL_BIN: lp->set_int(lp, j, TRUE);
                    break;
      default:      insist(kind != kind);
    }

    /* set column bounds */
    type = mpl_get_col_bnds(mpl, j, &lb, &ub);
    switch (type) {
      case MPL_FR: break;
      case MPL_LO: break;
      case MPL_UP: break;
      case MPL_DB: break;
      case MPL_FX: break;
      default:     insist(type != type);
    }
    if (kind == MPL_BIN) {
      if (type == MPL_FR || type == MPL_UP || lb < 0.0) lb = 0.0;
      if (type == MPL_FR || type == MPL_LO || ub > 1.0) ub = 1.0;
      type = MPL_DB;
    }
    if (type == MPL_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb))) {
      type = MPL_FX;
      if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
    }
    if(type == MPL_FR)
      lp->set_free(lp, j);
    else if(type == MPL_UP)
      lp->set_upbo(lp, j, ub);
    else if(type == MPL_LO)
      lp->set_lowbo(lp, j, lb);
    else {
      lp->set_upbo(lp, j, ub);
      lp->set_lowbo(lp, j, lb);
    }
  }

  /* allocate working arrays */
  ndx = malloc((1+n) * sizeof(*ndx));
  val = malloc((1+n) * sizeof(*val));

  /* build objective function (the first objective is used) */
  m = mpl_get_num_rows(mpl);
  for (i = 1; i <= m; i++) {

    kind = mpl_get_row_kind(mpl, i);
    if (kind == MPL_MIN || kind == MPL_MAX) {

      /* set optimization direction */
      lp->set_sense(lp, (MYBOOL) (kind == MPL_MAX) );

      /* set objective coefficients */
      len = mpl_get_mat_row(mpl, i, ndx, val);
      lp->set_obj_fnex(lp, len, val+1, ndx+1);

      /* set constant term */
      lp->set_rh(lp, 0, mpl_get_row_c0(mpl, i));

      /* set objective name */
      lp->set_row_name(lp, 0, mpl_get_row_name(mpl, i));

      break;
    }
  }

  /* build rows (constraints) */
  ii = 0;
  for (i = 1; i <= m; i++) {

    /* select constraint rows only */
    kind = mpl_get_row_kind(mpl, i);
    if (kind == MPL_MIN || kind == MPL_MAX)
      continue;

    /* set row bounds */
    type = mpl_get_row_bnds(mpl, i, &lb, &ub);
    switch (type) {
      case MPL_FR: break;
      case MPL_LO: break;
      case MPL_UP: break;
      case MPL_DB: break;
      case MPL_FX: break;
      default:     insist(type != type);
    }
    if(type == MPL_FR)
      continue;
    if (type == MPL_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb))) {
      type = MPL_FX;
      if (fabs(lb) <= fabs(ub))
        ub = lb;
      else
        lb = ub;
    }
    if(type == MPL_UP || type == MPL_DB) {
      type = ROWTYPE_LE;
      if(type == MPL_DB)
        lb = ub-lb;
      else
        lb = -1;
    }
    else if(type == MPL_LO) {
      type = ROWTYPE_GE;
      ub = lb;
      lb = -1;
    }
    else if(type == MPL_FX) {
      type = ROWTYPE_EQ;
      lb = -1;
    }

    /* adjust for non-zero constant term in the constraint */
    ub -= mpl_get_row_c0(mpl, i);

    /* set constraint coefficients */
    len = mpl_get_mat_row(mpl, i, ndx, val);
    lp->add_constraintex(lp, len, val+1, ndx+1, type, ub);
    ii++;
    if(lb > 0)
      lp->set_rh_range(lp, ii, lb);

    /* set constraint name */
    lp->set_row_name(lp, ii, mpl_get_row_name(mpl, i));

  }

  /* set status and free working arrays */
  status = TRUE;
  free(ndx);
  free(val);

  /* free resources used by the model translator */
done:
  lib_set_fault_hook((void *)lp, NULL);
  lib_set_print_hook((void *)lp, NULL);

  mpl_terminate(mpl);

  /* finally transpose the model A matrix */
  lp->set_add_rowmode(lp, FALSE);

  return( status );
}

/* CAN MODIFY */
MYBOOL XLI_CALLMODEL xli_writemodel(lprec *lp, char *filename, char *options, MYBOOL results)
{
  _MathProg_hook_fn(lp, "XLI_MathProg: xli_writemodel not implemented");
  return( FALSE );
}
