
/* Routines located in lp_XLI1.c; common for all language engines                     */
/* Cfr. lp_XLI.h for definitions                                                      */
/* ---------------------------------------------------------------------------------- */
/* Changes:                                                                           */
/* 13 June 2004      Initial version - KE                                             */
/* ---------------------------------------------------------------------------------- */

/* DON'T MODIFY */
MYBOOL XLI_CALLMODEL xli_compatible(lprec *lp, int xliversion, int lpversion)
{
  MYBOOL status = FALSE;
  
  if((lp != NULL) && (xliversion == XLIVERSION)) {
#if 0  
    if(lpversion == MAJORVERSION)  /* Forces XLI renewal at lp_solve major version changes */
#endif
      status = TRUE;
  }
  return( status );
}

