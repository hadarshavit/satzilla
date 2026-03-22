
/* Routines located in lp_XLI2.c; sample implementations                              */
/* Cfr. lp_XLI.h for definitions                                                      */
/* ---------------------------------------------------------------------------------- */


/* CAN MODIFY */
char * XLI_CALLMODEL xli_name(void)
{
  return( "XLI v1.0" );
}

/* CAN MODIFY */
MYBOOL XLI_CALLMODEL xli_readmodel(lprec *lp, char *filename, int verbose)
{
  return( NULL );
}

/* CAN MODIFY */
MYBOOL XLI_CALLMODEL xli_writemodel(lprec *lp, char *filename, MYBOOL results)
{
  return( TRUE );
}

