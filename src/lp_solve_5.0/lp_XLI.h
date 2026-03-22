
/* ---------------------------------------------------------------------------------- */
/* lp_solve v5+ headers for external language interface libraries (XLIs)              */
/* ---------------------------------------------------------------------------------- */

#ifndef XLI_CALLMODEL
  #ifdef WIN32
    #define XLI_CALLMODEL __stdcall   /* "Standard" call model */
  #else
    #define XLI_CALLMODEL
  #endif
#endif

#ifdef RoleIsExternalLanguageEngine
  #define __XLI_EXPORT_TYPE __EXPORT_TYPE
#else
  #define __XLI_EXPORT_TYPE
#endif
 

/* Routines with UNIQUE implementations for each XLI engine (sample in lp_XLI2.c)     */
/* ---------------------------------------------------------------------------------- */
char   __XLI_EXPORT_TYPE *(XLI_CALLMODEL xli_name)(void);
MYBOOL __XLI_EXPORT_TYPE (XLI_CALLMODEL xli_readmodel)(lprec *lp, char *modelname, char *dataname, char *options, int verbose);
MYBOOL __XLI_EXPORT_TYPE (XLI_CALLMODEL xli_writemodel)(lprec *lp, char *filename, char *options, MYBOOL results);


/* Routines SHARED for all XLI implementations; located in lp_XLI1.c                  */
/* ---------------------------------------------------------------------------------- */
MYBOOL __XLI_EXPORT_TYPE (XLI_CALLMODEL xli_compatible)(lprec *lp, int xliversion, int lpversion);


