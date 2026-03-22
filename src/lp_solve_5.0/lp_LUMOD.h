
#include "sparselib.h"

/* LU matrix defines */
#define LU_START_SIZE        10000  /* Start size of LU; realloc'ed if needed */
#define LU_DELTAROWS             1  /* Additional rows inserted at the top    */
#define DEF_MAXPIVOT            80  /* Maximum number of pivots before reinversion;
                                       40 appears optimal for a model like 25fv47.mps without MDO
                                       50 often gives better results with MDO */
#define MINTIMEPIVOT       5.0e-02  /* Minimum time per pivot for reinversion optimization
                                       purposes; use active monitoring only if a pivot
                                       takes more than MINTIMEPIVOT seconds.  5.0e-2 is
                                       roughly suitable for a 1GHz system.  */

/* typedef */ struct _INVrec
{
  int          status;          /* Last operation status code */
  int          dimcount;        /* The actual number of LU rows/columns */
  int          dimalloc;        /* The allocated LU rows/columns size */
  int          user_colcount;   /* The number of user LU columns */
  sparseMatrix *L, *U;
  double       *x, *y, *w;
  int          *colindex;       /* Structure containing data matrix column indexes of LU */
  int          colchange;       /* The column to be changed at the next update using data in value[.]*/
  REAL         *value;
  REAL         *pcol;           /* Reference to the elimination vector */
  REAL         *theta;          /* Structure containing the column theta values */

  int          max_Bsize;       /* The largest B matrix of user variables */
  int          max_colcount;    /* The maximum number of user columns in LU */
  int          max_LUsize;      /* The largest NZ-count of LU-files generated */
  int          num_refact;      /* Number of times the basis was refactored */
  int          num_pivots;      /* Number of pivots since last refactorization */
  int          extraP;          /* The primal RHS offset for the current inverse */
  REAL         extraD;          /* The dual objective function offset for the current inverse */
  MYBOOL       is_dirty;        /* Specifies if a column is incompletely processed */
  MYBOOL       set_Bidentity;   /* Force B to be the identity matrix at the next refactorization */
} /* INVrec */;


#ifdef __cplusplus
namespace LUMOD
extern "C" {
#endif

/* Put function headers here */

#ifdef __cplusplus
 }
#endif


