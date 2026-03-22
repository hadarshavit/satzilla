#ifndef __SATZILLA_GLOBAL_H__
#define __SATZILLA_GLOBAL_H__
#define MAX_FEATURES 2000

#define FEAT_OK 666
#define DONE 123456789
#define LP_ERROR 42
#define LP_UNSAT 314159
#define LP_TIMEOUT 271828
#define LP_SAT 6021023
#define TOTAL_TIMEOUT_CODE 111
#define VCG_TIMEOUT_CODE 111

#define DEFAULT_LP_TIME_LIMIT 30
#define DEFAULT_CG_TIME_LIMIT 20
#define DEFAULT_LOBJOIS_TIME_LIMIT 2
#define UBCSAT_NUM_RUNS "10000"
#define DEFAULT_UBCSAT_TIME_LIMIT 5
#define DEFAULT_CL_TIME_LIMIT 5
#define DEFAULT_DIA_TIME_LIMIT 10
#define DEFAULT_SP_TIME_LIMIT 5
#define DEFAULT_SBVA_TIME_LIMIT 30
#define DEFAULT_GROUP_TIME_LIMIT 180
#define DEFAULT_PREPROCESS_TIMEOUT 35

#define DEB 1

// == NOTE: Added after randsat paper
#define EPSILON 1e-10

#define DIAMETER_FIX( __a ) \
{\
    if( var_depth[ __a ] == 0 ) \
    {   var_depth[ __a ] = current_depth + 1; \
	*(dstackp++) = __a; } \
}

#include <cstdarg>
#include <cstdio>

inline void LogDebug(const char *fmt, ...)
{
    if (!DEB)
        return;
    va_list args;
    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
    va_end(args);
}

#include "stopwatch.h"


namespace ubcsat { int main(int, char**); }
namespace varsat{double**  main(char* input, int& spsize, int timeout);}

extern double preTime;
extern int  OrigNumVars, OrigNumClauses;
extern Stopwatch gSW;
extern int gTimeOut;
extern double myTime;
extern const char* mypath;
extern int gGroupTimeoutSeconds;
extern int gPreprocessTimeoutSeconds;

int remainingTotalTimeoutSeconds();
int resolveGroupTimeoutSeconds(int defaultSeconds);
int resolvePreprocessTimeoutSeconds();
bool totalTimeoutReached();
bool groupTimeoutReached(double startTimeSeconds, int defaultSeconds);

#endif
