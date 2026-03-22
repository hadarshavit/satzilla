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

// #define DIA_TIME_LIMIT 10
// #define SP_TIME_LIMIT 5

// == NOTE: Change these to control timout parameters

#ifdef FEATURES

#define LP_TIME_LIMIT 30
#define CG_TIME_LIMIT 20
#define LOBJOIS_TIME_LIMIT  2
#define UBCSAT_TIME_LIMIT "2"
#define UBCSAT_TIME_LIMIT_INT 5
#define UBCSAT_NUM_RUNS "10000"
#define CL_TIME_LIMIT 5
#define DIA_TIME_LIMIT 10
#define SP_TIME_LIMIT 5
#define SBVA_TIMEOUT "30"
#define TOTAL_TIMEOUT 180
#define VCG_TIMEOUT_CODE 111

#else


#define LP_TIME_LIMIT 30
#define CG_TIME_LIMIT 20
#define LOBJOIS_TIME_LIMIT  2
#define UBCSAT_TIME_LIMIT "2"
#define UBCSAT_TIME_LIMIT_INT 5
#define UBCSAT_NUM_RUNS "10000"
#define CL_TIME_LIMIT 5
#define DIA_TIME_LIMIT 10
#define SP_TIME_LIMIT 5
#define SBVA_TIMEOUT "30"
#define TOTAL_TIMEOUT 180

#endif

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

#endif
