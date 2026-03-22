#include <ctime>
#include <cstring>
#include <stdint.h>
#include <errno.h>
#include <signal.h>
#include <zlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


#include "Solver.h"

double** main(char* input, int& spsize);
inline double cpuTime(void);
inline double cpuTime(void);
static inline int memReadStat(int field);
static inline uint64_t memUsed() ;
static inline uint64_t memUsed(void);
static inline uint64_t memUsed() ;
void printStats(Solver& solver);
static void SIGINT_handler(int signum);
void printUsage(char** argv);
const char* hasPrefix(const char* str, const char* prefix);

