/******************************************************************************************[Main.C]
VarSat -- Copyright (c) 2008, Eric Hsu.
Built upon MiniSat code by Niklas Een and Niklas Sorensson.
Their original copyright notice is repeated below.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <ctime>
#include <cstring>
#include <stdint.h>
#include <errno.h>
#include <signal.h>
#include <zlib.h>

#include "Solver.h"

namespace varsat {

/*************************************************************************************/
#ifdef _MSC_VER
#include <ctime>

//static inline double cpuTime(void) {          EIH 22 FEB 08 removed static
inline double cpuTime(void) {
    return (double)clock() / CLOCKS_PER_SEC; }
#else

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

//static inline double cpuTime(void) {          EIH 22 FEB 08 removed static
inline double cpuTime(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; }
#endif


#if defined(__linux__)
static inline int memReadStat(int field)
{
    char    name[256];
    pid_t pid = getpid();
    sprintf(name, "/proc/%d/statm", pid);
    FILE*   in = fopen(name, "rb");
    if (in == NULL) return 0;
    int     value;
    for (; field >= 0; field--)
        fscanf(in, "%d", &value);
    fclose(in);
    return value;
}
static inline uint64_t memUsed() { return (uint64_t)memReadStat(0) * (uint64_t)getpagesize(); }


#elif defined(__FreeBSD__)
static inline uint64_t memUsed(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss*1024; }


#else
static inline uint64_t memUsed() { return 0; }
#endif

#if defined(__linux__)
#include <fpu_control.h>
#endif





//=================================================================================================
// DIMACS Parser:

#define CHUNK_LIMIT 1048576

class StreamBuffer {
    gzFile  in;
    char    buf[CHUNK_LIMIT];
    int     pos;
    int     size;

    void assureLookahead() {
        if (pos >= size) {
            pos  = 0;
            size = gzread(in, buf, sizeof(buf)); } }

public:
    StreamBuffer(gzFile i) : in(i), pos(0), size(0) {
        assureLookahead(); }

    int  operator *  () { return (pos >= size) ? EOF : buf[pos]; }
    void operator ++ () { pos++; assureLookahead(); }
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<class B>
static void skipWhitespace(B& in) {
    while ((*in >= 9 && *in <= 13) || *in == 32)
        ++in; }

template<class B>
static void skipLine(B& in) {
    for (;;){
        if (*in == EOF || *in == '\0') return;
        if (*in == '\n') { ++in; return; }
        ++in; } }

template<class B>
static int parseInt(B& in) {
    int     val = 0;
    bool    neg = false;
    skipWhitespace(in);
    if      (*in == '-') neg = true, ++in;
    else if (*in == '+') ++in;
    if (*in < '0' || *in > '9') reportf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
    while (*in >= '0' && *in <= '9')
        val = val*10 + (*in - '0'),
        ++in;
    return neg ? -val : val; }

template<class B>
static void readClause(B& in, Solver& S, vec<Lit>& lits) {
    int     parsed_lit, var;
    lits.clear();
    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        while (var >= S.nVars()) S.newVar();
        lits.push( (parsed_lit > 0) ? Lit(var) : ~Lit(var) );
    }
}

template<class B>
static bool match(B& in, char* str) {
    for (; *str != 0; ++str, ++in)
        if (*str != *in)
            return false;
    return true;
}


template<class B>
static void parse_DIMACS_main(B& in, Solver& S) {
    vec<Lit> lits;
    for (;;){
        skipWhitespace(in);
        if (*in == EOF)
            break;
        else if (*in == 'p'){
            if (match(in, "p cnf")){
                int vars    = parseInt(in);
                int clauses = parseInt(in);
				if (!S.cjm_mode) {
					reportf("|  Number of variables:  %-12d                                         |\n", vars);
					reportf("|  Number of clauses:    %-12d                                         |\n", clauses);
				}

				// Initialize data structures for computing surveys.  EIH 12 JAN 08
				S.pos_biases.growTo(vars);
				S.neg_biases.growTo(vars);
				S.star_biases.growTo(vars);
				S.neg_clauses.growTo(vars);
				S.pos_clauses.growTo(vars);

				// Initialize flag for whether to turn off surveys.  EIH 13 FEB 08
				S.heuristic_default_override = false;

	        }else{
                reportf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }
        } else if (*in == 'c' || *in == 'p')
            skipLine(in);
        else
            readClause(in, S, lits),
            S.addClause(lits);
    }
}

// Inserts problem into solver.
//
static void parse_DIMACS(gzFile input_stream, Solver& S) {
    StreamBuffer in(input_stream);
    parse_DIMACS_main(in, S); }


//=================================================================================================


void printStats(Solver& solver)
{
    //double   cpu_time = cpuTime();   EIH 22 FEB 08  divided up the time statistics
	double	total_time = cpuTime();
	double  surv_time = solver.surveys_time;
	double  solv_time = total_time - surv_time;
	double	cpu_time = solv_time;     // (alias)

    uint64_t mem_used = memUsed();
	reportf("surveys--CPU time (sec)   : %-12.4f   (%4.2f %% of total)\n", surv_time, 100.0 * surv_time / total_time);                                                              // EIH 22 FEB 08
	reportf("surveys--attempted        : %-12lld   (%.0f /sec)\n", solver.surveys_attempted, solver.surveys_attempted / surv_time);                                         // EIH 22 FEB 08
	reportf("surveys--converged        : %-12lld   (%4.2f %% successful)\n", solver.surveys_converged, 100.0 * (float)solver.surveys_converged / (float)solver.surveys_attempted);  // EIH 22 FEB 08
	reportf("restarts                  : %lld\n", solver.starts);
    reportf("conflicts                 : %-12lld   (%.0f /sec)\n", solver.conflicts   , solver.conflicts   /cpu_time);
    reportf("decisions                 : %-12lld   (%4.2f %% random) (%.0f /sec)\n", solver.decisions, (float)solver.rnd_decisions*100 / (float)solver.decisions, solver.decisions   /cpu_time);
    reportf("propagations              : %-12lld   (%.0f /sec)\n", solver.propagations, solver.propagations/cpu_time);
    reportf("conflict literals         : %-12lld   (%4.2f %% deleted)\n", solver.tot_literals, (solver.max_literals - solver.tot_literals)*100 / (double)solver.max_literals);
    if (mem_used != 0) reportf("Memory used               : %.2f MB\n", mem_used / 1048576.0);
    reportf("Total CPU time (sec)      : %-12.4f\n", total_time);  // EIH 22 FEB 08 changed to total_time, adjusted formatting
}

Solver* solver;
static void SIGINT_handler(int signum) {
    reportf("\n"); reportf("*** INTERRUPTED ***\n");
    printStats(*solver);
    reportf("\n"); reportf("*** INTERRUPTED ***\n");
    exit(1); }


//=================================================================================================
// Main:

void printUsage(char** argv)
{
    reportf("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n\n", argv[0]);
    reportf("OPTIONS:\n\n");
    reportf("  -heuristic     = {BP,EMBPL,EMBPG,EMBPGV2,SP,EMSPL,EMSPG,EMSPGV2,LC,CC,default}\n"); // EIH 10 APR 08
	reportf("  -threshold     = <threshold for survey deactivation, default = 0.6>\n");  // EIH 12 JAN 08
	reportf("  -branching     = {solution,conflict,mixed,streamline,default}\n");        // EIH 29 MAR 08
	reportf("  -seed          = <random seed, default = 1>\n");                          // EIH 22 FEB 08
	reportf("  -timeout       = <seconds before giving up, default = 43200>\n");         // EIH 10 APR 08
	reportf("  -polarity-mode = {true,false,rnd}\n");
    reportf("  -decay         = <num> [ 0 - 1 ]\n");
    reportf("  -rnd-freq      = <num> [ 0 - 1 ]\n");
    reportf("  -verbosity     = {0,1,2,3,4}\n");
	reportf("  -cjm (to get specially formatted one-survey output)\n");
    reportf("\n");
}


const char* hasPrefix(const char* str, const char* prefix)
{
    int len = strlen(prefix);
    if (strncmp(str, prefix, len) == 0)
        return str + len;
    else
        return NULL;
}


double** main(char* input, int & spsize, int timeout)
{
    Solver      S;
    S.verbosity = 1;

    int         i, j;
    const char* value;
  
    S.heuristic_mode = Solver::heuristic_EMSPG;
    S.cjm_mode = true;
    S.timeout = timeout;

	if (!S.cjm_mode) reportf("This is VarSAT 1.0 beta\n");
#if defined(__linux__)
	    fpu_control_t oldcw, newcw;
	    _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
	    if (!S.cjm_mode) reportf("WARNING: for repeatability, setting FPU to use double precision\n");
#endif
	

    double cpu_time = cpuTime();

    solver = &S;
    signal(SIGINT,SIGINT_handler);
    signal(SIGHUP,SIGINT_handler);

    gzFile in = gzopen(input, "rb");
    if (in == NULL)
        reportf("ERROR! Could not open file: %s\n", input), exit(1);

	// Bypass everything here if we are in cjm mode.
	// EIH 17 MAR 08
	//
	if (S.cjm_mode) {
		parse_DIMACS(in, S);
		gzclose(in);

		// Hey Christian I think you can uncommment this if you want to do unit propagation within minisat
		// instead of in your own code...
	    //if (!S.simplify()){
		//	reportf("Solved by unit propagation\n");
		//	if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
		//	printf("UNSATISFIABLE\n");
		//	exit(20);
		//}
		
		S.doCJM();
		double ** myresult;
		myresult = new double* [S.nVars()];           // Allocate the rows.
        for(i = 0; i < S.nVars(); i++)
            myresult[i] = new double [3];
            
		for (int v = 0; v < S.nVars(); v++) {
		myresult[v][0]=S.pos_biases[v];
		myresult[v][1]=S.neg_biases[v];
		myresult[v][2]=S.star_biases[v];
    }
        spsize=S.nVars();
  //      reportf("finish\n");
		return myresult;	
	}
}
}
