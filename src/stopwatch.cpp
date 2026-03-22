// Windows specific implementation
#include "stopwatch.h"

#ifdef _WIN32

#ifndef STOPWATCH_CPU_TIME
	#include <windows.h>
#endif

Stopwatch::Stopwatch() 
{
	#ifndef STOPWATCH_CPU_TIME
		LARGE_INTEGER li;
		mErr = !QueryPerformanceFrequency(&li);
		mFreq = (double)li.QuadPart;
	#endif
}

void Stopwatch::Start() 
{
	#ifdef STOPWATCH_CPU_TIME
		t_start = clock();
	#else
		LARGE_INTEGER li;
		mErr = QueryPerformanceCounter(&li);
		mStart = (double)li.QuadPart/mFreq;
	#endif
}

double Stopwatch::Lap() 
{
	#ifdef STOPWATCH_CPU_TIME
		return (double)(clock() - t_start) / CLOCKS_PER_SEC;
	#else
		LARGE_INTEGER li;
		mErr = QueryPerformanceCounter(&li);
		return (double)li.QuadPart/mFreq - mStart;
	#endif
}


#else
/* UNIX version */
#ifdef STOPWATCH_CPU_TIME
	#include <math.h>
	#include <ctype.h>
	#include <time.h>
	#include <sys/times.h>
	#include <limits.h> 
        #include <unistd.h>

	#include <string>
	#include <sstream>
	#include <fstream>
using namespace std;

#else
	#include <sys/time.h>
#endif

Stopwatch::Stopwatch() 
{
#ifdef STOPWATCH_CPU_TIME
  _SYS_CLK_TCK = sysconf(_SC_CLK_TCK);
#endif
}

void Stopwatch::Start() 
{
	#ifdef STOPWATCH_CPU_TIME
		struct tms buf;
		clock_t t;
		t = times(&buf);
		t_start_u = buf.tms_utime;
		t_start_s = buf.tms_stime;
	#else
		struct timeval tp;
		gettimeofday(&tp,(void *)NULL);
		start1 = tp.tv_sec;
		start2 = tp.tv_usec;
	#endif
}

double Stopwatch::Lap() 
{
	#ifdef STOPWATCH_CPU_TIME
		struct tms buf;
		clock_t t;
		t = times(&buf);
		return(1.0*((buf.tms_utime+buf.tms_stime)-(t_start_u+t_start_s))/_SYS_CLK_TCK);
		//		return(1.0*((buf.tms_utime)-(t_start_u))/_SYS_CLK_TCK);
	#else
		struct timeval tp;
		gettimeofday(&tp,(void *)NULL);
		return (double)(tp.tv_sec - start1) + 1e-6 * (double)(tp.tv_usec - start2);
	#endif
}

double Stopwatch::ChildrenLap()
{
#ifdef STOPWATCH_CPU_TIME
  struct tms buf;
  clock_t t = times(&buf);
  return double(buf.tms_cutime + buf.tms_cstime) / double(_SYS_CLK_TCK);
  
#else
  return -1.0; // -- not implemented
#endif
}

#endif
