#ifndef STOPWATCH_H
#define STOPWATCH_H

//#include "defines.h"

#define STOPWATCH_CPU_TIME

// includes
#ifdef STOPWATCH_CPU_TIME
	#include <time.h>
#endif
#ifndef _WIN32
	#ifdef STOPWATCH_CPU_TIME
		#include <sys/times.h>
		#include <limits.h>
		#ifndef CLK_TCK
			#define CLK_TCK 60
		#endif
	#endif
#endif

// Stopwatch class
class Stopwatch 
{
public:
	Stopwatch();
	void Start();
	double Lap();
	double ChildrenLap();
	double TotalLap() { return Lap() + ChildrenLap(); };
private:
#ifdef _WIN32
	#ifdef STOPWATCH_CPU_TIME
		clock_t t_start;				
	#else
		double mFreq;
		double mStart;
		double mStop;
		int mErr;
	#endif
#else
	#ifdef STOPWATCH_CPU_TIME
		clock_t t_start_u, t_start_s;
		long _SYS_CLK_TCK;
	#else
		int start1;
		int start2;
	#endif
#endif
};

#endif // !INCLUDE_STPWTCH_H
