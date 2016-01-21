/*
 * skyline_i.cpp
 *
 *  Created on: Mar 7, 2014
 *      Author: dariuss
 */

#include <cstdio>
#include "common/skyline_i.h"

/* Profiling stuff (for breakdown charts). */
#if PROFILER == 1

void SkylineI::initProfiler() {
  struct timeval tv;
  gettimeofday( &tv, NULL );
  prev_time_.tv_sec = tv.tv_sec;
  prev_time_.tv_usec = tv.tv_usec;

  // Initialize predefined sub-run-times:
  if ( profiler_.empty() ) {
    profiler_["01 pq-filter"] = 0;
    profiler_["02 select pivot"] = 0;
    profiler_["03 partition"] = 0;
    profiler_["11 phaseI"] = 0;
    profiler_["12 phaseII"] = 0;
    profiler_["13 compress"] = 0;
  }
}

void SkylineI::updateProfiler( std::string map_key ) {
  struct timeval tv;
  gettimeofday( &tv, NULL );
  double this_step_t = (tv.tv_sec - prev_time_.tv_sec) * 1000000 + tv.tv_usec
      - prev_time_.tv_usec;

  if ( profiler_.find( map_key ) == profiler_.end() )
    profiler_[map_key] = this_step_t;
  else
    profiler_[map_key] += this_step_t;

  prev_time_.tv_sec = tv.tv_sec;
  prev_time_.tv_usec = tv.tv_usec;
}

void SkylineI::printProfile() {
#if NVERBOSE
  // print only pre-defined values
//    printf( " %.0lf %.0lf %.0lf", profiler_["11 phaseI"] / 1000,
//        profiler_["12 phaseII"] / 1000, profiler_["13 sort/shift"] / 1000 );
  for (std::map<std::string, double>::iterator it = profiler_.begin();
      it != profiler_.end(); ++it) {
    printf( " %.0lf", it->second / 1000 );
  }
#else
  double sum = 0;
  printf( " profiling results:\n" );
  for (std::map<std::string, double>::iterator it = profiler_.begin();
      it != profiler_.end(); ++it) {
    printf( "  %s:\t%.0lfms\n", it->first.c_str(), it->second / 1000 );
    sum += it->second;
  }
  printf( "  total_t:\t%.0lfms\n", sum / 1000 );
#endif
}

#endif
