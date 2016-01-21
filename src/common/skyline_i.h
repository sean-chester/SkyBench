/*
 * skyline_i.h
 *
 *  Created on: Feb 18, 2014
 *      Author: dariuss
 *
 *  Interface for skyline algorithms.
 */

#ifndef SKYLINE_I_H_
#define SKYLINE_I_H_

#include <stdint.h>

#include <vector>
#include <string>
#include <map>
#include <sys/time.h>

// Use these MACROS to gather run-times at different
// algorithm stages (instead of function calls as MACROS
// are easy disabled with DPROFILER=0 during compilation).
#if PROFILER == 1
  #define INI_PROFILER() ( initProfiler() )
  #define UPD_PROFILER(STR_KEY) ( updateProfiler(STR_KEY) )
  #define PRI_PROFILER() ( printProfile() )
#else
  #define INI_PROFILER() ((void)0)
  #define UPD_PROFILER(STR) ((void)0)
  #define PRI_PROFILER() ((void)0)
#endif

class SkylineI {
public:
  SkylineI() { }
  virtual ~SkylineI() { }

  /* Pure virtual methods */
  virtual void Init(float** data) = 0;
  virtual std::vector<int> Execute() = 0;

  /* Profiling stuff (for breakdown charts). */
#if PROFILER == 1
  std::map<std::string, double> profiler_;
  struct timeval prev_time_;

  // Use the above MACROS instead of these function calls:
  void initProfiler();
  void printProfile();
  void updateProfiler(std::string map_key);
#endif
};

#endif /* SKYLINE_I_H_ */
