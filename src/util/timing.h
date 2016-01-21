/*
 * timing.h
 *
 *  Created on: Feb 4, 2014
 *      Author: dariuss
 */

#ifndef TIMING_H_
#define TIMING_H_

/*
 * Returns current milliseconds
 */
long GetTime();

void PrintTime(long milli_sec, const char *const prefix_header = "");

#endif /* TIMING_H_ */
