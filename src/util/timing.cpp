/*
 * timing.cpp
 *
 *  Created on: Feb 4, 2014
 *      Author: dariuss
 */

#include "util/timing.h"

#include <sys/time.h>
#include <cstdlib>
#include <cstdio>

long GetTime() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec * 1000 + t.tv_usec / 1000;
}

void PrintTime(long milli_sec, const char *const prefix_header) {
	long v = milli_sec;
	long hours = v / (1000 * 60 * 60);
	v %= (1000 * 60 * 60);
	long minutes = v / (1000 * 60);
	v %= (1000 * 60);
	long seconds = v / 1000;
	v %= 1000;
	long milli_seconds = v;
	long first = 1;
	printf("%s[", prefix_header);
	if (hours) {
		if (!first)
			printf(":");
		printf("%ldh", hours);
		first = 0;
	}
	if (minutes) {
		if (!first)
			printf(":");
		printf("%ldm", minutes);
		first = 0;
	}
	if (seconds) {
		if (!first)
			printf(":");
		printf("%lds", seconds);
		first = 0;
	}
	if (milli_seconds) {
		if (!first)
			printf(":");
		printf("%ldms", milli_seconds);
		first = 0;
	}
	printf("]\n");
}
