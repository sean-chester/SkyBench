/*
 * pskyline.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: dariuss
 */

#include "pskyline.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <strings.h>

#if defined(_OPENMP)
#include <omp.h>
#include <parallel/algorithm>
#else
#include <algorithm>
#define omp_get_thread_num() 0
#define omp_set_num_threads( t ) 0
#endif

PSkyline::PSkyline(uint32_t threads, uint32_t n, uint32_t d, float** data) :
    num_threads_( threads ), n_( n ), d_( d ), block_size_( n / threads ) {
  skyline_.reserve( 1024 );
  omp_set_num_threads( num_threads_ );
  data_ = NULL;
  input_ = NULL;
  flag_ = NULL;
}

PSkyline::~PSkyline() {
  delete[] data_;
  delete[] input_;
  delete[] flag_;
}

vector<int> PSkyline::Execute() {
  INI_PROFILER();
  Block* output = PMap( input_ );
  UPD_PROFILER("11 phaseI");
  Block result = SReduce( output );
  UPD_PROFILER("12 phaseII");

  for (uint32_t i = 0; i < result.size; ++i) {
    skyline_.push_back( data_[i].pid );
  }

  PRI_PROFILER();
  delete[] output;
  return skyline_;
}

void PSkyline::Init(float** data) {
  data_ = new TUPLE[n_];
#pragma omp parallel for
  for (uint32_t i = 0; i < n_; i++) {
    data_[i].pid = i;
    memcpy( data_[i].elems, data[i], sizeof(float) * NUM_DIMS );
  }

  input_ = new Block[num_threads_];
  flag_ = new int[n_];
  int start = 0, end = 0;
  uint32_t i;
  for (i = 0; i < num_threads_; i++) {
    end = start + block_size_ - 1;
    input_[i].start = start;
    input_[i].end = end;
    start = end + 1;
  }
  input_[i - 1].end = n_ - 1; // .end : inclusive
}

/*
 * Implementation of PMap using OMP:
 *  PMap(f, D) = { f(D[1]), ..., f(D[n]) }
 *
 * I.e., f (sskyline) is applied to each element of D in parallel.
 *
 */
Block* PSkyline::PMap(Block* input) {
  Block* output = new Block[num_threads_];
  #pragma omp parallel for default(shared)
  for (uint32_t i = 0; i < num_threads_; i++)
    output[i] = sskyline( input[i] );
  // END PARALLEL FOR
  return output;
}

/*
 * Simple Skyline
 */
Block PSkyline::sskyline(Block input) {
  int i;
  const int size = input.end - input.start + 1;
  int head = 0, tail = size - 1;

  TUPLE *const D = data_ + input.start;

  while ( head < tail ) {
    i = head + 1;
    while ( i <= tail ) {
      const int dtest = DominanceTest( D[head], D[i] );
      if ( dtest == DOM_LEFT )
        D[i] = D[tail--];
      else if ( dtest == DOM_RIGHT ) {
        D[head] = D[i];
        D[i] = D[tail--];
        i = head + 1;
      } else
        ++i;
    }
    head++;
  }

  // keep the output in memory
  input.size = tail + 1;

  return input;
}

/*
 * SReduce (Sequential Reduce): collapses several skylines into
 * one by sequentially calling parallel merge.
 */
Block PSkyline::SReduce(Block* input) {
  if ( num_threads_ == 1 )
    return input[0];

  Block* buf = new Block[num_threads_];
  memcpy( buf, input, num_threads_ * sizeof(Block) );
  Block ret = buf[0];
  for (uint32_t i = 1; i <= num_threads_ - 1; i++)
    ret = PMerge( ret, buf[i] );

  delete[] buf;
  return ret;
}

/*
 * Parallel Merge
 */
Block PSkyline::PMerge(Block left, Block right) {
  int* flag;
  int* left_flag;
  int* right_flag;

  TUPLE* left_skyline;
  TUPLE* right_skyline;

  // assumes that left[] appears to the left of right[]
  memmove( data_ + left.start + left.size, data_ + right.start,
      sizeof(TUPLE) * right.size );

  left_skyline = data_ + left.start;
  right_skyline = data_ + left.start + left.size;

  /* Set up flag */
  flag = flag_ + left.start;
  left_flag = flag;
  right_flag = flag + left.size;

  // assumption: LIVE == 0
  //bzero( flag, (left.size + right.size) * sizeof(int) );
  memset( flag, 0, (left.size + right.size) * sizeof(int) );

#pragma omp parallel for schedule(dynamic, 64)
  for (uint32_t i = 0; i < left.size; i++) {
    left_flag[i] = CheckSurvival( left_skyline[i], right_skyline, right_flag,
        right.size );
  } // END PARALLEL FOR

  /* Compact skylines */
  uint32_t cnt = 0;
  for (uint32_t i = 0; i < left.size + right.size; i++)
    if ( flag[i] == LIVE )
      left_skyline[cnt++] = left_skyline[i];

  left.size = cnt;

  return left;
}

