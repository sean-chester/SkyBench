/*
 * qflow.cpp
 *
 *  Created on: Feb 12, 2014
 *      Author: schester
 */

#include "qflow/qflow.h"

#include <cstdio>
#include <cassert>

#if defined(_OPENMP)
#include <omp.h>
#include <parallel/algorithm>
#else
#include <algorithm>
#define omp_get_thread_num() 0
#define omp_set_num_threads( t ) 0
#endif

QFlow::QFlow( uint32_t threads, uint32_t n, uint32_t d, float** data,
    uint32_t accum ) :
    num_threads_( threads ), n_( n ), accum_(accum) {

  omp_set_num_threads( threads );
  skyline_.reserve( 1024 );
  data_ = NULL;
}

QFlow::~QFlow() {
  delete[] data_;
}

void QFlow::Init( float** data ) {
  data_ = new STUPLE[n_];
  for (uint32_t i = 0; i < n_; i++) {
    data_[i].pid = i;
    for (uint32_t j = 0; j < NUM_DIMS; j++) {
      data_[i].elems[j] = data[i][j];
    }
  }
}

bool QFlow::STupleComp( STUPLE const &lhs, STUPLE const &rhs ) {
  if ( lhs.score < rhs.score )
    return true;
  if ( rhs.score < lhs.score )
    return false;
  // if a tie, sort by manhattan norm:
  float lhs_sum = lhs.elems[0];
  float rhs_sum = rhs.elems[0];
  for (uint32_t d = 1; d < NUM_DIMS; ++d) {
    lhs_sum += lhs.elems[d];
    rhs_sum += rhs.elems[d];
  }
  return lhs_sum < rhs_sum;
}

vector<int> QFlow::Execute() {
  INI_PROFILER();
  // sort:
  ComputeScores();
#if defined(_OPENMP)
  std::__parallel::sort( data_, data_ + n_, STupleComp );
#else
  std::sort( data_, data_ + n_, STupleComp );
#endif
  UPD_PROFILER("01 pq-filter");

  const int num_survive = skyline();
  for (uint32_t i = 0; i < num_survive; ++i) {
    skyline_.push_back( data_[i].pid );
  }
  PRI_PROFILER();

  return skyline_;
}

// return = number of surviving tuples
int QFlow::skyline() {
  int i, j;
  int head1, head2, start, stop;
  float stop_val, candidate_stop_val;
  bool* sky = new bool[n_]();

  // D[0...(head1 - 1)] = skyline tuples
  // D[head1...(head2 - 1)] = candidate tuples
  head1 = 0;
  head2 = 1;
  start = 1;
  sky[0] = true;

  // D[next] = tuple to be considered next
  while ( start < n_ ) {

    /* Check in parallel each of the next N_ACCUM
     * points to see if any are dominated by the
     * so-far-confirmed skyline points.
     */
    stop = start + accum_;
    if ( stop > n_ )
      stop = n_;
#pragma omp parallel for default(shared) private(i,j)
    for (i = start; i < stop; i++) {
      for (j = 0; j <= head1; j++) {
        if ( DominateLeft( data_[j], data_[i] ) ) {
          sky[i] = false;
          break;
        }
      }
      if ( j == head1 + 1 )
        sky[i] = true; /* Candidate to continue on. */
    }
    UPD_PROFILER("11 phaseI");

    /* Sequentially compress these points in advance
     * of comparing amongst themselves.
     */
    for (i = start, head2 = head1; i < stop; i++) {
      if ( sky[i] )
        data_[++head2] = data_[i];
    }
    UPD_PROFILER( "13 compress" );

    /* In parallel, confirm all new candidates against
     * each other to see if any are dominated.
     */
#pragma omp parallel for default(shared) private(i,j)
    for (i = head1 + 1; i <= head2; i++) {
      const uint32_t end = i;
      for (j = head1 + 1; j < end; j++) {
        if ( DominateLeft( data_[j], data_[i] ) ) {
          sky[i] = false;
          break;
        }
      }
      if ( j == end )
        sky[i] = true; /* Legitimately confirmed as skyline. */
    }
    UPD_PROFILER( "12 phaseII" );

    /* Finally, sequentially compress the confirmed
     * skyline points again and check if any update
     * the SaLSa stop condition.
     */
    for (i = head1 + 1; i <= head2; i++) {
      if ( sky[i] ) {
        data_[++head1] = data_[i];
      }
    }
    UPD_PROFILER( "13 compress" );
    start = stop;
  }

  delete[] sky;
  return head1 + 1;
}

void QFlow::ComputeScores() {
#pragma omp parallel for
  for (uint32_t i = 0; i < n_; i++) {
    data_[i].score = data_[i].elems[0];
    for (uint32_t j = 1; j < NUM_DIMS; j++) {
      data_[i].score += data_[i].elems[j];
    }
  } // END PARALLEL FOR
}
