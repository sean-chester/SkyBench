/*
 * pq_filter.h
 *
 *  Created on: Jul 16, 2014
 *      Author: dariuss
 *
 *
 *  Skyline filter based on priority queue
 */

#ifndef PQ_FILTER_H_
#define PQ_FILTER_H_

#include <queue>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#include <parallel/algorithm>
#else
#include <algorithm>
#define omp_get_thread_num() 0
#define omp_set_num_threads( t ) 0
#endif

#include <common/pq_filter.h>

#include "common/common.h"

using namespace std;

typedef std::pair<uint32_t, float> mn_w_idx;

struct PQComparator {
  bool operator()( const mn_w_idx &a, const mn_w_idx &b ) {
    return a.second < b.second;
  }
};

typedef priority_queue<mn_w_idx, vector<mn_w_idx>, PQComparator> PQ;

class PQFilter {
public:

  /*
   * Executes priority queue based filtering on data using num_threads
   * queues each of pq_size.
   *
   * Side affect: simultaneously computes Manhattan norm in TUPLE.score.
   */
  template<typename T>
  static uint32_t Execute( T* data, const uint32_t n,
      const uint32_t pq_size, const uint32_t num_threads );

};

// Templated static function has to be defined in a header file..

template<typename T>
uint32_t PQFilter::Execute( T* data, const uint32_t n, const uint32_t pq_size,
    const uint32_t num_threads ) {
  PQ * const PQs_ = new PQ[num_threads];

  /* Init all threads to first q_size points and score them. */
  for (uint32_t i = 0; i < pq_size; ++i) {
    data[i].score = 0;
    for (uint32_t j = 0; j < NUM_DIMS; ++j) {
      data[i].score += data[i].elems[j];
    }
    for (uint32_t j = 0; j < num_threads; ++j) {
      PQs_[j].push( mn_w_idx( i, data[i].score ) );
    }
  }

  /* Computing top man norm scores and remember best q_size ones. */
#pragma omp parallel num_threads(num_threads)
  {
    const uint32_t th_id = omp_get_thread_num();
    mn_w_idx worst_of_bests = PQs_[th_id].top();
#pragma omp for nowait
    for (uint32_t i = 0; i < n; ++i) {
      float sum = 0;
      for (uint32_t j = 0; j < NUM_DIMS; j++) {
        sum += data[i].elems[j];
      }
      data[i].score = sum;

      /* Compare to best found man norms for this thread. */
      if ( worst_of_bests.second > sum ) {
        PQs_[th_id].pop();
        PQs_[th_id].push( mn_w_idx( i, sum ) );
        worst_of_bests = PQs_[th_id].top();
      }
    }
  } // END PARALLEL FOR

  /* Take top pruners and merge them into one set. */
  vector<uint32_t> pruners;
  pruners.reserve( num_threads * pq_size );
  for (uint32_t i = 0; i < num_threads; ++i) {
    while ( !PQs_[i].empty() ) {
      mn_w_idx top = PQs_[i].top();
      pruners.push_back( top.first );
      PQs_[i].pop();
    }
  }

  //  UPD_PROFILER( "01 calc mns" );

  /* Pre-filter dataset using top pruners. */
#pragma omp parallel for
  for (uint32_t i = 0; i < n; ++i) {
    for (vector<uint32_t>::iterator it = pruners.begin(); it != pruners.end();
        ++it) {
      if ( DominateLeft( data[*it], data[i] ) ) {
        data[i].markPruned();
        break;
      }
    }
  } // END PARALLEL FOR

  /* Determine how many points were pruned. */
  uint32_t new_n = n;
  for (uint32_t i = 0; i < new_n; ++i) {
    if ( data[i].isPruned() ) {
      data[i--] = data[--new_n];
    }
  }

//#ifndef NVERBOSE
//  printf( " pq_filter: %0.2f %% pruned\n", (n - new_n) / (double) n * 100.0 );
//#endif

  return new_n;
}

#endif /* PQ_FILTER_H_ */
