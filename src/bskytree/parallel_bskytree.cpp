/*
 * parallel_bskytree.cpp
 *
 *  Created on: Jul 6, 2014
 *      Author: dariuss
 */

#include <bskytree/parallel_bskytree.h>

#if defined(_OPENMP)
#include <omp.h>
#include <parallel/algorithm>
#else
#include <algorithm>
#define omp_get_thread_num() 0
#define omp_set_num_threads( t ) 0
#endif

#include <cassert>

#include "common/common.h"
#include "util/timing.h"

ParallelBSkyTree::ParallelBSkyTree( const uint32_t num_threads,
    const uint32_t n, const uint32_t d, float** dataset ) :
    num_threads_( num_threads ), n_( n ), d_( d ) {

  omp_set_num_threads( num_threads_ );
  skyline_.reserve( 1024 );
  eqm_.reserve( 1024 );
}

void ParallelBSkyTree::Init( float** dataset ) {
  data_.reserve( n_ );
  for (uint32_t i = 0; i < n_; i++) {
    TUPLE t;
    t.pid = i;
    memcpy( t.elems, dataset[i], sizeof(float) * NUM_DIMS );
    data_.push_back( TUPLE_S( t, -1 ) );
  }
}

ParallelBSkyTree::~ParallelBSkyTree() {
  skyline_.clear();
  data_.clear();
  eqm_.clear();
}

vector<int> ParallelBSkyTree::Execute( void ) {

//  initProfiler();
  BSkyTreeS_ALGO();

//  printProfile();

// Add missing points from "equivalence matrix"
  skyline_.insert( skyline_.end(), eqm_.begin(), eqm_.end() );

  return skyline_;
}

void ParallelBSkyTree::BSkyTreeS_ALGO() {
//  initProfiler();
  SelectBalanced(); // pivot selection in the data_
//  updateProfiler( "pivot" );

  DoPartioning(); // mapping points to binary vectors representing sub-regions
//  updateProfiler( "partitioning" );

//  initProfiler();
  vector<TUPLE_S> &S = data_; // Alias
  volatile bool *const dead = new bool[S.size()];
  memset( const_cast<bool*>(dead), 0, sizeof(bool) * S.size() );
  int head = 1; // always points to the 1st tuple after confirmed heads
  int tail = S.size() - 1; // always points to the last (active) tuple
  while ( head < tail ) {
    int htail = head + BSKYTREE_ACCUM - 1 < tail ? head + BSKYTREE_ACCUM - 1 : tail;
    #pragma omp parallel num_threads(num_threads_)
    { // private constants:
      const uint32_t p_head = head, p_tail = tail, p_htail = htail;
      #pragma omp for nowait
      for (uint32_t th = p_head; th <= p_htail; ++th) { // th -> temporal head
        int cur = p_htail + 1;
        while ( cur <= p_tail ) {
          if ( dead[cur] ) {
            ++cur;
            continue;
          }
          // Check for partial dominance:
          if ( (S[th].partition & S[cur].partition) == S[th].partition
              || (S[th].partition & S[cur].partition) == S[cur].partition ) {

            const int dt_test = DT_dvc( S[th], S[cur] );
            if ( dt_test == DOM_LEFT ) {
              dead[cur] = true;
//            S[cur++] = S[tail--]; no compression because of multi-threading
              cur++;
            } else if ( dt_test == DOM_RIGHT ) {
              dead[cur] = true;
              S[th] = S[cur];
              cur = p_htail + 1;
//            S[cur] = S[tail--]; no compression because of multi-threading
            } else {
              cur++; // Point-level incomparability
            }
          } else {
            cur++; // Region-level incomparability
          }
        } // with cur we did one pass (till tail)
      } //  (all temporal heads were processed)
    } // END OF PARALLEL
//    updateProfiler( "parallel head processing" );

    // Single-thread execution:
    for (uint32_t th = head; th <= htail; ++th) { // th -> temporal head
      int c = th + 1;
      while ( c <= htail ) {
        if ( S[th].pid == S[c].pid ) {
          dead[htail] = true;
          S[c] = S[htail--];
        } else {
          const uint32_t dt_test = DT_dvc( S[th], S[c] );
          if ( dt_test == DOM_LEFT ) {
            dead[htail] = true;
            S[c] = S[htail--];
          } else if ( dt_test == DOM_RIGHT ) {
            S[th] = S[c];
            dead[htail] = true;
            S[c] = S[htail--];
            c = th + 1;
          } else {
            c++; // two heads are incomparable
          }
        }
      }
    }
    // heads are processed, so update it:
    head = htail + 1;
//    updateProfiler( "sequential head processing" );

    // Compress by removing dead tuples:
    int head_dead = head; // will store 1st dead tuple from head & upwards
    int tail_alive = tail; // will store 1st alive tuple from tail & downwards
    while ( head_dead < tail_alive ) {
      while ( !dead[head_dead] && head_dead < tail_alive ) {
        ++head_dead;
      }

      while ( dead[tail_alive] ) {
        --tail_alive;
      }
      if ( tail_alive > head_dead ) {
        dead[head_dead] = false;
        dead[tail_alive] = true;
        S[head_dead++] = S[tail_alive--];
      }
    }
    // Update tail after compression:
    tail = tail_alive;
    while ( dead[tail] )
      --tail;
//    updateProfiler( "sequential compression" );
  } // Skyline computed!

  for (uint32_t i = 0; i <= tail; ++i) {
    skyline_.push_back( S[i].pid );
  }
  delete [] dead;
}

/*
 * Partitions the data using the pivot point (data_[0]) by
 * assigning partition bitmap to each tuple. Also, removes
 * the points that are pruned (ALL_ONES partition).
 */
void ParallelBSkyTree::DoPartioning() {
  const uint32_t pruned = SHIFTS[NUM_DIMS] - 1;
  const TUPLE_S &pivot = data_[0];
  for (uint32_t i = 1; i < data_.size(); ++i) {
    if ( EqualityTest( pivot, data_[i] ) ) {
      eqm_.push_back( data_[i].pid );
      data_[i] = data_.back();
      data_.pop_back();
      continue;
    }
    const uint32_t lattice = DT_bitmap_dvc( data_[i], pivot );
    if ( lattice < pruned ) {
      assert( !DominateLeft( pivot, data_[i] ) );
      data_[i].partition = lattice;
    } else {
      data_[i] = data_.back();
      data_.pop_back();
    }
  }
}

/*
 * Chooses a pivot based on minimum range. The chosen pivot
 * is a skyline point.
 *
 * In addition to that, removes points from data_ that are
 * dominated by the (current) pivot point.
 *
 * The pivot point is stored in data_[0].
 */
void ParallelBSkyTree::SelectBalanced() {
  const vector<float> min_list( NUM_DIMS, 0.0 );
  const vector<float> max_list( NUM_DIMS, 1.0 );

  const uint32_t head = 0;
  uint32_t tail = data_.size() - 1, cur_pos = 1;
  float* hvalue = data_[head].elems;

  vector<float> range_list = SetRangeList( min_list, max_list );
  float min_dist = ComputeDistance( hvalue, min_list, range_list );

  while ( cur_pos <= tail ) {
    float* cvalue = data_[cur_pos].elems;

    const uint32_t dtest = DominanceTest( data_[head], data_[cur_pos] );
    if ( dtest == DOM_LEFT ) {
      data_[cur_pos] = data_[tail];
      data_.pop_back();
      tail--;
    } else if ( dtest == DOM_RIGHT ) {
      data_[head] = data_[cur_pos];
      data_[cur_pos] = data_[tail];
      data_.pop_back();
      tail--;

      hvalue = data_[head].elems;
      min_dist = ComputeDistance( hvalue, min_list, range_list );
      cur_pos = 1; // THIS IS THE SAME BUG AS IN QSkyCube: cur_pos is not reseted
    } else {
      assert( dtest == DOM_INCOMP );
      float cur_dist = ComputeDistance( cvalue, min_list, range_list );

      if ( cur_dist < min_dist ) {
        if ( EvaluatePoint( cur_pos ) ) {
          std::swap( data_[head], data_[cur_pos] );

          hvalue = data_[head].elems;
          min_dist = cur_dist;
          cur_pos++;
        } else {
          data_[cur_pos] = data_[tail];
          data_.pop_back();
          tail--;
        }
      } else
        cur_pos++;
    }
  }
}

vector<float> ParallelBSkyTree::SetRangeList( const vector<float>& min_list,
    const vector<float>& max_list ) {
  vector<float> range_list( NUM_DIMS, 0 );
  for (uint32_t d = 0; d < NUM_DIMS; d++)
    range_list[d] = max_list[d] - min_list[d];

  return range_list;
}

/*
 * Note that here we do not need to do normalization (we assume
 * the data is pre-normalized). Though, we do (TODO!)
 */
float ParallelBSkyTree::ComputeDistance( const float* value,
    const vector<float>& min_list, const vector<float>& range_list ) {
  float max_d, min_d;

  max_d = min_d = (value[0] - min_list[0]) / range_list[0];
  for (uint32_t d = 1; d < NUM_DIMS; d++) {
    float norm_value = (value[d] - min_list[d]) / range_list[d];
    if ( min_d > norm_value )
      min_d = norm_value;
    else if ( max_d < norm_value )
      max_d = norm_value;
  }

  return max_d - min_d;
}

/*
 * Checks if the point dataset[pos] is not dominated by any of points
 * before pos (dataset[0..pos-1]).
 *
 * Note that here we can remove additionally dominated points, but the
 * code does not do it (the paper suggests, though).
 */
bool ParallelBSkyTree::EvaluatePoint( const uint32_t pos ) {
  const TUPLE &cur_tuple = data_[pos];
  for (uint32_t i = 0; i < pos; ++i) {
    const TUPLE &prev_value = data_[i];
    if ( DominatedLeft( cur_tuple, prev_value ) )
      return false;
  }

  return true;
}
