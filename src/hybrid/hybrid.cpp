/**
 * The Hybrid multi-core algorithm for computing skyline queries.
 *
 * @author Sean Chester (schester)
 * @date Feb 12, 2014
 * @see S. Chester et al. (2015) "Scalable parallelization of skyline 
 * computation for multi-core processors." Proceedings of the 
 * 31st IEEE International Conference on Data Engineering (ICDE 
 * 2015). 12 pages. http://cs.au.dk/~schester/publications/chester_icde2015_mcsky.pdf
 *
 */

#include "hybrid/hybrid.h"

#include <cassert>
#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#include <parallel/algorithm>
#else
#include <algorithm>
#define omp_get_thread_num() 0
#define omp_set_num_threads( t ) 0
#endif

#include "common/pq_filter.h"
#include "util/timing.h"

/**
 * Constructs a new instance of a Hybrid skyline solver.
 *
 * @param threads the number of threads to launch.
 * @param n The number of input tuples in the dataset.
 * @param d The number of dimensions in the input dataset.
 * @param accum The blocksize, alpha, of points to process in each parallel batch.
 * @param pq_size Size of the priority queues to use in the pre-filter (i.e., the
 * maximum number of points that each thread should reserve for pre-pruning).
 * @note After instantiating, a Hybrid skyline solver still requires a call to
 * Init() to copy data locally.
 */
Hybrid::Hybrid( uint32_t threads, uint32_t n, uint32_t d,
    const uint32_t accum, const uint32_t pq_size ) :
    num_threads_( threads ), n_( n ), accum_( accum ), pq_size_( pq_size ) {

  omp_set_num_threads( threads );
  skyline_.reserve( 1024 );
  part_map_.reserve( 1024 );
  data_ = NULL;
}

/**
 * Destroys a Hybrid skyline solver and deletes the data associated with it.
 */
Hybrid::~Hybrid() {
  delete[] data_;
  part_map_.clear();
  skyline_.clear();
}

/**
 * Initializes the Hybrid skyline solver by copying the input set, running
 * the pre-filter, and partitioning the data.
 *
 * @param data An array of float arrays containing the actual data points
 * that will be memcopied into this Hybrid skyline solver.
 */
void Hybrid::Init( float** data ) {
  data_ = new EPTUPLE[n_];
  for (uint32_t i = 0; i < n_; i++) {
    data_[i].pid = i;
    data_[i].partition = 0;
    memcpy( data_[i].elems, data[i], sizeof(float) * NUM_DIMS );
  }

  /* Pre-filter */
  INI_PROFILER();
  n_ = PQFilter::Execute<EPTUPLE>( data_, n_, pq_size_, num_threads_ );
  UPD_PROFILER( "01 pq-filter" );

  partition();
#if defined(_OPENMP)
  std::__parallel::sort( data_, data_ + n_ );
#else
  std::sort( data_, data_ + n_ );
#endif
}

/**
 * Executes the Hybrid skyline solver to produce a skyline 
 * from the data that it has currently stored.
 */
vector<int> Hybrid::Execute() {

  /* Overwrite local dataset with skyline. */
  const int num_survive = skyline();
  
  /* Copy skyline into skyline_ result vector. */
  for (uint32_t i = 0; i < num_survive; ++i) {
    skyline_.push_back( data_[i].pid );
  }

  PRI_PROFILER();
  return skyline_;
}

/**
 * Phase II of the Hybrid algorithm, comparing a data point p to a block 
 * of other data points to see if they dominate p. Does not use sophisticated
 * data structures.
 *
 * @param me The index in the data array of the point p that should be tested
 * whether or not it's dominated. Also the upper bound on the index that should
 * be iterated to, on account of the Manhattan Norm sort.
 * @param start The index in the data array of the first point in the block against
 * which me should be tested.
 * @post The data point me is internally marked as a side-effect if it is 
 * determined to be dominated.
 */
void inline Hybrid::compare_to_peers( const uint32_t me, const uint32_t start ) {

  /* First, iterate points in partitions below me's, assuming 
   * distinct value condition.
   */
  uint32_t i, mylev = data_[me].getLevel();
  for (i = start; i < me; ++i) {
    if ( data_[i].isPruned() )
      continue;
    if ( data_[i].getLevel() == mylev )
      break;
    if ( !data_[me].canskip_partition( data_[i].getPartition() ) ) {
      if ( DominateLeftDVC( data_[i], data_[me] ) ) {
        data_[me].markPruned();
        return;
      }
    }
  }

  /* Skip all other partitions on the same level that
   * are not the same as me's: they clearly cannot contain
   * points that dominate me. Eventually will find my partition
   * (at position me, if not earlier).
   */
  for (; data_[i].getPartition() < data_[me].getPartition(); ++i)
    ;

  /* Finally, compare to points within same partition, 
   * up to me, since only those have a Manhattan Norm
   * <= to that of i. (equal Man Norm implies equal or 
   * incomparable points, neither of which dominate me).
   */
  for (; data_[i].score < data_[me].score; ++i) {
    if ( DominateLeftDVC( data_[i], data_[me] ) ) {
      data_[me].markPruned();
      return;
    }
  }
}

/**
 * Compares tuple t to all known skyline points, using the
 * two-level list of as-yet-created partitions.
 *
 * @pre Assumes that t comes from a partition that 
 * has not yet been added to part_maps_; therefore, 
 * distinct value can be assumed.
 */
void inline Hybrid::compare_to_skyline_points( EPTUPLE &t ) {

  /* Iterate through all partitions. */
  vector<pair<uint32_t, uint32_t> >::iterator it;
  for (it = part_map_.begin(); it != part_map_.end() - 1; ++it) {

    /* If tuple t cannot skip this partition, do work. */
    if ( !t.canskip_partition( it->first ) ) {

      /* Set boundaries [begin, end) of partition. */
      const uint32_t begin = it->second;
      const uint32_t end = (it + 1)->second;

      /* Compare to head/pivot of partition, constructing 
       * comparison bitmap. Return if it dominates t.
       */
      const uint32_t bitmap = DT_bitmap_dvc( t, data_[begin] );
      if ( bitmap == ALL_ONES && !EqualityTest( t, data_[begin] ) ) {
        t.markPruned();
        return;
      }

      /* Iterate rest of partition, looking for a 
       * point to dominate t, and aborting if found. 
       * Skips points based on mutual relationship to 
       * head/pivot of partition. Can skip if t has a clear
       * bit where point i has one set.
       */
      for (uint32_t i = begin + 1; i < end; ++i) {
        if ( !(~bitmap & data_[i].partition) || !data_[i].partition ) {
          if ( DominateLeft( data_[i], t ) ) {
            t.markPruned();
            return;
          }
        }
      }
    }
  }
}

/**
 * Updates the data structure of skyline points to reflect the newly added
 * points in the range [start, end).
 *
 * @param start The first index of newly added skyline points.
 * @param end One past the last index of newly added skyline points.
 */
void inline Hybrid::update_partition_map( const uint32_t start, const uint32_t end ) {
  /* Remove sentinel and recall id, start of last partition. */
  part_map_.pop_back();
  uint32_t last_val = part_map_.at( part_map_.size() - 1 ).first;
  uint32_t part_start = part_map_.at( part_map_.size() - 1 ).second;

  /* Iterate all new points to find partitions. */
  for (uint32_t i = start; i < end; ++i) {

    /* New partition if id doesn't match previous. */
    if ( data_[i].getPartition() != last_val ) {
      last_val = data_[i].getPartition();
      part_start = i;
      part_map_.push_back( pair<uint32_t, uint32_t>( last_val, i ) );
    }

    /* Otherwise, use the first point in partition to further partition
     * this point one level deeper. Can modify .partition member directly,
     * since these points will never again be sorted; no need to update
     * partition_level, since it will no longer be used.
     */
    else {
      const uint32_t bitcode = DT_bitmap_dvc( data_[i], data_[part_start] );
      data_[i].partition = bitcode;
    }
  }

  /* Replace sentinel at end. */
  part_map_.push_back( pair<uint32_t, uint32_t>( 0, end + 1 ) );
}

/**
 * Computes the skyline of data_ using the Hybrid algorithm.
 * Uses two (fixed) levels of partitioning on top of a parallel-
 * friendly traversal order adopted from the Q-Flow algorithm (similar 
 * to PSFS from the PSkyline paper and to GGS from the GPU skyline
 * paper).
 *
 * @note Modifies the data_ member so that the skyline tuples
 * appear at the front. May overwrite/delete other data.
 * @return The number of skyline tuples in data_.
 * @see KS Bøgh et al. "Efficient GPU-based skyline computation,"
 * Proc. 9th DaMoN workshop. 2013.
 * @see H Im et al. "Parallel skyline computation on multicore 
 * architectures." Information Systems: 36(4). 808--823. 2011.
 */
int Hybrid::skyline() {
  uint32_t i, head, start, stop; //cursors

  // D[0...(head - 1)] = skyline tuples
  // D[start...stop - 1] = current working window
  head = 0;
  start = 0;

  /* Init partition map. Consists of pairs: ( bitmap, start index in D ). */
  part_map_.push_back( pair<uint32_t, uint32_t>( data_[0].getPartition(), 0 ) ); //first part.
  part_map_.push_back( pair<uint32_t, uint32_t>( data_[0].getPartition(), 1 ) ); //sentinel

  // D[next] = tuple to be considered next
  while ( start < n_ ) {
    INI_PROFILER();
    /* Check in parallel each of the next N_ACCUM
     * points to see if any are dominated by the
     * so-far-confirmed skyline points.
     */
    stop = start + accum_;
    if ( stop > n_ )
      stop = n_;

#pragma omp parallel for schedule(dynamic, 16) default(shared) private(i)
    for (i = start; i < stop; ++i) {
      compare_to_skyline_points( data_[i] );
    } // END PARALLEL FOR
    UPD_PROFILER( "11 phaseI" );

    /* Sequentially compress these points in advance
     * of comparing amongst themselves.
     */
    sort( data_ + start, data_ + stop );
    for (i = start; i < stop && !data_[i].isPruned(); ++i)
      ;
    stop = i;
    UPD_PROFILER( "13 compress" );

    /* In parallel, confirm all new candidates against
     * each other to see if any are dominated.
     */
#pragma omp parallel for schedule(dynamic, 16) default(shared) private(i)
    for (i = start; i < stop; ++i) {
      compare_to_peers( i, start );
    } // END PARALLEL FOR
    UPD_PROFILER( "12 phaseII" );

    /* Finally, sequentially compress the confirmed
     * skyline points again.
     */
    const uint32_t head_old = head;
    sort( data_ + start, data_ + stop );
    for (i = start; i < stop && !data_[i].isPruned(); ++i, ++head) {
      data_[head] = data_[i];
    }
    /* Update partition map with new tuples and advance start pos. */
    update_partition_map( head_old, head );
    start += accum_;
    UPD_PROFILER( "13 compress" );
  }
  return head;
}

/**
 * Partitions the data by the median values on each dimension 
 * and sorts the data by their new partition.
 */
void inline Hybrid::partition() {
  /* transpose data for median calculation. */
  float *data = new float[NUM_DIMS * n_];
#pragma omp parallel num_threads(num_threads_)
  {
    const uint32_t th_id = omp_get_thread_num();
#pragma omp for nowait
    for (uint32_t i = 0; i < n_; i++) {
      for (uint32_t j = 0; j < NUM_DIMS; j++) {
        data[j * n_ + i] = data_[i].elems[j];
      }
    }
  } // END PARALLEL FOR

  /* Sort each dimension and retrieve median. Note that these
   * indices are not linked really with the indices in the data_ array
   * anymore (nor need be).
   */
  PTUPLE median;
  for (uint32_t i = 0; i < NUM_DIMS; i++) {
#if defined(_OPENMP)
    std::__parallel::sort( data + i * n_, data + (i + 1) * n_ );
#else
    std::sort( data + i * n_, data + (i + 1) * n_ );
#endif
    median.elems[i] = data[i * n_ + n_ / 2];
  }
  delete[] data;
  UPD_PROFILER( "02 select pivot" );

  /* Calc partition relative to median values. */
#pragma omp parallel for
  for (uint32_t i = 0; i < n_; i++) {
    data_[i].setPartition( DT_bitmap( data_[i], median ) );
  } // END PARALLEL FOR
  UPD_PROFILER( "03 partition" );

  return;
}
