/**
 * hybrid.h
 *
 * @date Feb 12, 2014
 * @author Sean Chester (schester)
 */

#ifndef HYBRID_H_
#define HYBRID_H_

#include <cstdio>

#include "common/common.h"
#include "common/skyline_i.h"

using namespace std;

class Hybrid: public SkylineI {
public:
  Hybrid(uint32_t threads, uint32_t tuples, uint32_t dims,
      const uint32_t accum, const uint32_t q_size );
  virtual ~Hybrid();

  vector<int> Execute();
  void Init(float** data);

  void printPartitionSizes() {
    printf( "Created %lu non-empty partitions:\n", part_map_.size() );
    for (uint32_t i = 1; i < part_map_.size(); i++) {
      printf( "%u\n", part_map_.at( i ).second - part_map_.at( i - 1 ).second );
    }
  }

private:
  int skyline();
  void inline partition();
  void inline compare_to_skyline_points( EPTUPLE &t );
  void inline compare_to_peers( const uint32_t i, const uint32_t start );
  void inline update_partition_map( const uint32_t start, const uint32_t end );

  // Data members:
  const uint32_t num_threads_; /**< Number of threads with which to execute */
  uint32_t n_; /**< Number of input tuples remaining */
  const uint32_t accum_; /**< Size of alpha block of points to concurrently process */
  const uint32_t pq_size_; /**< Number of points to use for each thread in the pre-filter */

  EPTUPLE* data_; /**< Array of input data points */
  vector<int> skyline_; /**< Vector in which the skyline result will be copied */
  vector<pair<uint32_t, uint32_t> > part_map_; /**< Data structure used in Phase I computation */
};

#endif /* HYBRID_H_ */
