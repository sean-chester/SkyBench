/*
 * parallel_bskytree.h
 *
 *  Created on: Jul 6, 2014
 *      Author: dariuss
 *
 *
 *  Our parallel variant of BSkyTree-S algorithm.
 */

#ifndef PARALLEL_BSKYTREE_H_
#define PARALLEL_BSKYTREE_H_

#include <map>
#include <vector>

#include <common/skyline_i.h>
#include <bskytree/node.h>

using namespace std;

class ParallelBSkyTree: public SkylineI {
public:
  ParallelBSkyTree( const uint32_t num_threads, const uint32_t n,
      const uint32_t d, float** dataset );
  virtual ~ParallelBSkyTree();

  void Init( float** dataset );
  vector<int> Execute( void );

private:
  void BSkyTreeS_ALGO();
  void DoPartioning();

  // PivotSelection methods
  void SelectBalanced();
  vector<float> SetRangeList( const vector<float>& min_list,
      const vector<float>& max_list );
  float ComputeDistance( const float* value, const vector<float>& min_list,
      const vector<float>& range_list );
  bool EvaluatePoint( const uint32_t pos );

  const uint32_t num_threads_;
  const uint32_t n_;
  const uint32_t d_;
  vector<TUPLE_S> data_;

  vector<int> skyline_;
  vector<int> eqm_; // "equivalence matrix"
};

#endif /* PARALLEL_BSKYTREE_H_ */
