#pragma once

//#include "bskytree/point.h"
#include "bskytree/node.h"
#include "bskytree/pivot_selection.h"
#include "common/skyline_i.h"
#include "common/common.h"

#include <map>
#include <vector>

using namespace std;

class SkyTree: public SkylineI {

public:
	SkyTree(const uint32_t n, const uint32_t d, float** dataset, 
    const bool useTree, const bool useDnC );
	~SkyTree(void);

	void Init(float** dataset);
	vector<int> Execute(void);

private:
	void ComputeSkyTree(const vector<float> min_list,
			const vector<float> max_list, vector<TUPLE> &dataset,
			Node& skytree );

	map<uint32_t, vector<TUPLE> > MapPointToRegion(vector<TUPLE>& dataset);

  void PartialDominance(const uint32_t lattice, vector<TUPLE>& dataset,
			Node& skytree );
  bool PartialDominance_with_trees(const uint32_t lattice, Node& left_tree,
      Node& right_tree );
	bool FilterPoint(const TUPLE &cur_value, Node& skytree);
  bool FilterPoint_without_skytree(const TUPLE &cur_value, Node& skytree);
	void TraverseSkyTree(const Node& skytree);

#ifndef NVERBOSE
	int MaxDepth(const Node& skytree, int d);
#endif

	const uint32_t n_;
	const uint32_t d_;
	vector<TUPLE> data_;

	vector<float> min_list_;
	vector<float> max_list_;

	Node skytree_;
	vector<int> skyline_;
	vector<int> eqm_; // "equivalence matrix"
  
  /* runtime params. */
  bool useTree_; //using SkyTree data structure in FilterPoints()
  bool useDnC_; //divide-and-conquer
  bool *dominated_; //for DnC variant

#ifndef NVERBOSE
  map<int, int> skytree_levels_;
#endif
};
