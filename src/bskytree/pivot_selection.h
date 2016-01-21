#pragma once

#include <limits.h>
#include <math.h>

#include <vector>

#include "common/common.h"

using namespace std;

class PivotSelection {
public:
	PivotSelection( const vector<float> &min_list, const vector<float> &max_list);
	~PivotSelection(void);

	void Execute( vector<TUPLE>& dataset );

private:

	vector<float> SetRangeList(const vector<float>& min_list,
			const vector<float>& max_list);
	float ComputeDistance(const float* value, const vector<float>& min_list,
			const vector<float>& range_list);

	bool EvaluatePoint(const unsigned pos, vector<TUPLE>& dataset);

	const vector<float> &min_list_;
	const vector<float> &max_list_;
};

