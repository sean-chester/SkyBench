#include "bskytree/pivot_selection.h"

#include <cstdio>
#include <cassert>

PivotSelection::PivotSelection(const vector<float> &min_list, const vector<float> &max_list) :
    min_list_( min_list ), max_list_( max_list ) {
    
}

PivotSelection::~PivotSelection(void) {

}


/**
 * Chooses a pivot based on minimum range. The chosen pivot
 * is a skyline point. In addition to that, removes points 
 * from dataset that are dominated by the (current) pivot point.
 */
void PivotSelection::Execute(vector<TUPLE>& dataset) {

  const uint32_t head = 0;
  uint32_t tail = dataset.size() - 1, cur_pos = 1;
  float* hvalue = dataset[head].elems;

  vector<float> range_list = SetRangeList( min_list_, max_list_ );
  float min_dist = ComputeDistance( hvalue, min_list_, range_list );

  while ( cur_pos <= tail ) {
    float* cvalue = dataset[cur_pos].elems;

    const uint32_t dtest = DominanceTest( dataset[head], dataset[cur_pos] );
    if ( dtest == DOM_LEFT ) {
      dataset[cur_pos] = dataset[tail];
      dataset.pop_back();
      tail--;
    } else if ( dtest == DOM_RIGHT ) {
      dataset[head] = dataset[cur_pos];
      dataset[cur_pos] = dataset[tail];
      dataset.pop_back();
      tail--;

      hvalue = dataset[head].elems;
      min_dist = ComputeDistance( hvalue, min_list_, range_list );
      cur_pos = 1; // THIS IS THE SAME BUG AS IN QSkyCube: cur_pos is not reseted
    } else {
      assert( dtest == DOM_INCOMP );
      float cur_dist = ComputeDistance( cvalue, min_list_, range_list );

      if ( cur_dist < min_dist ) {
        if ( EvaluatePoint( cur_pos, dataset ) ) {
          std::swap( dataset[head], dataset[cur_pos] );

          hvalue = dataset[head].elems;
          min_dist = cur_dist;
          cur_pos++;
        } else {
          dataset[cur_pos] = dataset[tail];
          dataset.pop_back();
          tail--;
        }
      } else
        cur_pos++;
    }
  }
}

vector<float> PivotSelection::SetRangeList(const vector<float>& min_list,
    const vector<float>& max_list) {
  vector<float> range_list( NUM_DIMS, 0 );
  for (uint32_t d = 0; d < NUM_DIMS; d++)
    range_list[d] = max_list[d] - min_list[d];

  return range_list;
}

/**
 * Note that here normalization must be done (even though we assume
 * the data is pre-normalized) because it spreads the values within
 * each (recursed) partition (where all values are within a range).
 */
float PivotSelection::ComputeDistance(const float* value,
    const vector<float>& min_list, const vector<float>& range_list) {
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

/**
 * Checks if the point dataset[pos] is not dominated by any of points
 * before pos (dataset[0..pos-1]).
 *
 * Note that here we can remove additionally dominated points, but the
 * code does not do it (the paper suggests, though).
 */
bool PivotSelection::EvaluatePoint(const uint32_t pos, vector<TUPLE>& dataset) {
  const TUPLE &cur_tuple = dataset[pos];
  for (uint32_t i = 0; i < pos; ++i) {
    const TUPLE &prev_value = dataset[i];
    if ( DominatedLeft( cur_tuple, prev_value ) )
      return false;
  }

  return true;
}
