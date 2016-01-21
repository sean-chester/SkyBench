/*
 * qflow.h
 *
 *  Created on: Feb 12, 2014
 *      Author: schester
 *
 *	QFLOW based:
 *		S. Chester et al., "Skyline in review,"
 *		PVLDB 8(2): 100--111. 2014.
 */

#ifndef QFLOW_H_
#define QFLOW_H_

#include "common/common.h"
#include "common/skyline_i.h"

using namespace std;

class QFlow: public SkylineI {
public:
  QFlow( uint32_t threads, uint32_t tuples, uint32_t dims, float** data,
      uint32_t accum );
  virtual ~QFlow();

  vector<int> Execute();

private:
  void Init( float** data );
  static bool STupleComp( STUPLE const &lhs, STUPLE const &rhs );
  int skyline();
  void ComputeScores();

  // Data members:
  const uint32_t num_threads_;
  const uint32_t n_;
  const uint32_t accum_;

  STUPLE* data_;
  vector<int> skyline_;

};

#endif /* QFLOW_H_ */
