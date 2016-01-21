/*
 * pskyline.h
 *
 *  Created on: Feb 5, 2014
 *      Author: dariuss
 *
 *  Based on acquired implementation from authors of:
 *    Sungwoo Park, Taekyung Kim, Jonghyun Park, Jinha Kim, Hyeonseung Im,
 *    "Parallel Skyline Computation on Multicore Architectures,"
 *    ICDE '09., pp.760,771, 2009.
 *
 *    Hyeonseung Im, Jonghyun Park, Sungwoo Park,
 *    "Parallel skyline computation on multicore architectures,"
 *    Information Systems, Volume 36, Issue 4, June 2011, Pages 808-823.
 */

#ifndef PSKYLINE_H_
#define PSKYLINE_H_

#include <stdint.h>

#include <vector>

#include "common/common.h"
#include "common/skyline_i.h"

using namespace std;

typedef struct Block {
  int start; // data[start --- end], inclusive
  int end; // flag[start --- end], inclusive
  uint32_t size; // skyline size after find_skyline
  // data[(start --- (start + size - 1)], inclusive after find_skyline
} Block;

class PSkyline: public SkylineI {
public:
  PSkyline(uint32_t threads, uint32_t tuples, uint32_t dims, float** data);
  virtual ~PSkyline();

  vector<int> Execute();

private:
  inline int CheckSurvival(TUPLE x, TUPLE* s, int* flag, int size) {
    for (int i = 0; i < size; i++) {
      if ( flag[i] == DEAD )
        continue;
      const int dtest = DominanceTest( x, s[i] );
      if ( dtest == DOM_LEFT )
        flag[i] = DEAD;
      else if ( dtest == DOM_RIGHT )
        return DEAD;

//      switch (DT( x, s[i] )) {
//      case DOM_LEFT:
//        flag[i] = DEAD;
//        break;
//      case DOM_RIGHT:
//        return DEAD;
//      }
    }
    return LIVE;
  }

  void Init(float** data);
  Block sskyline(Block input);
  Block PMerge(Block left, Block right);
  Block* PMap(Block* input);
  Block SReduce(Block* input);

  // Data members:
  const uint32_t num_threads_; // the same as num_blocks (one thread per block)
  const uint32_t n_; // #tuples
  const uint32_t d_; // #dims
  const uint32_t block_size_;

  TUPLE* data_;
  Block* input_;
  int* flag_;
  vector<int> skyline_;
};

#endif /* PSKYLINE_H_ */
