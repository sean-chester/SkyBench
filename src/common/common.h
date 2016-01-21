/*
 * common.h
 *
 *  Created on: Feb 11, 2014
 *      Author: dariuss
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <stdint.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cassert>

#define DOM_LEFT   0
#define DOM_RIGHT  1
#define DOM_INCOMP  2
#define P_ACCUM 256
#define BSKYTREE_ACCUM 256
#define DEFAULT_ALPHA 1024 // previous Q_ACCUM
#define DEFAULT_QP_SIZE 8

#define PRUNED (NUM_DIMS << 2)
#define ALL_ONES ((1<<NUM_DIMS) - 1)
#define DEAD 1
#define LIVE 0

#define PIVOT_RANDOM  0       // BSkyTree.Random
#define PIVOT_MEDIAN  1
#define PIVOT_BALANCED  2
#define PIVOT_BALSKY  3       // BSkyTree.Balanced
#define PIVOT_MANHATTAN 4     
#define PIVOT_VOLUME 5     // BSkyTree.MaxDom

static const uint32_t SHIFTS[] = { 1 << 0, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1
    << 5, 1 << 6, 1 << 7, 1 << 8, 1 << 9, 1 << 10, 1 << 11, 1 << 12, 1 << 13, 1
    << 14, 1 << 15, 1 << 16, 1 << 17, 1 << 18, 1 << 19, 1 << 20, 1 << 21, 1
    << 22, 1 << 23, 1 << 24, 1 << 25, 1 << 26, 1 << 27, 1 << 28, 1 << 29, 1
    << 30 };

typedef struct TUPLE {
  float elems[NUM_DIMS];
  int pid;
//  TUPLE(const int id, float* data): pid(id) {
//    for (uint32_t d = 0; d < NUM_DIMS; ++d) {
//      elems[d] = data[d];
//    }
//  }

  void printTuple() {
    printf( "[" );
    for (uint32_t i = 0; i < NUM_DIMS; ++i)
      printf( "%f ", elems[i] );
    printf( "]\n" );
  }

} TUPLE;

// Sort-based Tuple
typedef struct STUPLE: TUPLE {
  float score; // entropy, manhattan sum, or minC

  // By default, STUPLEs are sorted by score value
  bool operator<(const STUPLE &rhs) const {
    return score < rhs.score;
  }
} STUPLE;

typedef struct TUPLE_S: TUPLE {
  uint32_t partition; // bitset: 0 is <= pivot, 1 is > pivot
  TUPLE_S(const TUPLE t, const uint32_t p):
    TUPLE(t), partition(p) { }
} TUPLE_S;

// Partition-based Tuple
typedef struct PTUPLE: STUPLE {
  uint32_t partition; // bitset; 0 is <= pivot, 1 is >
  uint32_t partition_level; // level of this tuple's partition ends.
  uint32_t partition_end; // index where this tuple's partition ends.

  /* Natural order is first by partition level,
   * then by partition id, then
   * by score. */
  bool operator<(PTUPLE const &rhs) const {
    if ( partition_level < rhs.partition_level )
      return true;
    if ( rhs.partition_level < partition_level )
      return false;
    if ( partition < rhs.partition )
      return true;
    if ( rhs.partition < partition )
      return false;
    if ( score < rhs.score )
      return true;
    if ( rhs.score < score )
      return false;
    return false; // points are equal (by relation < ).
  }

  /* using partition == ALL_ONES to denote that a 
   * point is pruned (rather than adding a new bool
   * member) to keep the PTUPLE struct smaller.
   */
  inline void markPruned() {
    partition_level = NUM_DIMS;
  }
  inline bool isPruned() {
    return partition_level == NUM_DIMS;
  }

  /* Can skip other partition if there are bits that he has
   * and I don't (therefore, he cannot dominate me). */
  inline bool canskip_partition(const uint32_t other) {
    return (partition ^ other) & other;
  }
} PTUPLE;

// (Encoded) partition-based tuple with both partition level
// and partition bitmask encoded into uint32_t.
typedef struct EPTUPLE: STUPLE {
  uint32_t partition; // bitset; 0 is <= pivot, 1 is >

  /* Natural order is first by partition level,
   * then by partition id, then by score. */
  bool operator<(EPTUPLE const &rhs) const {
    if ( partition < rhs.partition )
      return true;
    if ( rhs.partition < partition )
      return false;
    if ( score < rhs.score )
      return true;
    if ( rhs.score < score )
      return false;
    return false; // points are equal (by relation < ).
  }

  /* (NUM_DIMS << NUM_DIMS) denotes that a
   * point is pruned (level with all bits set).
   */
  inline void markPruned() {
    partition = NUM_DIMS << NUM_DIMS;
  }
  inline bool isPruned() const {
    return partition == NUM_DIMS << NUM_DIMS;
  }

  /* Can skip other partition if there are bits that he has
   * and I don't (therefore, he cannot dominate me). */
  inline bool canskip_partition(const uint32_t other) const {
    return (getPartition() ^ other) & other;
  }

  inline uint32_t getLevel() const {
    return partition >> NUM_DIMS;
  }
  inline uint32_t getPartition() const {
    return partition & ALL_ONES;
  }
  inline void setPartition(const uint32_t p_bitmap) {
    partition = (__builtin_popcount(p_bitmap) << NUM_DIMS) | p_bitmap;
  }
} PTUPLE2;

extern uint64_t dt_count;
extern uint64_t dt_count_dom;
extern uint64_t dt_count_incomp;

// returns the maximum attribute value
inline float get_max(const STUPLE &p) {
  float maxc = p.elems[0];
  for (uint32_t d = 1; d < NUM_DIMS; d++) {
    maxc = std::max( maxc, p.elems[d] );
  }
  return maxc;
}



#if __AVX__

#include "common/dt_avx.h"

#else

/*
 * 2-way dominance test with NO assumption for distinct value condition.
 */
inline int DominanceTest(const TUPLE &t1, const TUPLE &t2) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  bool t1_better = false, t2_better = false;

  for (uint32_t i = 0; i < NUM_DIMS; i++) {
    if ( t1.elems[i] < t2.elems[i] )
      t1_better = true;
    else if ( t1.elems[i] > t2.elems[i] )
      t2_better = true;

    if ( t1_better && t2_better )
      return DOM_INCOMP;
  }
  if ( !t1_better && t2_better )
    return DOM_RIGHT;
  if ( !t2_better && t1_better )
    return DOM_LEFT;

//    if ( !t1_better && !t2_better )
//      return DOM_INCOMP; //equal
  return DOM_INCOMP;
}

/*
 * Dominance test returning result as a bitmap.
 * This is an original version (assuming distinct value
 * condition) used in in BSkyTree.
 *
 * In BSkyTree, it is by far the most frequent dominance test.
 */
inline uint32_t DT_bitmap_dvc(const TUPLE &cur_value, const TUPLE &sky_value) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif

  uint32_t lattice = 0;
  for (uint32_t dim = 0; dim < NUM_DIMS; dim++)
    if ( sky_value.elems[dim] <= cur_value.elems[dim] )
      lattice |= SHIFTS[dim];

#if COUNT_DT==1
  if ( lattice == ALL_ONES)
    __sync_fetch_and_add( &dt_count_dom, 1 );
  else
    __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif 
  return lattice;
}

/*
 * The same as above, but no assuming DVC.
 *
 * Note: is not called so frequently as DT_bitmap_dvc in BSkyTree,
 * so performance gain is negligible.
 */
inline uint32_t DT_bitmap(const TUPLE &cur_value, const TUPLE &sky_value) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif

  uint32_t lattice = 0;
  for (uint32_t dim = 0; dim < NUM_DIMS; dim++)
    if ( sky_value.elems[dim] < cur_value.elems[dim] )
      lattice |= SHIFTS[dim];

#if COUNT_DT==1
  if ( lattice == ALL_ONES)
    __sync_fetch_and_add( &dt_count_dom, 1 );
  else
    __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif 
  return lattice;
}


/*
 * One-way (optimized) dominance test.
 * No assumption for distinct value condition.
 */
inline bool DominateLeft(const TUPLE &t1, const TUPLE &t2) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  uint32_t i;
  for (i = 0; i < NUM_DIMS && t1.elems[i] <= t2.elems[i]; ++i)
    ;
  if ( i < NUM_DIMS ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
    return false; // Points are incomparable.
  }

  for (i = 0; i < NUM_DIMS; ++i) {
    if ( t1.elems[i] < t2.elems[i] ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_dom, 1 );
#endif
      return true; // t1 dominates t2
    }
  }
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
  return false; // Points are equal.
}

/*
 * Dominance test assuming distinct value condition.
 * DominateLeft(x, x) returns 1.
 */
inline int DominateLeftDVC(const TUPLE &t1, const TUPLE &t2) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  for (uint32_t i = 0; i < NUM_DIMS; i++) {
    if ( t1.elems[i] > t2.elems[i] ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
      return 0;
    }
  }
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_dom, 1 );
#endif
  return 1;
}

#endif

/*
 * Dominance test that computes a bitmap.
 * Produces side-effects in cur_value: sets score 
 * and partition.
 */
inline void DT_bitmap_withsum(EPTUPLE &cur_value, const EPTUPLE &sky_value) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  uint32_t partition = 0;
  cur_value.score = 0;
  for (uint32_t d = 0; d < NUM_DIMS; d++) {
    if ( sky_value.elems[d] < cur_value.elems[d] )
      partition |= SHIFTS[d];
    cur_value.score += cur_value.elems[d];
  }
#if COUNT_DT==1
  if ( partition == ALL_ONES)
    __sync_fetch_and_add( &dt_count_dom, 1 );
  else
    __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
  if ( partition == ALL_ONES ) {
    cur_value.markPruned();
  } else {
    cur_value.setPartition( partition );
  }
}

/*
 * Dominance test assuming distinct value condition.
 */
inline int DT_dvc(const TUPLE &t1, const TUPLE &t2) {
  bool t1_better = false, t2_better = false;
  for (uint32_t d = 0; d < NUM_DIMS; d++) {
//    if ( t1.elems[d] < t2.elems[d] )
//      t1_better = true;
//    else if ( t1.elems[d] > t2.elems[d] )
//      t2_better = true;
    t1_better = t1.elems[d] < t2.elems[d] || t1_better;
    t2_better = t1.elems[d] > t2.elems[d] || t2_better;

    if ( t1_better && t2_better ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
      return DOM_INCOMP;
    }
  }

  if ( !t1_better && t2_better ) {
    return DOM_RIGHT;
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_dom, 1 );
#endif
  }
  if ( !t2_better && t1_better ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_dom, 1 );
#endif
    return DOM_LEFT;
  }

  //    if ( !t1_better && !t2_better )
  //      return DOM_UNCOMP; //equal
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
  assert( false );
  return DOM_INCOMP;
}

/*
 * Returns true if left tuple is dominated.
 *
 * Used in EvaluatePoint. Assumes DVC. Identical to DominateRightDVC().
 */
inline bool DominatedLeft(const TUPLE &t1, const TUPLE &t2) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif

  for (uint32_t d = 0; d < NUM_DIMS; d++) {
    if ( t1.elems[d] < t2.elems[d] ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
      return false;
    }
  }
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_dom, 1 );
#endif
  return true;
}

/*
 * Dominance test assuming distinct value condition.
 * DominateLeft(x, x) returns 1.
 */
inline int DominateRightDVC(const TUPLE &t1, const TUPLE &t2) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  for (uint32_t i = 0; i < NUM_DIMS; i++) {
    if ( t1.elems[i] < t2.elems[i] ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
      return 0;
    }
  }
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_dom, 1 );
#endif
  return 1;
}

/*
 * Test for equality.
 */
inline bool EqualityTest(const TUPLE &t1, const TUPLE &t2) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  bool eq = true;
  for (uint32_t d = 0; d < NUM_DIMS; d++)
    if ( t1.elems[d] != t2.elems[d] ) {
      eq = false;
      break;
    }
  return eq;
}

inline float calc_norm_range ( const TUPLE &t, const float *mins, const float *ranges ) {
  float min, max;
  min = max = ( t.elems[0] - mins[0] ) / ranges[0];
  for (uint32_t j = 1; j < NUM_DIMS; j++) {
    const float v_norm = ( t.elems[j] - mins[j] ) / ranges[j];
    if ( min > v_norm )
      min = v_norm;
    else if ( max < v_norm )
      max = v_norm;
  }
  return max - min;
}

/* Below code is replicated, but these versions are needed 
 * for micro-benchmarking different types of DTs.
 */
inline bool DominateLeftNOAVX(const TUPLE &t1, const TUPLE &t2) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  uint32_t i;
  for (i = 0; i < NUM_DIMS && t1.elems[i] <= t2.elems[i]; ++i)
    ;
  if ( i < NUM_DIMS ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
    return false; // Points are incomparable.
  }

  for (i = 0; i < NUM_DIMS; ++i) {
    if ( t1.elems[i] < t2.elems[i] ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_dom, 1 );
#endif
      return true; // t1 dominates t2
    }
  }
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
  return false; // Points are equal.
}

inline int DominanceTestNOAVX(const TUPLE &t1, const TUPLE &t2) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  bool t1_better = false, t2_better = false;

  for (uint32_t i = 0; i < NUM_DIMS; i++) {
    if ( t1.elems[i] < t2.elems[i] )
      t1_better = true;
    else if ( t1.elems[i] > t2.elems[i] )
      t2_better = true;

    if ( t1_better && t2_better )
      return DOM_INCOMP;
  }
  if ( !t1_better && t2_better )
    return DOM_RIGHT;
  if ( !t2_better && t1_better )
    return DOM_LEFT;

//    if ( !t1_better && !t2_better )
//      return DOM_INCOMP; //equal
  return DOM_INCOMP;
}

inline uint32_t DT_bitmap_NOAVX(const TUPLE &cur_value, const TUPLE &sky_value) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif

  uint32_t lattice = 0;
  for (uint32_t dim = 0; dim < NUM_DIMS; dim++)
    if ( sky_value.elems[dim] < cur_value.elems[dim] )
      lattice |= SHIFTS[dim];

#if COUNT_DT==1
  if ( lattice == ALL_ONES)
    __sync_fetch_and_add( &dt_count_dom, 1 );
  else
    __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
  return lattice;
}

inline uint32_t DT_bitmap_dvc_NOAVX(const TUPLE &cur_value, const TUPLE &sky_value) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif

  uint32_t lattice = 0;
  for (uint32_t dim = 0; dim < NUM_DIMS; dim++)
    if ( sky_value.elems[dim] <= cur_value.elems[dim] )
      lattice |= SHIFTS[dim];

#if COUNT_DT==1
  if ( lattice == ALL_ONES)
    __sync_fetch_and_add( &dt_count_dom, 1 );
  else
    __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
  return lattice;
}

inline int DominateLeftDVC_NOAVX(const TUPLE &t1, const TUPLE &t2) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  for (uint32_t i = 0; i < NUM_DIMS; i++) {
    if ( t1.elems[i] > t2.elems[i] ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif
      return 0;
    }
  }
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count_dom, 1 );
#endif
  return 1;
}

#endif /* COMMON_H_ */
