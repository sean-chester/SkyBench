/*
 * dt_avx.h
 *
 *  Created on: Mar 14, 2014
 *      Author: dariuss
 *
 *  Dominance tests using AVX/SSE instructions.
 */

#ifndef DT_AVX_H_
#define DT_AVX_H_

#include <immintrin.h>  // AVX

inline uint32_t DT_bitmap_dvc( const TUPLE &cur, const TUPLE &sky ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif

  uint32_t lattice = 0;
  uint32_t dim = 0;

#if NUM_DIMS >= 8
  const TUPLE __attribute__ ((aligned(32))) cur_value = cur;
  const TUPLE __attribute__ ((aligned(32))) sky_value = sky;

  for (; dim+8 <= NUM_DIMS; dim+=8) {
    __m256 p_ymm = _mm256_load_ps( cur_value.elems + dim);
    __m256 sky_ymm = _mm256_load_ps( sky_value.elems + dim);
    __m256 comp_le = _mm256_cmp_ps(sky_ymm, p_ymm, 2);
    uint32_t le_mask = _mm256_movemask_ps(comp_le);
    lattice = lattice | (le_mask << dim);
  }

  for (; dim+4 <= NUM_DIMS; dim+=4) {
    __m128 p_xmm = _mm_load_ps(cur_value.elems + dim);
    __m128 sky_xmm = _mm_load_ps(sky_value.elems + dim);
    __m128 le128 = _mm_cmp_ps(sky_xmm, p_xmm, 2);
    uint32_t le_mask = _mm_movemask_ps(le128);
    lattice = lattice | (le_mask << dim);
  }

  for (; dim < NUM_DIMS; ++dim)
  if ( sky.elems[dim] <= cur.elems[dim] )
  lattice |= SHIFTS[dim];

#elif NUM_DIMS >= 4
  const TUPLE __attribute__ ((aligned(32))) cur_value = cur;
  const TUPLE __attribute__ ((aligned(32))) sky_value = sky;

  for (; dim + 4 <= NUM_DIMS; dim += 4) {
    __m128 p_xmm = _mm_load_ps( cur_value.elems + dim );
    __m128 sky_xmm = _mm_load_ps( sky_value.elems + dim );
    __m128 comp_le = _mm_cmp_ps( sky_xmm, p_xmm, 2 );
    uint32_t le_mask = _mm_movemask_ps( comp_le );
    lattice = lattice | (le_mask << dim);
  }

  for (; dim < NUM_DIMS; ++dim)
    if ( sky.elems[dim] <= cur.elems[dim] )
      lattice |= SHIFTS[dim];

#else

  for (; dim < NUM_DIMS; ++dim)
  if ( sky.elems[dim] <= cur.elems[dim] )
  lattice |= SHIFTS[dim];

#endif
#if COUNT_DT==1
  if ( lattice == ALL_ONES )
    __sync_fetch_and_add( &dt_count_dom, 1 );
  else
    __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif

  return lattice;
}

inline uint32_t DT_bitmap( const TUPLE cur, const TUPLE sky ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  uint32_t lattice = 0;
  uint32_t dim = 0;

#if NUM_DIMS >= 8
  const TUPLE __attribute__ ((aligned(32))) cur_value = cur;
  const TUPLE __attribute__ ((aligned(32))) sky_value = sky;

  for (; dim+8 <= NUM_DIMS; dim+=8) {
    __m256 p_ymm = _mm256_load_ps( cur_value.elems + dim);
    __m256 sky_ymm = _mm256_load_ps( sky_value.elems + dim);
    __m256 comp_lt = _mm256_cmp_ps(sky_ymm, p_ymm, 1);
    uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
    lattice = lattice | (lt_mask << dim);
  }

  for (; dim+4 <= NUM_DIMS; dim+=4) {
    __m128 p_xmm = _mm_load_ps(cur_value.elems + dim);
    __m128 sky_xmm = _mm_load_ps(sky_value.elems + dim);
    __m128 comp_lt = _mm_cmp_ps(sky_xmm, p_xmm, 1);
    uint32_t lt_mask = _mm_movemask_ps(comp_lt);
    lattice = lattice | (lt_mask << dim);
  }

  for (; dim < NUM_DIMS; dim++)
  if ( sky.elems[dim] < cur.elems[dim] )
  lattice |= SHIFTS[dim];

#elif NUM_DIMS >= 4
  const TUPLE __attribute__ ((aligned(32))) cur_value = cur;
  const TUPLE __attribute__ ((aligned(32))) sky_value = sky;

  for (; dim + 4 <= NUM_DIMS; dim += 4) {
    __m128 p_xmm = _mm_load_ps( cur_value.elems + dim );
    __m128 sky_xmm = _mm_load_ps( sky_value.elems + dim );
    __m128 comp_lt = _mm_cmp_ps( sky_xmm, p_xmm, 1 );
    uint32_t lt_mask = _mm_movemask_ps( comp_lt );
    lattice = lattice | (lt_mask << dim);
  }

  for (; dim < NUM_DIMS; dim++)
    if ( sky.elems[dim] < cur.elems[dim] )
      lattice |= SHIFTS[dim];

#else
  for (dim = 0; dim < NUM_DIMS; dim++)
  if ( sky.elems[dim] < cur.elems[dim] )
  lattice |= SHIFTS[dim];
#endif

#if COUNT_DT==1
  if ( lattice == ALL_ONES )
    __sync_fetch_and_add( &dt_count_dom, 1 );
  else
    __sync_fetch_and_add( &dt_count_incomp, 1 );
#endif

  return lattice;
}

/*
 * 2-way dominance test with NO assumption for distinct value condition.
 */
inline int DominanceTest( const TUPLE &left, const TUPLE &right ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  uint32_t dim = 0, left_better = 0, right_better = 0;

#if NUM_DIMS >= 8
  const TUPLE __attribute__ ((aligned(32))) right_value = right;
  const TUPLE __attribute__ ((aligned(32))) left_value = left;

  for (; dim+8 <= NUM_DIMS; dim+=8) {
    __m256 right_ymm = _mm256_load_ps( right_value.elems + dim);
    __m256 left_ymm = _mm256_load_ps( left_value.elems + dim);
    if ( !left_better ) {
      __m256 comp_lt = _mm256_cmp_ps(left_ymm, right_ymm, 1);
      uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
      left_better = lt_mask & (SHIFTS[8] - 1);
    }
    if ( !right_better ) {
      __m256 comp_lt = _mm256_cmp_ps(right_ymm, left_ymm, 1);
      uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
      right_better = lt_mask & (SHIFTS[8] - 1);
    }
    if( left_better && right_better ) return DOM_INCOMP;
  }

  for (; dim+4 <= NUM_DIMS; dim+=4) {
    __m128 right_ymm = _mm_load_ps( right_value.elems + dim);
    __m128 left_ymm = _mm_load_ps( left_value.elems + dim);
    if ( !left_better ) {
      __m128 comp_lt = _mm_cmp_ps(left_ymm, right_ymm, 1);
      uint32_t lt_mask = _mm_movemask_ps(comp_lt);
      left_better = lt_mask & (SHIFTS[4] - 1);
    }
    if ( !right_better ) {
      __m128 comp_lt = _mm_cmp_ps(right_ymm, left_ymm, 1);
      uint32_t lt_mask = _mm_movemask_ps(comp_lt);
      right_better = lt_mask & (SHIFTS[4] - 1);
    }
    if( left_better && right_better ) return DOM_INCOMP;
  }

  for (; dim < NUM_DIMS; dim++) {
    if ( !right_better && right.elems[dim] < left.elems[dim] )
    right_better = 1;
    if ( !left_better && left.elems[dim] < right.elems[dim] )
    left_better = 1;
    if( left_better && right_better ) return DOM_INCOMP;
  }

#elif NUM_DIMS >= 4
  const TUPLE __attribute__ ((aligned(32))) right_value = right;
  const TUPLE __attribute__ ((aligned(32))) left_value = left;

  for (; dim + 4 <= NUM_DIMS; dim += 4) {
    __m128 right_ymm = _mm_load_ps( right_value.elems + dim );
    __m128 left_ymm = _mm_load_ps( left_value.elems + dim );
    if ( !left_better ) {
      __m128 comp_lt = _mm_cmp_ps( left_ymm, right_ymm, 1 );
      uint32_t lt_mask = _mm_movemask_ps( comp_lt );
      left_better = lt_mask & (SHIFTS[4] - 1);
    }
    if ( !right_better ) {
      __m128 comp_lt = _mm_cmp_ps( right_ymm, left_ymm, 1 );
      uint32_t lt_mask = _mm_movemask_ps( comp_lt );
      right_better = lt_mask & (SHIFTS[4] - 1);
    }
    if ( left_better && right_better )
      return DOM_INCOMP;
  }

  for (; dim < NUM_DIMS; dim++) {
    if ( !right_better && right.elems[dim] < left.elems[dim] )
      right_better = 1;
    if ( !left_better && left.elems[dim] < right.elems[dim] )
      left_better = 1;
    if ( left_better && right_better )
      return DOM_INCOMP;
  }

#else
  for (dim = 0; dim < NUM_DIMS; dim++) {
    if ( !right_better ) right_better = ( right.elems[dim] < left.elems[dim] );
    if ( !left_better ) left_better = ( left.elems[dim] < right.elems[dim] );
    if( left_better && right_better ) return DOM_INCOMP;
  }
#endif

  if ( left_better && !right_better )
    return DOM_LEFT;
  else if ( right_better && !left_better )
    return DOM_RIGHT;
  else
    return DOM_INCOMP; //equal.
}

/*
 * One-way (optimized) dominance test.
 * No assumption for distinct value condition.
 */
inline bool DominateLeft( const TUPLE &left, const TUPLE &right ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  uint32_t dim = 0;

#if NUM_DIMS >= 8
  const TUPLE __attribute__ ((aligned(32))) right_value = right;
  const TUPLE __attribute__ ((aligned(32))) left_value = left;

  for (; dim+8 <= NUM_DIMS; dim+=8) {
    __m256 right_ymm = _mm256_load_ps( right_value.elems + dim);
    __m256 left_ymm = _mm256_load_ps( left_value.elems + dim);
    __m256 comp_lt = _mm256_cmp_ps(left_ymm, right_ymm, 2);
    uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
    if( lt_mask != 255 ) return false;
  }

  for (; dim+4 <= NUM_DIMS; dim+=4) {
    __m128 right_xmm = _mm_load_ps(right_value.elems + dim);
    __m128 left_xmm = _mm_load_ps(left_value.elems + dim);
    __m128 comp_lt = _mm_cmp_ps(left_xmm, right_xmm, 2);
    uint32_t lt_mask = _mm_movemask_ps(comp_lt);
    if( lt_mask != 15 ) return false;
  }

  for (; dim < NUM_DIMS; dim++)
  if ( right.elems[dim] < left.elems[dim] )
  return false;

#elif NUM_DIMS >= 4
  const TUPLE __attribute__ ((aligned(32))) right_value = right;
  const TUPLE __attribute__ ((aligned(32))) left_value = left;

  for (; dim + 4 <= NUM_DIMS; dim += 4) {
    __m128 right_xmm = _mm_load_ps( right_value.elems + dim );
    __m128 left_xmm = _mm_load_ps( left_value.elems + dim );
    __m128 comp_lt = _mm_cmp_ps( left_xmm, right_xmm, 2 );
    uint32_t lt_mask = _mm_movemask_ps( comp_lt );
    if ( lt_mask != 15 )
      return false;
  }

  for (; dim < NUM_DIMS; dim++)
    if ( right.elems[dim] < left.elems[dim] )
      return false;

#else
  for (dim = 0; dim < NUM_DIMS; dim++)
  if ( right.elems[dim] < left.elems[dim] )
  return false;
#endif

  //test equality.
  for (dim = 0; dim < NUM_DIMS; dim++)
    if ( right.elems[dim] != left.elems[dim] )
      return true;

  return false; //points are equal.
}

/*
 * One-way (optimized) dominance test.
 * With distinct value condition assumption.
 */
inline bool DominateLeftDVC( const TUPLE &left, const TUPLE &right ) {
#if COUNT_DT==1
  __sync_fetch_and_add( &dt_count, 1 );
#endif
  uint32_t dim = 0;

#if NUM_DIMS >= 8
  const TUPLE __attribute__ ((aligned(32))) right_value = right;
  const TUPLE __attribute__ ((aligned(32))) left_value = left;

  for (; dim+8 <= NUM_DIMS; dim+=8) {
    __m256 right_ymm = _mm256_load_ps( right_value.elems + dim);
    __m256 left_ymm = _mm256_load_ps( left_value.elems + dim);
    __m256 comp_lt = _mm256_cmp_ps(left_ymm, right_ymm, 2);
    uint32_t lt_mask = _mm256_movemask_ps(comp_lt);
    if( lt_mask != 255 ) return false;
  }

  for (; dim+4 <= NUM_DIMS; dim+=4) {
    __m128 right_xmm = _mm_load_ps(right_value.elems + dim);
    __m128 left_xmm = _mm_load_ps(left_value.elems + dim);
    __m128 comp_lt = _mm_cmp_ps(left_xmm, right_xmm, 2);
    uint32_t lt_mask = _mm_movemask_ps(comp_lt);
    if( lt_mask != 15 ) return false;
  }

  for (; dim < NUM_DIMS; dim++)
  if ( right.elems[dim] < left.elems[dim] )
  return false;

#elif NUM_DIMS >= 4
  const TUPLE __attribute__ ((aligned(32))) right_value = right;
  const TUPLE __attribute__ ((aligned(32))) left_value = left;

  for (; dim + 4 <= NUM_DIMS; dim += 4) {
    __m128 right_xmm = _mm_load_ps( right_value.elems + dim );
    __m128 left_xmm = _mm_load_ps( left_value.elems + dim );
    __m128 comp_lt = _mm_cmp_ps( left_xmm, right_xmm, 2 );
    uint32_t lt_mask = _mm_movemask_ps( comp_lt );
    if ( lt_mask != 15 )
      return false;
  }

  for (; dim < NUM_DIMS; dim++)
    if ( right.elems[dim] < left.elems[dim] )
      return false;

#else
  for (dim = 0; dim < NUM_DIMS; dim++)
  if ( right.elems[dim] < left.elems[dim] )
  return false;
#endif

  return true;
}

#endif /* DT_AVX_H_ */
