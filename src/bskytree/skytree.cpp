#include "bskytree/skytree.h"
#include <cstdio>
#include <cassert>

uint32_t countSkyTree( Node& skytree ) {
  uint32_t count = 0;
  uint32_t num_child = skytree.children.size();
  for (uint32_t c = 0; c < num_child; c++)
    count += countSkyTree( skytree.children[c] );
  return count + 1;
}

SkyTree::SkyTree( const uint32_t n, const uint32_t d, float** dataset,
    const bool useTree, const bool useDnC ) :
    n_( n ), d_( d ), useTree_( useTree ), useDnC_( useDnC ) {

  skytree_.lattice = 0;
  skyline_.reserve( 1024 );
  eqm_.reserve( 1024 );
  dominated_ = NULL;
}

void SkyTree::Init( float** dataset ) {
  data_.reserve( n_ );
  for (uint32_t i = 0; i < n_; i++) {
    TUPLE t;
    t.pid = i;
    for (uint32_t j = 0; j < NUM_DIMS; j++) {
      t.elems[j] = dataset[i][j];
    }
    data_.push_back( t );
  }
  if ( useDnC_ ) {
    const uint32_t n = n_;
    if ( n > 0 ) {
      dominated_ = new bool[n];
      for (uint32_t i = 0; i < n; ++i)
        dominated_[i] = false;
    } else
      useDnC_ = false; //so we don't try to delete[] dominated_.
  }
}

SkyTree::~SkyTree( void ) {
  min_list_.clear();
  max_list_.clear();
  skyline_.clear();
  ClearSkyTree( skytree_ );
  data_.clear();
  if ( useDnC_ )
    delete[] dominated_;
}

vector<int> SkyTree::Execute( void ) {
  const vector<float> min_list( NUM_DIMS, 0.0 );
  const vector<float> max_list( NUM_DIMS, 1.0 );

  ComputeSkyTree( min_list, max_list, data_, skytree_ );
  TraverseSkyTree( skytree_ );
//  printf( " %d\n", MaxDepth(skytree_, 0) );

#ifndef NVERBOSE
//  const uint32_t skytree_size = countSkyTree( skytree_ );
//  const int depth = MaxDepth( skytree_, 0 );
//  printf( " Skytree: total_children=%u, height=%d, skyline size=%lu\n",
//      skytree_size, depth, skyline_.size() );
//  uint32_t num_nodes = 0;
//  for (map<int, int>::const_iterator it = skytree_levels_.begin();
//      it != skytree_levels_.end(); ++it) {
//    printf( "  L%d: %d\n", it->first, it->second );
//    num_nodes += it->second;
//  }
//  assert( num_nodes == skyline_.size() );
//  assert( depth == skytree_levels_.size() );

  if ( useDnC_ ) {
    const uint32_t skytree_size = countSkyTree( skytree_ );
    printf( "Skytree size=%u, Skyline size=%lu\n", skytree_size,
        skyline_.size() );
  }
#endif

  // Add missing points from "equivalence matrix"
  skyline_.insert( skyline_.end(), eqm_.begin(), eqm_.end() );

  return skyline_;
}

void SkyTree::ComputeSkyTree( const vector<float> min_list,
    const vector<float> max_list, vector<TUPLE>& dataset, Node& skytree ) {
  // pivot selection in the dataset
  PivotSelection selection( min_list, max_list );
  selection.Execute( dataset );

  // mapping points to binary vectors representing subregions
  skytree.point = dataset[0];
  map<uint32_t, vector<TUPLE> > point_map = MapPointToRegion( dataset );

  for (map<uint32_t, vector<TUPLE> >::const_iterator it = point_map.begin();
      it != point_map.end(); it++) {
    uint32_t cur_lattice = (*it).first;
    vector<TUPLE> cur_dataset = (*it).second;

    if ( !useDnC_ && skytree.children.size() > 0 )
      PartialDominance( cur_lattice, cur_dataset, skytree ); // checking partial dominance relations

    if ( cur_dataset.size() > 0 ) {
      vector<float> min_list2( NUM_DIMS ), max_list2( NUM_DIMS );
      for (uint32_t d = 0; d < NUM_DIMS; d++) {
        const uint32_t bit = SHIFTS[d];
        if ( (cur_lattice & bit) == bit )
          min_list2[d] = dataset[0].elems[d], max_list2[d] = max_list[d];
        else
          min_list2[d] = min_list[d], max_list2[d] = dataset[0].elems[d];
      }

      Node child_node( cur_lattice );
      ComputeSkyTree( min_list2, max_list2, cur_dataset, child_node ); // recursive call

      if ( useDnC_ && skytree.children.size() > 0 )
        PartialDominance_with_trees( cur_lattice, skytree, child_node ); //pdom

      skytree.children.push_back( child_node );
    }
  }

//  delete point_map;
}

map<uint32_t, vector<TUPLE> > SkyTree::MapPointToRegion(
    vector<TUPLE>& dataset ) {
  const uint32_t pruned = SHIFTS[NUM_DIMS] - 1;

  map<uint32_t, vector<TUPLE> > data_map;

  const TUPLE &pivot = dataset[0];
  for (vector<TUPLE>::const_iterator it = dataset.begin() + 1;
      it != dataset.end(); it++) {

    if ( EqualityTest( pivot, *it ) ) {
      eqm_.push_back( it->pid );
      continue;
    }

    const uint32_t lattice = DT_bitmap_dvc( *it, pivot );
    if ( lattice < pruned ) { // <-- Same fix as below (same if condition). Also doubles dt's.
      assert( !DominateLeft( pivot, *it ) );
      data_map[lattice].push_back( *it );
    }
  }

  return data_map;
}

bool SkyTree::PartialDominance_with_trees( const uint32_t lattice,
    Node& left_tree, Node& right_tree ) {

  uint32_t num_child = right_tree.children.size();

  for (uint32_t c = 0; c < num_child; ++c) {
    if ( PartialDominance_with_trees( lattice, left_tree,
        right_tree.children[c] ) ) {
      if ( right_tree.children[c].children.size() == 0 ) {
        right_tree.children.erase( right_tree.children.begin() + c-- );
        --num_child;
      }
    }
  }

  num_child = left_tree.children.size();
  for (uint32_t c = 0; c < num_child; c++) {
    uint32_t cur_lattice = left_tree.children[c].lattice;
    if ( cur_lattice <= lattice ) {
      if ( (cur_lattice & lattice) == cur_lattice ) {
        if ( useTree_ ) {
          if ( FilterPoint( right_tree.point, left_tree.children[c] ) ) {
            dominated_[right_tree.point.pid] = true;
            return true;
          }
        } else {
          if ( FilterPoint_without_skytree( right_tree.point,
              left_tree.children[c] ) ) {
            dominated_[right_tree.point.pid] = true;
            return true;
          }
        }
      }
    } else
      break;
  }
  return false;
}

void SkyTree::PartialDominance( const uint32_t lattice, vector<TUPLE>& dataset,
    Node& skytree ) {
  const uint32_t num_child = skytree.children.size();

  for (uint32_t c = 0; c < num_child; c++) {
    uint32_t cur_lattice = skytree.children[c].lattice;
    if ( cur_lattice <= lattice ) {
      if ( (cur_lattice & lattice) == cur_lattice ) {
        // For each point, check whether the point is dominated by the existing skyline points.
        vector<TUPLE>::iterator it = dataset.begin();
        while ( it != dataset.end() ) {
          if ( useTree_ ) {
            if ( FilterPoint( *it, skytree.children[c] ) ) {
              *it = dataset.back();
              dataset.pop_back();
            } else
              ++it;
          } else {
            if ( FilterPoint_without_skytree( *it, skytree.children[c] ) ) {
              *it = dataset.back();
              dataset.pop_back();
            } else
              ++it;
          }
        }

        if ( dataset.empty() )
          break;
      }
    } else
      break;
  }
}

bool SkyTree::FilterPoint_without_skytree( const TUPLE &cur_value,
    Node& skytree ) {
  const uint32_t lattice = DT_bitmap_dvc( cur_value, skytree.point );
  const uint32_t pruned = SHIFTS[NUM_DIMS] - 1;

  if ( lattice < pruned ) {
    assert( !DominateLeft(skytree.point, cur_value) );
    if ( skytree.children.size() > 0 ) {
      const uint32_t num_child = skytree.children.size();
      for (uint32_t c = 0; c < num_child; c++) {
        if ( FilterPoint( cur_value, skytree.children[c] ) )
          return true;
      }
    }
    assert( !DominateLeft(skytree.point, cur_value) );
    return false;
  }
  assert( DominateLeft(skytree.point, cur_value) );
  return true;
}

bool SkyTree::FilterPoint( const TUPLE &cur_value, Node& skytree ) {
  const uint32_t lattice = DT_bitmap_dvc( cur_value, skytree.point );
  const uint32_t pruned = SHIFTS[NUM_DIMS] - 1;

  if ( lattice < pruned ) {
    assert( !DominateLeft(skytree.point, cur_value) );
    if ( skytree.children.size() > 0 ) {
      const uint32_t num_child = skytree.children.size();
      for (uint32_t c = 0; c < num_child; c++) {
        uint32_t cur_lattice = skytree.children[c].lattice;
        if ( cur_lattice <= lattice ) {
          if ( (cur_lattice & lattice) == cur_lattice ) {
            if ( FilterPoint( cur_value, skytree.children[c] ) )
              return true;
          }
        } else
          break;
      }
    }
    assert( !DominateLeft(skytree.point, cur_value) );
    return false;
  }
  assert( DominateLeft(skytree.point, cur_value) );
  return true;
}

void SkyTree::TraverseSkyTree( const Node& skytree ) {
  if ( !useDnC_ || !dominated_[skytree.point.pid] )
    skyline_.push_back( skytree.point.pid );

  uint32_t num_child = skytree.children.size();
  for (uint32_t c = 0; c < num_child; c++)
    TraverseSkyTree( skytree.children[c] );
}

#ifndef NVERBOSE
int SkyTree::MaxDepth( const Node& skytree, int d ) {
  skytree_levels_[d]++;

  if ( skytree.children.size() == 0 ) {
    return 1;
  } else {
    int depth = MaxDepth( skytree.children[0], d + 1 );
    for (uint32_t c = 1; c < skytree.children.size(); ++c) {
      int h = MaxDepth( skytree.children[c], d + 1 );
      if ( h > depth )
        depth = h;
    }

    return depth + 1;
  }
}
#endif
