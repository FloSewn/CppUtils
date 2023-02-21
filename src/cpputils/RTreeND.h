/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <memory>    // std::unique_ptr
#include <array>     // std::array
#include <list>      // std::list
#include <vector>    // std::vector
#include <stdexcept> // std::runtime_error
#include <iomanip>   // std::setprecision
#include <iostream>  // std::to_string
#include <algorithm> // std::sort
#include <numeric>   // std::iota

#include "VecND.h"
#include "BBoxND.h"
#include "MathUtility.h"
#include "Helpers.h"
#include "VtkIO.h"

namespace CppUtils {


// ObjectType.... Contained object
// M............. Max. element number 
// CoordType..... Coordinate type
// Dim........... Dimensions
// SortStrategy.. Sorting strategy
// SplitStrategy. Splitting strategy
  
#ifndef RTREE_NODE_DEF
#define RTREE_NODE_DEF       \
    typename    ObjectType,  \
    std::size_t M,           \
    typename    CoordType,   \
    std::size_t Dim,         \
    typename    SortStrategy  
#endif

#ifndef RTREE_TREE_DEF
#define RTREE_TREE_DEF        \
    typename    ObjectType,   \
    std::size_t M,            \
    typename    CoordType,    \
    std::size_t Dim,          \
    typename    SortStrategy, \
    typename    SplitStrategy  
#endif

#ifndef RTREE_NODE_ARG
#define RTREE_NODE_ARG \
  ObjectType, M, CoordType, Dim, SortStrategy
#endif

#ifndef RTREE_TREE_ARG
#define RTREE_TREE_ARG \
  ObjectType, M, CoordType, Dim, SortStrategy, SplitStrategy
#endif


/*********************************************************************
* Forward declarations 
*********************************************************************/
template<RTREE_NODE_DEF>
class RTreeNodeND;

template<RTREE_TREE_DEF>
class RTreeND;


/*********************************************************************
* This class defines a strategy for the sorting of RTree nodes
* within a layer
*********************************************************************/
class NearestXSort
{
public:

  /*------------------------------------------------------------------ 
  | This function takes a vector of objects and a corresponding 
  | vector of the object's bounding boxes and creates an array 
  | of sorted pointers to the given objects, where the sorting is 
  | based on the "choose_bbox()" function
  ------------------------------------------------------------------*/
  template <typename CoordType, std::size_t Dim>
  static inline std::vector<size_t> 
  sort(const std::vector<BBoxND<CoordType,Dim>>& bboxes)
  {
    // Create vector with indices for the object sorting
    std::vector<size_t> index( bboxes.size() );
    std::iota(index.begin(), index.end(), 0);

    // Use std::stable_sort() instead of sort to avoid unnecessary
    // index re-orderings when comparing equal values
    std::stable_sort(index.begin(), index.end(),
      [&bboxes](size_t i1, size_t i2)
    {
      return choose_bbox(bboxes[i1], bboxes[i2]);
    });

    return std::move( index );

  } // sort()

  /*------------------------------------------------------------------ 
  | This function returns true, if the given left-hand sided 
  | bounding box should be chosen, based on the nearest-X sort 
  ------------------------------------------------------------------*/
  template <typename CoordType, std::size_t Dim>
  static inline
  bool choose_bbox(const BBoxND<CoordType,Dim>& lhs, 
                   const BBoxND<CoordType,Dim>& rhs)
  {
    VecND<CoordType,Dim> xy_l = (lhs.lowleft() + lhs.upright());
    VecND<CoordType,Dim> xy_r = (rhs.lowleft() + rhs.upright());
    VecND<CoordType,Dim> delta = xy_l - xy_r;

    for (auto v : delta) 
    {
      if (v < CoordType{})
        return true;

      if ( !EQ0(v) )
        break;
    }

    return false;

  } // choose_bbox()

}; // NearestXSort



/*********************************************************************
* Entry to store the RTree data 
*********************************************************************/
template<RTREE_NODE_DEF>
class RTreeEntryND
{
public:
  using BBox     = BBoxND<CoordType,Dim>;
  using Node     = RTreeNodeND<ObjectType,M,CoordType,Dim,SortStrategy>;
  using Node_ptr = std::unique_ptr<Node>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  RTreeEntryND() { }

  /*------------------------------------------------------------------ 
  | Setter
  ------------------------------------------------------------------*/
  void bbox(const BBox& b) { bbox_ = b; }
  void child(Node_ptr& c) { child_ = std::move(c); }
  void object(const ObjectType* obj) { object_ = obj; }

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  BBox& bbox() { return bbox_; }
  const BBox& bbox() const { return bbox_; } 

  Node_ptr& child_ptr() { return child_; }
  const Node_ptr& child_ptr() const { return child_; }

  Node& child() { return *child_; }
  const Node& child() const { return *child_; }

  ObjectType* object() { return const_cast<ObjectType*>(object_); }
  const ObjectType* object() const { return object_; }

private:
  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  BBox               bbox_   {};
  Node_ptr           child_  { nullptr };
  const ObjectType*  object_ { nullptr };

}; // RTreeEntryND



/*********************************************************************
* Quadratic split strategy
*********************************************************************/
class QuadraticSplit
{
public:

  template<RTREE_NODE_DEF>
  struct SplitData
  {
    using Entry      = RTreeEntryND<RTREE_NODE_ARG>;
    using EntryArray = std::array<Entry,M>;
    using BBox       = BBoxND<CoordType,Dim>;

    SplitData(EntryArray& e) : entries { e } {}

    EntryArray&        entries;

    EntryArray         left  {};
    EntryArray         right {};

    BBox               bbox_left  {};
    BBox               bbox_right {};

    std::size_t        n_left  {};
    std::size_t        n_right {};

    std::array<bool,M> distributed { false };
    std::size_t        n_remaining { M };

    /*---------------------------------------------------------------- 
    | Add entry to the left side
    ----------------------------------------------------------------*/
    void add_to_left(std::size_t i)
    {
      ASSERT(!distributed[i],
      "QuadraticSplit::SplitData::add_to_left(): "
      "Can not redistribute entry " + std::to_string(i) + ".");

      ASSERT(n_remaining > 0,
      "QuadraticSplit::SplitData::add_to_left(): "
      "All entries are already distributed.");

      ASSERT(n_left < M-1,
      "QuadraticSplit::SplitData::add_to_left(): "
      "Array of left entries already full.");

      if ( n_left == 0 )
        bbox_left = entries[i].bbox();
      else
        bbox_left = bbox_left.bbox_cover( entries[i].bbox() );

      left[n_left] = std::move(entries[i]);
      ++n_left;
      distributed[i] = true;
      --n_remaining;

    } // add_to_left()

    /*---------------------------------------------------------------- 
    | Add entry to the right side
    ----------------------------------------------------------------*/
    void add_to_right(std::size_t i)
    {
      ASSERT(!distributed[i],
      "QuadraticSplit::SplitData::add_to_right(): "
      "Can not redistribute entry " + std::to_string(i) + ".");

      ASSERT(n_remaining > 0,
      "QuadraticSplit::SplitData::add_to_right(): "
      "All entries are already distributed.");

      ASSERT(n_right < M-1,
      "QuadraticSplit::SplitData::add_to_right(): "
      "Array of right entries already full.");

      if ( n_right == 0 )
        bbox_right = entries[i].bbox();
      else
        bbox_right = bbox_right.bbox_cover( entries[i].bbox() );

      right[n_right] = std::move(entries[i]);
      ++n_right;
      distributed[i] = true;
      --n_remaining;

    } // add_to_right()


    /*---------------------------------------------------------------- 
    | Put all remaining entries to the right side, so that it has at
    | least M/2 entries
    ----------------------------------------------------------------*/
    bool balance_left()
    {
      if ( n_right < M/2 )
        return false;

      if ( n_left + n_remaining != M/2 )
        return false;

      for (std::size_t i = 0; i < M; ++i)
        if ( !distributed[i] )
          add_to_left(i);

      return true;
    }
     
    /*---------------------------------------------------------------- 
    | Put all remaining entries to the right side, so that it has at
    | least M/2 entries
    ----------------------------------------------------------------*/
    bool balance_right()
    {
      if ( n_left < M/2 )
        return false;

      if ( n_right + n_remaining != M/2 )
        return false;

      for (std::size_t i = 0; i < M; ++i)
        if ( !distributed[i] )
          add_to_right(i);

      return true;
    }


    /*---------------------------------------------------------------- 
    | Compute the amount by which the covering rectangle of the 
    | left group has to be enlarged if entry "i" would be added
    ----------------------------------------------------------------*/
    CoordType rect_diff_left(std::size_t i) const 
    {
      ASSERT(!distributed[i],
      "QuadraticSplit::SplitData::rect_diff_left(): "
      "Entry " + std::to_string(i) + " already distributed.");
      const auto c_left = bbox_left.bbox_cover( entries[i].bbox() );
      return c_left.scale() - bbox_left.scale();
    }

    /*---------------------------------------------------------------- 
    | Compute the amount by which the covering rectangle of the 
    | right group has to be enlarged if entry "i" would be added
    ----------------------------------------------------------------*/
    CoordType rect_diff_right(std::size_t i) const
    {
      ASSERT(!distributed[i],
      "QuadraticSplit::SplitData::rect_diff_right(): "
      "Entry " + std::to_string(i) + " already distributed.");
      const auto c_right = bbox_right.bbox_cover( entries[i].bbox() );
      return c_right.scale() - bbox_right.scale();
    }

  }; // SplitData

  /*------------------------------------------------------------------ 
  | Split a given array of RTree entries
  ------------------------------------------------------------------*/
  template<RTREE_NODE_DEF>
  static inline std::unique_ptr<RTreeNodeND<RTREE_NODE_ARG>>
  split(RTreeNodeND<RTREE_NODE_ARG>& node, std::size_t new_node_id)
  {
    ASSERT( node.n_entries() == M, 
    "QuadraticSplit::split(): "
    "Can not split node that is not full.");

    SplitData<RTREE_NODE_ARG> split_data { node.entries() };

    // Pick first entry for each group
    pick_seeds( split_data );

    while ( split_data.n_remaining > 0 )
    {
      if ( split_data.balance_left() )
        break;

      if ( split_data.balance_right() )
        break;

      pick_next(split_data);
    }

    ASSERT( split_data.n_remaining == 0,
    "QuadraticSplit::split(): "
    "Entries to split not fully distributed.");

    ASSERT( split_data.n_left >= M/2,
    "QuadraticSplit::split(): "
    "Left group has only " + std::to_string(split_data.n_left) + 
    " entries.");

    ASSERT( split_data.n_right >= M/2,
    "QuadraticSplit::split(): "
    "Right group has only " + std::to_string(split_data.n_right) + 
    " entries.");

    return create_new_node( node, split_data, new_node_id );

  } // QuadraticSplit::split()


  /*------------------------------------------------------------------ 
  | Create a new node
  ------------------------------------------------------------------*/
  template<RTREE_NODE_DEF>
  static inline std::unique_ptr<RTreeNodeND<RTREE_NODE_ARG>>
  create_new_node(RTreeNodeND<RTREE_NODE_ARG>& node,
                  SplitData<RTREE_NODE_ARG>& split_data,
                  std::size_t new_node_id)
  {
    auto new_node_ptr 
      = std::make_unique<RTreeNodeND<RTREE_NODE_ARG>>(new_node_id);
    auto& new_node = *new_node_ptr;

    new_node.is_leaf( node.is_leaf() );
    node.n_entries( 0 );

    if ( node.is_leaf() )
    {
      for ( std::size_t i = 0; i < split_data.n_left; ++i )
      {
        auto& e = split_data.left[i];
        node.add_object( *e.object(), e.bbox() );
      }

      for ( std::size_t i = 0; i < split_data.n_right; ++i )
      {
        auto& e = split_data.right[i];
        new_node.add_object( *e.object(), e.bbox() );
      }
    }
    else
    {
      for ( std::size_t i = 0; i < split_data.n_left; ++i )
      {
        auto& e = split_data.left[i];
        node.add_child( e.child_ptr() );
      }

      for ( std::size_t i = 0; i < split_data.n_right; ++i )
      {
        auto& e = split_data.right[i];
        new_node.add_child( e.child_ptr() );
      }
    }

    return std::move(new_node_ptr);

  } // create_new_node()

  /*------------------------------------------------------------------ 
  | Compute inefficiency between two bounding boxes
  ------------------------------------------------------------------*/
  template<typename CoordType, std::size_t Dim>
  static inline CoordType 
  inefficiency(const BBoxND<CoordType,Dim>& b1,
               const BBoxND<CoordType,Dim>& b2)
  { return b1.bbox_cover(b2).scale() - b1.scale() - b2.scale(); } 


  /*------------------------------------------------------------------ 
  | Pick the two seeds from a given array of rtree entries that 
  | should be split. 
  ------------------------------------------------------------------*/
  template<RTREE_NODE_DEF>
  static inline void pick_seeds(SplitData<RTREE_NODE_ARG>& split_data)
  {
    std::size_t seed_left  = 0;
    std::size_t seed_right = 1;

    CoordType ineff_0 = CoordType{};

    // Calculate inefficiency of grouping entries together
    for (std::size_t i = 0; i < M; ++i)
    {
      const BBoxND<CoordType,Dim>& b_i = split_data.entries[i].bbox();

      for (std::size_t j = 0; j < M; ++j)
      {
        if ( i == j )
          continue;

        const BBoxND<CoordType,Dim>& b_j = split_data.entries[j].bbox();

        const CoordType ineff = inefficiency(b_i, b_j);

        // Choose the most wasteful pair
        if ( ineff > ineff_0 )
        {
          ineff_0    = ineff;
          seed_left  = i;
          seed_right = j;
        }
      }
    }

    ASSERT( seed_left != seed_right, 
    "QuadraticSplit::pick_seeds(): Invalid seeds for splitting.");

    // Distribute entries to split_data
    split_data.add_to_left( seed_left );
    split_data.add_to_right( seed_right );

  } // pick_seeds()

  /*------------------------------------------------------------------ 
  | Pick the next entries from a given array of rtree entries that 
  | should be split. 
  ------------------------------------------------------------------*/
  template<RTREE_NODE_DEF>
  static inline void pick_next(SplitData<RTREE_NODE_ARG>& split_data)
  {
    std::size_t i_picked = M;
    CoordType   max_diff = CoordType{};

    // Determine the cost of putting each entry in each group
    for (std::size_t i = 0; i < M; ++i)
    {
      if ( split_data.distributed[i] )
        continue;

      // Compute the size increases required in the covering rectangles
      // of the left and right group to include entry "i"
      const CoordType d_left  = split_data.rect_diff_left(i);
      const CoordType d_right = split_data.rect_diff_right(i);
      const CoordType diff    = ABS(d_left - d_right);

      // Choose the entry with the maximum difference 
      if ( diff >= max_diff )
      {
        i_picked = i;
        max_diff = diff;
      }
    }

    ASSERT( i_picked < M,
    "QuadraticSplit::SplitData::pick_next(): "
    "Did not find valid entry. "
    "Number of remaining entries: " 
    + std::to_string(split_data.n_remaining));

    ASSERT( !split_data.distributed[i_picked],
    "QuadraticSplit::SplitData::pick_next(): "
    "Entry " + std::to_string(i_picked) + " already distributed. ");

    // Add the picked entry to the group whose covering rectangle will
    // have to be enlarged least to accomodate it.
    const CoordType d_left  = split_data.rect_diff_left(i_picked);
    const CoordType d_right = split_data.rect_diff_right(i_picked);
    const CoordType diff    = d_left - d_right;
    bool tie                = EQ0(diff);

    if ( !tie && d_left < d_right )
    {
      split_data.add_to_left( i_picked );
      return;
    }

    if ( !tie && d_left > d_right )
    {
      split_data.add_to_right( i_picked );
      return;
    }

    // Tie: Add entry to set with smaller size
    if ( split_data.bbox_left.scale() < split_data.bbox_right.scale() )
    {
      split_data.add_to_left( i_picked );
      return;
    }

    if ( split_data.bbox_left.scale() > split_data.bbox_right.scale() )
    {
      split_data.add_to_right( i_picked );
      return;
    }

    // Further tie: Add entry to set with less entries
    if ( split_data.n_left < split_data.n_right )
    {
      split_data.add_to_left( i_picked );
      return;
    }

    split_data.add_to_right( i_picked );


  } // pick_next()


}; // QuadraticSplit



/*********************************************************************
* Some definitions on R-Trees:
* ----------------------------
* - Each leaf node (unless it is root) can host up to M entries
* - The minimum allowed number of entries is m = M/2
* - Each leaf node entry is of the form (bbox, oid) - where bbox is 
*   the minimum bounding rectangle that contains the object and oid
*   is the object id
* - Each internal node can store between m = M/2 and M entries
* - Each internal node entry is of the form (bbox, p), where p is a 
*   pointer to a child of the node and bbox is the minimum bounding 
*   rectangle that spatially contains the bboxs that are contained in 
*   the child
*
* References
* ----------
* - https://tildesites.bowdoin.edu/~ltoma/teaching/cs340/spring08/\
*   Papers/Rtree-chap1.pdf
*
*********************************************************************/
template 
<
  typename    ObjectType,                      // Contained object
  std::size_t M,                               // Max. element number 
  typename    CoordType,                       // Coordinate type
  std::size_t Dim,                             // Dimensions
  typename    SortStrategy = NearestXSort      // Sorting strategy
>
class RTreeNodeND
{
public:
  using BBox     = BBoxND<CoordType,Dim>;
  using Node     = RTreeNodeND<RTREE_NODE_ARG>;
  using Node_ptr = std::unique_ptr<Node>;
  using Entry    = RTreeEntryND<RTREE_NODE_ARG>;
  using Entries  = std::array<Entry,M>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  RTreeNodeND(std::size_t id) : id_ { id } {}

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  Entries& entries() { return entries_; }
  const Entries& entries() const { return entries_; }

  Node* parent() { return parent_; }
  const Node* parent() const { return parent_; }

  bool is_leaf() const { return is_leaf_; }
  std::size_t n_entries() const { return n_entries_; }
  std::size_t id() const { return id_; }
  std::size_t index() const { return index_; }

  /*------------------------------------------------------------------ 
  | Return the node, which is located in the same layer one place
  | to the left
  ------------------------------------------------------------------*/
  Node* left() const
  { 
    Node* parent_node = parent_;

    if ( !parent_node )
      return nullptr;

    std::size_t i = index_;

    if ( i > 0 ) 
      return (*parent_node).child_ptr(--i).get();

    parent_node = (*parent_node).left();

    if ( parent_node )
    {
      i = (*parent_node).n_entries();
      return (*parent_node).child_ptr(--i).get();
    }

    return nullptr; 
  }

  /*------------------------------------------------------------------ 
  | Return the node, which is located in the same layer one place
  | to the left
  ------------------------------------------------------------------*/
  Node* right() const
  { 
    Node* parent_node = parent_;

    if ( !parent_node )
      return nullptr;

    std::size_t i = index_;

    if ( i+2 <= (*parent_node).n_entries() ) 
      return (*parent_node).child_ptr(++i).get();

    parent_node = (*parent_node).right();

    if ( parent_node )
      return (*parent_node).child_ptr(0).get();

    return nullptr; 
  }

  /*------------------------------------------------------------------ 
  | Setter
  ------------------------------------------------------------------*/
  void parent(Node* p) { parent_ = p; }
  void is_leaf(bool t) { is_leaf_ = t; }
  void n_entries(std::size_t n) { n_entries_ = n; }
  void id(std::size_t i) { id_ = i; }
  void index(std::size_t index) { index_ = index; }

  /*------------------------------------------------------------------ 
  | Collect all entries of this node or all its following child nodes
  ------------------------------------------------------------------*/
  void collect_entry_bboxes(std::vector<const ObjectType*>& objects,
                            std::vector<BBox>& bboxes) const
  {
    if ( !is_leaf() )
    {
      for (std::size_t i = 0; i < n_entries(); ++i )
        child(i).collect_entry_bboxes(objects, bboxes);
    }
    else
    {
      for (std::size_t i = 0; i < n_entries(); ++i )
      {
        objects.push_back( &object(i) );
        bboxes.push_back( bbox(i) );
      }
    }

  } // collect_entry_bboxes()

  /*------------------------------------------------------------------ 
  | Function used to estimate the tree height
  ------------------------------------------------------------------*/
  std::size_t increment_tree_height(std::size_t height) const
  {
    if ( !is_leaf() )
    {
      ++height;
      height = child(0).increment_tree_height(height);
    }

    return height;
  }

  /*------------------------------------------------------------------ 
  | Get bounding box that encloses all objects stored in this node
  ------------------------------------------------------------------*/
  BBox bbox() const
  {
    if ( this->n_entries() < 1 )
      return {};

    BBox bbox = this->bbox(0);

    for ( size_t j = 1; j < this->n_entries(); ++j )
      bbox = bbox.bbox_cover(this->bbox(j));

    return bbox;
  };

  /*------------------------------------------------------------------ 
  | Access entries
  ------------------------------------------------------------------*/
  const Entry& entry(std::size_t i) const
  {
    ASSERT( i < n_entries(), 
    "RTreeNodeND::entry(): "
    "Unable to access entry[" + std::to_string(i) + "].");

    return entries_[i]; 
  }

  Entry& entry(std::size_t i)
  {
    ASSERT( i < n_entries(), 
    "RTreeNodeND::entry(): "
    "Unable to access entry[" + std::to_string(i) + "].");

    return entries_[i]; 
  }

  /*------------------------------------------------------------------ 
  | Access minimum bounding rectangles
  ------------------------------------------------------------------*/
  const BBox& bbox(std::size_t i) const
  { 
    ASSERT( i < n_entries(), 
    "RTreeNodeND::bbox(): "
    "Unable to access bbox[" + std::to_string(i) + "].");

    return entries_[i].bbox(); 
  }

  /*------------------------------------------------------------------ 
  | Access objects
  ------------------------------------------------------------------*/
  ObjectType& object(std::size_t i) 
  { 
    ASSERT( is_leaf(),
    "RTreeNodeND::object(): "
    "Unable to access object of non-leaf node.");

    ASSERT( i < n_entries(), 
    "RTreeNodeND::object(): " 
    "Unable to access object[" + std::to_string(i) + "].");

    return *(const_cast<ObjectType*>(entries_[i].object())); 
  }

  const ObjectType& object(std::size_t i) const 
  { 
    ASSERT( is_leaf(),
    "RTreeNodeND::object(): "
    "Unable to access object of non-leaf node.");

    ASSERT( i < n_entries(), 
    "RTreeNodeND::object(): " 
    "Unable to access object[" + std::to_string(i) + "].");

    return *(entries_[i].object());
  }

  /*------------------------------------------------------------------ 
  | Access children pointers
  ------------------------------------------------------------------*/
  Node_ptr& child_ptr(std::size_t i) 
  { 
    ASSERT( !is_leaf(),
    "RTreeNodeND::child_ptr(): "
    "Unable to access child of leaf node.");

    ASSERT( i < n_entries(), 
    "RTreeNodeND::child_ptr(): " 
    "Unable to access child[" + std::to_string(i) + "].");

    return entries_[i].child_ptr();
  }

  const Node_ptr& child_ptr(std::size_t i) const
  { 
    ASSERT( !is_leaf(),
    "RTreeNodeND::child_ptr(): "
    "Unable to access child of leaf node.");

    ASSERT( i < n_entries(), 
    "RTreeNodeND::child_ptr(): " 
    "Unable to access child[" + std::to_string(i) + "].");

    return entries_[i].child_ptr(); 
  }

  /*------------------------------------------------------------------ 
  | Access children as references
  ------------------------------------------------------------------*/
  Node& child(std::size_t i) 
  { 
    ASSERT( !is_leaf(),
    "RTreeNodeND::child(): "
    "Unable to access child of leaf node.");

    ASSERT( i < n_entries(), 
    "RTreeNodeND::child(): " 
    "Unable to access child[" + std::to_string(i) + "].");

    return entries_[i].child();
  }

  const Node& child(std::size_t i) const 
  { 
    ASSERT( !is_leaf(),
    "RTreeNodeND::child(): "
    "Unable to access child of leaf node.");

    ASSERT( i < n_entries(), 
    "RTreeNodeND::child(): " 
    "Unable to access child[" + std::to_string(i) + "].");

    return entries_[i].child();
  }

  /*------------------------------------------------------------------ 
  | Set bounding boxes of children / objects
  ------------------------------------------------------------------*/
  void bbox(std::size_t i, const BBox& b) 
  { 
    ASSERT( i < n_entries(), 
    "RTreeNodeND::bbox(): "
    "Unable to set bbox[" + std::to_string(i) + "]." );

    entries_[i].bbox(b);
  }

  /*------------------------------------------------------------------ 
  | Set objects
  ------------------------------------------------------------------*/
  void object(std::size_t i, const ObjectType& obj)
  { 
    ASSERT( is_leaf(),
    "RTreeNodeND::object(): "
    "Unable to set object of non-leaf node.");

    ASSERT( i < n_entries(), 
    "RTreeNodeND::object(): " 
    "Unable to set object[" + std::to_string(i) + "].");

    entries_[i].object(&obj);
  }

  /*------------------------------------------------------------------ 
  | Set children
  ------------------------------------------------------------------*/
  void child(std::size_t i, std::unique_ptr<Node>& c)
  { 
    ASSERT( !is_leaf(),
    "RTreeNodeND::child(): "
    "Unable to set child of leaf node.");

    ASSERT( i < n_entries(), 
    "RTreeNodeND::child(): " 
    "Unable to set child[" + std::to_string(i) + "].");

    entries_[i].child(c);
  }

  /*------------------------------------------------------------------ 
  | Add new object to the node
  ------------------------------------------------------------------*/
  void add_object(const ObjectType& object, const BBox& obj_bbox)
  {
    const std::size_t n_entries = this->n_entries();

    ASSERT( n_entries < M, 
    "RTreeNodeND::add_object(): "
    "Can not add more than " + std::to_string(M) + " objects.");

    ASSERT( is_leaf(), 
    "RTreeNodeND:add_object(): "
    "Can not add object to non-leaf node.");

    std::size_t i = n_entries;

    // Obtain index i at which all former object entries are located
    // closer to the origin. Start searching from the back (faster)
    while ( i > 0 )
    {
      if ( SortStrategy::choose_bbox(this->bbox(i-1), obj_bbox) )
        break;
      --i;
    }

    this->n_entries( n_entries + 1 );

    // Shift remaining entries to the right
    for (std::size_t j = n_entries; j > i; --j)
    {
      this->bbox(j, this->bbox(j-1));
      this->object(j, this->object(j-1));
    }

    // Finally add the new object
    this->bbox(i, obj_bbox); 
    this->object(i, object);

  } // add_object()

  /*------------------------------------------------------------------ 
  | Add new child to the node
  ------------------------------------------------------------------*/
  void add_child(Node_ptr& child_ptr)
  {
    Node& child = *child_ptr;
    const std::size_t n_entries = this->n_entries();

    ASSERT( n_entries < M, 
    "RTreeNodeND::add_child(): "
    "Can not add more than " + std::to_string(M) + " children.");

    ASSERT( !is_leaf(), 
    "RTreeNodeND:add_child(): "
    "Can not add child to leaf node.");

    std::size_t i = n_entries;

    // Obtain index i at which all former object entries are located
    // closer to the origin. Start searching from the back (faster)
    while ( i > 0 )
    {
      if ( SortStrategy::choose_bbox(this->bbox(i-1), child.bbox()) )
        break;
      --i;
    }

    this->n_entries( n_entries + 1 );

    // Shift remaining entries to the right
    for (std::size_t j = n_entries; j > i; --j)
    {
      this->bbox(j, this->bbox(j-1));
      this->child(j, this->child_ptr(j-1));
      this->child(j).index(j);
    }

    // Finally add the new child
    this->bbox(i, child.bbox()); 
    this->child(i, child_ptr);
    child.parent(this);
    child.index(i);

  } // add_child()

  /*------------------------------------------------------------------ 
  | Remove the object at a given index 
  ------------------------------------------------------------------*/
  void remove_object(std::size_t index)
  {
    const std::size_t n_entries = this->n_entries();

    ASSERT( n_entries > 0, 
    "RTreeNodeND::remove_object(): "
    "Can not remove object from empty node.");

    ASSERT( is_leaf(), 
    "RTreeNodeND::remove_object(): "
    "Can not remove object from non-leaf node.");

    // Shift remaining entries to the left
    for (std::size_t j = index; j < n_entries-1; ++j)
    {
      this->bbox(j, this->bbox(j+1));
      this->object(j, this->object(j+1));
    }

    this->n_entries( n_entries - 1 );

  } // remove_object()

  /*------------------------------------------------------------------ 
  | Remove the child at a given index 
  ------------------------------------------------------------------*/
  void remove_child(std::size_t index)
  {
    const std::size_t n_entries = this->n_entries();
     
    ASSERT( n_entries > 0, 
    "RTreeNodeND::remove_child(): "
    "Can not remove child from empty node.");

    ASSERT( !is_leaf(), 
    "RTreeNodeND::remove_child(): "
    "Can not remove child from leaf node.");

    // Shift remaining entries to the left
    for (std::size_t j = index; j < n_entries-1; ++j)
    {
      this->bbox(j, this->bbox(j+1));
      this->child(j, this->child_ptr(j+1));

      this->child(j).index(j);
    }

    this->n_entries( n_entries - 1 );

  } // remove_child()

private:
  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  Entries             entries_   { };
  Node*               parent_    { nullptr };

  bool                is_leaf_   { true };
  std::size_t         n_entries_ { 0 };
  std::size_t         id_        { 0 };
  std::size_t         index_     { 0 };

}; // RTreeNodeND


/*********************************************************************
* This class defines the interface to an R-tree structure 
* for N-dimensional simplices
*
* References
* ----------
* - Custom iterator implementation: https://www.internalpointers.\
*   com/post/writing-custom-iterators-modern-cpp
*********************************************************************/
template 
<
  typename    ObjectType,                      // Contained object
  std::size_t M,                               // Max. element number 
  typename    CoordType,                       // Coordinate type
  std::size_t Dim,                             // Dimensions
  typename    SortStrategy = NearestXSort,     // Sorting strategy
  typename    SplitStrategy = QuadraticSplit   // Splitting strategy
>
class RTreeND
{
  static inline long unsigned node_id_ = 0;

public:

  using BBox       = BBoxND<CoordType,Dim>;
  using Node       = RTreeNodeND<RTREE_NODE_ARG>;
  using NodeVector = std::vector<std::unique_ptr<Node>>;
  using Entry      = RTreeEntryND<RTREE_NODE_ARG>;
  using Entries    = std::array<Entry,M>;

  /*------------------------------------------------------------------ 
  | 
  ------------------------------------------------------------------*/
  template <typename QueryObject = VecND<CoordType,Dim>>
  using ObjectDistFunction 
    = CoordType (*)(const QueryObject&, const ObjectType&);

  template <typename QueryObject = VecND<CoordType,Dim>>
  using BBoxDistFunction
    = CoordType (*)(const QueryObject&, const BBox&);

  /*------------------------------------------------------------------ 
  | Iterator implementation
  ------------------------------------------------------------------*/
  struct Iterator
  {
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = ObjectType;
    using pointer           = ObjectType*;
    using reference         = ObjectType&;

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Constructor
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    Iterator() {}
    Iterator(Node& node) : cur_node_ {&node} 
    { 
      if ( (*cur_node_).n_entries() == 0 )
        cur_node_ = nullptr;
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Operators
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    reference operator*() const 
    { return cur_node_->object(cur_index_); }

    pointer operator->() 
    { return &(cur_node_->object(cur_index_)); }

    Iterator& operator++() 
    { 
      ++cur_index_;

      if ( cur_index_ == (*cur_node_).n_entries() )
      {
        cur_node_  = cur_node_->right();
        cur_index_ = 0;
      }

      return *this;
    }

    Iterator operator++(int) 
    { Iterator tmp = *this; ++(*this); return tmp; }

    friend bool operator== (const Iterator& a, const Iterator& b) 
    { return (   a.cur_node_ == b.cur_node_ 
              && a.cur_index_ == b.cur_index_); }

    friend bool operator!= (const Iterator& a, const Iterator& b) 
    { return !(a == b); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Attributes
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
  private:
    Node* cur_node_ { nullptr };
    std::size_t cur_index_ { 0 };

  }; // Iterator 

  /*------------------------------------------------------------------ 
  | ConstantIterator implementation
  ------------------------------------------------------------------*/
  struct ConstantIterator
  {
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = ObjectType;
    using pointer           = ObjectType*;
    using reference         = ObjectType&;

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Constructor
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    ConstantIterator() {}
    ConstantIterator(Node& node) : cur_node_ {&node} 
    { 
      if ( (*cur_node_).n_entries() == 0 )
        cur_node_ = nullptr;
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Operators
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    const reference operator*() const 
    { return cur_node_->object(cur_index_); }

    const pointer operator->() 
    { return &(cur_node_->object(cur_index_)); }

    ConstantIterator& operator++() 
    { 
      ++cur_index_;

      if ( cur_index_ == (*cur_node_).n_entries() )
      {
        cur_node_  = cur_node_->right();
        cur_index_ = 0;
      }

      return *this;
    }

    ConstantIterator operator++(int) 
    { ConstantIterator tmp = *this; ++(*this); return tmp; }

    friend bool operator== 
    (const ConstantIterator& a, const ConstantIterator& b) 
    { return (   a.cur_node_ == b.cur_node_ 
              && a.cur_index_ == b.cur_index_); }

    friend bool operator!= 
    (const ConstantIterator& a, const ConstantIterator& b) 
    { return !(a == b); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Attributes
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
  private:
    Node* cur_node_ { nullptr };
    std::size_t cur_index_ { 0 };

  }; // ConstantIterator 

  Iterator begin() { return Iterator( leaf_leftmost() ); }
  Iterator end() { return Iterator(); }

  ConstantIterator cbegin() const 
  { return ConstantIterator( leaf_leftmost() ); }
  ConstantIterator cend() const 
  { return ConstantIterator(); }


  /*------------------------------------------------------------------ 
  | A container for the nearest neighbor query 
  ------------------------------------------------------------------*/
  struct NearestNeighbor
  {
    const ObjectType* object_ptr {nullptr};
    CoordType         dist_sqr   {std::numeric_limits<CoordType>::max()};
  }; 

  /*------------------------------------------------------------------ 
  | A sub-class that is returned after the Query for k nearst 
  | neighbors of the tree
  ------------------------------------------------------------------*/
  class KNearestList
  {
  public:
    using List           = std::list<const ObjectType*>;
    using iterator       = typename List::iterator;
    using const_iterator = typename List::const_iterator;

    iterator begin() { return nbrs_.begin(); }
    iterator end() { return nbrs_.end(); }

    const_iterator cbegin() const { return nbrs_.cbegin(); }
    const_iterator cend() const { return nbrs_.cend(); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Constructor
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    KNearestList(std::size_t k) 
    : k_            { k } 
    , max_dist_sqr_ { std::numeric_limits<CoordType>::max() }
    {}

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Getters
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    std::size_t k() const { return k_; }
    CoordType max_dist_sqr() const { return max_dist_sqr_; }

    List& values() { return nbrs_; }
    const List& values() const { return nbrs_; }

    std::list<CoordType>& squared_distances() 
    { return dists_sqr_; }
    const std::list<CoordType>& squared_distances() const 
    { return dists_sqr_; }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Insert an element to the list
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
    void insert(const ObjectType& obj, CoordType obj_dist)
    {
      // Find place to add new entry into list 
      // and eventually remove last entry
      auto dist_iter = dists_sqr_.begin();
      auto nbrs_iter = nbrs_.begin();

      for (; nbrs_iter != nbrs_.end(); ++nbrs_iter, ++dist_iter)
      {
        if ( obj_dist >= *dist_iter )
          continue;

        // Add new entry to list
        dists_sqr_.insert(dist_iter, obj_dist);
        nbrs_.insert(nbrs_iter, &obj);

        // Remove old entry
        if ( nbrs_.size() > k_ )
        {
          dists_sqr_.erase( --( dists_sqr_.end() ) );
          nbrs_.erase( --( nbrs_.end() ) );
          max_dist_sqr_ = *(--( dists_sqr_.end() ));
        }

        ASSERT( dists_sqr_.size() <= k_,
        "RTreeND::KNearestList::insert(): "
        "Invalid number of distance list entries.");

        ASSERT( nbrs_.size() <= k_,
        "RTreeND::KNearestList::insert(): "
        "Invalid number of neighbor list entries.");
        
        return;
      }

      // Append entry to the end 
      nbrs_.push_back( &obj );
      dists_sqr_.push_back( obj_dist );

      if ( nbrs_.size() == k_ )
        max_dist_sqr_ = obj_dist;

    } // insert()

  private:
    std::size_t                   k_            {};
    std::list<const ObjectType*>  nbrs_         {};
    std::list<CoordType>          dists_sqr_    {};
    CoordType                     max_dist_sqr_ {};

  }; // KNearestList


  /*------------------------------------------------------------------ 
  | RTreeND Constructor
  ------------------------------------------------------------------*/
  RTreeND()
  {
    root_ = std::make_unique<Node>(node_id_++);
  }

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  Node& root() { return *root_; }
  const Node& root() const { return *root_; }

  /*------------------------------------------------------------------ 
  | Get the height of the entire tree
  ------------------------------------------------------------------*/
  std::size_t height() const 
  {
    std::size_t height = 0;
    height = (*root_).increment_tree_height(height);
    return height;
  } 

  /*------------------------------------------------------------------ 
  | This function returns the leftmost leaf in the tree
  ------------------------------------------------------------------*/
  Node& leaf_leftmost() const
  {
    Node* node = root_.get();
    while ( !(*node).is_leaf() )
      node = (*node).child_ptr(0).get();
    return *node;
  }

  /*------------------------------------------------------------------ 
  | This function returns the rightmost leaf in the tree
  ------------------------------------------------------------------*/
  Node& leaf_rightmost() const
  {
    Node* node = root_.get();
    while ( !(*node).is_leaf() )
    {
      std::size_t n = (*node).n_entries();
      node = (*node).child_ptr(n-1).get();
    }
    return *node;
  }

  /*------------------------------------------------------------------ 
  | Query all objects in a given bounding box
  ------------------------------------------------------------------*
  std::vector<const ObjectType*> 
  query(const BBox& query_bbox,
        const QueryFunction& query_fun) const
  {

  } // query() */

  /*------------------------------------------------------------------ 
  | Query the k nearest entries to a given input position
  | The input function "sqr_dist_fun" must return the SQUARED distance
  | between the query point and a given object in the tree.
  |
  | Reference:
  | ----------
  | - Roussopoulos, Nick, Stephen Kelley, and Frederic Vincent: 
  |   Nearest neighbor queries, Proceedings of the 1995 ACM SIGMOD 
  |   international conference on Management of data, 1995 
  | 
  | - Hjaltason and Samet: Distance browsing in spatial database,
  |   ACM Transactions on Database Systems (TODS) 24.2 (1999): 265-318
  |
  ------------------------------------------------------------------*/
  template <typename QueryObject = VecND<CoordType,Dim>>
  KNearestList nearest(const QueryObject& query, 
                       std::size_t k,
                       const ObjectDistFunction<QueryObject>& object_dist_fun, 
                       const BBoxDistFunction<QueryObject>& bbox_dist_fun = 
                       [](const QueryObject& p, const BBox& b)
                       { return b.point_dist_sqr(p); }) const
  {
    KNearestList nearest_list { k };

    k_nearest_traversal(nearest_list, query, *root_, 
                        object_dist_fun, bbox_dist_fun);

    return std::move( nearest_list );

  } // RTreeND::nearest()

  /*------------------------------------------------------------------ 
  | Query the nearest entry to a given input position
  ------------------------------------------------------------------*/
  template <typename QueryObject = VecND<CoordType,Dim>>
  NearestNeighbor nearest(const QueryObject& query, 
                          const ObjectDistFunction<QueryObject>& object_dist_fun, 
                          const BBoxDistFunction<QueryObject>& bbox_dist_fun = 
                          [](const QueryObject& p, const BBox& b)
                          { return b.point_dist_sqr(p); }) const
  {
    NearestNeighbor nearest_neighbor {};
    
    nearest_traversal(nearest_neighbor, query, *root_,
                      object_dist_fun, bbox_dist_fun);

    return std::move(nearest_neighbor); 

  } // RTreeND::nearest()

  /*------------------------------------------------------------------ 
  | Remove a given object from the tree 
  ------------------------------------------------------------------*/
  bool remove(const ObjectType& obj, const BBox& obj_bbox)
  {
    if ( root_->n_entries() < 1 )
      return false;

    // Check if the object is contained in the R-Tree region
    if ( !root_->bbox().bbox_inside_touch( obj_bbox ) )
      return false;

    // Find the node & index of the respective entry to remove
    std::size_t i_entry {};
    Node* node_ptr = find_object_to_remove(obj, obj_bbox, 
                                           *root_, i_entry);

    // Return if no node is found 
    if ( node_ptr == nullptr )
      return false;

    // Remove the object
    Node& node = *node_ptr;
    node.remove_object( i_entry );

    // Re-balance the tree
    condense_tree( node );

    // Shorten the tree
    if ( !(*root_).is_leaf() && (*root_).n_entries() == 1 )
    {
      ASSERT( (*root_).child_ptr(0) != nullptr,
      "RTreeND::remove(): "
      "Invalid root node.");

      root_ = std::move( (*root_).child_ptr(0) );
      (*root_).parent( nullptr );
    }

    return true;

  } // RTreeND::remove(ObjectType)

  /*------------------------------------------------------------------ 
  | Insert a new object into the RTree structure
  ------------------------------------------------------------------*/
  void insert(const ObjectType& object, const BBox& obj_bbox)
  {
    // Handle the case where the root node is full
    if ( root_->n_entries() == M )
    {
      add_root_node();
      split_child(*root_, 0);
    }

    insert_nonfull(*root_, object, obj_bbox);

  } // RTreeND::insert()

  /*------------------------------------------------------------------ 
  | Insert a bulk of objects into the RTree structure
  ------------------------------------------------------------------*/
  void insert(const std::vector<ObjectType>& objects,
              const std::vector<BBox>& obj_bboxes)
  {
    // Put all objects in a temporary container, which will be 
    // sorted
    std::vector<std::size_t> sorted_indices 
      = SortStrategy::sort( obj_bboxes );

    // Group objects into (N / M) leaf nodes
    NodeVector node_layer {};
    node_layer.push_back( std::make_unique<Node>(node_id_++) );

    for ( std::size_t i = 0; i < sorted_indices.size(); ++i)
    {
      std::size_t index = sorted_indices[i];

      Node& cur_leaf = *node_layer.back().get();

      const ObjectType& cur_obj  = objects[index];
      const BBox&       cur_bbox = obj_bboxes[index];

      cur_leaf.add_object( cur_obj, cur_bbox );

      if ( cur_leaf.n_entries() >= M && i+1 < sorted_indices.size() )
      {
        node_layer.push_back( std::make_unique<Node>(node_id_++) );
        (*node_layer.back()).index( node_layer.size() - 1 );
      }
    }

    // Recursively pack nodes into a node layer at the 
    // next level until root is reached
    while ( node_layer.size() > M )
      node_layer = build_tree_bulk_insertion(node_layer);

    // Place remaining nodes into root node
    Node& root = *root_;

    for ( std::size_t i = 0; i < node_layer.size(); ++i )
    {
      Node& cur_child = *node_layer[i];

      cur_child.parent( &root );
      root.n_entries( i + 1 );
      root.is_leaf( false );
      root.bbox(i, cur_child.bbox() );
      root.child(i, node_layer[i] );
    }
      
  } // RTreeND::insert()


private:

  /*------------------------------------------------------------------ 
  | Traverse the tree to find the single nearest neighbor to a given 
  | query point in the tree
  ------------------------------------------------------------------*/
  template <typename QueryObject = VecND<CoordType,Dim>>
  void nearest_traversal(NearestNeighbor& nearest_neighbor,
                         const QueryObject& query,
                         const Node& node,
                         const ObjectDistFunction<QueryObject>& object_dist_fun,
                         const BBoxDistFunction<QueryObject>& bbox_dist_fun) const 
  {
    // Obtain the k nearest neighbors among all objects of 
    // this leaf node
    if ( node.is_leaf() )
    {
      std::vector<CoordType> dists_sqr (node.n_entries(), {});

      for (std::size_t i = 0; i < node.n_entries(); ++i)
        dists_sqr[i] = object_dist_fun(query, node.object(i));

      std::size_t winner = std::distance(
        std::begin(dists_sqr), 
        std::min_element(std::begin(dists_sqr), std::end(dists_sqr))
      );

      nearest_neighbor.object_ptr = &node.object(winner);
      nearest_neighbor.dist_sqr = dists_sqr[winner];
    }
    else
    {
      // Sort entries of current node in ascending order by their 
      // distance to query position
      std::vector<size_t> index( node.n_entries() );
      std::iota(index.begin(), index.end(), 0);

      std::stable_sort(index.begin(), index.end(),
        [&node, &query, &bbox_dist_fun](size_t i1, size_t i2)
      {
        CoordType d1 = bbox_dist_fun(query, node.bbox(i1));
        CoordType d2 = bbox_dist_fun(query, node.bbox(i2));
        return d1 < d2;
      });

      // Recursively call this function on children
      for ( std::size_t i = 0; i < node.n_entries(); ++i )
      {
        CoordType dist_sqr = bbox_dist_fun(query, node.bbox(index[i]));

        if ( dist_sqr >= nearest_neighbor.dist_sqr )
          break;
        
        nearest_traversal(nearest_neighbor, query, node.child(index[i]), 
                          object_dist_fun, bbox_dist_fun);
      }
    }

  } // nearest_traversal()

  /*------------------------------------------------------------------ 
  | Traverse the tree to find k-nearest neighbors to a given 
  | query point in the tree
  ------------------------------------------------------------------*/
  template <typename QueryObject = VecND<CoordType,Dim>>
  void k_nearest_traversal(KNearestList& nearest_list,
                           const QueryObject& query,
                           const Node& node,
                           const ObjectDistFunction<QueryObject>& object_dist_fun,
                           const BBoxDistFunction<QueryObject>& bbox_dist_fun) const 
  {
    // Obtain the k nearest neighbors among all objects of 
    // this leaf node
    if ( node.is_leaf() )
    {
      for (std::size_t i = 0; i < node.n_entries(); ++i)
      {
        CoordType dist_sqr = object_dist_fun(query, node.object(i));

        if ( dist_sqr < nearest_list.max_dist_sqr() )
          nearest_list.insert( node.object(i), dist_sqr );
      }
    }
    else
    {
      // Sort entries of current node in ascending order by their 
      // distance to query position
      std::vector<size_t> index( node.n_entries() );
      std::iota(index.begin(), index.end(), 0);

      std::stable_sort(index.begin(), index.end(),
        [&node, &query, &bbox_dist_fun](size_t i1, size_t i2)
      {
        CoordType d1 = bbox_dist_fun(query, node.bbox(i1));
        CoordType d2 = bbox_dist_fun(query, node.bbox(i2));
        return d1 < d2;
      });

      // Recursively call this function on children
      for ( std::size_t i = 0; i < node.n_entries(); ++i )
      {
        CoordType dist_sqr = bbox_dist_fun(query, node.bbox(index[i]));

        if ( dist_sqr >= nearest_list.max_dist_sqr() )
          break;
        
        k_nearest_traversal(nearest_list, query, node.child(index[i]), 
                            object_dist_fun, bbox_dist_fun);
      }
    }

  } // k_nearest_traversal()

  /*------------------------------------------------------------------ 
  | This function splits a given node. If the parent of this given
  | node is full, it alos invokes the split on the parent node.
  | Thus, it can only be called upon nodes that have parent nodes.
  ------------------------------------------------------------------*/
  void split_node(Node& node)
  {
    ASSERT( node.n_entries() == M, 
    "RTreeND::split_node(): "
    "Can not split non-full RTree child node.");

    ASSERT( node.parent(),
    "RTreeND::split_node(): "
    "Can not call function on node without parent.");

    if ( (*node.parent()).n_entries() == M )
      split_node( *node.parent() );

    split_child( *node.parent(), node.index() );

  } // split_node()

  /*------------------------------------------------------------------ 
  | This function splits the i-th child of a given "parent_node"
  ------------------------------------------------------------------*/
  void split_child(Node& parent_node, std::size_t i)
  {
    ASSERT(parent_node.n_entries() != M, 
    "RTreeND::split_child(): "
    "Can not split child of full RTree node.");

    // This is the child node, whose entries will be splitted 
    Node& child_node = parent_node.child(i);

    ASSERT(child_node.n_entries() == M, 
    "RTreeND::split_child(): "
    "Can not split non-full RTree child node.");

    // This is the new node, which will get half of the entries 
    // of "child_node"
    auto new_node_ptr = SplitStrategy::split( child_node, node_id_++ );

    // Add "new_node" and its bounding boxx to "parent_node"
    parent_node.add_child( new_node_ptr );

    // Compute the covering bboxes for all entries of "parent node"
    for (std::size_t j = 0; j < parent_node.n_entries(); ++j)
    {
      const Node& child = parent_node.child(j);

      BBox cover = child.bbox(0);

      for ( std::size_t k = 1; k < child.n_entries(); ++k )
        cover = cover.bbox_cover( child.bbox(k) );

      parent_node.bbox( j, cover );
    }

  } // RTreeND::split_child()

  /*------------------------------------------------------------------ 
  | Insert a new object into the RTree structure
  ------------------------------------------------------------------*/
  void insert_nonfull(Node&             node, 
                      const ObjectType& object,
                      const BBox&       obj_bbox)
  {
    // Choose an appropriate leaf to insert the object
    Node& leaf = choose_leaf_insertion(node, obj_bbox);

    ASSERT( leaf.n_entries() < M,
    "RTreeND::insert_nonfull(): "
    "Can not insert object in full leaf node.");

    leaf.add_object( object, obj_bbox );

    // Update all BBoxes in the path from root to this leaf, 
    // so that all of them cover the object's bbox
    update_parent_bbox(leaf);

  } // RTreeND::insert_nonfull()

  /*------------------------------------------------------------------ 
  | Choose appropriate leaf node for RTree-insertion
  ------------------------------------------------------------------*/
  Node& choose_leaf_insertion(Node& node,
                              const BBox& object_bbox)
  {
    // Choos only leaf nodes
    if ( node.is_leaf() )
      return node;


    // Compute the enlargements of all elements 
    std::vector<CoordType> enlargements ( node.n_entries() );
    std::vector<CoordType> covers ( node.n_entries() );
    std::size_t cur_id = 0;

    std::transform(node.entries().cbegin(), 
                   node.entries().cbegin() + node.n_entries(),
                   enlargements.begin(), 
    [&object_bbox, &covers, &cur_id](const Entry& cur_entry)
    {
      const BBox& cur_bbox      = cur_entry.bbox();
      const CoordType cur_union = cur_bbox.bbox_union(object_bbox);
      const CoordType cur_cover = cur_bbox.bbox_cover(object_bbox).scale();
      covers[cur_id]            = cur_cover;
      ++cur_id;
      return cur_cover - cur_union;
    });


    // Choose the entry, that needs the least enlargement to 
    // include the object's bbox
    // In case of ties, use the entry with the smalles size
    CoordType max_enlarge = std::numeric_limits<CoordType>::max();
    CoordType winner_cover = node.bbox(0).bbox_cover(object_bbox).scale();
    std::size_t winner_id = 0;
    cur_id = 0;

    std::for_each(enlargements.begin(), enlargements.end(), 
    [&covers, &max_enlarge, &winner_cover, &winner_id, &cur_id]
    (const CoordType& e)
    {
      if ( ( e < max_enlarge ) ||
           ( EQ(e, max_enlarge) && covers[cur_id] < winner_cover ) )
      {
        max_enlarge  = e;
        winner_id    = cur_id;
        winner_cover = covers[cur_id];
      }
      ++cur_id;
    });


    // In case that the found child is full, split it in two
    // nodes & and run search for leaf in this layer again
    if ( node.child(winner_id).n_entries() == M )
    {
      split_node( node.child(winner_id) );

      return choose_leaf_insertion( node, object_bbox );
    }

    return choose_leaf_insertion( node.child(winner_id), object_bbox );

  } // RTreeND::choose_leaf_insertion()

  /*------------------------------------------------------------------ 
  | Update bounding box
  ------------------------------------------------------------------*/
  void update_parent_bbox(Node& node)
  {
    // Stop at root node
    if ( node.parent() == nullptr )
      return;

    Node& parent_node = (*node.parent());

    ASSERT( !parent_node.is_leaf(), 
    "RTreeND::update_parent_bbox(): "
    "Invalid parent node.");

    // If "node" is the "i"-th child of "parent_node", then 
    // "parent_node.bbox(i)" is the bbox, that covers all entries
    // of "node"
    //
    std::size_t i = 0; 

    while ( &parent_node.child(i) != &node )
    {
      ++i;
      ASSERT( i < M, "Invalid RTree data structure.");
    }

    BBox cover = node.bbox(0);

    for ( std::size_t j = 1; j < node.n_entries(); ++j )
      cover = cover.bbox_cover( node.bbox(j) );

    parent_node.bbox( i, cover );

    return update_parent_bbox(parent_node);

  } // RTreeND::update_parent_bbox()

  /*------------------------------------------------------------------ 
  | This functin creates a new root node and adds it to a new top 
  | level of the tree
  ------------------------------------------------------------------*/
  void add_root_node()
  {
    auto new_root = std::make_unique<Node>(node_id_++);

    (*new_root).is_leaf(false);
    (*new_root).add_child( root_ );

    root_ = std::move(new_root);

    ASSERT( (*root_).n_entries() == 1, 
    "RTReeND::add_root_node(): "
    "Generation of new root node failed.");

  } // RTreeND::add_root_node()

  /*------------------------------------------------------------------ 
  | This function is called during the packing insertion of the tree.
  | It distributes a given set of nodes to a new layer
  | of nodes in the tree. The distribution is based on the nearest-X
  | algorithm.
  ------------------------------------------------------------------*/
  NodeVector build_tree_bulk_insertion(NodeVector& children)
  {
    // For each child, compute the bounding box that contains all
    // its objects 
    std::vector<BBox> bboxes {};

    for ( auto& child_ptr : children )
    {
      const Node& child = *(child_ptr);

      ASSERT( child.n_entries() > 0, "Invalid bulk insertion node");

      BBox bbox = child.bbox(0);

      for ( size_t j = 1; j < child.n_entries(); ++j )
        bbox = bbox.bbox_cover( child.bbox(j) );

      bboxes.push_back( bbox );
    }

    // Create vector of indices, which sort the given children by
    // means of their bounding boxes
    std::vector<size_t> index( bboxes.size() );
    std::iota(index.begin(), index.end(), 0);

    std::stable_sort(index.begin(), index.end(),
      [&bboxes](size_t i1, size_t i2)
    {
      return SortStrategy::choose_bbox(bboxes[i1], bboxes[i2]);
    });

    // Create new layer of parent nodes and fill distribute the given
    // child nodes to them
    NodeVector parent_nodes {};

    parent_nodes.push_back(std::make_unique<Node>(node_id_++));

    for ( size_t i = 0; i < children.size(); ++i )
    {
      if ( parent_nodes.back().get()->n_entries() >= M-1 )
      {
        parent_nodes.push_back(
          std::make_unique<Node>(node_id_++)
        );
        (*parent_nodes.back()).index( parent_nodes.size() - 1 );
      }

      Node& cur_node = *parent_nodes.back().get();

      cur_node.is_leaf( false );

      cur_node.add_child( children[index[i]] );
    }

    return std::move(parent_nodes);

  } // RTreeND::build_tree_bulk_insertion()

  /*------------------------------------------------------------------ 
  | This function is used to re-balance the entire tree after an
  | object's removal.
  ------------------------------------------------------------------*/
  void condense_tree(Node& node)
  {
    // Just update the parent bboxes if the given node still contains
    // enough other objects
    // --> Do this also if root is leaf
    if ( node.n_entries() >= M/2 || !node.parent() )
    {
      update_parent_bbox(node);
      return;
    }

    // Remove all nodes from the tree that contain too less elements
    // after the removal and collect the in a vector
    std::vector<std::unique_ptr<Node>> eliminated {}; 

    Node *cur_node_ptr = &node;

    while ( cur_node_ptr )
    {
      Node& cur_node = *cur_node_ptr;

      if ( !cur_node.parent() )
        break;

      Node& parent_node = *cur_node.parent();

      if ( cur_node.n_entries() < M/2 )
      {
        std::size_t cur_i = cur_node.index();

        eliminated.push_back( 
          std::move( parent_node.child_ptr( cur_i ) ) 
        );

        parent_node.remove_child(cur_i);
      }

      cur_node_ptr = &parent_node;
    }

    // Gather the entries from all removed nodes
    std::vector<const ObjectType*> objects { };
    std::vector<BBox> obj_bboxes {};

    for (auto& cur_node : eliminated)
      (*cur_node).collect_entry_bboxes(objects, obj_bboxes);

    // Re-insert remaining entries in tree
    for ( std::size_t i = 0; i < objects.size(); ++i )
      this->insert(*objects[i], obj_bboxes[i]);

  } // condense_tree()

  /*------------------------------------------------------------------ 
  | Find a given object in the tree which should be removed 
  | and also find its position in the node that contains the object
  ------------------------------------------------------------------*/
  Node* find_object_to_remove(const ObjectType& obj, 
                              const BBox& obj_bbox,
                              Node& node,
                              std::size_t &i_entry) 
  {
    // Node is a leaf node: Check if the given object is one of 
    // the node's entries
    if ( node.is_leaf() )
    {
      for (std::size_t i = 0; i < node.n_entries(); ++i)
        if ( &obj == &(node.object(i)) )
        {
          i_entry = i;
          return &node;
        }
      return nullptr;
    }

    // Otherwise: Check if the given object is contained in one
    // of the node's children
    for (std::size_t i = 0; i < node.n_entries(); ++i)
      if ( node.bbox(i).bbox_inside_touch( obj_bbox ) )
        return find_object_to_remove(obj, obj_bbox, 
                                     node.child(i), i_entry);

    return nullptr;

  } // RTreeND::find_object_to_remove()


  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  std::unique_ptr<Node> root_ {nullptr};

}; // RTreeND


/*********************************************************************
* This class is used to export the RTreeND structure
*********************************************************************/
template 
<
  typename    ObjectType,                      // Contained object
  std::size_t M,                               // Max. element number 
  typename    CoordType,                       // Coordinate type
  std::size_t Dim,                             // Dimensions
  typename    SortStrategy = NearestXSort,     // Sorting strategy
  typename    SplitStrategy = QuadraticSplit   // Splitting strategy
>
class RTreeNDWriter
{

public:

  using Node  = RTreeNodeND<RTREE_NODE_ARG>;
  using RTree = RTreeND<RTREE_TREE_ARG>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  RTreeNDWriter(const RTree& rtree)
  : rtree_ { &rtree }
  { }

  /*------------------------------------------------------------------ 
  | Write the RTree data to a VTU file
  ------------------------------------------------------------------*/
  void write_to_vtu(const std::string& export_prefix) const
  {
    std::string file_name = export_prefix;

    if (file_name.substr(file_name.find_last_of(".") + 1) != "vtu")
      file_name += ".vtu";

    std::vector<CoordType>   points {};
    std::vector<std::size_t> connectivity {};
    std::vector<std::size_t> offsets {};
    std::vector<std::size_t> types {};
    std::vector<int>         heights {};

    std::size_t tree_height = (*rtree_).height();

    write_vtu_data((*rtree_).root(), points, connectivity, 
                   offsets, types, heights, tree_height);

    VtuWriter writer { points, connectivity, offsets, types };
    writer.add_cell_data( heights, "height", 1);
    writer.write( file_name );

  } // write_to_vtu()

  /*------------------------------------------------------------------ 
  | Write the RTree data to a TXT file
  ------------------------------------------------------------------*/
  void write_to_txt(const std::string& export_prefix) const
  {
    std::ofstream outfile;

    std::string file_name = export_prefix;

    if (file_name.substr(file_name.find_last_of(".") + 1) != "txt")
      file_name += ".txt";

    outfile.open( file_name );

    outfile << "# Node entries, tree-level, node-index, parent-index\n";
    outfile << "# x_low, y-low, x-up, y-up\n";

    std::size_t height = (*rtree_).height();

    write_txt_data((*rtree_).root(), outfile, height);

    outfile.close();

  } // write_to_txt()

  /*------------------------------------------------------------------ 
  | Print the RTree data to the command line
  ------------------------------------------------------------------*/
  std::ostream& print(std::ostream& os) const
  {
    return print_node((*rtree_).root(), os);
  } 

private:

  /*------------------------------------------------------------------ 
  | Write the data of a rtree node to a vtu file
  ------------------------------------------------------------------*/
  std::size_t write_vtu_data(const Node&               node,
                             std::vector<CoordType>&   points,
                             std::vector<std::size_t>& connectivity,
                             std::vector<std::size_t>& offsets,
                             std::vector<std::size_t>& types,
                             std::vector<int>&         heights,
                             std::size_t               cur_height,
                             std::size_t               cur_offset=(1<<Dim))
  const
  {
    // Call method for child nodes
    if ( !node.is_leaf() )
      for (std::size_t i = 0; i < node.n_entries(); ++i)
        cur_offset = write_vtu_data(node.child(i), points, connectivity, 
                                    offsets, types, heights, 
                                    cur_height-1, cur_offset);

    std::size_t n_verts = points.size() / 3;

    for (std::size_t i = 0; i < node.n_entries(); ++i)
    {
      auto vertices = node.bbox(i).vertices();

      std::size_t v_start = n_verts;

      for (std::size_t j = 0; j < vertices.size(); ++j)
      {
        for (std::size_t k = 0; k < Dim; ++k)
          points.push_back(vertices[j][k]);

        for (std::size_t k = Dim; k < 3; ++k)
          points.push_back( -1.0f * cur_height );

        connectivity.push_back( v_start + vtk_conn_map_[ j ]);
        ++n_verts;
      }

      offsets.push_back( cur_offset );
      cur_offset += vertices.size();

      if (Dim == 3)
        types.push_back( 12 ); // VTK_HEXAHEDRON

      if (Dim == 2)
        types.push_back( 9 ); // VTK_QUAD

      if (Dim == 1)
        types.push_back( 3 ); // VTK_LINE

      heights.push_back( cur_height );
    }

    return cur_offset;

  } // write_vtu_data()

  /*------------------------------------------------------------------ 
  | Write the data of an rtree node to a txt file
  ------------------------------------------------------------------*/
  void write_txt_data(const Node&   node, 
                      std::ostream& os,
                      std::size_t   level=0, 
                      std::size_t   index=0,
                      std::size_t   parent_index=0) const
  {
    if ( !node.is_leaf() )
    {
      for (std::size_t i = 0; i < node.n_entries(); ++i)
      {
        write_txt_data(node.child(i), os, level-1, i, index);
        os << "\n";
      }
    }

    os << node.n_entries() << ", " 
       << level << ", " 
       << index << ", "
       << parent_index << "\n";

    for (std::size_t i = 0; i < node.n_entries(); ++i)
    {
      const VecND<CoordType,Dim>& ll = node.bbox(i).lowleft();
      const VecND<CoordType,Dim>& ur = node.bbox(i).upright();

      os << std::setprecision(5) << std::fixed 
         << ll.x << ", " << ll.y << ", " 
         << ur.x << ", " << ur.y;

      if ( i < node.n_entries() - 1 )
        os << "\n";
    }

  } // write_txt_data() 

  /*------------------------------------------------------------------ 
  | Print the data of an rtree-node to the command line
  ------------------------------------------------------------------*/
  std::ostream& print_node(const Node&   node, 
                           std::ostream& os, 
                           std::size_t   level=0) const
  {
    for ( std::size_t i = 0; i < node.n_entries(); ++i )
    {
      if ( level > 0 )
      {
        for ( std::size_t j = 0; j < level; ++j )
          os << "   ";
        os << "|\n";
        for ( std::size_t j = 0; j < level; ++j )
          os << "   ";
      }

      os << "*-[" << node.id() << " | " << level 
         << " - " << i+1 << "/" << node.n_entries() << "]: ";
      os << node.bbox(i);

      os << " - Local index: " << node.index();

      if ( node.left() )
        os << " - Left: " << (*node.left()).id();
      else
        os << " - Left: NULL";

      if ( node.right() )
        os << " - Right: " << (*node.right()).id();
      else
        os << " - Right: NULL";
      os << "\n";

      if ( !node.is_leaf() )
        print_node(node.child(i), os, level+1);
    }

    if ( level == 1 )
      os << "\n";

    return os;

  } // print_node()

  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  const RTree* rtree_;

  // Array to map from BBoxND vertex coordinates
  // to VTK hexahedral coordinates
  const std::array<std::size_t,8> vtk_conn_map_ 
  { 0, 1, 2, 3, 7, 6, 5, 4 };

}; // RTreeNDWriter

} // namespace CppUtils
