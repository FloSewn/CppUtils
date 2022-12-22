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
#include "Log.h"
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
    std::size_t seed_right = 0;

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
  void parent(Node& p) { parent_ = &p; }
  void parent(const Node& p) { parent_ = const_cast<Node*>(&p); }

  void is_leaf(bool t) { is_leaf_ = t; }
  void n_entries(std::size_t n) { n_entries_ = n; }
  void id(std::size_t i) { id_ = i; }
  void index(std::size_t index) { index_ = index; }

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
        "RTreeNodeND: Unable to access entry at position " 
        + std::to_string(i) );
    return entries_[i]; 
  }

  Entry& entry(std::size_t i)
  {
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to access entry at position " 
        + std::to_string(i) );
    return entries_[i]; 
  }

  /*------------------------------------------------------------------ 
  | Access minimum bounding rectangles
  ------------------------------------------------------------------*/
  const BBox& bbox(std::size_t i) const
  { 
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to access key at position " 
        + std::to_string(i) );
    return entries_[i].bbox(); 
  }

  /*------------------------------------------------------------------ 
  | Access objects
  ------------------------------------------------------------------*/
  ObjectType& object(std::size_t i) 
  { 
    ASSERT( i < n_entries() && is_leaf(), 
        "RTreeNodeND: Unable to access key at position " 
        + std::to_string(i) );
    return *(const_cast<ObjectType*>(entries_[i].object())); 
  }

  const ObjectType& object(std::size_t i) const 
  { 
    ASSERT( i < n_entries() && is_leaf(), 
        "RTreeNodeND: Unable to access key at position " 
        + std::to_string(i) );
    return *(entries_[i].object());
  }

  /*------------------------------------------------------------------ 
  | Access children pointers
  ------------------------------------------------------------------*/
  Node_ptr& child_ptr(std::size_t i) 
  { 
    ASSERT( i < n_entries() && !is_leaf(), 
        "RTreeNodeND: Unable to access child at position " 
        + std::to_string(i) );
    return entries_[i].child_ptr();
  }

  const Node_ptr& child_ptr(std::size_t i) const
  { 
    ASSERT( i < n_entries() && !is_leaf(), 
        "RTreeNodeND: Unable to access child at position " 
        + std::to_string(i) );
    return entries_[i].child_ptr(); 
  }

  /*------------------------------------------------------------------ 
  | Access children as references
  ------------------------------------------------------------------*/
  Node& child(std::size_t i) 
  { 
    ASSERT( i < n_entries() && !is_leaf(), 
        "RTreeNodeND: Unable to access child at position " 
        + std::to_string(i) );
    return entries_[i].child();
  }

  const Node& child(std::size_t i) const 
  { 
    ASSERT( i < n_entries() && !is_leaf(), 
        "RTreeNodeND: Unable to access child at position " 
        + std::to_string(i) );
    return entries_[i].child();
  }

  /*------------------------------------------------------------------ 
  | Set bounding boxes of children / objects
  ------------------------------------------------------------------*/
  void bbox(std::size_t i, const BBox& b) 
  { 
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to access key at position " 
        + std::to_string(i) );
    entries_[i].bbox(b);
  }

  /*------------------------------------------------------------------ 
  | Set objects
  ------------------------------------------------------------------*/
  void object(std::size_t i, const ObjectType& obj)
  { 
    ASSERT( i < n_entries() && is_leaf(), 
        "RTreeNodeND: Unable to set object at position " 
        + std::to_string(i) );
    entries_[i].object(&obj);
  }

  /*------------------------------------------------------------------ 
  | Set children
  ------------------------------------------------------------------*/
  void child(std::size_t i, std::unique_ptr<Node>& c)
  { 
    ASSERT( i < n_entries() && !is_leaf(), 
        "RTreeNodeND: Unable to set child at position " 
        + std::to_string(i) );
    entries_[i].child(c);
  }

  /*------------------------------------------------------------------ 
  | Add new object to the node
  ------------------------------------------------------------------*/
  void add_object(const ObjectType& object, const BBox& obj_bbox)
  {
    const std::size_t n_entries = this->n_entries();

    ASSERT( n_entries < M, "RTreeNodeND: Can not add more than " 
        + std::to_string(M) + " objects.");
    ASSERT( is_leaf(), "RTreeNodeND: Can not add object to "
        "non-leaf node.");

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

    ASSERT( n_entries < M, "RTreeNodeND: Can not add more than " 
        + std::to_string(M) + " children.");
    ASSERT( !is_leaf(), "RTreeNodeND: Can not add child to "
        "non-leaf node.");

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
    child.parent(*this);
    child.index(i);

  } // add_child()

  /*------------------------------------------------------------------ 
  | Remove the object at a given index 
  ------------------------------------------------------------------*/
  void remove_object(std::size_t index)
  {
    const std::size_t n_entries = this->n_entries();

    ASSERT( n_entries > 0, "RTreeNodeND: Invalid object removal." );
    ASSERT( is_leaf(), "RTreeNodeND: Can not remove object from "
        "non-leaf node.");

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
     
    ASSERT( n_entries > 0, "RTreeNodeND: Invalid child node removal." );
    ASSERT( !is_leaf(), "RTreeNodeND: Can not remove child from "
        "leaf node.");

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
  | Iterator implementation
  ------------------------------------------------------------------*/
  struct Iterator
  {
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = ObjectType;
    using pointer           = ObjectType*;
    using reference         = ObjectType&;

    Iterator() {}
    Iterator(Node& node) : cur_node_ {&node} 
    { 
      if ( (*cur_node_).n_entries() == 0 )
        cur_node_ = nullptr;
    }

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
    }

    Iterator operator++(int) 
    { Iterator tmp = *this; ++(*this); return tmp; }

    friend bool operator== (const Iterator& a, const Iterator& b) 
    { return (   a.cur_node_ == b.cur_node_ 
              && a.cur_index_ == b.cur_index_); }

    friend bool operator!= (const Iterator& a, const Iterator& b) 
    { return !(a == b); }

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

    ConstantIterator() {}
    ConstantIterator(Node& node) : cur_node_ {&node} 
    { 
      if ( (*cur_node_).n_entries() == 0 )
        cur_node_ = nullptr;
    }

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
  | Constructor
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
  | 
  | 
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

    return true;

  } // RTreeND::remove(ObjectType)


  /*------------------------------------------------------------------ 
  | 
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

    // Pick all objects of children that will be removed
    std::vector<const ObjectType*> objects { };
    std::vector<BBox> obj_bboxes {};

    Node* cur_node_ptr = &node;

    // Traverse until root 
    while ( (*cur_node_ptr).parent() )
    {
      Node& cur_node    = *cur_node_ptr;
      Node& parent_node = *cur_node.parent();

      if ( cur_node.n_entries() < M/2 )
      {
        // Store remaining objects of leaf nodes
        if ( cur_node.is_leaf() )
        {
          for (std::size_t i = 0; i < cur_node.n_entries(); ++i )
          {
            objects.push_back( &cur_node.object(i) );
            obj_bboxes.push_back( cur_node.bbox(i) );
          }
        }

        // Remove node
        std::size_t i = cur_node.index();
        parent_node.remove_child(i);
      }

      cur_node_ptr = &parent_node;
    }

    // Re-insert remaining objects in tree
    for ( std::size_t i = 0; i < objects.size(); ++i )
      this->insert(*objects[i], obj_bboxes[i]);

  } // condense_tree()

  /*------------------------------------------------------------------ 
  | 
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

      cur_child.parent( root );
      root.n_entries( i + 1 );
      root.is_leaf( false );
      root.bbox(i, cur_child.bbox() );
      root.child(i, node_layer[i] );
    }
      
  } // RTreeND::insert()


private:

  /*------------------------------------------------------------------ 
  | This function splits the i-th child of a given "parent_node"
  ------------------------------------------------------------------*/
  void split_child(Node& parent_node, std::size_t i)
  {
    // Check that parent node is not full
    ASSERT(parent_node.n_entries() != M, 
    "RTreeND::split_child(): "
    "Can not split child of full RTree node.");

    // This is the child node, whose entries will be splitted 
    Node& child_node = parent_node.child(i);

    // Check that child node is full
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

    // If the leaf has enough space to store the object, add it
    if ( leaf.n_entries() < M )
    {
      leaf.add_object( object, obj_bbox );

      // Update all BBoxes in the path from root to this leaf, 
      // so that all of them cover the object's bbox
      update_parent_bbox(leaf);
    }
    // Otherwise, split the leaf in two nodes. This might lead its 
    // parent to overflow, thus leading it to be splitted recursively.
    // If the root node must be splitted, a new root node will be 
    // created, which then keeps the splitted "old" root node 
    // as children.
    else
    {
      ASSERT( false, "IMPLEMENTATION ERROR");
    } 

  } // RTreeND::insert_nonfull()

  /*------------------------------------------------------------------ 
  | Choose appropriate leaf node for RTree-insertion
  ------------------------------------------------------------------*/
  Node& choose_leaf_insertion(Node&       node,
                              const BBox& object_bbox)
  {
    // Choos only leaf nodes
    if ( node.is_leaf() )
      return node;

    CoordType max_enlarge = std::numeric_limits<CoordType>::max();
    std::size_t j = 0;

    CoordType cover_j = node.bbox(j).bbox_cover(object_bbox).scale();

    // Find the child-node that has the least enlargement with the 
    // object's bbox
    for ( std::size_t i = 0; i < node.n_entries(); ++i )
    {
      const BBox& child_bbox = node.bbox(i);

      // Compute enlargement
      const CoordType union_i = child_bbox.bbox_union(object_bbox);
      const CoordType cover_i 
        = node.bbox(i).bbox_cover(object_bbox).scale();

      const CoordType enlarge = cover_i - union_i;

      // Choose the entry, that needs the least enlargement to 
      // include the object's bbox
      // In case of ties, use the entry with the smalles size
      if ( ( enlarge < max_enlarge ) ||
           ( EQ(enlarge, max_enlarge) && cover_i < cover_j ) )
      {
        max_enlarge = enlarge;
        j = i;
        cover_j = cover_i;
      }
    }

    // In case that the found child is full, split it in two
    // nodes & and run search for leaf in this layer again
    if ( node.child(j).n_entries() == M )
    {
      split_child(node, j);

      if ( node.parent() )
        return choose_leaf_insertion( *node.parent(), object_bbox );

      return choose_leaf_insertion( node, object_bbox );
    }

    return choose_leaf_insertion( node.child(j), object_bbox );

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

    ASSERT( !parent_node.is_leaf(), "Invalid RTree data structure.");

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

    ASSERT( (*root_).n_entries() == 1, "Invalid root node.");
    //root_->child(0).parent(*root_);   // maybe needed?

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
