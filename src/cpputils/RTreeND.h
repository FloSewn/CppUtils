/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
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

// Array to map from BBoxND vertex coordinates
// to VTK hexahedral coordinates
static std::array<std::size_t,8> BBOXND_VTK_CONN_MAP 
{ 0, 1, 2, 3, 7, 6, 5, 4 };

static inline long unsigned RTreeNodeID = 0;

/*********************************************************************
* Forward declarations 
*********************************************************************/
template 
<
  typename    ObjectType,                  // Contained object
  std::size_t M,                           // Max. element number 
  typename    CoordType,                   // Coordinate type
  std::size_t Dim,                         // Dimensions
  typename    SortStrategy                 // Sorting strategy
>
class RTreeNodeND;

template 
<
  typename    ObjectType,                  // Contained object
  std::size_t M,                           // Max. element number 
  typename    CoordType,                   // Coordinate type
  std::size_t Dim,                         // Dimensions
  typename    SortStrategy,                // Sorting strategy
  typename    SplitStrategy                // Splitting strategy
>
class RTreeND;



/*********************************************************************
* This class defines a strategy for the sorting of RTree nodes
* within a layer
*********************************************************************/
class NearestXSort
{
public:

  /*------------------------------------------------------------------ 
  | 
  ------------------------------------------------------------------*/
  template <typename ObjectType>
  static inline
  std::vector<const ObjectType*> 
  sort(const std::vector<ObjectType>& objects)
  {
    std::vector<const ObjectType*> sorted_objects;

    for ( const ObjectType& obj : objects )
      sorted_objects.push_back( &obj );

    std::sort(sorted_objects.begin(), sorted_objects.end(),
              [](const ObjectType* lhs, const ObjectType* rhs)
    {
      return choose_bbox((*lhs).bbox(), (*rhs).bbox());
    });

    return std::move( sorted_objects );

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
    VecND<CoordType,Dim> delta = lhs.lowleft() - rhs.lowleft();

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
* Quadratic split strategy
*********************************************************************/
class QuadraticSplit
{
public:
  /*------------------------------------------------------------------ 
  | This function is called upon the splitting of the entries of 
  | a given "node" in two sets "n1" and "n2".
  | The output array "add_to_n2" marks all entries of the input node, 
  | wether they will be distributed to the set "n1" (false) or 
  | to the set "n2" (true).
  ------------------------------------------------------------------*/
  template<
    typename    ObjectType,                  // Contained object
    std::size_t M,                           // Max. element number 
    typename    CoordType,                   // Coordinate type
    std::size_t Dim,                         // Dimensions
    typename    SortStrategy                 // Sorting strategy
  >
  static inline
  std::array<bool,M> 
  split(const RTreeNodeND<ObjectType,M,CoordType,Dim,SortStrategy>& node)
  {
    std::array<bool,M> add_to_n2 { false };

    // This array marks all distributed entries
    std::array<bool,M> distributed { false };

    // This variable referes to the remaining number of entries, which 
    // have not yet been distributed
    std::size_t n_remaining = node.n_entries();

    // These are the bounding boxes that enclose all objects in the 
    // two sets
    BBoxND<CoordType,Dim> bbox_1 {};
    BBoxND<CoordType,Dim> bbox_2 {};

    // This function picks the first entry for each group
    // The entry for n1 is added indirectly through its marker 
    // in the array "distributed"
    std::size_t seed = pick_seeds(node, distributed, n_remaining, 
                                  bbox_1, bbox_2);
    add_to_n2[seed]  = true;
    ASSERT( n_remaining == node.n_entries()-2, 
        "Error in function RTree::pick_seeds() failed.");

    // These are the number of entries that are already distributed 
    // to both sets "n1" and "n2" 
    std::size_t n_n1 = 1;
    std::size_t n_n2 = 1;

    while ( true )
    {
      // Stop if all entries have been assigned
      if ( n_remaining == 0 )
        return add_to_n2;

      // First node has so few entries, that all remaining entries must
      // be assigned to it in order to have the minimum number
      if ( n_n1 + n_remaining == M/2 )
      {
        ASSERT( n_n2+1 >= M/2, "Error in function RTree::quadratic_split");
        return add_to_n2;
      }

      // Second node has so few entries, that all remaining entries must
      // be assigned to it in order to have the minimum number
      if ( n_n2 + n_remaining == M/2 )
      {
        for ( std::size_t i = 0; i < M; ++i )
          add_to_n2[i] = (distributed[i] == false) ? true : add_to_n2[i];

        ASSERT( n_n1+1 >= M/2, "Error in function RTree::quadratic_split");
        return add_to_n2;
      }

      // Pick the next entry to assign. This function also 
      // decrements "n_remaining" and handles the markation
      // of entry "i" in the array "distributed"
      std::size_t i = pick_next(node, distributed, n_remaining, 
                                bbox_1, bbox_2);

      // Add entry "i" to the group whose covering rectangle will
      // have to be enlarged least to accomodate it.
      const BBoxND<CoordType,Dim>& E_i = node.bbox(i);

      const BBoxND<CoordType,Dim> cover_1 = bbox_1.bbox_cover(E_i);
      const BBoxND<CoordType,Dim> cover_2 = bbox_2.bbox_cover(E_i);
      
      const double d1 = cover_1.scale() - bbox_1.scale(); 
      const double d2 = cover_2.scale() - bbox_2.scale(); 

      bool add_entry_to_n2 = false;

      // Add entry to set "n2", if its covering rectangle will be 
      // enlarged less than for set "n1"
      if ( d2 < d1 )
        add_entry_to_n2 = true;

      // In case of ties 
      if ( EQ(d1, d2) ) 
      {
        // Add entry to set with smaller size
        if ( bbox_2.scale() < bbox_1.scale() )
          add_entry_to_n2 = true;

        // Then to the one with fewer entries
        if ( EQ(bbox_1.scale(), bbox_2.scale()) && (n_n2 < n_n1) )
          add_entry_to_n2 = true;
      }

      // Otherwise, add entry to set "n1"
      if ( add_entry_to_n2 )
      {
        add_to_n2[i] = true;
        bbox_2       = cover_2;
      }
      else
      {
        bbox_1 = cover_1;
      }
    }

  } // QuadraticSplit::split()


private:

  /*------------------------------------------------------------------ 
  | Select two entries of the given full "node" which are used as 
  | seeds for the quadratic_split() algorithm during a node splitting
  | operation.
  ------------------------------------------------------------------*/
  template<
    typename    ObjectType,                  // Contained object
    std::size_t M,                           // Max. element number 
    typename    CoordType,                   // Coordinate type
    std::size_t Dim,                         // Dimensions
    typename    SortStrategy                 // Sorting strategy
  >
  static inline
  std::size_t 
  pick_seeds(const RTreeNodeND<ObjectType,M,CoordType,Dim,SortStrategy>& node, 
             std::array<bool,M>     distributed,
             std::size_t&           n_remaining,
             BBoxND<CoordType,Dim>& bbox_1,
             BBoxND<CoordType,Dim>& bbox_2) 
  {
    ASSERT(node.n_entries() == M, "Invalid data structure passed to "
        "function RTree::pick_seeds()");

    std::size_t seed_1 = -1;
    std::size_t seed_2 = -1;

    double ineff = 0.0;

    // Calculate inefficiency of grouping entries together
    for (std::size_t i = 0; i < M; ++i)
    {
      const BBoxND<CoordType,Dim>& E_i = node.bbox(i);

      for (std::size_t j = 0; j < M; ++j)
      {
        if ( i == j )
          continue;

        // Compose a rectangle J including E_i and E_j 
        const BBoxND<CoordType,Dim>& E_j = node.bbox(j);

        const BBoxND<CoordType,Dim> J = E_i.bbox_cover(E_j);

        // Compute the inefficiency size
        const double diff = J.scale() - E_i.scale() - E_j.scale();

        // Choose the most wasteful pair
        if ( diff > ineff )
        {
          ineff  = diff;
          seed_1 = i;
          seed_2 = j;
        }
      }
    }

    ASSERT( seed_1 != seed_2, "Function RTree::pick_seeds() failed.");
    
    // Mark chosen elements
    bbox_1 = node.bbox(seed_1);
    distributed[seed_1] = true;
    --n_remaining;

    bbox_2 = node.bbox(seed_2);
    distributed[seed_2] = true;
    --n_remaining;

    return seed_2;

  } // QuadraticSplit::pick_seeds()

  /*------------------------------------------------------------------ 
  | 
  ------------------------------------------------------------------*/
  template<
    typename    ObjectType,                  // Contained object
    std::size_t M,                           // Max. element number 
    typename    CoordType,                   // Coordinate type
    std::size_t Dim,                         // Dimensions
    typename    SortStrategy                 // Sorting strategy
  >
  static inline
  std::size_t 
  pick_next(const RTreeNodeND<ObjectType,M,CoordType,Dim,SortStrategy>& node, 
            std::array<bool,M>     distributed,
            std::size_t&           n_remaining,
            BBoxND<CoordType,Dim>& bbox_1,
            BBoxND<CoordType,Dim>& bbox_2) 
  {
    std::size_t entry = 0;

    double max_diff = 0.0;

    // Determine the cost of putting each entry in each group
    for (std::size_t i = 0; i < M; ++i)
    {
      if ( distributed[i] )
        continue;

      const BBoxND<CoordType,Dim>& E_i = node.bbox(i);

      // Compute the size increases required in the covering rectangles
      // "bbox_1" and "bbox_2" to include entry "i"
      BBoxND<CoordType,Dim> C_1 = bbox_1.bbox_cover(E_i);
      BBoxND<CoordType,Dim> C_2 = bbox_2.bbox_cover(E_i);

      const double d1 = C_1.scale() - bbox_1.scale();
      const double d2 = C_2.scale() - bbox_2.scale();

      const double diff = ABS(d1 - d2);

      // Choose the entry with the maximum difference between d1 and d2
      if ( diff > max_diff )
      {
        entry    = i;
        max_diff = diff;
      }
    }

    distributed[entry] = true;
    --n_remaining;

    return entry;

  } // QuadraticSplit::pick_next()

}; // QuadraticSplit


/*********************************************************************
* Entry to store the rtree data 
*********************************************************************/
template 
<
  typename    ObjectType,                  // Contained object
  std::size_t M,                           // Max. element number 
  typename    CoordType,                   // Coordinate type
  std::size_t Dim,                         // Dimensions
  typename    SortStrategy                 // Sorting strategy
>
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
  void parent(Node* p) { parent_ = p; }

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  BBox& bbox() { return bbox_; }
  const BBox& bbox() const { return bbox_; } 

  Node_ptr& child_ptr() { return child_; }
  const Node_ptr& child_ptr() const { return child_; }

  Node& child() { return *child_; }
  const Node& child() const { return *child_; }

  const ObjectType* object() const { return object_; }

  Node* parent() { return parent_; }
  const Node* parent() const { return parent_; }


private:
  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  BBox               bbox_   {};
  Node_ptr           child_  { nullptr };
  const ObjectType*  object_ { nullptr };
  Node*              parent_ { nullptr };

}; // RTreeEntryND


/*********************************************************************
* References
* ----------
* - https://tildesites.bowdoin.edu/~ltoma/teaching/cs340/spring08/\
*   Papers/Rtree-chap1.pdf
*
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
  using Node     = RTreeNodeND<ObjectType,M,CoordType,Dim,SortStrategy>;
  using Node_ptr = std::unique_ptr<Node>;
  using Entry    = RTreeEntryND<ObjectType,M,CoordType,Dim,SortStrategy>;
  using Entries  = std::array<Entry,M>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  RTreeNodeND(int id) { id_ = id; }

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

  Node* left() { return left_; }
  const Node* left() const { return left_; }

  Node* right() { return right_; }
  const Node* right() const { return right_; }

  /*------------------------------------------------------------------ 
  | Setter
  ------------------------------------------------------------------*/
  void parent(Node& p) { parent_ = &p; }
  void parent(const Node& p) { parent_ = const_cast<Node*>(&p); }

  void is_leaf(bool t) { is_leaf_ = t; }
  void n_entries(std::size_t n) { n_entries_ = n; }
  void id(std::size_t i) { id_ = i; }

  void left(Node& n) { left_ = &n; }
  void right(Node& n) { right_ = &n; }

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
  | Export the node data to the VTK format
  ------------------------------------------------------------------*/
  std::size_t export_to_vtu(std::vector<CoordType>&   points,
                            std::vector<std::size_t>& connectivity,
                            std::vector<std::size_t>& offsets,
                            std::vector<std::size_t>& types,
                            std::vector<int>&         heights,
                            std::size_t               cur_height,
                            std::size_t               cur_offset=(1<<Dim))
  {
    // Call method for child nodes
    if ( !is_leaf() )
      for (std::size_t i = 0; i < n_entries(); ++i)
        cur_offset = child(i).export_to_vtu(points, connectivity, 
                                            offsets,types, heights, 
                                            cur_height-1, cur_offset);

    std::size_t n_verts = points.size() / 3;

    for (std::size_t i = 0; i < n_entries(); ++i)
    {
      auto vertices = bbox(i).vertices();

      std::size_t v_start = n_verts;

      for (std::size_t j = 0; j < vertices.size(); ++j)
      {
        for (std::size_t k = 0; k < Dim; ++k)
          points.push_back(vertices[j][k]);

        for (std::size_t k = Dim; k < 3; ++k)
          points.push_back( -1.0f * cur_height );

        connectivity.push_back( v_start + BBOXND_VTK_CONN_MAP[ j ]);
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

  } // RTreeNodeND::export_to_vtu()

  /*------------------------------------------------------------------ 
  | Export the node data 
  ------------------------------------------------------------------*/
  std::ostream& export_to_txt(std::ostream& os,  
                              std::size_t   level=0, 
                              std::size_t   index=0,
                              std::size_t   parent_index=0) const
  {
    if ( !is_leaf() )
    {
      for (std::size_t i = 0; i < n_entries(); ++i)
      {
        child(i).export_to_txt(os, level-1, i, index);
        os << "\n";
      }
    }

    os << n_entries() << ", " 
       << level << ", " 
       << index << ", "
       << parent_index << "\n";

    for (std::size_t i = 0; i < n_entries(); ++i)
    {
      const VecND<CoordType,Dim>& ll = bbox(i).lowleft();
      const VecND<CoordType,Dim>& ur = bbox(i).upright();

      os << std::setprecision(5) << std::fixed 
         << ll.x << ", " << ll.y << ", " 
         << ur.x << ", " << ur.y;

      if ( i < n_entries() - 1 )
        os << "\n";
    }

    return os;

  } // RTreeNodeND::export_to_txt()

  /*------------------------------------------------------------------ 
  | Print out the tree structure to the command line
  ------------------------------------------------------------------*/
  std::ostream& print(std::ostream& os,  
                      std::size_t level=0) const
  {
    for ( std::size_t i = 0; i < n_entries(); ++i )
    {
      if ( level > 0 )
      {
        for ( std::size_t j = 0; j < level; ++j )
          os << "   ";
        os << "|\n";
        for ( std::size_t j = 0; j < level; ++j )
          os << "   ";
      }

      os << "*-[" << id_ << " | " << level 
         << " - " << i+1 << "/" << n_entries() << "]: ";
      os << bbox(i) << "\n";

      if ( !is_leaf() )
        child(i).print(os, level+1);
    }

    if ( level == 1 )
      os << "\n";

  } // print()

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
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to access key at position " 
        + std::to_string(i) );
    return *(const_cast<ObjectType*>(entries_[i].object())); 
  }

  const ObjectType& object(std::size_t i) const 
  { 
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to access key at position " 
        + std::to_string(i) );
    return *(entries_[i].object());
  }

  /*------------------------------------------------------------------ 
  | Access children pointers
  ------------------------------------------------------------------*/
  Node_ptr& child_ptr(std::size_t i) 
  { 
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to access child at position " 
        + std::to_string(i) );
    return entries_[i].child_ptr();
  }

  const Node_ptr& child_ptr(std::size_t i) const
  { 
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to access child at position " 
        + std::to_string(i) );
    return entries_[i].child_ptr(); 
  }

  /*------------------------------------------------------------------ 
  | Access children as references
  ------------------------------------------------------------------*/
  Node& child(std::size_t i) 
  { 
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to access child at position " 
        + std::to_string(i) );
    return entries_[i].child();
  }

  const Node& child(std::size_t i) const 
  { 
    ASSERT( i < n_entries(), 
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
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to set object at position " 
        + std::to_string(i) );
    entries_[i].object(&obj);
  }

  /*------------------------------------------------------------------ 
  | Set children
  ------------------------------------------------------------------*/
  void child(std::size_t i, std::unique_ptr<Node>& c)
  { 
    ASSERT( i < n_entries(), 
        "RTreeNodeND: Unable to set child at position " 
        + std::to_string(i) );
    entries_[i].child(c);
  }

  /*------------------------------------------------------------------ 
  | Add new object to the node
  ------------------------------------------------------------------*/
  void add_object(const ObjectType& object)
  {
    const std::size_t n_entries = this->n_entries();
    ASSERT( n_entries < M, "RTreeNodeND: Can not add more than " 
        + std::to_string(M) + " objects.");

    //std::size_t i = 0;
    std::size_t i = n_entries;

    // Obtain index i at which all former object entries are located
    // closer to the origin. Start searching from the back (faster)
    while ( i > 0 )
    {
      if ( SortStrategy::choose_bbox(this->bbox(i-1), object.bbox()) )
        break;
      --i;
    }

    this->n_entries( n_entries + 1 );

    // Shift remaining entries to the right
    for (std::size_t j = n_entries; j > i; --j)
    {
      this->bbox(j, this->bbox(j-1));
      this->child(i, this->child_ptr(j-1));
    }

    // Finally add the new object
    this->bbox(i, object.bbox()); 
    this->object(i, object);

  } // add_object()

  /*------------------------------------------------------------------ 
  | Add new child to the node
  ------------------------------------------------------------------*/
  void add_child(Node_ptr& child_ptr, bool sort=true)
  {
    Node& child = *child_ptr;
    const std::size_t n_entries = this->n_entries();

    ASSERT( n_entries < M, "RTreeNodeND: Can not add more than " 
        + std::to_string(M) + " children.");

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
      this->child(i, this->child_ptr(j-1));
    }

    // Finally add the new child
    this->bbox(i, child.bbox()); 
    this->child(i, child_ptr);
    this->is_leaf( false );
    child.parent( *this );

    // Set connectivity between child nodes
    if ( i > 0 )
    {
      child.left( this->child(i-1) );
      this->child(i-1).right( child );
    }
    
    if ( i < this->n_entries() - 1 )
    {
      child.right( this->child(i+1) );
      this->child(i+1).left( child );
    }


  } // add_child()

private:
  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  Entries             entries_   { };
  Node*               parent_    { nullptr };

  bool                is_leaf_   { true };
  std::size_t         n_entries_ { 0 };
  std::size_t         id_        { 0 };
  Node*               left_      { nullptr };
  Node*               right_     { nullptr };

}; // RTreeNodeND


/*********************************************************************
* This class defines the interface to an R-tree structure 
* for N-dimensional simplices
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
public:

  using BBox       = BBoxND<CoordType,Dim>;
  using Node       = RTreeNodeND<ObjectType,M,CoordType,Dim,SortStrategy>;
  using NodeVector = std::vector<std::unique_ptr<Node>>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  RTreeND()
  {
    root_ = std::make_unique<Node>(RTreeNodeID++);
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

  } // RTreeND::height()

  /*------------------------------------------------------------------ 
  | Print out the r-tree structure to the command line
  ------------------------------------------------------------------*/
  std::ostream& print(std::ostream& os)
  { (*root_).print(os); }

  /*------------------------------------------------------------------ 
  | 
  ------------------------------------------------------------------*/
  void write_to_vtu(const std::string& path) const
  {
    std::string file_name = path;

    if (file_name.substr(file_name.find_last_of(".") + 1) != "vtu")
      file_name += ".vtu";
    
    std::vector<CoordType> points {};
    std::vector<std::size_t> connectivity {};
    std::vector<std::size_t> offsets {};
    std::vector<std::size_t> types {};
    std::vector<int> heights {};

    std::size_t tree_height = height();

    (*root_).export_to_vtu(points, connectivity, offsets, 
                           types, heights, tree_height);

    VtuWriter writer { points, connectivity, offsets, types };
    writer.add_cell_data( heights, "height", 1);
    writer.write( file_name );


  } // RTreeND::write_to_vtu()

  /*------------------------------------------------------------------ 
  | Write the R-Tree structure to a text file
  ------------------------------------------------------------------*/
  void write_to_file(const std::string& path) const
  {
    std::ofstream outfile;

    std::string file_name = path;

    if (file_name.substr(file_name.find_last_of(".") + 1) != "txt")
      file_name += ".txt";

    outfile.open( file_name );

    outfile << "# Node entries, tree-level, node-index, parent-index\n";
    outfile << "# x_low, y-low, x-up, y-up\n";

    outfile << (*this);

    outfile.close();

  } // RTreeND::write_to_file()



  /*------------------------------------------------------------------ 
  | 
  | 
  ------------------------------------------------------------------*/
  bool remove(const ObjectType& obj)
  {
    if ( root_->n_entries() < 1 )
      return false;

    // Check if the object is contained in the R-Tree region
    if ( !root_->bbox().bbox_inside_touch( obj.bbox() ) )
      return false;

    // Find the node & index of the respective entry to remove
    std::size_t i_entry {};
    Node* node = find_object_to_remove(obj, node, i_entry);

    // Return if no node is found 
    if ( node == nullptr )
      return false;

    // Remove entry from node
    //
    // If node.n_entries() >= M/2 + 1
    //   Simply remove the entry and maybe shift following 
    //   entries to the front
    // Else
    //   Compute number N of entries among all nodes that belong
    //   to the curren node's parent 
    //
    //   Either take an entry from the left node in the current layer
    //   -> take_from_left_node()
    //
    //   Otherwise if not enough entries, remove layer and distribute
    //   remaining entries among nodes in next layer
    //
    //remove_node_entry(node, i_entry);



  } // RTreeND::remove(ObjectType)


  /*------------------------------------------------------------------ 
  | 
  ------------------------------------------------------------------*/
  Node* find_object_to_remove(const ObjectType& obj, 
                              const Node& node,
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
      if ( node.bbox(i).bbox_inside_touch( obj.bbox() ) )
        return find_object_to_remove(obj, node.child(i), i_entry);

    return nullptr;

  } // RTreeND::find_object_to_remove()



  /*------------------------------------------------------------------ 
  | Insert a new object into the RTree structure
  ------------------------------------------------------------------*/
  void insert(const ObjectType& object)
  {
    // Handle the case where the root node is full
    if ( root_->n_entries() == M )
    {
      add_root_node();
      split_child(*root_, 0);
    }
     
    insert_nonfull(*root_, object);

  } // RTreeND::insert()

  /*------------------------------------------------------------------ 
  | Insert a bulk of objects into the RTree structure
  ------------------------------------------------------------------*/
  void insert(const std::vector<ObjectType>& objects)
  {
    // Put all objects in a temporary container, which will be 
    // sorted
    std::vector<const ObjectType*> sorted_objects 
      = SortStrategy::sort( objects );

    // Group objects into (N / M) leaf nodes
    NodeVector node_layer {};
    node_layer.push_back( std::make_unique<Node>(RTreeNodeID++) );

    for ( const ObjectType* obj : sorted_objects )
    {
      Node& cur_leaf = *node_layer.back().get();

      cur_leaf.add_object( *obj );

      if ( cur_leaf.n_entries() >= M-1 )
        node_layer.push_back( std::make_unique<Node>(RTreeNodeID++) );
    }

    // Recursively pack nodes into a node layer at the 
    // next level until root is reached
    while ( node_layer.size() > M )
      node_layer = build_tree_bulk_insertion(node_layer);

    // Place remaining nodes into root node
    for ( std::size_t i = 0; i < node_layer.size(); ++i )
    {
      Node& root = *root_;
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
    ASSERT(parent_node.n_entries() != M, "Invalid R-Tree structure.");

    // Check for valid child pointers
#ifndef NDEBUG
    for (std::size_t j = 0; j <= i; ++j)
      ASSERT(&parent_node.child(j) != nullptr,
          "Invalid child pointer at position " + std::to_string(j));
#endif

    // This is the child node, whose entries will be splitted 
    Node& child_node = parent_node.child(i);

    // Check that child node is full
    ASSERT(child_node.n_entries() == M, "Invalid R-Tree structure.");

    // This is the new node, which will get half of the entries 
    // of "child_node"
    auto new_node_ptr = std::make_unique<Node>(RTreeNodeID++);
    Node& new_node = *new_node_ptr;

    new_node.is_leaf( child_node.is_leaf() );

    // This array contains the information, which entries of the 
    // child node "child_node" will be added to "new_node"
    std::array<bool,M> add_to_new = SplitStrategy::split( child_node );

    // Distribute entries from "child_node" to "new_node"
    for ( std::size_t j = 0; j < M; ++j )
    {
      if ( !add_to_new[j] )
        continue;

      if ( new_node.is_leaf() )
        new_node.add_object( child_node.object(j) );
      else
        new_node.add_child( child_node.child_ptr(j) );
    }

    // Re-distribute remining entries in "child_node"
    for ( std::size_t j = M; j > 0; --j )
    {
      if ( !add_to_new[j-1] )
        continue;

      // Put all elements right of j one to the left
      for ( std::size_t k = j; k < M; ++k )
      {
        child_node.object( k-1, child_node.object(k) );
        child_node.bbox( k-1, child_node.bbox(k) );
        child_node.child( k-1, child_node.child_ptr(k) );
      }

      child_node.n_entries( child_node.n_entries() - 1 );
    }

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
  void insert_nonfull(Node& node, const ObjectType& object)
  {
    const BBox& bb_obj = object.bbox();

    // Choose an appropriate leaf to insert the object
    Node& leaf = choose_leaf_insertion(node, bb_obj);

    // If the leaf has enough space to store the object, add it
    if ( leaf.n_entries() < M )
    {
      leaf.add_object( object );

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

    double m = CPPUTILS_MAX;
    std::size_t j = 0;

    double cover_j = node.bbox(j).bbox_cover(object_bbox).scale();

    // Find the child-node that has the least enlargement with the 
    // object's bbox
    for ( std::size_t i = 0; i < node.n_entries(); ++i )
    {
      const BBox& child_bbox = node.bbox(i);

      // Compute enlargement
      const double c = child_bbox.bbox_union(object_bbox)
                     - child_bbox.scale();

      const double cover_i = node.bbox(i).bbox_cover(object_bbox).scale();

      // Choose the entry, that needs the least enlargement to 
      // include the object's bbox
      // In case of ties, use the entry with the smalles size
      if ( ( c < m ) ||
           ( EQ(c, m) && cover_i < cover_j ) )
      {
        m = c;
        j = i;
      }
    }

    // In case that the found child is full, split it in two
    // nodes
    if ( node.child(j).n_entries() == M )
    {
      split_child(node, j);
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

    ASSERT( !parent_node.is_leaf(),
        "Invalid data structure of rtree.");

    // If "node" is the "i"-th child of "parent_node", then 
    // "parent_node.bbox(i)" is the bbox, that covers all entries
    // of "node"
    //
    std::size_t i = 0; 

    while ( &parent_node.child(i) != &node )
    {
      ++i;
      ASSERT( i < M, "Invalid data structure of R-tree");
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
    auto new_root = std::make_unique<Node>(RTreeNodeID++);

    (*new_root).add_child( root_ );
    root_ = std::move(new_root);
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

    parent_nodes.push_back(std::make_unique<Node>(RTreeNodeID++));

    // Performance might be better if we add from the back...
    for ( size_t i = 0; i < children.size(); ++i )
    {
      if ( parent_nodes.back().get()->n_entries() >= M-1 )
      {
        parent_nodes.push_back(
          std::make_unique<Node>(RTreeNodeID++)
        );
      }

      Node& cur_node = *parent_nodes.back().get();

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
* ostream
*********************************************************************/
template 
<
  typename    ObjectType, 
  std::size_t M, 
  typename    CoordType, 
  std::size_t Dim
>
std::ostream& operator<<(std::ostream& os, 
                         const RTreeND<ObjectType,M,CoordType,Dim>& tree)
{
  std::size_t height = tree.height();

  return tree.root().export_to_txt(os, height);
}


} // namespace CppUtils
