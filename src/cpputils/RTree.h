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
#include <stdexcept> // std::runtime_error
#include <iomanip>   // std::setprecision
#include <iostream>  // std::to_string

#include "Vec2.h"
#include "Geometry.h"
#include "MathUtility.h"
#include "Helpers.h"
#include "Log.h"

namespace CppUtils {


template <typename T, long unsigned int M>
class RTree;


/*********************************************************************
* Minimum bounding rectangle 
*********************************************************************/
class RTreeBBox
{
public:

  /*------------------------------------------------------------------ 
  | Constructors
  ------------------------------------------------------------------*/
  RTreeBBox() {}

  RTreeBBox(const Vec2d& ll, const Vec2d& ur)
  : lowleft_ { ll }
  , upright_ { ur }
  {
    area_ = (ur.x - ll.x) * (ur.y - ll.y);
    ASSERT( area_ > 0, "Invalid RTreeBBox definition." );
  }

  /*------------------------------------------------------------------ 
  | Getters
  ------------------------------------------------------------------*/
  double area() const { return area_; }

  Vec2d& lowleft() { return lowleft_; }
  const Vec2d& lowleft() const { return lowleft_; }

  Vec2d& upright() { return upright_; }
  const Vec2d& upright() const { return upright_; }

  /*------------------------------------------------------------------ 
  | Compute the intersection area between two RTreeBBoxes
  ------------------------------------------------------------------*/
  double intersection_area(const RTreeBBox& bbox) const
  {
    const Vec2d& a = lowleft_;
    const Vec2d& b = upright_;

    const Vec2d& p = bbox.lowleft();
    const Vec2d& q = bbox.upright();

    // BBox of intersection
    const Vec2d ll { MAX(a.x, p.x), MAX(a.y, p.y) };
    const Vec2d ur { MIN(b.x, q.x), MIN(b.y, q.y) };

    const double dx = (ur.x - ll.x);
    const double dy = (ur.y - ll.y);

    if (dx >= 0.0 && dy >= 0.0)
      return dx*dy;

    return 0.0;

  } // intersection_area()

  /*------------------------------------------------------------------ 
  | Compute the union area between two RTreeBBoxes
  ------------------------------------------------------------------*/
  double union_area(const RTreeBBox& bbox) const
  {
    return area_ + bbox.area() - intersection_area(bbox);
  }


  /*------------------------------------------------------------------ 
  | Compose an RTreeBBox that includes both this and another BBox
  ------------------------------------------------------------------*/
  RTreeBBox bounding_box(const RTreeBBox& bbox) const
  {
    const double xmin = MIN(lowleft_.x, bbox.lowleft_.x);
    const double ymin = MIN(lowleft_.y, bbox.lowleft_.y);

    const double xmax = MAX(upright_.x, bbox.upright_.x);
    const double ymax = MAX(upright_.y, bbox.upright_.y);

    return { {xmin,ymin}, {xmax,ymax} }; 
  }



private:

  Vec2d  lowleft_ { 0.0, 0.0 };
  Vec2d  upright_ { 0.0, 0.0 };
  double area_    { 0.0 };

}; // RTreeBBox


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
*
* The RTree node implementation
* -> M... minimum degree of the R-Tree
* -> T... object-type
*
*********************************************************************/
template <typename T, long unsigned int M>
class RTreeNode
{
public:

  friend RTree<T,M>;

  using BBoxes   = std::array<RTreeBBox,M>;
  using Children = std::array<std::unique_ptr<RTreeNode<T,M>>,M>;
  using Objects  = std::array<T*,M>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  RTreeNode() {}

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  BBoxes& bboxes() { return bboxes_; }
  const BBoxes& bboxes() const { return bboxes_; }

  Children& children() { return children_; }
  const Children& children() const { return children_; }

  Objects& objects() { return objects_; }
  const Objects& objects() const { return objects_; }

  RTreeNode<T,M>* parent() { return parent_; }
  const RTreeNode<T,M>* parent() const { return parent_; }

  bool is_leaf() const { return is_leaf_; }
  long unsigned int n_entries() const { return n_entries_; }

  /*------------------------------------------------------------------ 
  | Function used to estimate the tree height
  ------------------------------------------------------------------*/
  std::size_t add_height(std::size_t height) const
  {
    if ( !is_leaf_ )
    {
      ++height;
      height = (*child(0)).add_height(height);
    }

    return height;
  }

  /*------------------------------------------------------------------ 
  | Print out node data
  ------------------------------------------------------------------*/
  std::ostream& print(std::ostream& os,  
                      std::size_t   level=0, 
                      std::size_t   index=0,
                      std::size_t   parent_index=0) const
  {
    if ( !is_leaf_ )
    {
      for (std::size_t i = 0; i < n_entries_; ++i)
      {
        (*child(i)).print(os, level-1, i, index);
        os << "\n";
      }
    }

    os << n_entries_ << ", " 
       << level << ", " 
       << index << ", "
       << parent_index << "\n";

    for (std::size_t i = 0; i < n_entries_; ++i)
    {
      const Vec2d& ll = bbox(i).lowleft();
      const Vec2d& ur = bbox(i).upright();

      os << std::setprecision(5) << std::fixed 
         << ll.x << ", " << ll.y << ", " 
         << ur.x << ", " << ur.y;

      if ( i < n_entries_ - 1 )
        os << "\n";
    }

    return os;

  } // print()




  /*------------------------------------------------------------------ 
  | Access minimum bounding rectangles
  ------------------------------------------------------------------*/
  const RTreeBBox& bbox(std::size_t i) const
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to access key at position " 
        + std::to_string(i) );
    return bboxes_[i]; 
  }

  /*------------------------------------------------------------------ 
  | Access objects
  ------------------------------------------------------------------*/
  T* object(std::size_t i) 
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to access key at position " 
        + std::to_string(i) );
    return objects_[i]; 
  }
  const T* object(std::size_t i) const 
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to access key at position " 
        + std::to_string(i) );
    return objects_[i]; 
  }

  /*------------------------------------------------------------------ 
  | Access children
  ------------------------------------------------------------------*/
  std::unique_ptr<RTreeNode<T,M>>& child_ptr(std::size_t i) 
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to access child at position " 
        + std::to_string(i) );
    return children_[i]; 
  }

  const std::unique_ptr<RTreeNode<T,M>>& child_ptr(std::size_t i) const
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to access child at position " 
        + std::to_string(i) );
    return children_[i]; 
  }

  RTreeNode<T,M>* child(std::size_t i) 
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to access child at position " 
        + std::to_string(i) );
    return children_[i].get(); 
  }
  const RTreeNode<T,M>* child(std::size_t i) const 
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to access child at position " 
        + std::to_string(i) );
    return children_[i].get(); 
  }


private:

  /*------------------------------------------------------------------ 
  | Set bbox
  ------------------------------------------------------------------*/
  void bbox(std::size_t i, const RTreeBBox& b) 
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to access key at position " 
        + std::to_string(i) );
    bboxes_[i] = b; 
  }

  /*------------------------------------------------------------------ 
  | Set objects
  ------------------------------------------------------------------*/
  void object(std::size_t i, T* obj)
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to set object at position " 
        + std::to_string(i) );
    objects_[i] = obj;
  }

  /*------------------------------------------------------------------ 
  | Set children
  ------------------------------------------------------------------*/
  void child(std::size_t i, std::unique_ptr<RTreeNode<T,M>>& c)
  { 
    ASSERT( i < n_entries_, 
        "RTreeNode: Unable to set child at position " 
        + std::to_string(i) );
    children_[i] = std::move(c);
  }




  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  BBoxes            bboxes_    { };
  Children          children_  { nullptr };
  Objects           objects_   { nullptr };
  RTreeNode<T,M>*   parent_    { nullptr };

  bool              is_leaf_   { true };
  long unsigned int n_entries_ { 0 };

}; // RTreeNode


/*********************************************************************
* This class defines the interface to an R-tree structure 
* for 2D simplices
*********************************************************************/
template <typename T, long unsigned int M>
class RTree
{
public:

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  RTree()
  {
    root_ = std::make_unique<RTreeNode<T,M>>();
  }

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  RTreeNode<T,M>& root() { return *root_; }
  const RTreeNode<T,M>& root() const { return *root_; }

  /*------------------------------------------------------------------ 
  | Get the height of the entire tree
  ------------------------------------------------------------------*/
  std::size_t height() const 
  {
    std::size_t height = 0;

    height = (*root_).add_height(height);

    return height;

  } // height()

  /*------------------------------------------------------------------ 
  | Search the leaf node that contains a given object in the tree
  ------------------------------------------------------------------*/
  const RTreeNode<T,M>* search_node(const T& obj,
                                    const RTreeNode<T,M>& node) const
  {
  } // RTree::search_node()

  /*------------------------------------------------------------------ 
  | Insert a new object into the RTree structure
  ------------------------------------------------------------------*/
  void insert(T* object)
  {
    // Handle the case where the root node is full
    if ( root_->n_entries() == M )
    {
      auto new_root = std::make_unique<RTreeNode<T,M>>();
      new_root->is_leaf_   = false;
      new_root->n_entries_ = 1;
      new_root->child(0, root_);
      root_ = std::move(new_root);

      split_child(*root_, 0);
    }
     
    insert_nonfull(*root_, object);

  } // insert()

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

  } // write_to_file()

private:

  /*------------------------------------------------------------------ 
  | This function splits the i-th child of a given "parent_node"
  ------------------------------------------------------------------*/
  void split_child(RTreeNode<T,M>& parent_node, std::size_t i)
  {
    LOG(DEBUG) << "ENTRY FUNCTION \"split_child()\"";

    // Check that parent node is not full
    ASSERT(parent_node.n_entries() != M, "Invalid R-Tree structure.");

    // Check for valid child pointers
#ifndef NDEBUG
    for (std::size_t j = 0; j <= i; ++j)
      ASSERT(parent_node.child(j) != nullptr,
          "Invalid child pointer at position " + std::to_string(j));
#endif

    // This is the child node, whose entries will be splitted 
    RTreeNode<T,M>& child_node = *(parent_node.child(i));

    // Check that child node is full
    ASSERT(child_node.n_entries() == M, "Invalid R-Tree structure.");

    // This is the new node, which will get half of the entries 
    // of "child_node"
    auto new_node = std::make_unique<RTreeNode<T,M>>();
    (*new_node).is_leaf_ = child_node.is_leaf_;
    (*new_node).parent_  = &parent_node;

    // This array contains the information, which entries of the 
    // child node "child_node" will be added to "new_node"
    std::array<bool,M> add_to_new = quadratic_split( child_node );

#ifndef NDEBUG
    std::string log_str = "QUADRATIC SPLIT ARRAY: [";
    for ( std::size_t j = 0; j < M; ++j )
    {
      log_str += std::to_string(add_to_new[j]);
      if ( j < M-1 )
        log_str += ", ";
    }
    LOG(DEBUG) << log_str << "]";
#endif

    // Distribute entries from "child_node" to "new_node"
    for ( std::size_t j = 0; j < M; ++j )
    {
      if ( !add_to_new[j] )
        continue;

      size_t n = (*new_node).n_entries_;
      ++(*new_node).n_entries_;
      (*new_node).object( n, child_node.object(j) );
      (*new_node).bbox( n, child_node.bbox(j) );
      (*new_node).child( n, child_node.child_ptr(j) );
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

      --child_node.n_entries_;
    }

    // Add "new_node" and its bounding boxx to "parent_node"
    size_t n = parent_node.n_entries_;
    ++parent_node.n_entries_;
    parent_node.child( n, new_node );

    // Compute the covering bboxes for all entries of "parent node"
    for (std::size_t j = 0; j < parent_node.n_entries(); ++j)
    {
      const RTreeNode<T,M>& child = *parent_node.child(j);

      RTreeBBox cover = child.bbox(0);

      for ( std::size_t k = 1; k < child.n_entries(); ++k )
        cover = cover.bounding_box( child.bbox(k) );

      parent_node.bbox( j, cover );
    }

#ifndef NDEBUG
    LOG(DEBUG) << "CURRENT TREE STRUCTURE: ";
    LOG(DEBUG) << "\n" << (*this) << "\n";
#endif

  } // split_child()

  /*------------------------------------------------------------------ 
  | This function is called upon the splitting of the entries of 
  | a given "node" in two sets "n1" and "n2".
  | The output array "add_to_n2" marks all entries of the input node, 
  | wether they will be distributed to the set "n1" (false) or 
  | to the set "n2" (true).
  ------------------------------------------------------------------*/
  std::array<bool,M> quadratic_split(const RTreeNode<T,M>& node) const
  {
    std::array<bool,M> add_to_n2 { false };

    // This array marks all distributed entries
    std::array<bool,M> distributed { false };

    // This variable referes to the remaining number of entries, which 
    // have not yet been distributed
    std::size_t n_remaining = node.n_entries();

    // These are the bounding boxes that enclose all objects in the 
    // two sets
    RTreeBBox bbox_1 {};
    RTreeBBox bbox_2 {};

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
        ASSERT( n_n2 >= M/2, "Error in function RTree::quadratic_split");
        return add_to_n2;
      }

      // Second node has so few entries, that all remaining entries must
      // be assigned to it in order to have the minimum number
      if ( n_n2 + n_remaining == M/2 )
      {
        for ( std::size_t i = 0; i < M; ++i )
          add_to_n2[i] = (distributed[i] == false) ? true : add_to_n2[i];

        ASSERT( n_n1 >= M/2, "Error in function RTree::quadratic_split");
        return add_to_n2;
      }

      // Pick the next entry to assign. This function also 
      // decrements "n_remaining" and handles the markation
      // of entry "i" in the array "distributed"
      std::size_t i = pick_next(node, distributed, n_remaining, 
                                bbox_1, bbox_2);

      // Add entry "i" to the group whose covering rectangle will
      // have to be enlarged least to accomodate it.
      const RTreeBBox& E_i = node.bbox(i);

      const RTreeBBox cover_1 = bbox_1.bounding_box(E_i);
      const RTreeBBox cover_2 = bbox_2.bounding_box(E_i);
      
      const double d1 = cover_1.area() - bbox_1.area(); 
      const double d2 = cover_2.area() - bbox_2.area(); 

      bool add_entry_to_n2 = false;

      // Add entry to set "n2", if its covering rectangle will be 
      // enlarged less than for set "n1"
      if ( d2 < d1 )
        add_entry_to_n2 = true;

      // In case of ties 
      if ( EQ(d1, d2) ) 
      {
        // Add entry to set with smaller area
        if ( bbox_2.area() < bbox_1.area() )
          add_entry_to_n2 = true;

        // Then to the one with fewer entries
        if ( EQ(bbox_1.area(), bbox_2.area()) && (n_n2 < n_n1) )
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

  } // quadratic_split()

  /*------------------------------------------------------------------ 
  | Select two entries of the given full "node" which are used as 
  | seeds for the quadratic_split() algorithm during a node splitting
  | operation.
  ------------------------------------------------------------------*/
  std::size_t pick_seeds(const RTreeNode<T,M>& node, 
                         std::array<bool,M>    distributed,
                         std::size_t&          n_remaining,
                         RTreeBBox&            bbox_1,
                         RTreeBBox&            bbox_2) const 
  {
    ASSERT(node.n_entries() == M, "Invalid data structure passed to "
        "function RTree::pick_seeds()");

    std::size_t seed_1 = -1;
    std::size_t seed_2 = -1;

    double ineff = 0.0;

    // Calculate inefficiency of grouping entries together
    for (std::size_t i = 0; i < M; ++i)
    {
      const RTreeBBox& E_i = node.bbox(i);

      for (std::size_t j = 0; j < M; ++j)
      {
        if ( i == j )
          continue;

        // Compose a rectangle J including E_i and E_j 
        const RTreeBBox& E_j = node.bbox(j);

        const RTreeBBox J = E_i.bounding_box(E_j);

        // Compute the inefficiency area
        const double diff = J.area() - E_i.area() - E_j.area();

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

  } // pick_seeds()

  /*------------------------------------------------------------------ 
  | 
  ------------------------------------------------------------------*/
  std::size_t pick_next(const RTreeNode<T,M>& node, 
                        std::array<bool,M>    distributed,
                        std::size_t&          n_remaining,
                        RTreeBBox&            bbox_1,
                        RTreeBBox&            bbox_2) const 
  {
    std::size_t entry = 0;

    double max_diff = 0.0;

    // Determine the cost of putting each entry in each group
    for (std::size_t i = 0; i < M; ++i)
    {
      if ( distributed[i] )
        continue;

      const RTreeBBox& E_i = node.bbox(i);

      // Compute the area increases required in the covering rectangles
      // "bbox_1" and "bbox_2" to include entry "i"
      RTreeBBox C_1 = bbox_1.bounding_box(E_i);
      RTreeBBox C_2 = bbox_2.bounding_box(E_i);

      const double d1 = C_1.area() - bbox_1.area();
      const double d2 = C_2.area() - bbox_2.area();

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

  } // pick_next()


  /*------------------------------------------------------------------ 
  | Insert a new object into the RTree structure
  ------------------------------------------------------------------*/
  void insert_nonfull(RTreeNode<T,M>& node, T* object)
  {
    LOG(DEBUG) << "ENTRY FUNCTION \"insert_nonfull()\"";

    RTreeBBox& bb_obj = object->bbox();

    // Choose an appropriate leaf to insert the object
    RTreeNode<T,M>& leaf = choose_leaf_insertion(node, bb_obj);

    // If the leaf has enough space to store the object, add it
    if ( leaf.n_entries() < M )
    {
      std::size_t i = leaf.n_entries();
      
      ++leaf.n_entries_;
      leaf.bbox(i, bb_obj); 
      leaf.object(i, object);

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

  } // insert_nonfull()

  /*------------------------------------------------------------------ 
  | Choose appropriate leaf node for RTree-insertion
  ------------------------------------------------------------------*/
  RTreeNode<T,M>& choose_leaf_insertion(RTreeNode<T,M>&  node,
                                        const RTreeBBox& object_bbox)
  {
    LOG(DEBUG) << "ENTRY FUNCTION \"choose_leaf_insertion()\"";

    // Choos only leaf nodes
    if ( node.is_leaf() )
      return node;

    double m = CPPUTILS_MAX;
    std::size_t j = 0;

    double cover_j = node.bbox(j).bounding_box(object_bbox).area();

    // Find the child-node that has the least enlargement with the 
    // object's bbox
    for ( std::size_t i = 0; i < node.n_entries(); ++i )
    {
      const RTreeBBox& child_bbox = node.bbox(i);

      // Compute enlargement
      const double c = child_bbox.union_area(object_bbox)
                     - child_bbox.area();

      const double cover_i = node.bbox(i).bounding_box(object_bbox).area();

      // Choose the entry, that needs the least enlargement to 
      // include the object's bbox
      // In case of ties, use the entry with the smalles area
      if ( ( c < m ) ||
           ( EQ(c, m) && cover_i < cover_j ) )
      {
        m = c;
        j = i;
      }
    }

    // In case that the found child is full, split it in two
    // nodes
    if ( (*node.child(j)).n_entries() == M )
    {
      split_child(node, j);
    }

    return choose_leaf_insertion( *node.child(j), object_bbox );

  } // choose_leaf_insertion()

  /*------------------------------------------------------------------ 
  | Update bounding box
  ------------------------------------------------------------------*/
  void update_parent_bbox(RTreeNode<T,M>& node)
  {
    // Stop at root node
    if ( node.parent() == nullptr )
      return;

    RTreeNode<T,M>& parent_node = (*node.parent());

    ASSERT( !parent_node.is_leaf(),
        "Invalid data structure of rtree.");

    // Loop over all children of "parent_node", in order to estimate
    // each bounding box
    for (std::size_t i = 0; i < parent_node.n_entries(); ++i)
    {
      const RTreeNode<T,M>& child = *parent_node.child(i);
      RTreeBBox cover = child.bbox(0);

      for ( std::size_t j = 1; j < child.n_entries(); ++j )
        cover = cover.bounding_box( child.bbox(j) );

      parent_node.bbox( i, cover );
    }

    return update_parent_bbox(parent_node);

  } // update_parent_bbox()


  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  std::unique_ptr<RTreeNode<T,M>> root_ {nullptr};

}; // RTree


/*********************************************************************
* Stream to std::cout
*********************************************************************/
template <typename T, long unsigned int M>
std::ostream& operator<<(std::ostream& os, 
                         const RTree<T,M>& tree)
{
  std::size_t height = tree.height();

  return tree.root().print(os, height);
}

} // CppUtils
