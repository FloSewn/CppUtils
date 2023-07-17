/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <climits>      // CHAR_BIT
#include <list>      
#include <map>      
#include <vector>    
#include <array>     
#include <stdexcept> 
#include <iomanip>   
#include <iostream>  

#include "VecND.h"
#include "BBoxND.h"
#include "Helpers.h"
#include "Log.h"
#include "VtkIO.h"

namespace CppUtils {


// ObjectType.... Contained object
// M............. Max. element number 
// CoordType..... Coordinate type
// Dim........... Dimensions
    
#ifndef OCTREE_DEF
#define OCTREE_DEF          \
    typename    ObjectType, \
    std::size_t M,          \
    typename    CoordType,  \
    std::size_t Dim           
#endif

#ifndef OCTREE_ARG
#define OCTREE_ARG \
  ObjectType, M, CoordType, Dim 
#endif

/*********************************************************************
* Forward declarations 
*********************************************************************/
template<OCTREE_DEF>
class OcTreeNodeND;

template<OCTREE_DEF>
class OcTreeND;

/*********************************************************************
* 
*********************************************************************/
template<OCTREE_DEF>
class OcTreeNodeND
{
  // Number of bits used to store the tree in a space-filling curve
  static constexpr std::size_t bit_size_ 
  {sizeof(std::size_t) * CHAR_BIT};

public:

  /*------------------------------------------------------------------ 
  | Container for node entries
  ------------------------------------------------------------------*/
  struct Entry 
  {
    Entry(const ObjectType& obj, const VecND<CoordType,Dim>& pos)
    : object   { &obj }
    , position { pos }
    {}

    const ObjectType*    object   { nullptr };
    VecND<CoordType,Dim> position { };

  }; // Entry

  /*------------------------------------------------------------------ 
  | Typedefs
  ------------------------------------------------------------------*/
  using OcTree       = OcTreeND<OCTREE_ARG>;
  using BBox         = BBoxND<CoordType,Dim>;
  using Vec          = VecND<CoordType,Dim>;
  using Node         = OcTreeNodeND<OCTREE_ARG>;
  using Node_ptr     = std::unique_ptr<Node>;
  using Entries      = std::list<Entry>;
  using Children     = std::array<Node_ptr, (1<<Dim)>;
  using EntryVector  = std::vector<Entry>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  OcTreeNodeND(OcTree& tree, const BBox& bbox, CoordType min_scale, 
               std::size_t height=0, std::size_t local_index=0,
               Node* parent=nullptr)
  : tree_      { &tree }
  , bbox_      { bbox }
  , min_scale_ { min_scale }
  , height_    { height }
  , parent_    { parent }
  {
    center_ = ( bbox_.lowleft() + bbox_.upright() ) / two_;

    if ( height == 0 || !parent )
      return;

    // Obtain string identifier that defines the 
    // position on the space-filling curve
    std::bitset<bit_size_> bit_id 
      = std::bitset<bit_size_>( (*parent).curve_id() );
    bit_id <<= Dim;
    bit_id |= local_index;

    curve_id_ = bit_id.to_string( );

    // Remove leading bits
    std::size_t id_rem = curve_id_.length() - Dim * height_;
    curve_id_.erase(curve_id_.begin(), curve_id_.begin()+id_rem);

    // Add curve-position in tree map
    (*tree_).add_to_space_filling_curve( *this );
  }

  /*------------------------------------------------------------------ 
  | Destructor
  ------------------------------------------------------------------*/
  ~OcTreeNodeND()
  { 
    if ( height_ > 0 )
      (*tree_).remove_from_space_filling_curve( *this ); 
  }

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  const BBox& bbox() const { return bbox_; }
  std::size_t height() const { return height_; }
  Node* parent() const { return parent_; }

  const Vec& center() const { return center_; }
  bool is_split() const { return is_split_; }
  std::size_t n_total_entries() const { return n_tot_entries_; }
  std::size_t n_entries() const { return entries_.size(); }

  Entries& entries() { return entries_; }
  const Entries& entries() const { return entries_; }

  Children& children() { return children_; }
  const Children& children() const { return children_; }

  std::string curve_id() const { return curve_id_; }

  OcTree& tree()  { return *tree_; }
  const OcTree& tree() const { return *tree_; }

  Node& child(std::size_t i)
  {
    ASSERT( is_split(),
    "OcTreeNodeND::child(): "
    "Unable to access child of leaf node.");

    ASSERT( i < (1<<Dim), 
    "OcTreeNodeND::child(): " 
    "Unable to access child[" + std::to_string(i) + "].");

    return children_[i].child();
  }

  const Node& child(std::size_t i) const
  {
    ASSERT( is_split(),
    "OcTreeNodeND::child(): "
    "Unable to access child of leaf node.");

    ASSERT( i < (1<<Dim), 
    "OcTreeNodeND::child(): " 
    "Unable to access child[" + std::to_string(i) + "].");

    return *children_[i];
  }

  /*------------------------------------------------------------------ 
  | Get all items inside a given bounding box
  ------------------------------------------------------------------*/
  void query(const BBox& bbox, EntryVector& found_objects) const
  {
    // Skip empty leaf nodes
    if ( !is_split() && n_entries() < 1 )
      return;

    // Skip if node does not overlap with query bbox
    if ( !bbox_.bbox_intersect_touch( bbox ) )
      return;

    // Query child nodes
    if ( is_split() )
    {
      for ( auto& child : children_ )
        (*child).query(bbox, found_objects);

      return;
    }

    // Check for objects in this node
    for ( auto& entry : entries_ )
      if ( bbox.point_inside_touch( entry.position ) )
        found_objects.push_back( entry );

    return;
    
  } // query()

  /*------------------------------------------------------------------ 
  | Get all items inside a sphere / circle of given center and radius
  ------------------------------------------------------------------*/
  void query(const Vec& center, CoordType radius, 
             EntryVector& found_objects) const
  {
    // Skip empty leaf nodes
    if ( !is_split() && n_entries() < 1 )
      return;

    // Skip if node does not overlap with query bbox
    if ( !bbox_.sphere_intersect( center, radius ) )
      return;

    // Query child nodes
    if ( is_split() )
    {
      for ( auto& child : children_ )
        (*child).query(center, radius, found_objects);

      return;
    }

    // Check for objects in this node
    for ( auto& entry : entries_ )
      if ( bbox.point_inside_touch( entry.position ) )
        found_objects.push_back( entry );

  } // query()

  /*------------------------------------------------------------------ 
  | Function used to estimate the tree height
  ------------------------------------------------------------------*/
  std::size_t tree_height() const
  {
    if ( is_split() )
    {
      std::size_t h = height();

      for ( auto& child : children_ )
        h = std::max(h, (*child).tree_height());

      return h;
    }

    return height();
  }

  /*------------------------------------------------------------------ 
  | Add a new item to the tree
  ------------------------------------------------------------------*/
  bool insert(const ObjectType& object, const Vec& obj_pos)
  {
    if ( !bbox_.point_inside_touch( obj_pos ) )
      return false;

    if ( is_split() )
      return pass_to_children( object, obj_pos );

    entries_.push_back( { object, obj_pos } );
    change_total_number_of_entries( 1 );

    const Vec delta = bbox_.upright() - bbox_.lowleft();
    const CoordType scale_sqr = delta.norm_sqr();
    const CoordType child_scale_sqr = scale_sqr / four_;

    bool splitable = child_scale_sqr < scale_sqr;
    splitable     &= child_scale_sqr > min_scale_ * min_scale_;
    splitable     &= height_ < (bit_size_ / Dim);

    if ( splitable && n_total_entries() > M )
      split();

    return true;

  } // insert()

  /*------------------------------------------------------------------ 
  | Remove an item from the tree
  ------------------------------------------------------------------*/
  bool remove(const ObjectType& object, const Vec& obj_pos)
  {
    if ( !bbox_.point_inside_touch( obj_pos ) )
      return false;

    // Remove object from children
    if ( is_split_ )
    {
      std::size_t i_child = locate_child(obj_pos);

      bool removed = (*children_[i_child]).remove(object, obj_pos);

      ASSERT( removed,
      "OcTreeNodeND::remove(): "
      "Failed to remove object from child node");
      (void) removed;

      if (n_total_entries() <= M)
      {
        bool merged = merge();

        ASSERT( merged,
        "OcTreeNodeND::remove(): "
        "Failed to merge children");
        (void) merged;
      }

      return true;
    }

    // Remove object from this node
    std::size_t old_size = entries_.size();
    auto it = entries_.begin();

    while ( it != entries_.end() )
    {
      if ( (*it).object == &object )
      {
        entries_.erase( it );
        break;
      }

      ++it;
    }

    std::size_t new_size = entries_.size();

    change_total_number_of_entries( -1 );

    return (new_size < old_size);

  } // remove()

private:

  /*------------------------------------------------------------------ 
  | Locate the quadrant of a given input position and return the 
  | respective child node index of this quadrant
  ------------------------------------------------------------------*/
  std::size_t locate_child(const Vec& pos)
  {
    ASSERT( is_split(),
    "OcTreeNodeND::get_child(): "
    "Can not locate any child nodes of non-split octree node.");

    const Vec orient_pos = pos - center_;

    const CoordType zero = {};

    std::size_t i_child = (orient_pos[0] >= zero);

    for ( std::size_t i = 1; i < Dim; ++i )
      i_child |= (orient_pos[i] >= zero) << i; 

    return i_child;

  } // locate_child()

  /*------------------------------------------------------------------ 
  | Pass object to child nodes
  ------------------------------------------------------------------*/
  bool pass_to_children(const ObjectType& object, const Vec& obj_pos)
  {
    std::size_t i_child = locate_child(obj_pos);

    if ( (*children_[i_child]).insert(object, obj_pos) )
      return true;

    return false;

  } // pass_to_children()

  /*------------------------------------------------------------------ 
  | Merge child nodes 
  ------------------------------------------------------------------*/
  bool merge()
  {
#ifndef NDEBUG 
    for ( std::size_t i = 0; i < (1<<Dim); ++i )
      ASSERT( !(*children_[i]).is_split(),
      "OcTreeNodeND::merge(): "
      "Failed to merge children - child " + std::to_string(i) + " "
      "is not a leaf node.");
#endif 

    for ( auto& child : children_ )
    {
      for ( auto& entry : (*child).entries() )
        entries_.push_back( entry );
      child = nullptr;
    }

    is_split_ = false;

    return true;

  } // merge()


  /*------------------------------------------------------------------ 
  | Split this node
  ------------------------------------------------------------------*/
  void split()
  {
    std::size_t child_height { height_ + 1 };

    const Vec delta = (bbox_.upright() - bbox_.lowleft()) / two_;
    std::size_t nv  = (1<<Dim);

    for (std::size_t i=0; i < nv; ++i)
    {
      // Create a bitset to access all bbox vertices, 
      // e.g. for Dim==2: { (0,0), (0,1), (1,0), (1,1) }
      std::bitset<Dim> bits = std::bitset<Dim>(i);

      // Fill bitset into a Vec
      Vec dir {};
      for (std::size_t j=0; j < Dim; ++j)
        dir[j] = bits[j];

      // Bounding box of new child node
      const Vec lowleft = bbox_.lowleft() + dir * delta;
      const Vec upright = lowleft + delta;

      const BBox bb { lowleft, upright };

      children_[i] = std::make_unique<Node>(*tree_, bb, min_scale_, 
                                            child_height, i, this);
    }

    // Distribute entries among children
    is_split_ = true;
    distribute_entries();
    
    ASSERT( entries_.size() == 0,
    "OcTreeNodeND::split(): "
    "Distribution of node entries failed.");

  } // split()

  /*------------------------------------------------------------------ 
  | Distribute the entries of the node to its children
  ------------------------------------------------------------------*/
  void distribute_entries()
  {
    while ( entries_.size() > 0 )
    {
      Entry& entry = entries_.back();

      std::size_t i_child = locate_child(entry.position);

      (*children_[i_child]).insert(*entry.object, entry.position);

      change_total_number_of_entries( -1 );
      entries_.pop_back();
    }

  } // distribute_entries() 

  /*------------------------------------------------------------------ 
  | Increment the total number of entries in this node and all
  | its parent nodes
  ------------------------------------------------------------------*/
  void change_total_number_of_entries(std::size_t i)
  {
    n_tot_entries_ += i;
    if ( parent_ )
      (*parent_).change_total_number_of_entries( i );
  }

  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  OcTree*     tree_          { nullptr };
  BBox        bbox_          {};
  CoordType   min_scale_     {};
  std::size_t height_        {};
  Node*       parent_        { nullptr };

  Vec         center_        {};
  bool        is_split_      { false };
  std::size_t n_tot_entries_ { 0 };

  Entries     entries_       { };     
  Children    children_      { nullptr };

  std::string curve_id_      {""};

  CoordType   two_           { static_cast<CoordType>(2) };
  CoordType   four_          { static_cast<CoordType>(4) };

}; // OcTreeNodeND



/*********************************************************************
* 
*********************************************************************/
template<OCTREE_DEF>
class OcTreeND
{
  friend OcTreeNodeND<OCTREE_ARG>;

public:

  using BBox        = BBoxND<CoordType,Dim>;
  using Vec         = VecND<CoordType,Dim>;
  using Node        = OcTreeNodeND<OCTREE_ARG>;
  using Node_ptr    = std::unique_ptr<Node>;
  using NodeMap     = std::map<std::string,Node*>;

  using Entry       = typename Node::Entry;
  using EntryVector = std::vector<Entry>;

  static inline CoordType minimum_scale 
    = std::numeric_limits<CoordType>::epsilon();

  /*------------------------------------------------------------------ 
  | Iterator implementation
  | -> This iterator loops over all child nodes that carry obects
  ------------------------------------------------------------------*/
  struct Iterator
  {
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = Node;
    using pointer           = Node*;
    using reference         = Node&;

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Constructor
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Iterator() {}
    Iterator(const NodeMap& curve)
    {
      auto  iter  = curve.begin();
      Node* start = iter->second;
      while ( iter != curve.end() && (*start).is_split() )
      {
        ++iter; 
        start = iter->second;
      }
      cur_node_ = start;
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Operators
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    reference operator*() const { return *cur_node_; }
    pointer operator->() { return cur_node_; }

    Iterator& operator++() 
    { 
      const auto& curve = (*cur_node_).tree().curve();
      auto iter = curve.find( (*cur_node_).curve_id() );

      // Traverse only nodes that contain objects
      do 
      {
        ++iter;

        if ( iter == curve.end() )
          cur_node_ = nullptr;
        else
          cur_node_ = iter->second;

      } while( cur_node_ && (*cur_node_).is_split() );

      return *this;
    }

    Iterator operator++(int) 
    { Iterator tmp = *this; ++(*this); return tmp; }

    friend bool operator== (const Iterator& a, const Iterator& b) 
    { return (   a.cur_node_ == b.cur_node_ ); }

    friend bool operator!= (const Iterator& a, const Iterator& b) 
    { return !(a == b); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Attributes
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  private:
    Node* cur_node_ { nullptr };

  }; // Iterator

  /*------------------------------------------------------------------ 
  | ConstantIterator implementation
  ------------------------------------------------------------------*/
  struct ConstantIterator
  {
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = Node;
    using pointer           = Node*;
    using reference         = Node&;

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Constructor
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ConstantIterator() {}
    ConstantIterator(const NodeMap& curve)
    {
      auto  iter  = curve.begin();
      Node* start = iter->second;
      while ( iter != curve.end() && (*start).is_split() )
      {
        ++iter; 
        start = iter->second;
      }
      cur_node_ = start;
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Operators
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    const reference operator*() const { return *cur_node_; }
    const pointer operator->() { return cur_node_; }

    ConstantIterator& operator++() 
    { 
      const auto& curve = (*cur_node_).tree().curve();
      auto iter = curve.find( (*cur_node_).curve_id() );

      // Traverse only nodes that contain objects
      do 
      {
        ++iter;

        if ( iter == curve.end() )
          cur_node_ = nullptr;
        else
          cur_node_ = iter->second;

      } while( cur_node_ && (*cur_node_).is_split() );

      return *this;
    }

    ConstantIterator operator++(int) 
    { ConstantIterator tmp = *this; ++(*this); return tmp; }

    friend bool operator== 
    (const ConstantIterator& a, const ConstantIterator& b) 
    { return (   a.cur_node_ == b.cur_node_ ); }

    friend bool operator!= 
    (const ConstantIterator& a, const ConstantIterator& b) 
    { return !(a == b); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    | Attributes
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  private:
    Node* cur_node_ { nullptr };

  }; // ConstantIterator

  Iterator begin() const { return Iterator( curve_ ); }
  Iterator end() const { return Iterator(); } 

  ConstantIterator cbegin() const { return ConstantIterator( curve_ ); }
  ConstantIterator cend() const { return ConstantIterator(); } 

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  OcTreeND(const BBox& bbox, 
           CoordType s=std::numeric_limits<CoordType>::epsilon())
  : min_scale_ { s }
  { root_ = std::make_unique<Node>(*this, bbox, min_scale_); }

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  Node& root() { return *root_; }
  const Node& root() const { return *root_; }

  NodeMap& curve() { return curve_; }
  const NodeMap& curve() const { return curve_; }

  /*------------------------------------------------------------------ 
  | Query for objects inside a given bounding box
  ------------------------------------------------------------------*/
  EntryVector query(const BBox& bbox) const
  {
    EntryVector output {};

    (*root_).query( bbox, output );

    return std::move( output );
  }

  /*------------------------------------------------------------------ 
  | Query for objects inside a given sphere
  ------------------------------------------------------------------*/
  EntryVector query(const Vec& center, CoordType radius) const
  {
    EntryVector output {};

    (*root_).query( center, radius );

    return std::move( output );
  }

  /*------------------------------------------------------------------ 
  | Insert element to OcTreeND
  ------------------------------------------------------------------*/
  bool insert(const ObjectType& object, const Vec& obj_pos)
  { return (*root_).insert( object, obj_pos ); } 

  /*------------------------------------------------------------------ 
  | Remove element from OcTreeND
  ------------------------------------------------------------------*/
  bool remove(const ObjectType& object, const Vec& obj_pos)
  { return (*root_).remove( object, obj_pos ); }

  /*------------------------------------------------------------------ 
  | Get the height of the entire tree
  ------------------------------------------------------------------*/
  std::size_t height() const 
  { return (*root_).tree_height(); }

protected:
  /*------------------------------------------------------------------ 
  | Add a new node to the space-filling curve of the tree
  ------------------------------------------------------------------*/
  void add_to_space_filling_curve(Node& node)
  { 
    ASSERT( curve_.count(node.curve_id()) == 0,
    "OcTreeND::add_to_space_filling_curve(): "
    "Node-ID already contained in space-filling curve");
    this->curve_[node.curve_id()] = &node;
  }

  /*------------------------------------------------------------------ 
  | Remove node from the space-filling curve of the tree
  ------------------------------------------------------------------*/
  void remove_from_space_filling_curve(Node& node)
  {
    ASSERT( curve_.count(node.curve_id()) > 0,
    "OcTreeND::remove_from_space_filling_curve(): "
    "Node-ID not contained in space-filling curve");
    auto it = curve_.find( node.curve_id() );
    curve_.erase( it );
  }

private:
  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  CoordType min_scale_ { };
  NodeMap   curve_ {};
  Node_ptr  root_ {nullptr};

}; // OcTreeND


/*********************************************************************
* This class is used to export the OcTreeND structure
*********************************************************************/
template<OCTREE_DEF>
class OcTreeNDWriter
{
public:
  using Node   = OcTreeNodeND<OCTREE_ARG>;
  using OcTree = OcTreeND<OCTREE_ARG>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  OcTreeNDWriter(const OcTree& octree)
  : octree_ { &octree }
  { }

  /*------------------------------------------------------------------ 
  | Write the OcTree data to a VTU file
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

    std::size_t tree_height = (*octree_).height();

    write_vtu_data((*octree_).root(), points, connectivity, 
                   offsets, types, heights, tree_height);

    VtuWriter writer { points, connectivity, offsets, types };
    writer.add_cell_data( heights, "height", 1);
    writer.write( file_name );

  } // write_to_vtu()

private:
  /*------------------------------------------------------------------ 
  | Write the data of a octree node to a vtu file
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
    if ( node.is_split() )
    {
      for (std::size_t i = 0; i < (1<<Dim); ++i)
        cur_offset = write_vtu_data(node.child(i), points, connectivity, 
                                    offsets, types, heights, 
                                    cur_height-1, cur_offset);
    }

    std::size_t n_verts = points.size() / 3;

    auto vertices = node.bbox().vertices();

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

    return cur_offset;

  } // write_vtu_data()

  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  const OcTree* octree_;

  // Array to map from BBoxND vertex coordinates
  // to VTK hexahedral coordinates
  const std::array<std::size_t,8> vtk_conn_map_ 
  { 0, 1, 2, 3, 7, 6, 5, 4 };


}; // OcTreeNDWriter


} // namespace CppUtils
