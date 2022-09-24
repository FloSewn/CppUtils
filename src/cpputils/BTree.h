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

#include "Helpers.h"
#include "Log.h"

namespace CppUtils {

template <typename T, long unsigned int M>
class BTree;

/*********************************************************************
* References
* ----------
* - https://www.cs.yale.edu/homes/aspnes/pinewiki/BTrees.html
* - https://github.com/solangii/b-plus-tree
*
* Some definitions on B-Trees:
* ----------------------------
* - M is the minimum degree of the B-Tree
* - Every node must have at least (M-1) keys
* - Every internal node other than the root has at least M children
* - If the tree is nonempty, the root must have at least one key
* - Every node may contain at most (2M-1) keys
* - An internal node may have at most 2M children
* - A node is full if it contains exactly (2M-1) keys
*
*        Keys:       14   38   57
*    Children:     A    B    C    D
*
*    A... holding values < 14
*    B... holding values >= 14 && < 38
*    C... holding values >= 38 && < 57
*    D... holding values >= 57 
*
*   
*********************************************************************/

/*********************************************************************
* The BTree node implementation
* -> T... key-type
* -> M... order of the B-Tree
*********************************************************************/
template <typename T, long unsigned int M>
class BTreeNode
{
public:

  friend BTree<T,M>;

  using Keys     = std::array<T*,2*M-1>;
  using Children = std::array<std::unique_ptr<BTreeNode<T,M>>,2*M>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  BTreeNode() {}

  /*------------------------------------------------------------------ 
  | Copy constructor
  ------------------------------------------------------------------*/
  BTreeNode(const BTreeNode& b)
  : keys_     { b.keys_ }
  , children_ { b.children_ }
  , parent_   { b.parent_ }
  , is_leaf_  { b.is_leaf_ }
  , n_keys_   { b.n_keys_ }
  {}

  /*------------------------------------------------------------------ 
  | Move constructor
  ------------------------------------------------------------------*/
  BTreeNode(const BTreeNode&& b)
  : keys_     { std::move(b.keys_) }
  , children_ { std::move(b.children_) }
  , parent_   { b.parent_ }
  , is_leaf_  { b.is_leaf_ }
  , n_keys_   { b.n_keys_ }
  {}

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  Keys& keys() { return keys_; }
  const Keys& keys() const { return keys_; }

  T* key(std::size_t i) { return keys_[i]; }
  const T* key(std::size_t i) const { return keys_[i]; }

  Children& children() { return children_; }
  const Children& children() const { return children_; }

  BTreeNode<T,M>* child(std::size_t i) { return children_[i].get(); }
  const BTreeNode<T,M>* child(std::size_t i) const 
  { return children_[i].get(); }

  BTreeNode<T,M>* parent() { return parent_; }
  const BTreeNode<T,M>* parent() const { return parent_; }

  bool is_leaf() const { return is_leaf_; }

  long unsigned int n_keys() const { return n_keys_; }


private:
  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  Keys              keys_     { nullptr };
  Children          children_ { nullptr };
  BTreeNode<T,M>*   parent_   { nullptr };

  bool              is_leaf_  { true };
  long unsigned int n_keys_   { 0 };

}; // BTreeNode



/*********************************************************************
* This class defines the interface to an R-tree structure 
* for 2D simplices
*********************************************************************/
template <typename T, long unsigned int M>
class BTree
{
public:

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  BTree()
  {
    root_ = std::make_unique<BTreeNode<T,M>>();
  }

  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  BTreeNode<T,M>& root() { return *root_; }
  const BTreeNode<T,M>& root() const { return *root_; }

  /*------------------------------------------------------------------ 
  | Search
  ------------------------------------------------------------------*/
  const BTreeNode<T,M>* search(const T& key, const BTreeNode<T,M>& node)
  const
  {
    const BTreeNode<T,M>* cur = &node;

    std::size_t i = 1;

    while ( i < node.n_keys()  &&  key > *node.key(i-1) )
      ++i;

    if ( i <= node.n_keys()  &&  key == *node.key(i-1) )
      return cur;
    else if ( node.is_leaf() )
      return nullptr;
    else
    {
      const BTreeNode<T,M>* child = node.child(i-1);

      if ( !child )
        throw std::runtime_error(
          "Invalid BTree access. Data structure seems to be broken.");

      return search(key, *child);
    }

  } // search() 


  /*------------------------------------------------------------------ 
  | Insert a key into the BTree structure
  ------------------------------------------------------------------*/
  void insert(T* key)
  {
    // Handle the case where root node r is full
    // --> The root splits and a new node becomes the root 
    if ( root_->n_keys() == 2*M-1 )
    {
      auto new_node = std::make_unique<BTreeNode<T,M>>();

      new_node->is_leaf_ = false;
      new_node->n_keys_  = 0;

      new_node->children_[0] = std::move(root_);
      root_                  = std::move(new_node);

      split_child(*root_, 0);
      insert_nonfull(*root_, key);
    }
    else
      insert_nonfull(*root_, key);

  } // insert() */

  /*------------------------------------------------------------------ 
  | Insert nonfull
  ------------------------------------------------------------------*/
  void insert_nonfull(BTreeNode<T,M>& node, T* key)
  {
    auto i = node.n_keys();

    // A leaf node is given -> insert key 
    if ( node.is_leaf() )
    {
      // Shift larger keys to the right
      while ( i > 0 && (*key < *node.keys_[i-1]) )
      {
        node.keys_[i] = node.keys_[i-1];
        --i;
      }

      // Place key
      node.keys_[i] = key;
      ++node.n_keys_;
    }
    else 
    {
      // Got to the right key
      while ( i > 0 && (*key < *node.keys_[i-1]) )
        --i;

      // Point to child
      ++i;

      if ( node.children_[i-1]->n_keys() == 2*M-1 )
      {
        split_child( node, i-1 );

        if ( *key > *node.keys_[i-1] )
          ++i;
      }

      insert_nonfull(*node.children_[i-1], key);
    }

  } // insert_nonfull() 


  /*------------------------------------------------------------------ 
  | Split a full node (with 2*M-1 keys) around its median key 
  | into two nodes with (M-1) keys each.
  | The median key moves up into the node's parent. 
  |
  | Arguments:
  | ----------
  | parent_node... nonfull internal node (parent node)
  | i............. index, such that node.child(i) is a full child 
  |                of node
  |
  |  Example: M = 4,  i = 3  
  |  -------
  |  * Maximum number of keys: (2M-1) = 7        
  |  * Minimum number of keys:  (M-1) = 3    median key is put     
  |                                          into parent node       
  |           parent node                        v
  |        .---------------.             .---------------.
  |        | A B C D E . . |    --->     | A B C S D E . | 
  |        '---------------'             '---------------'
  |               |                             | |     
  |               v   Full child node           v |     
  |               .---------------.             .-------. old 
  |               | P Q R S T U V |             | P Q R | node
  |               '---------------'             '-------' 
  |                       ^                       |
  |                   median key                  v
  |                                               .-------. new  
  |                                               | T U V | node 
  |                                               '-------' 
  |
  ------------------------------------------------------------------*/
  void split_child(BTreeNode<T,M>& parent_node, std::size_t i)
  {
    // This is the new node. It gets the keys of the full node,
    // that are greated than the median
    auto new_node = std::make_unique<BTreeNode<T,M>>();

    // Child must be a full node
    BTreeNode<T,M>* child_node = parent_node.children_[i].get();
    ASSERT( child_node->n_keys() == (2*M-1), "Invalid BTree structure");

    new_node->is_leaf_ = child_node->is_leaf();
    new_node->n_keys_  = M-1;

    // Put all keys of the full child into the new node,
    // which are greater than the median (=M)
    for (std::size_t j = 0; j < M-1; ++j)
      new_node->keys_[j] = child_node->keys_[j+M];

    // If the child node is not a leaf, pass also its 
    // children to the new node
    if ( !child_node->is_leaf() )
      for (std::size_t j = 0; j < M; ++j) 
        new_node->children_[j] = std::move(child_node->children_[j+M]);

    // Update number of keys 
    child_node->n_keys_ = M-1;

    // new_node will be placed at parent_node.child(i+1),
    // so all other children at j > i must be shifted to the right
    for (std::size_t j = parent_node.n_keys()+1; j > i; --j)
      parent_node.children_[j] = std::move(parent_node.children_[j-1]);

    parent_node.children_[i] = std::move( new_node );

    // All following keys in parent_node must be shifted to the right
    for (std::size_t j = parent_node.n_keys()-1; j > i-1; --j)
      parent_node.keys_[j+1] = parent_node.keys_[j];

    parent_node.keys_[i] = child_node->keys_[M-1];

  } // split_child()


private:
  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  std::unique_ptr<BTreeNode<T,M>> root_ {nullptr};

}; // BTree


} // CppUtils

