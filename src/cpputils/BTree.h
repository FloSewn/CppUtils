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
* The BTree node implementation
* -> T... key-type
* -> M... minimum degree of the B-Tree
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

  Children& children() { return children_; }
  const Children& children() const { return children_; }

  BTreeNode<T,M>* parent() { return parent_; }
  const BTreeNode<T,M>* parent() const { return parent_; }

  bool is_leaf() const { return is_leaf_; }

  long unsigned int n_keys() const { return n_keys_; }


  /*------------------------------------------------------------------ 
  | Access keys
  ------------------------------------------------------------------*/
  T* key(std::size_t i) 
  { 
    ASSERT( i < n_keys_, "Invalid BTree access.");
    return keys_[i]; 
  }
  const T* key(std::size_t i) const 
  { 
    ASSERT( i < n_keys_, "Invalid BTree access.");
    return keys_[i]; 
  }

  /*------------------------------------------------------------------ 
  | Access children
  ------------------------------------------------------------------*/
  std::unique_ptr<BTreeNode<T,M>>& child_ptr(std::size_t i) 
  { 
    ASSERT( i <= n_keys_, "Invalid BTree access.");
    return children_[i]; 
  }

  const std::unique_ptr<BTreeNode<T,M>>& child_ptr(std::size_t i) const
  { 
    ASSERT( i <= n_keys_, "Invalid BTree access.");
    return children_[i]; 
  }

  BTreeNode<T,M>* child(std::size_t i) 
  { 
    ASSERT( i <= n_keys_, "Invalid BTree access.");
    return children_[i].get(); 
  }
  const BTreeNode<T,M>* child(std::size_t i) const 
  { 
    ASSERT( i <= n_keys_, "Invalid BTree access.");
    return children_[i].get(); 
  }


private:
  /*------------------------------------------------------------------ 
  | Set keys
  ------------------------------------------------------------------*/
  void key(std::size_t i, T* k) 
  { 
    if ( i >= n_keys_ )
      throw std::runtime_error("Invalid BTree access.");
    keys_[i] = k;
  }

  /*------------------------------------------------------------------ 
  | Set children
  ------------------------------------------------------------------*/
  void child(std::size_t i, std::unique_ptr<BTreeNode<T,M>>& c) 
  { 
    if ( i > n_keys_ )
      throw std::runtime_error("Invalid BTree access.");
    children_[i] = std::move(c);
  }



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
  | Search for a key in the tree
  ------------------------------------------------------------------*/
  const BTreeNode<T,M>* search(const T& key, 
                               const BTreeNode<T,M>& node) const
  {
    const BTreeNode<T,M>* cur = &node;

    if ( node.n_keys() == 0 ) 
      return nullptr;

    std::size_t i = 0;

    while ( i < node.n_keys()  &&  key > *node.key(i) )
      ++i;

    if ( i < node.n_keys()  &&  key == *node.key(i) )
    {
      return cur;
    }
    else if ( node.is_leaf() )
    {
      return nullptr;
    }
    else
    {
      const BTreeNode<T,M>* child = node.child(i);

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
      // This node will be the new root
      auto new_node = std::make_unique<BTreeNode<T,M>>();

      // new_node gets its first key in split_child()
      new_node->is_leaf_ = false;
      new_node->n_keys_  = 0;

      new_node->child(0, root_);
      root_ = std::move(new_node);

      // The old root node will now be divided in two nodes with 
      // (M-1) keys each and move the remaining median key into 
      // the new root node
      split_child(*root_, 0);

      ASSERT(root_->n_keys_ == 1, 
             "Invalid number of keys in root node." );

      // Next, insert the key 
      insert_nonfull(*root_, key);
    }
    else
      insert_nonfull(*root_, key);

  } // insert() */

  /*------------------------------------------------------------------ 
  | Insert a key into a nonfull node
  |
  | Arguments:
  | ----------
  | node.... nonfull node (branch or leaf)
  | key..... the key that will be inserted to the tree
  |
  | 
  |  Example: M = 4,  node is leaf,  input key = 3
  |  -------
  |  
  |  :                       :
  |  |                       |
  |  .---------------.       .---------------.
  |  | 1 4 5 7 . . . |  -->  | 1 3 4 5 7 . . |
  |  '---------------'       '---------------'
  |  
  |  
  |  Example: M = 4,  node is branch,  input key = 15
  |  -------
  |                            1) Pick child C (>=10 & <23)
  |  :                  
  |  |                  
  |  .----------------------.  
  |  |  3 10 23 43 60  .  . |  
  |  '----------------------'  
  |   |  |  |  |  |
  |   |  |  |  |  |
  |   A  B  |  D  E
  |         | 
  |         .----------------------.
  |       C | 11 13 14 17 19 20 22 |
  |         '----------------------'
  |
  |  Step 1) Split C & add new node N
  |  --------------------------------
  |  :                  
  |  |                  
  |  .----------------------.  
  |  |  3 10 17 23 43 60  . |  
  |  '----------------------'  
  |   |  |  |  |  |  |
  |   |  |  |  |  |  |
  |   A  B  |  |  D  E
  |         |  |
  |         .----------------------.
  |       C | 11 13 14  .  .  .  . |
  |         '----------------------'
  |            |
  |            |
  |            .----------------------.
  |          N | 19 20 22  .  .  .  . |
  |            '----------------------'
  |
  |  Step 2) Add input key 15 to node C
  |  ----------------------------------
  |  :                  
  |  |                  
  |  .----------------------.  
  |  |  3 10 17 23 43 60  . |  
  |  '----------------------'  
  |   |  |  |  |  |  |
  |   |  |  |  |  |  |
  |   A  B  |  |  D  E
  |         |  |
  |         .----------------------.
  |       C | 11 13 14 15  .  .  . |
  |         '----------------------'
  |            |
  |            |
  |            .----------------------.
  |          N | 19 20 22  .  .  .  . |
  |            '----------------------'
  |                    
  ------------------------------------------------------------------*/
  void insert_nonfull(BTreeNode<T,M>& node, T* key)
  {
    auto i = node.n_keys();

    // Check that the given node is nonfull
    ASSERT(node.n_keys() != (2*M-1), "Invalid BTree structure");

    // ---- Leaf node ----
    if ( node.is_leaf() )
    {
      ++node.n_keys_;

      // Traverse all node keys until input key is greatest
      while ( i > 0 && (*key < *node.key(i-1)) )
      {
        node.key(i, node.key(i-1));
        --i;
      }

      // Place key
      node.key(i, key);
    }
    // ---- Branch node ----
    else 
    {
      // Traverse all node keys until input key is greatest
      while ( i > 0 && (*key < *node.key(i-1)) )
        --i;

      // In case that the found child is full, split it in two nodes
      if ( node.child(i)->n_keys() == 2*M-1 )
      {
        // Split the node's i-th child in two nodes
        // -> This inserts a new key to the current node 
        //    at position i
        split_child( node, i );

        // If input key is larger than new key of current node,
        // insert input key in newly generated node
        if ( *key > *node.key(i) )
          ++i;
      }

      insert_nonfull(*node.child(i), key);
    }

  } // insert_nonfull() 


  /*------------------------------------------------------------------ 
  | Split a full node (with 2*M-1 keys) around its median key 
  | into two nodes with (M-1) keys each.
  | The median key moves up into the nonfull node's parent. 
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
    // Check that the parent node is nonfull
    ASSERT(parent_node.n_keys() != (2*M-1), "Invalid BTree structure");

    // This is the new node. It gets the keys of the full node,
    // that are greated than the median
    auto new_node = std::make_unique<BTreeNode<T,M>>();

    // Child must be a full node
    BTreeNode<T,M>* child_node = parent_node.child(i);
    ASSERT(child_node->n_keys() == (2*M-1), "Invalid BTree structure");

    new_node->is_leaf_ = child_node->is_leaf();
    new_node->n_keys_  = M-1;

    // Put all keys of the full child into the new node,
    // which are greater than the median (=M)
    for (std::size_t j = 0; j < M-1; ++j)
      new_node->key(j, child_node->key(j+M));

    // If the child node is not a leaf, pass also its 
    // children to the new node
    if ( !child_node->is_leaf() )
      for (std::size_t j = 0; j < M; ++j) 
        new_node->child(j, child_node->child_ptr(j+M));

    // Update number of keys 
    child_node->n_keys_ = M-1;

    // new_node will be placed at parent node's (i+1)-th child,
    // so all other children at j > i must be shifted to the right
    for (std::size_t j = parent_node.n_keys()+1; j > i+1; --j)
      parent_node.child(j, parent_node.child_ptr(j-1));

    ++parent_node.n_keys_;
    parent_node.child(i+1, new_node);

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

