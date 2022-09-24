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
#include "Helpers.h"
#include "Log.h"

namespace CppUtils {


/*********************************************************************
* Resources:
* - https://tildesites.bowdoin.edu/~ltoma/teaching/cs340/spring08/Papers/Rtree-chap1.pdf
*********************************************************************/

/*********************************************************************
* Minimum bounding rectangle 
*********************************************************************/
struct MBR
{
  MBR() {}

  MBR(const Vec2d& ll, const Vec2d& ur)
  : lowleft_ { ll }
  , upright_ { ur }
  {}

  Vec2d lowleft_ { 0.0, 0.0 };
  Vec2d upright_ { 0.0, 0.0 };

}; // MBR


/*********************************************************************
* The RTree node implementation
*********************************************************************/
template <typename T, long unsigned int M>
class RTreeNode
{
public:
  using Values   = std::array<T*,M>;
  using Children = std::array<std::unique_ptr<RTreeNode<T,M>>,M>;
  using MBRs     = std::array<MBR,M>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  RTreeNode() {}


  /*------------------------------------------------------------------ 
  | Getter
  ------------------------------------------------------------------*/
  Values& values() { return values_; }
  const Values& values() const { return values_; }
  const MBRs& mbrs() const { return mbrs_; }
  bool is_leaf() const { return is_leaf; }

protected:
  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  Values   values_   { nullptr };
  Children children_ { nullptr };
  MBRs     mbrs_     { };

  bool     is_leaf_  { false };

}; // RTreeNode



/*********************************************************************
* This class defines the interface to an R-tree structure 
* for 2D simplices
*********************************************************************/
template <typename T, long unsigned int M>
class RTree
{
public:

  using Roots = std::array<RTreeNode<T,M>,M>;
  using MBRs  = std::array<MBR,M>;

  friend class RTreeNode<T,M>;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  RTree(const Vec2d& lowleft, const Vec2d& upright)
  : mbr_ { lowleft, upright } 
  {}

  /*------------------------------------------------------------------ 
  | Insert
  ------------------------------------------------------------------*
  bool insert(T* item)
  {
    if ( !item ) 
      return false;

    // The mbr of the item
    MBR mbr_item { item->lowleft(), item->upright() }; 

    // Traverse the tree from its root to the appropriate leaf
    // At each level, select the node L, whose MBR will require 
    // the minimum area enlargement to cover T's MBR
    RTreeNode* leaf = nullptr;

    return true;

  } // insert()*/

private:
  /*------------------------------------------------------------------ 
  | Attributes
  ------------------------------------------------------------------*/
  MBR   mbr_;

  Roots roots_ {};
  MBRs  mbrs_  {};

}; // RTree


} // CppUtils
