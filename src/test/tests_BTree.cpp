/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <memory>
#include <vector>
#include <cassert>

#include "tests.h"

#include "Testing.h"
#include "BTree.h"

namespace BTreeTests 
{
using namespace CppUtils;


/*--------------------------------------------------------------------
| Test constructor
--------------------------------------------------------------------*/
void constructor()
{
  BTree<int,2> btree;
  std::vector<int> data {};

  data.push_back(0);

  const BTreeNode<int,2>* found = nullptr;

  found = btree.search(data[0], btree.root());

  CHECK( found == nullptr );

} // constructor()

/*--------------------------------------------------------------------
| Test first branch of insert_nonfull(), where a large value is 
| put in first, followed by a smaller value
--------------------------------------------------------------------*/
void insert_nonfull_first_descending()
{
  BTree<int,2> btree;
  std::vector<int> data {};

  data.push_back(0);
  data.push_back(1);

  const BTreeNode<int,2>* found = nullptr;

  btree.insert( &data[1] );
  found = btree.search(data[1], btree.root());
  CHECK( found->key(0) == &data[1] );

  btree.insert( &data[0] );
  found = btree.search(data[0], btree.root());
  CHECK( found->key(0) == &data[0] );


} // insert_nonfull_first_descending()

/*--------------------------------------------------------------------
| Test first branch of insert_nonfull(), where a small value is 
| put in first, followed by a larger value
--------------------------------------------------------------------*/
void insert_nonfull_first_ascending()
{
  BTree<int,2> btree;
  std::vector<int> data {};

  data.push_back(0);
  data.push_back(1);

  const BTreeNode<int,2>* found = nullptr;

  btree.insert( &data[0] );
  found = btree.search(data[0], btree.root());
  CHECK( found->key(0) == &data[0] );

  btree.insert( &data[1] );
  found = btree.search(data[1], btree.root());
  CHECK( found->key(1) == &data[1] );


} // insert_nonfull_first_ascending()

/*--------------------------------------------------------------------
| Test first branch of insert_nonfull(), where a equal values are 
| added
--------------------------------------------------------------------*/
void insert_nonfull_first_equal()
{
  BTree<int,2> btree;
  std::vector<int> data {};

  data.push_back(0);
  data.push_back(0);

  const BTreeNode<int,2>* found = nullptr;

  btree.insert( &data[0] );
  found = btree.search(data[0], btree.root());
  CHECK( found->key(0) == &data[0] );

  btree.insert( &data[1] );
  found = btree.search(data[1], btree.root());
  CHECK( found->key(1) == &data[1] );

} // insert_nonfull_first_equal()

/*--------------------------------------------------------------------
| Test second branch of insert_nonfull()
--------------------------------------------------------------------*/
void insert_nonfull_second()
{
  BTree<int,2> btree;
  std::vector<int> data {};

  data.push_back(0);
  data.push_back(1);
  data.push_back(2);
  data.push_back(3);
  
  data.push_back(4);
  data.push_back(5);

  data.push_back(6);
  data.push_back(7);
  data.push_back(8);
  data.push_back(9);
  data.push_back(10);

  const BTreeNode<int,2>* found = nullptr;

  btree.insert( &data[0] );
  found = btree.search(data[0], btree.root());
  CHECK( found->key(0) == &data[0] );

  btree.insert( &data[1] );
  found = btree.search(data[1], btree.root());
  CHECK( found->key(1) == &data[1] );

  btree.insert( &data[2] );
  found = btree.search(data[2], btree.root());
  CHECK( found->key(2) == &data[2] );

  // BTree root should split now
  btree.insert( &data[3] );

  CHECK( !btree.root().is_leaf() );
  CHECK( btree.root().n_keys() == 1);
  CHECK( btree.root().child(0) != nullptr);
  CHECK( btree.root().child(1) != nullptr);

  found = btree.search(data[3], btree.root());
  CHECK( found != nullptr );
  CHECK( found->key(1) == &data[3] );

  btree.insert( &data[4] );
  found = btree.search(data[4], btree.root());
  CHECK( found != nullptr );
  CHECK( found->key(2) == &data[4] );

  btree.insert( &data[5] );
  found = btree.search(data[5], btree.root());
  CHECK( found != nullptr );
  CHECK( found->key(1) == &data[5] );

  btree.insert( &data[6] );
  btree.insert( &data[7] );
  btree.insert( &data[8] );
  btree.insert( &data[9] );
  btree.insert( &data[10] );

  found = btree.search(data[10], btree.root());
  CHECK( found != nullptr );
  CHECK( found->key(2) == &data[10] );



  //std::cout << btree << std::endl;


} // insert_nonfull_second()


/*--------------------------------------------------------------------
| Test the BTree insertion
--------------------------------------------------------------------*/
void insertion()
{
  BTree<int,2> btree;
  std::vector<int> data {};

  // Vector acts as container - data must be stored prior to insertion
  // into BTree, since element addresses change when vector needs to
  // grow!
  data.push_back(1);
  data.push_back(2);
  data.push_back(3);
  data.push_back(4);
  data.push_back(5);
  data.push_back(6);
  data.push_back(7);
  data.push_back(8);
  data.push_back(9);

  // Add 1st value
  // --------------------------------
  //
  //      .-----------.
  //      | 1 |   |   |
  //      '-----------'
  //
  btree.insert( &data[0] );

  {
    const auto& root = btree.root();
    CHECK( root.n_keys() == 1 );
    CHECK( *root.key(0) == data[0] );

    const BTreeNode<int,2>* found = nullptr;
    found = btree.search(data[0], root);
    CHECK( *found->key(0) == data[0] );
  }

  // Add 2nd value
  // --------------------------------
  //
  //      .-----------.
  //      | 1 | 2 |   |
  //      '-----------'
  //
  btree.insert( &data[1] );

  {
    const auto& root = btree.root();
    CHECK( root.n_keys() == 2 );
    CHECK( *root.key(0) == data[0] );
    CHECK( *root.key(1) == data[1] );

    const BTreeNode<int,2>* found = nullptr;
    found = btree.search(data[1], root);
    CHECK( *found->key(1) == data[1] );

    found = btree.search(data[0], root);
    CHECK( *found->key(0) == data[0] );
  }

  // Add 3rd value
  // --------------------------------
  //
  //      .-----------.
  //      | 1 | 2 | 3 |
  //      '-----------'
  //
  btree.insert( &data[2] );

  {
    const auto& root = btree.root();
    CHECK( root.n_keys() == 3 );
    CHECK( *root.key(0) == data[0] );
    CHECK( *root.key(1) == data[1] );
    CHECK( *root.key(2) == data[2] );

    const BTreeNode<int,2>* found = nullptr;
    found = btree.search(data[2], root);
    CHECK( *found->key(2) == data[2] );

    found = btree.search(data[1], root);
    CHECK( *found->key(1) == data[1] );

    found = btree.search(data[0], root);
    CHECK( *found->key(0) == data[0] );
  }


  // Add 4th value
  // --------------------------------
  // -> Root node will be split in two nodes
  // -> Median key (2) remains in root
  // -> Key (1) is placed in first child
  // -> Keys (3) & (4) is placed in second child
  //
  //      .-----------.
  //      | 2 |   |   |
  //      '-----------'
  //      :   :
  //      :   .-----------.
  //      :   | 3 | 4 |   |
  //      :   '-----------'
  //      :
  //      .-----------.
  //      | 1 |   |   |
  //      '-----------'
  //
  //
  btree.insert( &data[3] );

  {
    const auto& root = btree.root();
    CHECK( root.n_keys() == 1 );
    CHECK( *root.key(0) == data[1] );

    const auto& child_0 = *root.child(0);
    CHECK( child_0.n_keys() == 1 );
    CHECK( *child_0.key(0) == data[0] );

    const auto& child_1 = *root.child(1);
    CHECK( child_1.n_keys() == 2 );
    CHECK( *child_1.key(0) == data[2] );
    CHECK( *child_1.key(1) == data[3] );


    const BTreeNode<int,2>* found = nullptr;
    found = btree.search(data[0], root);
    CHECK( *found->key(0) == data[0] );

    found = btree.search(data[1], root);
    CHECK( *found->key(0) == data[1] );

    found = btree.search(data[2], root);
    CHECK( *found->key(0) == data[2] );

    found = btree.search(data[3], root);
    CHECK( *found->key(1) == data[3] );
  }


  // Add 5th value
  // --------------------------------
  //
  //      .-----------.
  //      | 2 |   |   |
  //      '-----------'
  //      :   :
  //      :   .-----------.
  //      :   | 3 | 4 | 5 |
  //      :   '-----------'
  //      :
  //      .-----------.
  //      | 1 |   |   |
  //      '-----------'
  //
  //
  btree.insert( &data[4] );

  {
    const auto& root = btree.root();
    CHECK( root.n_keys() == 1 );
    CHECK( *root.key(0) == data[1] );

    const auto& child_0 = *root.child(0);
    CHECK( child_0.n_keys() == 1 );
    CHECK( *child_0.key(0) == data[0] );

    const auto& child_1 = *root.child(1);
    CHECK( child_1.n_keys() == 3 );
    CHECK( *child_1.key(0) == data[2] );
    CHECK( *child_1.key(1) == data[3] );
    CHECK( *child_1.key(2) == data[4] );


    const BTreeNode<int,2>* found = nullptr;
    found = btree.search(data[0], root);
    CHECK( *found->key(0) == data[0] );

    found = btree.search(data[1], root);
    CHECK( *found->key(0) == data[1] );

    found = btree.search(data[2], root);
    CHECK( *found->key(0) == data[2] );

    found = btree.search(data[3], root);
    CHECK( *found->key(1) == data[3] );

    found = btree.search(data[4], root);
    CHECK( *found->key(2) == data[4] );
  }

  // Add 6th value
  // --------------------------------
  // -> New child node will be added
  // -> Median key (4) of full second child is placed into root node
  // -> Remaining keys (5) and (6) are placed in new child node
  //
  //      .-----------.
  //      | 2 | 4 |   |
  //      '-----------'
  //      :   :   : 
  //      :   :   .-----------.
  //      :   :   | 5 | 6 |   |
  //      :   :   '-----------'
  //      :   :
  //      :   .-----------.
  //      :   | 3 |   |   |
  //      :   '-----------'
  //      :
  //      .-----------.
  //      | 1 |   |   |
  //      '-----------'
  //
  //
  btree.insert( &data[5] );

  {
    const auto& root = btree.root();
    CHECK( root.n_keys() == 2 );
    CHECK( *root.key(0) == data[1] );

    const auto& child_0 = *root.child(0);
    CHECK( child_0.n_keys() == 1 );
    CHECK( *child_0.key(0) == data[0] );

    const auto& child_1 = *root.child(1);
    CHECK( child_1.n_keys() == 1 );
    CHECK( *child_1.key(0) == data[2] );

    const auto& child_2 = *root.child(2);
    CHECK( child_2.n_keys() == 2 );
    CHECK( *child_2.key(0) == data[4] );
    CHECK( *child_2.key(1) == data[5] );


    const BTreeNode<int,2>* found = nullptr;
    found = btree.search(data[0], root);
    CHECK( *found->key(0) == data[0] );

    found = btree.search(data[1], root);
    CHECK( *found->key(0) == data[1] );

    found = btree.search(data[2], root);
    CHECK( *found->key(0) == data[2] );

    found = btree.search(data[3], root);
    CHECK( *found->key(1) == data[3] );

    found = btree.search(data[4], root);
    CHECK( *found->key(0) == data[4] );

    found = btree.search(data[5], root);
    CHECK( *found->key(1) == data[5] );
  }


  // Add 7th value
  // --------------------------------
  //
  //      .-----------.
  //      | 2 | 4 |   |
  //      '-----------'
  //      :   :   : 
  //      :   :   .-----------.
  //      :   :   | 5 | 6 | 7 |
  //      :   :   '-----------'
  //      :   :
  //      :   .-----------.
  //      :   | 3 |   |   |
  //      :   '-----------'
  //      :
  //      .-----------.
  //      | 1 |   |   |
  //      '-----------'
  //
  //
  btree.insert( &data[6] );

  {
    const auto& root = btree.root();
    CHECK( root.n_keys() == 2 );
    CHECK( *root.key(0) == data[1] );

    const auto& child_0 = *root.child(0);
    CHECK( child_0.n_keys() == 1 );
    CHECK( *child_0.key(0) == data[0] );

    const auto& child_1 = *root.child(1);
    CHECK( child_1.n_keys() == 1 );
    CHECK( *child_1.key(0) == data[2] );

    const auto& child_2 = *root.child(2);
    CHECK( child_2.n_keys() == 3 );
    CHECK( *child_2.key(0) == data[4] );
    CHECK( *child_2.key(1) == data[5] );
    CHECK( *child_2.key(2) == data[6] );
  }

  // Add 8th value
  // --------------------------------
  //
  //      .-----------.
  //      | 2 | 4 | 6 |
  //      '-----------'
  //      :   :   :   :
  //      :   :   :   .-----------. 
  //      :   :   :   | 7 | 8 |   |
  //      :   :   :   '-----------'
  //      :   :   : 
  //      :   :   .-----------.
  //      :   :   | 5 |   |   |
  //      :   :   '-----------'
  //      :   :
  //      :   .-----------.
  //      :   | 3 |   |   |
  //      :   '-----------'
  //      :
  //      .-----------.
  //      | 1 |   |   |
  //      '-----------'
  //
  //
  btree.insert( &data[7] );

  {
    const auto& root = btree.root();
    CHECK( root.n_keys() == 3 );
    CHECK( *root.key(0) == data[1] );
    CHECK( *root.key(1) == data[3] );
    CHECK( *root.key(2) == data[5] );

    const auto& child_0 = *root.child(0);
    CHECK( child_0.n_keys() == 1 );
    CHECK( *child_0.key(0) == data[0] );

    const auto& child_1 = *root.child(1);
    CHECK( child_1.n_keys() == 1 );
    CHECK( *child_1.key(0) == data[2] );

    const auto& child_2 = *root.child(2);
    CHECK( child_2.n_keys() == 1 );
    CHECK( *child_2.key(0) == data[4] );

    const auto& child_3 = *root.child(3);
    CHECK( child_3.n_keys() == 2 );
    CHECK( *child_3.key(0) == data[6] );
    CHECK( *child_3.key(1) == data[7] );
  }


  // Add 9th value
  // --------------------------------
  // -> root is full and thus splits
  // -> introduce new root node, which gets median key of old root
  //
  //      .-----------.
  //      | 4 |   |   |
  //      '-----------'
  //      :   :
  //      :   .-----------.
  //      :   | 6 |   |   |
  //      :   '-----------'
  //      :   :   : 
  //      :   :   .-----------.
  //      :   :   | 7 | 8 | 9 |
  //      :   :   '-----------'
  //      :   .-----------.   
  //      :   | 5 |   |   |
  //      :   '-----------'
  //      :
  //      .-----------.
  //      | 2 |   |   |
  //      '-----------'
  //      :   : 
  //      :   .-----------.
  //      :   | 3 |   |   |
  //      :   '-----------'
  //      .-----------.   
  //      | 1 |   |   |
  //      '-----------'
  //
  //
  btree.insert( &data[8] );

  {
    const auto& root = btree.root();
    CHECK( root.n_keys() == 1 );
    CHECK( *root.key(0) == data[3] );

    const auto& child_0 = *root.child(0);
    CHECK( child_0.n_keys() == 1 );
    CHECK( *child_0.key(0) == data[1] );

    const auto& child_1 = *root.child(1);
    CHECK( child_1.n_keys() == 1 );
    CHECK( *child_1.key(0) == data[5] );

    const auto& child_0_0 = *child_0.child(0);
    CHECK( child_0_0.n_keys() == 1 );
    CHECK( *child_0_0.key(0) == data[0] );

    const auto& child_0_1 = *child_0.child(1);
    CHECK( child_0_1.n_keys() == 1 );
    CHECK( *child_0_1.key(0) == data[2] );

    const auto& child_1_0 = *child_1.child(0);
    CHECK( child_1_0.n_keys() == 1 );
    CHECK( *child_1_0.key(0) == data[4] );

    const auto& child_1_1 = *child_1.child(1);
    CHECK( child_1_1.n_keys() == 3 );
    CHECK( *child_1_1.key(0) == data[6] );
    CHECK( *child_1_1.key(1) == data[7] );
    CHECK( *child_1_1.key(2) == data[8] );


    const BTreeNode<int,2>* found = nullptr;
    found = btree.search(data[6], root);
    CHECK(  found->n_keys() == 3 );
    CHECK( *found->key(0) == data[6] );
    CHECK( *found->key(1) == data[7] );
    CHECK( *found->key(2) == data[8] );

  }


  
  //LOG(INFO) << "\n\n------- BTree insertion() --------\n";
  //LOG(INFO) << "\n" << btree;






} // insertion()


} // namespace BTreeTests


/*********************************************************************
* Run tests for: BTree.h
*********************************************************************/
void run_tests_BTree()
{
  BTreeTests::constructor();
  BTreeTests::insertion();
  BTreeTests::insert_nonfull_first_descending();
  BTreeTests::insert_nonfull_first_ascending();
  BTreeTests::insert_nonfull_first_equal();
  BTreeTests::insert_nonfull_second();

} // run_tests_BTree()
