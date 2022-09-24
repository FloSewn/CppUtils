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
  found = btree.search(data[3], btree.root());
  CHECK( found->key(1) == &data[3] );


} // insert_nonfull_second()


} // namespace BTreeTests


/*********************************************************************
* Run tests for: BTree.h
*********************************************************************/
void run_tests_BTree()
{
  BTreeTests::constructor();
  BTreeTests::insert_nonfull_first_descending();
  BTreeTests::insert_nonfull_first_ascending();
  BTreeTests::insert_nonfull_first_equal();
  BTreeTests::insert_nonfull_second();

} // run_tests_BTree()
