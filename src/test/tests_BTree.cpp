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
  DEBUG_LOG("\n-------------------------------------------");
  DEBUG_LOG(" Test function: insert_nonfull_second() ");

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

  DEBUG_LOG("ROOT");
  for (std::size_t i = 0; i < btree.root().n_keys(); ++i)
    DEBUG_LOG(*btree.root().key(i));


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


  /*
  DEBUG_LOG("ROOT");
  for (std::size_t i = 0; i < btree.root().n_keys(); ++i)
    DEBUG_LOG(*btree.root().key(i));

  DEBUG_LOG("CHILD 1");
  for (std::size_t i = 0; i < btree.root().child(0)->n_keys(); ++i)
    DEBUG_LOG(*btree.root().child(0)->key(i));

  DEBUG_LOG("CHILD 2");
  for (std::size_t i = 0; i < btree.root().child(1)->n_keys(); ++i)
    DEBUG_LOG(*btree.root().child(1)->key(i));

  //DEBUG_LOG("CHILD 3");
  for (std::size_t i = 0; i < btree.root().child(2)->n_keys(); ++i)
    DEBUG_LOG(*btree.root().child(2)->key(i));
  */


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
