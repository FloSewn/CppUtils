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
#include <algorithm>
#include <random>

#include <CppUtilsConfig.h>

#include "tests.h"

#include "VecND.h"
#include "Testing.h"
#include "RTreeND.h"

namespace RTreeNDTests 
{
using namespace CppUtils;

/*--------------------------------------------------------------------
| Test constructor
--------------------------------------------------------------------*/
void constructor()
{
  RTreeND<int,3,double,2> tree {};

  std::vector<int> values;
  std::vector<BBoxND<double,2>> bboxes;

  values.push_back( 1 );
  bboxes.push_back( { {0.0,0.0}, {1.0,1.0} } );

  values.push_back( 2 );
  bboxes.push_back( { {0.5,1.0}, {1.5,2.0} } );

  values.push_back( 3 );
  bboxes.push_back( { {2.0,0.5}, {3.0,1.5} } );

  values.push_back( 4 );
  bboxes.push_back( { {2.0,2.0}, {3.0,3.0} } );

  values.push_back( 5 );
  bboxes.push_back( { {4.0,4.0}, {5.0,5.0} } );

  values.push_back( 6 );
  bboxes.push_back( { {4.0,1.0}, {4.5,1.5} } );

  values.push_back( 7 );
  bboxes.push_back( { {3.0,4.5}, {5.0,5.5} } );

  for (std::size_t i = 0; i < values.size(); ++i)
    tree.insert( values[i], bboxes[i] );


  // Export the tree structure
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };

  std::string file_name 
  { source_dir + "/auxiliary/test_data/RTree_constructor" };

  RTreeNDWriter writer { tree };

  writer.write_to_vtu( file_name );
  writer.write_to_txt( file_name );
  writer.print(std::cout);


  // Check internal connectivity
  auto& root = tree.root();
  auto& c_0 = root.child(0);
  auto& c_1 = root.child(1);

  CHECK( c_0.left()  == nullptr );
  CHECK( c_0.right() == &c_1 );
  CHECK( c_1.left()  == &c_0 );
  CHECK( c_1.right() == nullptr );

  auto& c_0_0 = c_0.child(0);
  auto& c_0_1 = c_0.child(1);

  auto& c_1_0 = c_1.child(0);

  CHECK( c_0_0.right() == &c_0_1 );
  CHECK( c_0_1.right() == &c_1_0 );
  CHECK( c_1_0.right() == nullptr );

  CHECK( c_1_0.left() == &c_0_1 );
  CHECK( c_0_1.left() == &c_0_0 );
  CHECK( c_0_0.left()  == nullptr );


  // Removal
  tree.remove(values[0], bboxes[0]);
  tree.remove(values[1], bboxes[1]);
  tree.remove(values[2], bboxes[2]);
  tree.remove(values[3], bboxes[3]);
  tree.remove(values[4], bboxes[4]);
  tree.remove(values[5], bboxes[5]);
  tree.remove(values[6], bboxes[6]);

} // constructor()

/*--------------------------------------------------------------------
| Test 1D tree
--------------------------------------------------------------------*/
void insertion_1d()
{
  RTreeND<int,3,int,1> tree {};
  RTreeNDWriter writer { tree };

  std::vector<int> values;
  std::vector<BBoxND<int,1>> bboxes;

  int k = 13;
  int n = 13;

  for ( std::size_t i = 0; i < k; ++i )
    values.push_back( i );

  for ( auto v : values )
    bboxes.push_back( {v, v} );

  // Permutate objects
  std::vector<size_t> perm_ids( values.size() );
  std::iota(perm_ids.begin(), perm_ids.end(), 0);

  auto rng = std::default_random_engine {};
  std::shuffle(perm_ids.begin(), perm_ids.end(), rng);

  for (std::size_t i = 0; i < n; ++i)
    tree.insert( values[perm_ids[i]], bboxes[perm_ids[i]] );

  //tree.insert( values, bboxes );

  writer.print(std::cout);


  // Check for correct leaf order and check tree iterator behaviour
  std::size_t i = 0;
  for ( auto& it : tree )
    CHECK( &it == &(values[i++]) );

  i = 0;
  for ( auto it : tree )
  {
    CHECK(  it ==   values[i]  );
    CHECK( &it != &(values[i++]) );
  }

  i = 0;
  for ( const auto& it : tree )
    CHECK( &it == &(values[i++]) );

  i = 0;
  for ( const auto it : tree )
  {
    CHECK(  it ==   values[i]  );
    CHECK( &it != &(values[i++]) );
  }


  //int s1 = std::accumulate(tree.cbegin(), tree.cend(), 0);
  //int s2 = std::accumulate(values.cbegin(), values.cend(), 0);

  //CHECK( s1 == s2 );


} // insertion_1d()

/*--------------------------------------------------------------------
| Test bulk insertion of many objects
--------------------------------------------------------------------*/
void bulk_insertion_2d()
{
  RTreeND<int,3,double,2> tree {};

  std::vector<int> values;
  std::vector<BBoxND<double,2>> bboxes;


  // Fill the tree 
  values.push_back( 1 );
  bboxes.push_back( { {3.0,4.5}, {5.0,5.5} } );

  values.push_back( 2 );
  bboxes.push_back( { {0.0,0.0}, {1.0,1.0} } );

  values.push_back( 3 );
  bboxes.push_back( { {4.0,4.0}, {5.0,5.0} } );

  values.push_back( 4 );
  bboxes.push_back( { {2.0,0.5}, {3.0,1.5} } );

  values.push_back( 5 );
  bboxes.push_back( { {0.5,1.0}, {1.5,2.0} } );

  values.push_back( 6 );
  bboxes.push_back( { {2.0,2.0}, {3.0,3.0} } );

  values.push_back( 7 );
  bboxes.push_back( { {4.0,1.0}, {4.5,1.5} } );

  tree.insert( values, bboxes );


  // Export the tree structure
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };

  std::string file_name 
  { source_dir + "/auxiliary/test_data/RTree_bulk_insertion_2d" };

  RTreeNDWriter writer { tree };

  writer.write_to_vtu( file_name );
  writer.write_to_txt( file_name );
  writer.print(std::cout);


  // Check internal connectivity
  auto& root = tree.root();
  auto& c_0 = root.child(0);
  auto& c_1 = root.child(1);
  auto& c_2 = root.child(2);

  CHECK( c_0.right() == &c_1 );
  CHECK( c_1.right() == &c_2 );
  CHECK( c_2.right() == nullptr );

  CHECK( c_0.left() == nullptr );
  CHECK( c_1.left() == &c_0 );
  CHECK( c_2.left() == &c_1 );


} // bulk_insertion_2d()

/*--------------------------------------------------------------------
| Test bulk insertion of many objects
--------------------------------------------------------------------*/
void bulk_insertion_3d()
{
  RTreeND<int, 3, double, 3> tree {};

  std::vector<int> values;
  std::vector<BBoxND<double,3>> bboxes;

  values.push_back( 1 );
  bboxes.push_back( { {0.0,0.0,0.0}, {2.0,2.0,2.0} } );

  values.push_back( 2 );
  bboxes.push_back( { {3.0,3.0,3.0}, {5.0,5.0,5.0} } );

  values.push_back( 3 );
  bboxes.push_back( { {6.0,7.0,3.0}, {8.0,8.0,5.0} } );

  values.push_back( 4 );
  bboxes.push_back( { {1.0,3.0,1.0}, {8.0,8.0,5.0} } );

  values.push_back( 5 );
  bboxes.push_back( { {2.0,2.0,2.0}, {2.5,2.5,2.5} } );

  tree.insert( values, bboxes );

  //for (std::size_t i = 0; i < values.size(); ++i)
  //  tree.insert( values[i], bboxes[i] );

  // Export tree structure
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };

  std::string file_name 
  { source_dir + "/auxiliary/test_data/RTree_bulk_insertion_3d" };

  RTreeNDWriter writer { tree };

  writer.write_to_vtu( file_name );
  writer.write_to_txt( file_name );
  writer.print(std::cout);

} // bulk_insertion_3d()

} // namespace RTreeNDTests


/*********************************************************************
* Run tests for: RTreeND.h
*********************************************************************/
void run_tests_RTreeND()
{
  RTreeNDTests::constructor();

  std::cout << "-----------------------------------------" << std::endl;

  RTreeNDTests::insertion_1d();

  std::cout << "-----------------------------------------" << std::endl;

  RTreeNDTests::bulk_insertion_2d();

  std::cout << "-----------------------------------------" << std::endl;

  RTreeNDTests::bulk_insertion_3d();

} // run_tests_RTreeND()
