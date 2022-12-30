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
#include "BBoxND.h"
#include "OcTreeND.h"
#include "Testing.h"
#include "Geometry.h"

#include "Timer.h"

namespace OcTreeNDTests 
{
using namespace CppUtils;

/*--------------------------------------------------------------------
| Test constructor
--------------------------------------------------------------------*/
void constructor()
{
  /*
  OcTreeND<int, 1, int, 1> tree { {-100,100} };

  std::vector<int> values {};

  values.push_back( 1 );
  values.push_back( 2 );
  values.push_back( 3 );
  values.push_back( 4 );
  values.push_back( 5 );
  values.push_back( 6 );

  for ( auto& v : values )
    tree.insert( v, {v} );

  */


  OcTreeND<int, 1, double, 2> tree { {{-10.0,-10.0}, {10.0,10.0}} };

  int value = 1;

  std::vector<VecND<double,2>> positions {};

  positions.push_back( { 0.0,  0.0} );
  positions.push_back( { 1.0,  1.0} );
  positions.push_back( { 2.0,  2.0} );
  positions.push_back( {-3.0, -3.0} );
  positions.push_back( {-4.0, -4.0} );
  positions.push_back( {-5.0, -5.0} );

  for ( auto& pos : positions )
    CHECK( tree.insert( value, pos ) );

  for ( auto& it : tree )
  {
    LOG(INFO) << it.curve_id() << " => " << it.bbox();

    for ( auto& e : it.entries() )
      LOG(INFO) << "  -> " << e.position;

    if ( it.n_entries() < 1 )
      LOG(INFO) << "  -> HAS NO ENTRIES";
  }

  for ( auto& pos : positions )
    CHECK( tree.remove( value, pos ) );

  // Export the tree structure
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };

  std::string file_name 
  { source_dir + "/auxiliary/test_data/OcTree_constructor" };

  OcTreeNDWriter writer { tree };

  writer.write_to_vtu( file_name );

} // constructor() 

/*--------------------------------------------------------------------
| Test 3D octree
--------------------------------------------------------------------*/
void test_3d()
{
  OcTreeND<int, 1, double, 3> tree { {{0.0,0.0,0.0}, {10.0,10.0,10.0}} };

  int value = 1;

  std::vector<VecND<double,3>> positions {};

  positions.push_back( { 0.5,  0.8, 1.2 } );
  positions.push_back( { 1.5,  2.8, 3.2 } );
  positions.push_back( { 7.5,  5.8, 3.1 } );
  positions.push_back( { 0.2,  0.2, 0.8 } );
  positions.push_back( { 0.1,  0.3, 1.2 } );

  for ( auto& pos : positions )
    CHECK( tree.insert( value, pos ) );

  // Export the tree structure
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };

  std::string file_name 
  { source_dir + "/auxiliary/test_data/OcTree_3d" };

  OcTreeNDWriter writer { tree };

  writer.write_to_vtu( file_name );

} // test_3d()

} // namespace OcTreeNDTests


/*********************************************************************
* Run tests for: OcTreeND.h
*********************************************************************/
void run_tests_OcTreeND()
{
  OcTreeNDTests::constructor();
  OcTreeNDTests::test_3d();

} // run_tests_OcTreeND()
