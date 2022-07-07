/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <cassert>

#include "CppUtilsConfig.h"

#include "tests.h"

#include "Testing.h"
#include "VtkIO.h"

namespace VtkIOTests 
{

using namespace CppUtils;

/*--------------------------------------------------------------------
| Test writing
--------------------------------------------------------------------*/
void write()
{
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };
  std::string file_name {source_dir + "/aux/VTU_File.vtu"};


  std::vector<double> points {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    2.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    1.0, 1.0, 0.0,
    2.0, 1.0, 0.0,
    0.0, 2.0, 0.0,
    1.0, 2.0, 0.0,
  };

  std::vector<size_t> connectivity {
    0, 1, 3,
    1, 4, 3,
    1, 2, 4,
    2, 5, 4,
    3, 4, 6,
    4, 7, 6,
  };

  std::vector<size_t> offsets {
    3, 6, 9, 12, 15, 18,
  };

  std::vector<size_t> types {
    5, 5, 5, 5, 5, 5,
  };

  VtuWriter writer { points, connectivity, offsets, types };



  std::vector<double> pressure {
    0.0, 1.0, 2.0, 1.0, 2.0, 3.0, 2.0, 3.0,
  };

  std::vector<float> velocity {
    1.0, 1.0, 0.0,
    2.0, 1.0, 0.0,
    3.0, 1.0, 0.0,
    1.0, 2.0, 0.0,
    2.0, 2.0, 0.0,
    3.0, 2.0, 0.0,
    1.0, 3.0, 0.0,
    2.0, 3.0, 0.0,
  };

  std::vector<int> node_max {
    3, 4, 4, 5, 6, 7,
  };

  writer.add_point_data( pressure, "pressure", 1 );
  writer.add_point_data( velocity, "velocity", 3 );
  writer.add_cell_data( node_max, "node_max", 1 );


  writer.write( file_name );


  CHECK( true );

} // write()


} // namespace VtkIOTests

/*********************************************************************
* Run tests for: VtkIO.h
*********************************************************************/
void run_tests_VtkIO()
{
  VtkIOTests::write();

} // run_tests_VtkIO()
