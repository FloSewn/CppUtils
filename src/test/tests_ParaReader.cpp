/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <cassert>
#include <string>

#include "CppUtilsConfig.h"

#include "tests.h"
#include "Helpers.h"
#include "MathUtility.h"

#include "ParaReader.h"
#include "Testing.h"

namespace ParaReaderTests 
{

using namespace CppUtils;


/*--------------------------------------------------------------------
| Test reader constructor
--------------------------------------------------------------------*/
static void constructor()
{
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };
  std::string file_name {source_dir + "/auxiliary/ParameterFile.txt"};

  // Handle invalid input file path
  try 
  {
    ParaReader reader { file_name + "INVALID" };
    CHECK(false);
  }
  catch (...)
  {
    CHECK( true );
  }

  // Handle correct input file path
  ParaReader reader { file_name };
  CHECK( true );

} // constructor

/*--------------------------------------------------------------------
| Test scalar parameters
--------------------------------------------------------------------*/
static void scalar_parameters()
{
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };
  std::string file_name {source_dir + "/auxiliary/ParameterFile.txt"};

  ParaReader reader { file_name };

  // Create working parameter definitions
  //------------------------------------------------------------------
  reader.new_scalar_parameter<bool>( 
      "bool_param",  "Scalar bool parameter:" );
  reader.new_scalar_parameter<int>( 
      "int_param", "Scalar int parameter:" );
  reader.new_scalar_parameter<float>( 
      "float_param", "Scalar float parameter:" );
  reader.new_scalar_parameter<double>( 
      "double_param", "Scalar double parameter:" );
  reader.new_scalar_parameter<int>( 
      "invalid_int_param", "Invalid int parameter:" );
  reader.new_scalar_parameter<std::string>(
      "str_param", "String parameter:");

  // Check that invalid parameter definition does not succeed
  //------------------------------------------------------------------
  try
  {
    reader.new_scalar_parameter<double>( 
        "double_param", "Scalar double parameter:" );
    CHECK(false);
  }
  catch (...)
  {
    CHECK(true);
  }

  // Check that empty parameter queries do not work
  //------------------------------------------------------------------
  try
  {
    reader.new_scalar_parameter<double>( 
        "next_double_param", "" );
    CHECK(false);
  }
  catch (...)
  {
    CHECK(true);
  }


  // Check queries
  //------------------------------------------------------------------
  reader.query<bool>( "bool_param" );
  reader.query<int>( "int_param" );
  reader.query<std::string>( "str_param" );

  bool bool_param = reader.get_value<bool>("bool_param");
  CHECK( (bool_param == true) );
  CHECK( reader.found("bool_param") );

  int int_param = reader.get_value<int>("int_param");
  CHECK( (int_param == 1) );
  CHECK( reader.found("int_param") );

  CHECK( !reader.query<int>( "invalid_int_param" ) );
  CHECK( !reader.found("invalid_int_param") );

  std::string str_param = reader.get_value<std::string>("str_param");

  CHECK( reader.found("str_param") );
  CHECK( str_param == "Test" );


  // Check second queries
  //------------------------------------------------------------------
  CHECK( reader.query<int>( "int_param" ) );
  CHECK( reader.found("int_param") );
  int_param = reader.get_value<int>("int_param");
  CHECK( (int_param == 2) );
  CHECK( reader.found("int_param") );


} // scalar_parameters() 


/*--------------------------------------------------------------------
| Test list parameters
--------------------------------------------------------------------*/
static void list_parameters()
{
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };
  std::string file_name {source_dir + "/auxiliary/ParameterFile.txt"};

  ParaReader reader { file_name };

  // Create working parameter definitions
  //------------------------------------------------------------------
  reader.new_matrix_parameter<int>( 
      "para_list_1",  "Parameter list start:", "Parameter list end", 3 );
  reader.new_vector_parameter<double>( 
      "para_list_2",  "List parameter:", 4 );


  reader.new_matrix_parameter<int>( 
      "invalid_list_1",  "Invalid list start:", "Invalid list end", 3 );
  reader.new_vector_parameter<double>( 
      "invalid_list_2",  "Invalid list parameter:", 4 );
  reader.new_matrix_parameter<int>( 
      "valid_list",  "Invalid list start:", "Invalid list end", 3 );


  // Check that invalid parameter definition does not succeed
  //------------------------------------------------------------------
  try
  {
    reader.new_matrix_parameter<int>( 
      "para_list_1",  "Parameter list start:", "Parameter list end:", 3);
    CHECK(false);
  }
  catch (...)
  {
    CHECK(true);
  }

  // Check that empty parameter queries do not work
  //------------------------------------------------------------------
  try
  {
    reader.new_matrix_parameter<int>( 
      "para_list_3",  "", "Parameter list end:", 4 );
    CHECK(false);
  }
  catch (...)
  {
    CHECK(true);
  }
  try
  {
    reader.new_matrix_parameter<int>( 
      "para_list_4",  "Parameter list start:", "", 4 );
    CHECK(false);
  }
  catch (...)
  {
    CHECK(true);
  }
  try
  {
    reader.new_vector_parameter<int>( 
      "para_list_5",  "", 5 );
    CHECK(false);
  }
  catch (...)
  {
    CHECK(true);
  }

  // Check queries of valid lists
  //------------------------------------------------------------------
  reader.query<double>( "para_list_2" );
  CHECK( reader.found("para_list_2") );
  CHECK( EQ(reader.get_value<double>(0, "para_list_2"), 1.0) );
  CHECK( EQ(reader.get_value<double>(1, "para_list_2"), 2.0) );
  CHECK( EQ(reader.get_value<double>(2, "para_list_2"), 3.0) );
  CHECK( EQ(reader.get_value<double>(3, "para_list_2"), 4.0) );
  CHECK( EQ(reader.get_value<double>(4, "para_list_2"), 1.0) );

  reader.query<int>( "para_list_1" );
  CHECK( reader.found("para_list_1") );
  CHECK( reader.get_parameter<int>("para_list_1").columns() == 3 );
  CHECK( reader.get_parameter<int>("para_list_1").rows() == 3 );
  CHECK( reader.get_value<int>(0, "para_list_1") == 1 );
  CHECK( reader.get_value<int>(1, "para_list_1") == 2 );

  CHECK( reader.get_value<int>(1, 1, "para_list_1") == 5 );
  CHECK( reader.get_value<int>(1, 2, "para_list_1") == 8 );


  // Check queries of invalid lists
  //------------------------------------------------------------------
  reader.query<double>( "invalid_list_2" );
  CHECK( reader.found("invalid_list_2") );
  CHECK( EQ(reader.get_value<double>(0, "invalid_list_2"), 1.0) );
  CHECK( EQ(reader.get_value<double>(1, "invalid_list_2"), 2.0) );
  CHECK( EQ(reader.get_value<double>(2, "invalid_list_2"), 4.5) );
  CHECK( EQ(reader.get_value<double>(3, "invalid_list_2"), 5.0) );

  // Check that default values have been added to invalid list
  CHECK( !reader.query<int>( "invalid_list_1" ) );

  CHECK( reader.query<int>( "valid_list", true, -1 ) );
  CHECK( reader.found("valid_list") );
  CHECK( EQ(reader.get_value<int>(0, "valid_list"), 1) );
  CHECK( EQ(reader.get_value<int>(1, "valid_list"), 2) );
  CHECK( EQ(reader.get_value<int>(2, "valid_list"), 3) );
  CHECK( EQ(reader.get_value<int>(3, "valid_list"), -1) );
  CHECK( EQ(reader.get_value<int>(4, "valid_list"), -1) );
  CHECK( EQ(reader.get_value<int>(5, "valid_list"), 3) );
  CHECK( EQ(reader.get_value<int>(6, "valid_list"), -1) );
  CHECK( EQ(reader.get_value<int>(7, "valid_list"), 3) );
  CHECK( EQ(reader.get_value<int>(8, "valid_list"), -1) );
  CHECK( EQ(reader.get_value<int>(9, "valid_list"), 3) );
  CHECK( EQ(reader.get_value<int>(10, "valid_list"), 4) );
  CHECK( EQ(reader.get_value<int>(11, "valid_list"), 5) );
  CHECK( EQ(reader.get_value<int>(12, "valid_list"), -1) );
  CHECK( EQ(reader.get_value<int>(13, "valid_list"), -1) );
  CHECK( EQ(reader.get_value<int>(14, "valid_list"), -1) );


} // list_parameters() 


/*--------------------------------------------------------------------
| Test block parameters
--------------------------------------------------------------------*/
static void block_parameters()
{
  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };
  std::string file_name {source_dir + "/auxiliary/ParameterFile.txt"};

  ParaReader reader { file_name };


  // Create valid parameter definitions
  //------------------------------------------------------------------
  reader.new_block_parameter( 
      "block_para_1",  "Parameter block start:", "Parameter block end" );

  // Create parameters within block
  //------------------------------------------------------------------
  ParaBlock& block_1 = reader.get_block( "block_para_1" );
  block_1.new_scalar_parameter<int>(
      "int_param",  "Int parameter in block:" );
  block_1.new_scalar_parameter<std::string>(
      "str_param", "String parameter in block:");
  block_1.new_vector_parameter<float>(
      "vec_param", "List parameter in block:", 4);
  block_1.new_matrix_parameter<int>( 
      "mat_param", "Parameter list in block start:", 
      "Parameter list in block end", 3 );

  // Check that valid parameter of first block succeeds
  //------------------------------------------------------------------
  CHECK( reader.query( "block_para_1" ) );
  CHECK( reader.get_block( "block_para_1" ).block_start() > 0 );
  CHECK( reader.get_block( "block_para_1" ).block_end() > 0 );
  CHECK( reader.get_block( "block_para_1" ).block_end() > 
         reader.get_block( "block_para_1" ).block_start()  );

  CHECK( block_1.query<int>("int_param") );
  CHECK( block_1.query<std::string>("str_param") );
  CHECK( block_1.query<float>("vec_param") );
  CHECK( block_1.query<int>("mat_param") );


  CHECK( block_1.found("int_param") );
  int int_param = block_1.get_value<int>("int_param");
  CHECK( (int_param == 2) );

  CHECK( block_1.found("str_param") );
  std::string str_param = block_1.get_value<std::string>("str_param");
  CHECK( (str_param == "Test-1") );

  CHECK( block_1.found("vec_param") );
  CHECK( EQ(block_1.get_value<float>(0, "vec_param"), 1.0f) );
  CHECK( EQ(block_1.get_value<float>(1, "vec_param"), 2.0f) );
  CHECK( EQ(block_1.get_value<float>(2, "vec_param"), 3.0f) );
  CHECK( EQ(block_1.get_value<float>(3, "vec_param"), 4.0f) );

  CHECK( block_1.found("mat_param") );
  CHECK( EQ(block_1.get_value<int>(0, 0, "mat_param"), 1) );
  CHECK( EQ(block_1.get_value<int>(1, 0, "mat_param"), 2) );
  CHECK( EQ(block_1.get_value<int>(2, 0, "mat_param"), 3) );
  CHECK( EQ(block_1.get_value<int>(2, 2, "mat_param"), 9) );

  // Check that valid parameter of second block succeeds
  //------------------------------------------------------------------
  std::cout << "\n \n " << std::endl;
  std::cout << "Query next block \n\n";
  CHECK( reader.query( "block_para_1" ) );

  CHECK( block_1.query<int>("int_param") );
  CHECK( block_1.query<std::string>("str_param") );

  str_param = block_1.get_value<std::string>("str_param");
  CHECK( (str_param == "Test-2") );


} // block_parameters()



} // namespace ParaReaderTests  


/*********************************************************************
* Run tests for: ParaReader.h
*********************************************************************/
void run_tests_ParaReader()
{
  ParaReaderTests::constructor();
  ParaReaderTests::scalar_parameters();
  ParaReaderTests::list_parameters();
  ParaReaderTests::block_parameters();

} // run_tests_ParaReader()
