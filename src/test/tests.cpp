/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <vector>
#include <string>

#include "tests.h"
#include "Helpers.h"
#include "Log.h"
#include "Testing.h"

/*********************************************************************
* Log utils
*********************************************************************/
using CppUtils::LOG;
using CppUtils::LogLevel::INFO;
using CppUtils::LogColor::GREEN;
using CppUtils::LogColor::RED;

/*********************************************************************
* The main test function
*********************************************************************/
int run_tests(const std::string& library)
{
  /*------------------------------------------------------------------
  | Print header
  ------------------------------------------------------------------*/
  LOG(INFO) << "";
  LOG(INFO) << "   -------------------------   ";
  LOG(INFO) << "   | CppUtils - Test suite |   ";
  LOG(INFO) << "   -------------------------   ";
  LOG(INFO) << "";

  /*------------------------------------------------------------------
  | Run all tests
  ------------------------------------------------------------------*/
  if ( !library.compare("MathUtility") )
  {
    LOG(INFO) << "  Running tests for \"MathUtility\" library...";
    run_tests_MathUtility();
  }
  else if ( !library.compare("StringOps") )
  {
    LOG(INFO) << "  Running tests for \"StringOps\" library...";
    run_tests_StringOps();
  }
  else if ( !library.compare("Vec2") )
  {
    LOG(INFO) << "  Running tests for \"Vec2\" library...";
    run_tests_Vec2();
  }
  else if ( !library.compare("Geometry") )
  {
    LOG(INFO) << "  Running tests for \"Geometry\" library...";
    run_tests_Geometry();
  }
  else if ( !library.compare("QuadTree") )
  {
    LOG(INFO) << "  Running tests for \"QuadTree\" library...";
    run_tests_QuadTree();
  }
  else if ( !library.compare("Container") )
  {
    LOG(INFO) << "  Running tests for \"Container\" library...";
    run_tests_Container();
  }
  else if ( !library.compare("ParaReader") )
  {
    LOG(INFO) << "  Running tests for \"ParaReader\" library...";
    run_tests_ParaReader();
  }
  else if ( !library.compare("VtkIO") )
  {
    LOG(INFO) << "  Running tests for \"VtkIO\" library...";
    run_tests_VtkIO();
  }
  else if ( !library.compare("Log") )
  {
    LOG(INFO) << "  Running tests for \"Log\" library...";
    run_tests_Log();
  }
  else if ( !library.compare("Matrix") )
  {
    LOG(INFO) << "  Running tests for \"Matrix\" library...";
    run_tests_Matrix();
  }
  else if ( !library.compare("BTree") )
  {
    LOG(INFO) << "  Running tests for \"BTree\" library...";
    run_tests_BTree();
  }
  else if ( !library.compare("RTree") )
  {
    LOG(INFO) << "  Running tests for \"RTree\" library...";
    run_tests_RTree();
  }
  else if ( !library.compare("VecND") )
  {
    LOG(INFO) << "  Running tests for \"VecND\" library...";
    run_tests_VecND();
  }
  else if ( !library.compare("BBoxND") )
  {
    LOG(INFO) << "  Running tests for \"BBoxND\" library...";
    run_tests_BBoxND();
  }
  else
  {
    LOG(INFO) << "";
    LOG(INFO, RED) << "  No library \"" << library 
                   << "\" found to test";
    LOG(INFO) << "";
    return EXIT_FAILURE;
  }

  /*------------------------------------------------------------------
  | Check for failed tests
  ------------------------------------------------------------------*/
  std::vector<CppUtils::TestData>& test_data 
    = CppUtils::TestDataSingleton::instance();

  bool   state = true;
  size_t error_count = 0;
  size_t total_tests = test_data.size();

  for (auto data : test_data )
  {
    if ( !data.state() )
    {
      ++error_count;
      LOG(INFO, RED) << "[ERROR] Test (" << error_count 
              << "/" << total_tests << ") failed.";
      LOG(INFO) << "        --> " << data;
    }
    state &= data.state();
  }

  /*------------------------------------------------------------------
  | Succeess / fail
  ------------------------------------------------------------------*/
  LOG(INFO) << "";
  if (!state)
  {
    LOG(INFO, RED) << "  --> (" << error_count << "/" 
            << total_tests << ") tests failed.";
  }
  else
  {
    LOG(INFO, GREEN) << "  --> (" << total_tests-error_count << "/" 
            << total_tests << ") tests succeeded.";
  }
  LOG(INFO) << "";
  LOG(INFO) << "";

  if (!state)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

} // run_tests()
