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

#include "benchmarks.h"
#include "Helpers.h"
#include "Log.h"

/*********************************************************************
* Log utils
*********************************************************************/
using CppUtils::LOG;
using CppUtils::LogLevel::INFO;
using CppUtils::LogColor::GREEN;
using CppUtils::LogColor::RED;

using CppUtils::LOG;
using CppUtils::LogLevel::INFO;

/*********************************************************************
* The main benchmark function
*********************************************************************/
int run_benchmarks(const std::string& library)
{
  /*------------------------------------------------------------------
  | Print header
  ------------------------------------------------------------------*/
  LOG(INFO) << "";
  LOG(INFO) << "   -------------------------   ";
  LOG(INFO) << "   | CppUtils - Benchmarks |   ";
  LOG(INFO) << "   -------------------------   ";
  LOG(INFO) << "";

  /*------------------------------------------------------------------
  | Run all benchmarks
  ------------------------------------------------------------------*/
  if ( !library.compare("Container") )
  {
    LOG(INFO) << "  Running benchmarks for \"Container\" library...";
    run_benchmarks_Container();
  }
  else if ( !library.compare("QuadTree") )
  {
    LOG(INFO) << "  Running benchmarks for \"QuadTree\" library...";
    run_benchmarks_QuadTree();
  }
  else if ( !library.compare("RTreeND") )
  {
    LOG(INFO) << "  Running benchmarks for \"RTreeND\" library...";
    run_benchmarks_RTreeND();
  }
  else if ( !library.compare("OcTreeND") )
  {
    LOG(INFO) << "  Running benchmarks for \"OcTreeND\" library...";
    run_benchmarks_OcTreeND();
  }
  else
  {
    LOG(INFO) << "";
    LOG(INFO, RED) << "  No library \"" << library 
                   << "\" found to benchmark";
    LOG(INFO) << "";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

} // run_tests()
