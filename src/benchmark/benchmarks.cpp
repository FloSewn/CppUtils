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

/*********************************************************************
* Color text
*********************************************************************/
#define NC "\e[0m"
#define RED "\e[0;31m"
#define GRN "\e[0;32m"
#define CYN "\e[0;36m"
#define REDB "\e[41m"

/*********************************************************************
* The main benchmark function
*********************************************************************/
int run_benchmarks(const std::string& library)
{
  CppUtils::SimpleLogger MSG(std::clog, "  ");

  /*------------------------------------------------------------------
  | Print header
  ------------------------------------------------------------------*/
  MSG << "" << std::endl;
  MSG << "   -------------------------   " << std::endl;
  MSG << "   | CppUtils - Benchmarks |   " << std::endl;
  MSG << "   -------------------------   " << std::endl;
  MSG << "" << std::endl;

  /*------------------------------------------------------------------
  | Run all benchmarks
  ------------------------------------------------------------------*/
  if ( !library.compare("Container") )
  {
    MSG << "  Running benchmarks for \"Container\" library..." 
            << std::endl;
    run_benchmarks_Container();
  }
  else if ( !library.compare("QuadTree") )
  {
    MSG << "  Running benchmarks for \"QuadTree\" library..." 
            << std::endl;
    run_benchmarks_QuadTree();
  }
  else
  {
    MSG << std::endl;
    MSG << RED "  No library \"" << library 
            << "\" found to benchmark" NC << std::endl;
    MSG << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

} // run_tests()
