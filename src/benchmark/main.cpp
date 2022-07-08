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

#include "Helpers.h"

#include "benchmarks.h"

/*********************************************************************
* Main function
*********************************************************************/
int main(int argc, char *argv[])
{
  CppUtils::SimpleLogger MSG(std::clog, "  ");
  
  if ( argc < 2 )
  {
    MSG << "" << std::endl;
    MSG << "   -------------------------   " << std::endl;
    MSG << "   | CppUtils - Benchmarks |   " << std::endl;
    MSG << "   -------------------------   " << std::endl;
    MSG << "" << std::endl;
    MSG << "Usage: " << argv[0] << " <library-name-to-benchmark>" 
              << std::endl;
    MSG << "" << std::endl;
    MSG << "" << std::endl;
    return EXIT_FAILURE;
  }

  std::string input { argv[1] };

  return run_benchmarks( input );
}
