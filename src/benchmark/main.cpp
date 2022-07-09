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
#include "Log.h"

#include "benchmarks.h"

using CppUtils::LOG_PROPERTIES;
using CppUtils::LOG;
using CppUtils::LogLevel::INFO;

/*********************************************************************
* Main function
*********************************************************************/
int main(int argc, char *argv[])
{
  LOG_PROPERTIES.set_level( INFO );
  LOG_PROPERTIES.show_header( true );
  LOG_PROPERTIES.set_info_header( "  " );
  
  if ( argc < 2 )
  {
    LOG(INFO) << "";
    LOG(INFO) << "   -------------------------   ";
    LOG(INFO) << "   | CppUtils - Benchmarks |   ";
    LOG(INFO) << "   -------------------------   ";
    LOG(INFO) << "";
    LOG(INFO) << "Usage: " << argv[0] << " <library-name-to-benchmark>";
    LOG(INFO) << "";
    LOG(INFO) << "";
    return EXIT_FAILURE;
  }

  std::string input { argv[1] };

  return run_benchmarks( input );
}
