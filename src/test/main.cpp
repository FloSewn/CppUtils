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

#include "tests.h"

//using CppUtils::LOG_PROPERTIES;
//using CppUtils::LOG;
//using CppUtils::LogLevel::INFO;

using namespace CppUtils;

/*********************************************************************
* Main function
*********************************************************************/
int main(int argc, char *argv[])
{
  LOG_PROPERTIES.set_level( DEBUG );
  LOG_PROPERTIES.show_header( true );
  LOG_PROPERTIES.set_info_header( "  " );
  LOG_PROPERTIES.set_debug_header( "  " );

  if ( argc < 2 )
  {
    LOG(INFO) << "";
    LOG(INFO) << "   -------------------------   ";
    LOG(INFO) << "   | CppUtils - Test suite |   ";
    LOG(INFO) << "   -------------------------   ";
    LOG(INFO) << "";
    LOG(INFO) << "Usage: " << argv[0] << " <library-name-to-test>" ;
    LOG(INFO) << "";
    LOG(INFO) << "";
    return EXIT_FAILURE;
  }

  std::string input { argv[1] };

  return run_tests( input );
}
