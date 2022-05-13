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

#include "tests.h"

/*********************************************************************
* Main function
*********************************************************************/
int main(int argc, char *argv[])
{
  CppUtils::SimpleLogger TESTMSG(std::clog, "  ");

  if ( argc < 2 )
  {
    TESTMSG << std::endl;
    TESTMSG << "   -------------------------   " << std::endl;
    TESTMSG << "   | CppUtils - Test suite |   " << std::endl;
    TESTMSG << "   -------------------------   " << std::endl;
    TESTMSG << std::endl;
    TESTMSG << "Usage: " << argv[0] << " <library-name-to-test>" 
              << std::endl;
    TESTMSG << std::endl;
    TESTMSG << std::endl;
    return EXIT_FAILURE;
  }

  std::string input { argv[1] };

  return run_tests( input );
}
