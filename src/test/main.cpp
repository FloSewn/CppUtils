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

  if ( argc < 2 )
  {
    std::cout << std::endl;
    std::cout << "   -------------------------   " << std::endl;
    std::cout << "   | CppUtils - Test suite |   " << std::endl;
    std::cout << "   -------------------------   " << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: " << argv[0] << " <library-name-to-test>" 
              << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    return EXIT_FAILURE;
  }

  std::string input { argv[1] };

  return run_tests( input );
}
