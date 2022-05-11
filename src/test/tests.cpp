#include <iostream>
#include <vector>
#include <string>

#include "tests.h"


/*********************************************************************
* Test output format
*********************************************************************/
#define TESTMSG(str)                                                 \
  do { std::clog << str << std::endl; } while(false)
  
/*********************************************************************
* Color text
*********************************************************************/
#define NC "\e[0m"
#define RED "\e[0;31m"
#define GRN "\e[0;32m"
#define CYN "\e[0;36m"
#define REDB "\e[41m"


/*********************************************************************
* The main test function
*********************************************************************/
int run_tests()
{
  /*------------------------------------------------------------------
  | Print header
  ------------------------------------------------------------------*/
  TESTMSG("");
  TESTMSG("   -------------------------   " );
  TESTMSG("   | CppUtils - Test suite |   " );
  TESTMSG("   -------------------------   " );
  TESTMSG("");

  /*------------------------------------------------------------------
  | Run all tests
  ------------------------------------------------------------------*/
  run_tests_Vec2();

  /*------------------------------------------------------------------
  | Check for failed tests
  ------------------------------------------------------------------*/
  std::vector<TestData>& test_data = TestDataSingleton::instance();

  bool   state = true;
  size_t error_count = 0;
  size_t total_tests = test_data.size();

  for (auto data : test_data )
  {
    if ( !data.state() )
    {
      ++error_count;
      TESTMSG( RED "[ERROR] Test (" << error_count 
                << "/" << total_tests << ") failed.\n" NC
                << "        --> " << data );
    }
    state &= data.state();
  }

  /*------------------------------------------------------------------
  | Succeess / fail
  ------------------------------------------------------------------*/
  TESTMSG("");
  if (!state)
  {
    TESTMSG( RED "  --> (" << error_count << "/" 
             << total_tests << ") tests failed." NC );
  }
  else
  {
    TESTMSG( GRN "  --> (" << total_tests-error_count << "/" 
             << total_tests << ") tests succeeded." NC );
  }
  TESTMSG("\n");

  if (!state)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

} // run_tests()
