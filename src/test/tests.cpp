#include <iostream>
#include <vector>
#include <string>

#include "tests.h"
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
* The main test function
*********************************************************************/
int run_tests(const std::string& library)
{
  CppUtils::SimpleLogger TESTMSG(std::clog, "  ");

  /*------------------------------------------------------------------
  | Print header
  ------------------------------------------------------------------*/
  TESTMSG << "" << std::endl;
  TESTMSG << "   -------------------------   " << std::endl;
  TESTMSG << "   | CppUtils - Test suite |   " << std::endl;
  TESTMSG << "   -------------------------   " << std::endl;
  TESTMSG << "" << std::endl;

  /*------------------------------------------------------------------
  | Run all tests
  ------------------------------------------------------------------*/
  if ( !library.compare("MathUtility") )
  {
    TESTMSG << "  Running tests for \"MathUtility\" library..."
            << std::endl;
    run_tests_MathUtility();
  }
  else if ( !library.compare("StringOps") )
  {
    TESTMSG << "  Running tests for \"StringOps\" library..."
            << std::endl;
    run_tests_StringOps();
  }
  else if ( !library.compare("Vec2") )
  {
    TESTMSG << "  Running tests for \"Vec2\" library..."
            << std::endl;
    run_tests_Vec2();
  }
  else if ( !library.compare("Geometry") )
  {
    TESTMSG << "  Running tests for \"Geometry\" library..."
            << std::endl;
    run_tests_Geometry();
  }
  else if ( !library.compare("QuadTree") )
  {
    TESTMSG << "  Running tests for \"QuadTree\" library..." 
            << std::endl;
    run_tests_QuadTree();
  }
  else
  {
    TESTMSG << std::endl;
    TESTMSG << RED "  No library \"" << library 
            << "\" found to test" NC << std::endl;
    TESTMSG << std::endl;
    return EXIT_FAILURE;
  }

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
      TESTMSG << RED "[ERROR] Test (" << error_count 
              << "/" << total_tests << ") failed." NC << std::endl;
      TESTMSG << "        --> " << data << std::endl;
    }
    state &= data.state();
  }

  /*------------------------------------------------------------------
  | Succeess / fail
  ------------------------------------------------------------------*/
  TESTMSG << "" << std::endl;
  if (!state)
  {
    TESTMSG << RED "  --> (" << error_count << "/" 
            << total_tests << ") tests failed." NC  << std::endl;
  }
  else
  {
    TESTMSG << GRN "  --> (" << total_tests-error_count << "/" 
            << total_tests << ") tests succeeded." NC << std::endl;
  }
  TESTMSG << std::endl;
  TESTMSG << std::endl;

  if (!state)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

} // run_tests()
