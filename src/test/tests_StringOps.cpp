#include <iostream>
#include <string>

#include "tests.h"

#include "StringOps.h"
#include "MathUtility.h"

namespace StringOpsTests 
{


/*--------------------------------------------------------------------
| Test split()
--------------------------------------------------------------------*/
static void split()
{
  std::string test_1 { "This is a string" };

  std::vector<std::string> split_test_1 = CppUtils::split(test_1, ' ');

  CHECK( (split_test_1.size() == 4) );
  CHECK( (split_test_1[0].compare("This") == 0) );
  CHECK( (split_test_1[1].compare("is") == 0) );
  CHECK( (split_test_1[2].compare("a") == 0) );
  CHECK( (split_test_1[3].compare("string") == 0) );



  std::string test_2 { "Ends with whitespace " };

  auto split_test_2 = CppUtils::split(test_2, ' ');

  CHECK( (split_test_2.size() == 3) );
  CHECK( (split_test_2[2].compare("whitespace") == 0) );



  std::string test_3 { " Starts with whitespace" };

  auto split_test_3 = CppUtils::split(test_3, ' ');

  CHECK( (split_test_3.size() == 3) );



  std::string test_4 { "Has__two__delimiters" };

  auto split_test_4 = CppUtils::split(test_4, '_');

  CHECK( (split_test_4.size() == 3) );

  // Create empty elements between delimiters
  auto split_test_5 = CppUtils::split(test_4, '_', false);

  CHECK( (split_test_5.size() == 5) );

} // split()

/*--------------------------------------------------------------------
| Test converter
--------------------------------------------------------------------*/
static void convert()
{
  std::string s_i {"1"};
  int i = CppUtils::sto(s_i);
  CHECK( (i == 1) );

  std::string s_f {"1.2"};
  float f = CppUtils::sto(s_f);
  CHECK( CppUtils::EQ(f, 1.2f) );

  std::string s_d {"1.2"};
  double d = CppUtils::sto(s_d);
  CHECK( CppUtils::EQ(d, 1.2) );

} // convert()


} // namespace StringOps 


/*********************************************************************
* Run tests for: StringOps.h
*********************************************************************/
void run_tests_StringOps()
{
  StringOpsTests::split();
  StringOpsTests::convert();

} // run_tests_StringOps()
