/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <iostream>
#include <cassert>

#include "CppUtilsConfig.h"

#include "tests.h"

#include "Testing.h"
#include "Log.h"

namespace LogTests 
{

using namespace CppUtils;

/*--------------------------------------------------------------------
| Test writing
--------------------------------------------------------------------*/
void fiddle()
{
  LOG_PROPERTIES.set_level( DEBUG );

  LOG_PROPERTIES.show_header( true );
  LOG_PROPERTIES.set_info_header( "# " );
  LOG_PROPERTIES.set_debug_header( "# " );

  std::string source_dir { CPPUTILSCONFIG__SOURCE_DIR };
  std::string file_name {source_dir + "/auxiliary/Log_Messages.log"};
  LOG_PROPERTIES.set_debug_ostream( TO_FILE, file_name );

  LOG(ERROR) << "This is an error message.";

  LOG(WARNING) << "This is a warn message.";

  LOG(INFO) << "This is an info message.";

  LOG(DEBUG) << "This is a debug message.";

  DEBUG_LOG( "ANOTHER DEBUG TEST " << 1 << " WORKS" );

  LOG_PROPERTIES.set_info_header( "  " );
  LOG_PROPERTIES.set_debug_header( "  " );

} // fiddle()


} // namespace LogTests

/*********************************************************************
* Run tests for: Log.h
*********************************************************************/
void run_tests_Log()
{
  //LogTests::fiddle();

} // run_tests_Log()
