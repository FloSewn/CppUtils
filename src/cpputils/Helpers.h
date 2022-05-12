/*
* This source file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <cassert>
#include <iostream>
#include <string>

namespace CppUtils {

/*********************************************************************
* Asserts messages
*********************************************************************/
#ifndef NDEBUG
static inline void ASSERT(bool cond, const std::string& msg)
{
  if (!cond)
  {
    std::cerr << "[ERROR] " << msg << std::endl;
    assert(cond);
  }
}
#else
static inline void ASSERT(bool cond, const std::string& msg)
{}
#endif

/*********************************************************************
* A simple logging device
*********************************************************************/
class SimpleLogger
{
public:
  SimpleLogger(std::ostream& o, std::string s)
  : out_(o)
  , msg_(std::move(s))
  {}
  
  template <typename T>
  std::ostream& operator<<(T const& x)
  { return out_ << msg_ << x; }

private:
  std::ostream& out_;
  std::string const msg_;

}; // SimpleLogger

//SimpleLogger MSG(std::clog, "");
//SimpleLogger TESTMSG(std::clog, "[TEST] ");

} // namespace CppUtils
