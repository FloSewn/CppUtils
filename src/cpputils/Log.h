/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once


#include <fstream>
#include <string>
#include <memory>

namespace CppUtils {

/*********************************************************************
* Log levels
*********************************************************************/
enum LogLevel 
{ ERROR, WARNING, INFO, DEBUG };

enum OStreamType
{ TO_COUT, TO_CERR, TO_CLOG, TO_FILE };

/*********************************************************************
* Interface to create ostream unique_ptr, which gets properly 
* deleted.
*
* Reference:
* ----------
* -https://stackoverflow.com/questions/56521318/ostream-class-\
*  that-outputs-either-on-cout-or-on-a-file
*
*********************************************************************/
struct ConditionalDeleter
{
    bool must_delete;
    void operator()(std::ostream* os) const 
    { if (must_delete) delete os; }
};

using OStreamPtr = std::unique_ptr<std::ostream, ConditionalDeleter>;

OStreamPtr create_stream(OStreamType type, const std::string& path="")
{
  switch( type ) {
    case TO_COUT:
      return OStreamPtr { &std::cout, ConditionalDeleter {false} };
    case TO_CERR: 
      return OStreamPtr { &std::cerr, ConditionalDeleter {false} };
    case TO_CLOG:
      return OStreamPtr { &std::clog, ConditionalDeleter {false} };
    case TO_FILE: 
      return OStreamPtr { new std::ofstream {path}, 
                          ConditionalDeleter {true} };
  }
  return OStreamPtr { &std::cout, ConditionalDeleter {false} };
}


/*********************************************************************
* The global logging properties
*********************************************************************/
class LogProperties
{
public:
  LogProperties() 
  {
    error_os_ = create_stream( TO_COUT );
    warn_os_  = create_stream( TO_COUT );
    info_os_  = create_stream( TO_COUT );
    debug_os_ = create_stream( TO_COUT );
  }

  void set_level(LogLevel level) { level_ = level; }
  void show_header(bool show) { show_header_ = show; }
  void set_error_header(const std::string& msg) { error_header_ = msg; }
  void set_warn_header(const std::string& msg) { warn_header_ = msg; }
  void set_info_header(const std::string& msg) { info_header_ = msg; }
  void set_debug_header(const std::string& msg) { debug_header_ = msg; }

  void set_error_ostream(OStreamType type, const std::string& f="")
  { error_os_ = create_stream( type, f ); }
  void set_warn_ostream(OStreamType type, const std::string& f="")
  { warn_os_ = create_stream( type, f ); }
  void set_info_ostream(OStreamType type, const std::string& f="")
  { info_os_ = create_stream( type, f ); }
  void set_debug_ostream(OStreamType type, const std::string& f="")
  { debug_os_ = create_stream( type, f ); }

  const LogLevel& level() const { return level_; }
  bool show_header() const { return show_header_; }

  const std::string& get_header(LogLevel level) const 
  {
    switch( level ) {
      case ERROR:   return error_header_; 
      case WARNING: return warn_header_; 
      case INFO:    return info_header_; 
      case DEBUG:   return debug_header_; 
    }
    return info_header_;
  }

  std::ostream& get_ostream(LogLevel level)
  {
    switch( level ) {
      case ERROR:   return *error_os_; 
      case WARNING: return *warn_os_; 
      case INFO:    return *info_os_; 
      case DEBUG:   return *debug_os_; 
    }
    return *info_os_;
  }


private:

  LogLevel    level_         = INFO;
  bool        show_header_   = true;
  std::string error_header_  = "[ERROR] ";
  std::string warn_header_   = "[WARNING] ";
  std::string info_header_   = "[INFO] ";
  std::string debug_header_  = "[DEBUG] ";

  OStreamPtr error_os_ { nullptr };
  OStreamPtr warn_os_  { nullptr };
  OStreamPtr info_os_  { nullptr };
  OStreamPtr debug_os_ { nullptr };

};

static inline LogProperties LOG_PROPERTIES;

/*********************************************************************
* The interface for the actual SimpleLogger
*
* Reference:
* ----------
* -https://stackoverflow.com/questions/5028302/small-logger-class
*********************************************************************/
class LOG
{
public:
  LOG() {}

  LOG(LogLevel level)
  {
    level_ = level;
    if ( LOG_PROPERTIES.show_header() )
      operator << ( LOG_PROPERTIES.get_header( level_ ) );
  }

  ~LOG() 
  {
    if ( opened_ )
      LOG_PROPERTIES.get_ostream( level_ ) << std::endl;
    opened_ = false;
  }

  template<class T>
  LOG& operator<<(const T& msg)
  {
    if ( level_ <= LOG_PROPERTIES.level() )
    {
      LOG_PROPERTIES.get_ostream( level_ ) << msg;
      opened_ = true;
    }
    return *this;
  }

private:

  bool     opened_ = false;
  LogLevel level_  = DEBUG;

}; // LOG

} // namespace CppUtils
