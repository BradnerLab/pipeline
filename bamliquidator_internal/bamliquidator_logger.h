#ifndef PIPELINE_BAMLIQUIDATORINTERNAL_BAMLIQUIDATOR_LOGGER_H
#define PIPELINE_BAMLIQUIDATORINTERNAL_BAMLIQUIDATOR_LOGGER_H

#include <ostream>
#include <sstream>
#include <string>

// these functions are intended to match bamliquidator_batch.py logging output style

class Logger
{
public:
  // Configures logging.  If not configured, logs are just written to stderr.
  static void configure(const std::string& log_file_path,
                        bool include_warnings_in_stderr);
  
  /* e.g. Logger::warn() << "oops " << 123 results in a logged line like the following written to the log file:
   *
   *  2014-08-05 13:25:06 WARNING	oops 123
   * 
   * and if configured to log to stderr, the following line written there as well:
   *
   *   WARNING	oops 123
   *
   * Since copy constructor is private, returned value must be used as an anonymous temporary as in the example.
   */ 
  static Logger warn();

  /* e.g. Logger::error() << "oops " << 123 results in a logged line like the following written to the log file:
   *
   *  2014-08-05 13:25:06 ERROR	oops 123
   * 
   * and the following written to stderr as well:
   *
   *   ERROR	oops 123
   *
   * Since copy constructor is private, returned value must be used as an anonymous temporary as in the example.
   */
  static Logger error();

  template<typename T>
  Logger& operator<<(const T& v)
  {
    ss << v;
    return *this;
  }

  ~Logger();

private:
  Logger(const std::string& level, bool write_to_stderr);

  // Prevent copy construction so callers of warn() and error() can only use returned object as an anonymous
  // temporary with <<, which is appropriate since the actual logging occurs on destruction.
  Logger(const Logger& logger);

  Logger& operator=( const Logger& ) = delete;

  const std::string level;
  const bool write_to_stderr;
  std::stringstream ss;
  mutable bool copied;
};

/* The MIT License (MIT) 

   Copyright (c) 2014 John DiMatteo (jdimatteo@gmail.com)

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE. 
 */

#endif
