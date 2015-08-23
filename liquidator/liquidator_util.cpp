#include "liquidator_util.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <unistd.h>

namespace
{
  std::ofstream log_file;
  
  bool include_warnings(true);
}

namespace liquidator
{

void Logger::configure(const std::string& log_file_path, bool include_warnings_in_stderr)
{
  log_file.open(log_file_path.c_str(), std::ios::app);
  include_warnings = include_warnings_in_stderr;
}

Logger Logger::warn()
{
  return Logger("WARNING", include_warnings);
}

Logger Logger::error()
{
  return Logger("ERROR", true);
}

Logger::Logger(const std::string& a_level, bool a_write_to_stderr):
  level(a_level),
  write_to_stderr(a_write_to_stderr),
  copied(false)
{}

Logger::~Logger()
{
  // rvo should make the copied check unnecessary, but just in case...
  if (copied) return;

  try
  {
    if (write_to_stderr)
    {
      std::cerr << level << '\t' << ss.str() << std::endl;
    }

    std::time_t t = std::time(NULL);
    const char datefmt[] = "%Y-%m-%d %H:%M:%S ";
    char buffer[sizeof(datefmt)*4];
    const size_t written_bytes = std::strftime(buffer, sizeof(buffer), datefmt, std::localtime(&t));

    log_file << (written_bytes == 0 ? "" : buffer) << level << '\t' << ss.str() << std::endl;
  }
  catch(...) {} // don't let destructor throw
}

Logger::Logger(const Logger& logger):
  level(logger.level),
  write_to_stderr(logger.write_to_stderr)
{
  ss << logger.ss.str();
  logger.copied = true;
}

}

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
