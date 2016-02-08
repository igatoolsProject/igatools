//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
//
// This file is part of the igatools library.
//
// The igatools library is free software: you can use it, redistribute
// it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-+--------------------------------------------------------------------
#ifndef RAPIDJSON_FILESTREAM_H_
#define RAPIDJSON_FILESTREAM_H_

#include <cstdio>

namespace rapidjson
{

//! Wrapper of C file stream for input or output.
/*!
  This simple wrapper does not check the validity of the stream.
  \implements Stream
*/
class FileStream
{
public:
  typedef char Ch;  //!< Character type. Only support char.

  FileStream(FILE *fp) : fp_(fp), count_(0)
  {
    Read();
  }

  char Peek() const
  {
    return current_;
  }
  char Take()
  {
    char c = current_;
    Read();
    return c;
  }
  size_t Tell() const
  {
    return count_;
  }
  void Put(char c)
  {
    fputc(c, fp_);
  }

  // Not implemented
  char *PutBegin()
  {
    return 0;
  }
  size_t PutEnd(char *)
  {
    return 0;
  }

private:
  void Read()
  {
    RAPIDJSON_ASSERT(fp_ != 0);
    int c = fgetc(fp_);
    if (c != EOF)
    {
      current_ = (char)c;
      count_++;
    }
    else
      current_ = '\0';
  }

  FILE *fp_;
  char current_;
  size_t count_;
};

} // namespace rapidjson

#endif // RAPIDJSON_FILESTREAM_H_
