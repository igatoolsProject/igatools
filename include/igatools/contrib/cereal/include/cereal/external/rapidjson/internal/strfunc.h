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
#ifndef RAPIDJSON_INTERNAL_STRFUNC_H_
#define RAPIDJSON_INTERNAL_STRFUNC_H_

namespace rapidjson
{
namespace internal
{

//! Custom strlen() which works on different character types.
/*! \tparam Ch Character type (e.g. char, wchar_t, short)
  \param s Null-terminated input string.
  \return Number of characters in the string.
  \note This has the same semantics as strlen(), the return value is not number of Unicode codepoints.
*/
template <typename Ch>
inline SizeType StrLen(const Ch *s)
{
  const Ch *p = s;
  while (*p != '\0')
    ++p;
  return SizeType(p - s);
}

} // namespace internal
} // namespace rapidjson

#endif // RAPIDJSON_INTERNAL_STRFUNC_H_
