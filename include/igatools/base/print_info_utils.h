//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

#ifndef __PRINT_INFO_UTILS_H_
#define __PRINT_INFO_UTILS_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <type_traits>

IGA_NAMESPACE_OPEN
/**
 * Type traits to determine if a class provides
 * a print_info function
 */
template<class T>
using print_info_type =
    decltype(std::declval<T>().print_info(std::declval<LogStream &>()));


template<class T>
constexpr
EnableIf<std::is_void<print_info_type<T>>::value, bool >
                                       has_print_info(int)
{
    return true;
}

template<class T>
constexpr bool
has_print_info(long)
{
    return false;
}

///**
// * Ouput for vector onto a LogStream.
// * Mostly use for debugging.
// *
// * @relates LogStream
// */
//template <class T>
//EnableIf<has_print_info<T>(0), LogStream &>
//operator<<(LogStream &out, const vector<T> &vec)
//{
//    out << "Vector with: " << vec.size() << " entries." << std::endl;
//    for (auto &entry : vec)
//    {
//        entry.print_info(out);
//        out<<std::endl;
//    }
//    return out;
//}

//template <class T>
//EnableIf<!has_print_info<T>(0), LogStream &>
//operator<<(LogStream &out, const vector<T> &vector)
//{
//    out << "[ ";
//    for (auto &i:vector)
//        out << i << " ";
//    out << "]";
//    return out;
//}
//


IGA_NAMESPACE_CLOSE

#endif
