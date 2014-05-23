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

#ifndef VECTOR_TOOLS_INLINE_H_
#define VECTOR_TOOLS_INLINE_H_


#include <igatools/utils/vector_tools.h>

#include <algorithm>

using std::vector ;


IGA_NAMESPACE_OPEN

namespace vector_tools
{

template< class T >
inline
void
count_and_remove_duplicates(
    const vector<T> &vec_with_duplicates,
    vector<T> &vec_without_duplicates,
    vector<int> &multiplicities)
{
    //------------------------------------------------------------------------------------------
    const auto vec_with_duplicates_begin = vec_with_duplicates.cbegin() ;
    const auto vec_with_duplicates_end   = vec_with_duplicates.cend() ;


    vec_without_duplicates.resize(vec_with_duplicates.size()) ;
    auto vec_without_duplicates_begin = vec_without_duplicates.begin() ;

    // here we remove the duplicate values
    auto it = std::unique_copy(vec_with_duplicates_begin,
                               vec_with_duplicates_end,
                               vec_without_duplicates_begin) ;

    // resizing the vector of values without repetition
    vec_without_duplicates.resize(it - vec_without_duplicates_begin) ;
    //------------------------------------------------------------------------------------------


    //------------------------------------------------------------------------------------------
    // now we count how many entries with the same values (==>multiplicity)
    multiplicities.clear() ;
    for (const auto & v : vec_without_duplicates)
        multiplicities.push_back(count(vec_with_duplicates_begin, vec_with_duplicates_end, v)) ;
    //------------------------------------------------------------------------------------------

}

} ;


IGA_NAMESPACE_CLOSE

#endif //#ifndef VECTOR_TOOLS_INLINE_H_
