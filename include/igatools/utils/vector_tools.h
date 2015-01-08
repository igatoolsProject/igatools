//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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
#ifndef VECTOR_TOOLS_H_
#define VECTOR_TOOLS_H_

#include <igatools/base/config.h>


#include <vector>

IGA_NAMESPACE_OPEN


/**
 * @brief Utility functions for vector
 * @author M.Martinelli
 * @date 2012
 */
namespace vector_tools
{

/**
 * This function counts and removes the duplicates values from a vector.
 * @pre Before calling this function, the entries in the vector @vec_with_duplicates must be sorted.
 * @param[in] vec_with_duplicates Vector with duplicates
 * @param[out] vec_without_duplicates Vector without duplicates
 * @param[out] multiplicities Multiplicities of the values in the
 * <tt>vec_without_duplicates</tt> vector.
 * This vector has the same size of <tt>vec_without_duplicates</tt>.
 * <tt>multiplicities[i]</tt> is the multiplicity of
 * the value <tt>vec_without_duplicates[i]</tt>.
 * If a given value in the input <tt>vec_with_duplicates</tt> vector has no duplicates,
 * the correspondent multiplicity is set to 1.
 */
template< class T >
void
count_and_remove_duplicates(
    const vector<T> &vec_with_duplicates,
    vector<T> &vec_without_duplicates,
    vector<int> &multiplicities) ;

} ;


IGA_NAMESPACE_CLOSE



#endif // VECTOR_TOOLS_H_


#include <igatools/utils/vector_tools-inline.h>
