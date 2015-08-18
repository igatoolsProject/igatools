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

#ifndef __BASE_ELEMENT_H_
#define __BASE_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_index.h>
#include <igatools/utils/safe_stl_vector.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Cartesian Product Element iterator (without properties)
 *
 * @author pauletti, 2015
 *
 * @ingroup serializable
 */
template <int dim>
using BaseElement = TensorIndex<dim>;

/**
 * Generates a vector with the tensor indices of the given
 * rectangular range.
 *
 *  @relates TensorIndex
 *  @author pauletti 2015
 */
template<int k>
SafeSTLSet<BaseElement<k>> el_tensor_range(TensorIndex<k> first, TensorIndex<k> last)
{
    Assert(first <= last, ExcMessage("first bigger than last"));
    SafeSTLSet<BaseElement<k>> result;
    TensorIndex<k-1> ind(sequence<k-1>());
    auto vec = el_tensor_range<k-1>(first.get_sub_tensor(ind), last.get_sub_tensor(ind));

    for (int i=first[k-1]; i<last[k-1]; ++i)
    {
        for (auto &t_k_1 : vec)
        {
            BaseElement<k> t_k;
            for (int j=0; j<k-1; ++j)
                t_k[j] = t_k_1[j];
            t_k[k-1] = i;
            result.insert(t_k);
        }
    }
    return result;
}

template<>
SafeSTLSet<BaseElement<1> >
el_tensor_range(TensorIndex<1> first, TensorIndex<1> last);
template<>
SafeSTLSet<BaseElement<0> >
el_tensor_range(TensorIndex<0> first, TensorIndex<0> last);



IGA_NAMESPACE_CLOSE

#endif
