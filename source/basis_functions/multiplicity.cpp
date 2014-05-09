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


#include <igatools/basis_functions/multiplicity.h>

using std::vector;
using std::array;

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
Multiplicity<dim, range, rank>::
Multiplicity(std::shared_ptr<const Grid> knots,
             const DegreeTable &deg,
             const bool max_reg)
    :
    parent_t::StaticMultiArray(T(knots->get_num_knots_dim())),
    deg_(deg)
{
    fill_max_regularity();
}

//template <int dim>
//auto
//Multiplicity<dim>::
//accumulate() -> parent_t
//{
//    const TensorSize<dim> size = this->tensor_size();
//    parent_t result(size);
//
//    for (int i = 0; i < dim; ++i)
//    {
//        result.entry(i, 0) =  this->data_[i][0];
//
//        const Size size_i = size(i);
//        for (int k = 1 ; k < size_i ; ++k)
//            result.entry(i, k) = result.entry(i, k-1) + this->data_[i][k];
//    }
//
//    return result;
//}
//



template<int dim, int range, int rank>
auto Multiplicity<dim, range, rank>::
compute_index_space_offset() -> parent_t
{
    parent_t res;
    auto res_it = res.begin();
    auto mult_it = this->begin();
    auto end = this->end();
    for (; mult_it != end; ++mult_it, ++res_it)
    {
        auto size =mult_it->tensor_size();
        T comp(size);
        for (int i = 0; i < dim; ++i)
        {
            comp.entry(i, 0) =  0;
            const Size size_i = size(i);
            for (int k=1; k < size_i ; ++k)
                comp.entry(i, k) = comp.entry(i, k-1) + mult_it->entry(i,k);
        }
        *(res_it) = comp;
    }

    return res;
//    const TensorSize<dim> size = this->tensor_size();
//    parent_t result(size);
//
//    for (int i = 0; i < dim; ++i)
//    {
//        result.entry(i, 0) =  this->data_[i][0] - degree[i] - 1;
//
//        const Size size_i = size(i);
//        for (int k=1; k < size_i ; ++k)
//            result.entry(i, k) = result.entry(i, k-1) + this->data_[i][k];
//    }
//
//    return result;
}




template<int dim, int range, int rank>
void Multiplicity<dim, range, rank>::
fill_max_regularity()
{
    auto deg_it  = deg_.begin();
    auto end     = deg_.end();
    auto mult_it = this->begin();

    for (; deg_it != end; ++deg_it, ++mult_it)
    {
        for (int i = 0; i < dim; ++i)
        {
            std::vector<Size> vec(mult_it->get_data_direction(i));
            fill(vec.begin(), vec.end(),1);
            vec.front() += (*deg_it)[i];
            vec.back()  += (*deg_it)[i];
            mult_it->copy_data_direction(i, vec);
        }
    }
}








IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/multiplicity.inst>

