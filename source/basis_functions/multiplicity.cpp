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

//template <int dim>
//Multiplicity<dim>::
//Multiplicity(std::shared_ptr<const Grid> knots,
//             Degrees &deg,
//             bool max_reg)
//:
//parent_t::CartesianProductArray(knots->get_num_knots_dim())
//{
//    fill_max_regularity(deg);
//}

template <int dim>
Multiplicity<dim>
Multiplicity<dim>::
accumulate()
{
    const TensorSize<dim> size = this->tensor_size();
    Multiplicity<dim> result(size);

    for (int i = 0; i < dim; ++i)
    {
        result.data_[i][0] =  this->data_[i][0];

        const Size size_i = size(i);
        for (int k = 1 ; k < size_i ; ++k)
            result.data_[i][k] = result.data_[i][k-1] + this->data_[i][k];
    }

    return result;
}

template <int dim>
Multiplicity<dim>
Multiplicity<dim>::
compute_index_space_offset(const std::array<int,dim> &degree)
{
    const TensorSize<dim> size = this->tensor_size();
    Multiplicity<dim> result(size);

    for (int i = 0; i < dim; ++i)
    {
        result.data_[i][0] =  this->data_[i][0] - degree[i] - 1;

        const Size size_i = size(i);
        for (int k=1; k < size_i ; ++k)
            result.data_[i][k] = result.data_[i][k-1] + this->data_[i][k];
    }

    return result;
}



template <int dim>
void Multiplicity<dim>::
fill_max_regularity(const int degree)
{
    int_array<dim> deg;
    deg.fill(degree);
    fill_max_regularity(deg);
}



template <int dim>
void  Multiplicity<dim>::fill_max_regularity(const int_array<dim> degree)
{
    for (int i = 0; i < dim; ++i)
    {
        auto &vector = this->data_[i];
        fill(vector.begin(), vector.end(),1);
        vector.front() += degree[i];
        vector.back()  += degree[i];
    }
}








IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/multiplicity.inst>

