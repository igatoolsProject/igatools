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

/*
 *  Test for the
 *
 *  author: pauletti
 *  date: Sep 11, 2014
 *
 */

#include "../tests.h"

#include <igatools/utils/safe_stl_vector.h>
#include <igatools/utils/safe_stl_array.h>
#include <igatools/utils/tensor_index.h>
#include <igatools/base/array_utils.h>

#include <algorithm>




template<int size>
SafeSTLVector<TensorIndex<size>>
                              partition(const int n)
{
    SafeSTLVector<TensorIndex<size>> v;
    TensorIndex<size> arr(0);

    arr[0] = n;
    v.push_back(arr);

    for (int j=1; j<n+1; ++j)
    {
        auto w = partition<size-1>(j);
        for (auto a : w)
        {
            arr[0] = n-j;
            std::copy(a.begin(), a.end(), arr.begin()+1);
            v.push_back(arr);
        }
    }
    return v;
}

template<>
SafeSTLVector<TensorIndex<1>>
                           partition<1>(const int n)
{
    TensorIndex<1> arr(n);
    return SafeSTLVector<TensorIndex<1>>(1,arr);
}

//template<>
//SafeSTLVector<TensorIndex<0>>
//partition<0>(const int n)
//{
//    return SafeSTLVector<TensorIndex<0>>();
//}



template<int dim, int order>
class TensorFunctionDerivativesSymmetry
{
public:
//    static const int num_entries_total = pow(dim,order);
    static const int num_entries_eval = constexpr_binomial_coefficient(dim-1+order,order);

    TensorFunctionDerivativesSymmetry()
    {
        auto uni_indices = partition<dim>(order);
        std::copy(uni_indices.begin(), uni_indices.end(), univariate_order.begin());



        for (int j=0; j<num_entries_eval; ++j)
        {
            auto &der_ind = eval_indices[j];
            int s=0;
            for (int dir=0; dir<dim; ++dir)
            {
                for (int l=0; l<uni_indices[j][dir]; ++l)
                    der_ind[s+l] = dir;
                s += uni_indices[j][dir];
            }

            auto ind = sequence<order>();
            SafeSTLVector<TensorIndex<order>> v;
            do
            {
                TensorIndex<order> ti;
                for (int i=0; i<order; ++i)
                    ti[i] = eval_indices[j][ind[i]];
                v.push_back(ti);
            }
            while (std::next_permutation(ind.begin(),ind.end()));

            auto it = std::unique(v.begin(), v.end());
            v.resize(std::distance(v.begin(),it));

            copy_indices[j] = v;
        }
    }

    void print_info(LogStream &out) const
    {
        out.begin_item("univariate derivative orders:");
        univariate_order.print_info(out);
        out.end_item();

        out.begin_item("Assigment indices:");
        eval_indices.print_info(out);
        out.end_item();

        out.begin_item("all equal indices indices:");
        copy_indices.print_info(out);
        out.end_item();
    }
    SafeSTLArray<TensorIndex<dim>, num_entries_eval> univariate_order;

    SafeSTLArray<TensorIndex<order>, num_entries_eval> eval_indices;

    SafeSTLArray<SafeSTLVector<TensorIndex<order>>, num_entries_eval> copy_indices;

};




template <int order, int dim>
void indices()
{
    auto v = partition<dim>(order);
    v.print_info(out);
    out << endl;
}



int main()
{
    // indices<3,2>();
    //indices<2,3>();

    //indices<0,2>();

    TensorFunctionDerivativesSymmetry<2,0> a;
    a.print_info(out);


    return  0;
}
