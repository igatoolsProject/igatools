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

//TODO: inline, check constructor, fix headers
#include <igatools/utils/value_table.h>
#include <igatools/base/tensor.h>
#include <igatools/base/exceptions.h>

#include <iostream>


IGA_NAMESPACE_OPEN



template <class T>
ValueTable<T>::ValueTable()
    :
    ValueTable<T>(0,0)
{}

template <class T>
ValueTable<T>::ValueTable(
    const Size num_functions,
    const Size num_points)
    :
    ValueContainer<T>(num_functions,num_points)
{}








template <class T>
auto
ValueTable<T>::
get_function_view(const int i) -> view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;

    const auto n_pts = this->get_num_points();
    Assert(i >= 0 && i < this->get_num_functions(), ExcIndexRange(i,0,this->get_num_functions()));
    return view(
        iterator(*this, i    * n_pts, 1),
        iterator(*this,(i+1) * n_pts, 1));
}


template <class T>
auto
ValueTable<T>::
get_function_view(const int i) const -> const_view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;

    const auto n_pts = this->get_num_points();
    Assert(i >= 0 && i < this->get_num_functions(), ExcIndexRange(i,0,this->get_num_functions()));
    return const_view(
               const_iterator(*this, i    * n_pts, 1),
               const_iterator(*this,(i+1) * n_pts, 1));
}

template <class T>
auto
ValueTable<T>::
get_point_view(const int i) -> view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;

    const auto n_pts = this->get_num_points();
    Assert(i >= 0 && i < n_pts, ExcIndexRange(i,0,n_pts));
    return view(
        iterator(*this,i,n_pts),
        iterator(*this,IteratorState::pass_the_end,n_pts));
}


template <class T>
auto
ValueTable<T>::
get_point_view(const int i) const -> const_view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;

    const auto n_pts = this->get_num_points();
    Assert(i >= 0 && i < n_pts, ExcIndexRange(i,0,n_pts));
    return const_view(
               const_iterator(*this,i,n_pts),
               const_iterator(*this,IteratorState::pass_the_end,n_pts));
}


template <class T>
ValueVector<T>
ValueTable<T>::
evaluate_linear_combination(const SafeSTLVector<Real> &coefficients) const
{
    const auto n_funcs = this->get_num_functions();
    const auto n_pts = this->get_num_points();

    Assert(n_pts > 0, ExcLowerRange(n_pts,0));
    Assert(n_funcs > 0, ExcLowerRange(n_funcs,0));
    Assert(n_funcs == static_cast<int>(coefficients.size()),
           ExcDimensionMismatch(n_funcs,static_cast<int>(coefficients.size())));


    ValueVector<T> linear_combination(n_pts) ;

    for (int fn = 0 ; fn < n_funcs ; ++fn)
    {
        const auto func = this->get_function_view(fn) ;

        const Real coeff_fn = coefficients[fn] ;

        for (int pt = 0 ; pt < n_pts ; ++pt)
            linear_combination[pt] += coeff_fn * func[pt] ;
    }

    return linear_combination ;
}

template <class T>
void
ValueTable<T>::
resize(const Size num_functions, const Size num_points)
{
    ValueContainer<T>::resize(num_functions,num_points);
}

template <class T>
void
ValueTable<T>::
clear() noexcept
{
    ValueContainer<T>::clear();
}


template <class T>
void
ValueTable<T>::
print_info(LogStream &out) const
{
    // TODO (pauletti, Sep 12, 2014): should just called the parent print_info
    const auto n_funcs = this->get_num_functions();
    out.begin_item("ValueTable (num_functions=" + std::to_string(n_funcs) + ",num_points=" +
                   std::to_string(this->get_num_points()) + ") :");

    for (int iFunc = 0 ; iFunc < n_funcs ; ++iFunc)
    {
        out.begin_item("Function " + std::to_string(iFunc));

        auto value_func = this->get_function_view(iFunc) ;

        for (auto &value : value_func)
            out << value << " " ;
        out.end_item();
    }
    out.end_item();
}




IGA_NAMESPACE_CLOSE



#include <igatools/utils/value_table.inst>



