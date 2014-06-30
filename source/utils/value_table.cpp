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

//TODO: inline, check constructor, fix headers
#include <igatools/utils/value_table.h>
#include <igatools/base/tensor.h>
#include <igatools/base/exceptions.h>

#include <iostream>
using std::vector ;

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
    DynamicMultiArray<T,2>(TensorSize<2>( {num_points,num_functions})),
                  num_functions_ {num_functions},
num_points_ {num_points}
{
    Assert(num_functions >= 0, ExcLowerRange(num_functions,0));
    Assert(num_points >= 0, ExcLowerRange(num_points,0));
}



template <class T>
void
ValueTable<T>::
resize(const Size num_functions, const Size num_points)
{
    Assert(num_functions >= 0, ExcLowerRange(num_functions,0));
    Assert(num_points >= 0, ExcLowerRange(num_points,0));

    if (num_functions_ != num_functions ||
        num_points_ != num_points)
    {
        num_functions_ = num_functions ;
        num_points_ = num_points ;

        DynamicMultiArray<T,2>::resize(TensorSize<2>( {num_points_,num_functions_}));
    }
}

template <class T>
Size
ValueTable<T>::
size() const
{
    Assert(this->flat_size() == num_functions_ * num_points_,
           ExcDimensionMismatch(this->flat_size(), num_functions_ * num_points_)) ;

    return this->flat_size();
}

template <class T>
void
ValueTable<T>::
clear() noexcept
{
    num_functions_ = 0;
    num_points_ = 0;
    DynamicMultiArray<T,2>::resize(TensorSize<2>({num_points_,num_functions_}));
}


template <class T>
Size
ValueTable<T>::
get_num_functions() const noexcept
{
    Assert(num_functions_ == this->tensor_size()(1),
           ExcDimensionMismatch(num_functions_,this->tensor_size()(1)));
    return this->tensor_size()(1);
}


template <class T>
Size
ValueTable<T>::
get_num_points() const noexcept
{
    Assert(num_points_ == this->tensor_size()(0),
           ExcDimensionMismatch(num_points_,this->tensor_size()(0)));
    return this->tensor_size()(0);
}


template <class T>
auto
ValueTable<T>::
get_function_view(const int i) -> view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < num_functions_, ExcIndexRange(i,0,num_functions_));
    return view(
        iterator(*this, i    * num_points_, 1),
        iterator(*this,(i+1) * num_points_, 1));
}


template <class T>
auto
ValueTable<T>::
get_function_view(const int i) const -> const_view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < num_functions_, ExcIndexRange(i,0,num_functions_));
    return const_view(
               const_iterator(*this, i    * num_points_, 1),
               const_iterator(*this,(i+1) * num_points_, 1));
}

template <class T>
auto
ValueTable<T>::
get_point_view(const int i) -> view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < num_points_, ExcIndexRange(i,0,num_points_));
    return view(
        iterator(*this, i   ,num_points_),
        iterator(*this,IteratorState::pass_the_end,num_points_));
//    iterator(*this,(num_functions_-1) * num_points_ + i+1,num_points_));
}


template <class T>
auto
ValueTable<T>::
get_point_view(const int i) const -> const_view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < num_points_, ExcIndexRange(i,0,num_points_));
    return const_view(
               const_iterator(*this, i   ,num_points_),
               const_iterator(*this,IteratorState::pass_the_end,num_points_));
//               const_iterator(*this,(num_functions_-1) * num_points_ + i+1,num_points_));
}


template <class T>
ValueVector<T>
ValueTable<T>::
evaluate_linear_combination(const std::vector<Real> &coefficients) const
{
    Assert(num_points_ > 0, ExcLowerRange(num_points_,0));
    Assert(num_functions_ > 0, ExcLowerRange(num_functions_,0));
    Assert(num_functions_ == static_cast<int>(coefficients.size()),
           ExcDimensionMismatch(num_functions_,static_cast<int>(coefficients.size())));

    ValueVector<T> linear_combination(num_points_) ;

    for (int iFn = 0 ; iFn < num_functions_ ; ++iFn)
    {
        const auto func = this->get_function_view(iFn) ;

        Real coeff_iFn = coefficients[iFn] ;

        for (int jPt = 0 ; jPt < num_points_ ; ++jPt)
            linear_combination[jPt] += coeff_iFn * func[jPt] ;
    }

    return linear_combination ;
}



template <class T>
void
ValueTable<T>::
print_info(LogStream &out) const
{
    out << "ValueTable (num_functions=" << num_functions_ << ",num_points=" << num_points_ << ") :" << std::endl ;

    for (int iFunc = 0 ; iFunc < num_functions_ ; iFunc++)
    {
        out.push("\t");

        out << "Function[" << iFunc << "] : " ;

        auto value_func = this->get_function_view(iFunc) ;

        for (auto &value : value_func)
            out << value << " " ;

        out << std::endl ;
        out.pop();
    }
}


template <class T>
void
ValueTable<T>::
zero()
{
    for (auto & value : (*this))
        value = T() ;
}


IGA_NAMESPACE_CLOSE



#include <igatools/utils/value_table.inst>



