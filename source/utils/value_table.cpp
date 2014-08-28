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
    Assert(i >= 0 && i < this->num_functions_, ExcIndexRange(i,0,this->num_functions_));
    return view(
        iterator(*this, i    * this->num_points_, 1),
        iterator(*this,(i+1) * this->num_points_, 1));
}


template <class T>
auto
ValueTable<T>::
get_function_view(const int i) const -> const_view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < this->num_functions_, ExcIndexRange(i,0,this->num_functions_));
    return const_view(
               const_iterator(*this, i    * this->num_points_, 1),
               const_iterator(*this,(i+1) * this->num_points_, 1));
}

template <class T>
auto
ValueTable<T>::
get_point_view(const int i) -> view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < this->num_points_, ExcIndexRange(i,0,this->num_points_));
    return view(
        iterator(*this,i,this->num_points_),
        iterator(*this,IteratorState::pass_the_end,this->num_points_));
}


template <class T>
auto
ValueTable<T>::
get_point_view(const int i) const -> const_view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < this->num_points_, ExcIndexRange(i,0,this->num_points_));
    return const_view(
               const_iterator(*this,i,this->num_points_),
               const_iterator(*this,IteratorState::pass_the_end,this->num_points_));
}


template <class T>
ValueVector<T>
ValueTable<T>::
evaluate_linear_combination(const vector<Real> &coefficients) const
{
    Assert(this->num_points_ > 0, ExcLowerRange(this->num_points_,0));
    Assert(this->num_functions_ > 0, ExcLowerRange(this->num_functions_,0));
    Assert(this->num_functions_ == static_cast<int>(coefficients.size()),
           ExcDimensionMismatch(this->num_functions_,static_cast<int>(coefficients.size())));

    ValueVector<T> linear_combination(this->num_points_) ;

    for (int iFn = 0 ; iFn < this->num_functions_ ; ++iFn)
    {
        const auto func = this->get_function_view(iFn) ;

        Real coeff_iFn = coefficients[iFn] ;

        for (int jPt = 0 ; jPt < this->num_points_ ; ++jPt)
            linear_combination[jPt] += coeff_iFn * func[jPt] ;
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
    out << "ValueTable (num_functions=" << this->num_functions_ << ",num_points=" << this->num_points_ << ") :" << std::endl ;

    for (int iFunc = 0 ; iFunc < this->num_functions_ ; iFunc++)
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




IGA_NAMESPACE_CLOSE



#include <igatools/utils/value_table.inst>



