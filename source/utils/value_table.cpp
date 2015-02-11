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
    Assert(i >= 0 && i < this->get_num_functions(), ExcIndexRange(i,0,this->get_num_functions()));
    return view(
        iterator(*this, i    * this->get_num_points(), 1),
        iterator(*this,(i+1) * this->get_num_points(), 1));
}


template <class T>
auto
ValueTable<T>::
get_function_view(const int i) const -> const_view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < this->get_num_functions(), ExcIndexRange(i,0,this->get_num_functions()));
    return const_view(
               const_iterator(*this, i    * this->get_num_points(), 1),
               const_iterator(*this,(i+1) * this->get_num_points(), 1));
}

template <class T>
auto
ValueTable<T>::
get_point_view(const int i) -> view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < this->get_num_points(), ExcIndexRange(i,0,this->get_num_points()));
    return view(
        iterator(*this,i,this->get_num_points()),
        iterator(*this,IteratorState::pass_the_end,this->get_num_points()));
}


template <class T>
auto
ValueTable<T>::
get_point_view(const int i) const -> const_view
{
    Assert(this->size() > 0, ExcEmptyObject()) ;
    Assert(i >= 0 && i < this->get_num_points(), ExcIndexRange(i,0,this->get_num_points()));
    return const_view(
               const_iterator(*this,i,this->get_num_points()),
               const_iterator(*this,IteratorState::pass_the_end,this->get_num_points()));
}


template <class T>
ValueVector<T>
ValueTable<T>::
evaluate_linear_combination(const vector<Real> &coefficients) const
{
    Assert(this->get_num_points() > 0, ExcLowerRange(this->get_num_points(),0));
    Assert(this->get_num_functions() > 0, ExcLowerRange(this->get_num_functions(),0));
    Assert(this->get_num_functions() == static_cast<int>(coefficients.size()),
           ExcDimensionMismatch(this->get_num_functions(),static_cast<int>(coefficients.size())));

    ValueVector<T> linear_combination(this->get_num_points()) ;

    for (int iFn = 0 ; iFn < this->get_num_functions() ; ++iFn)
    {
        const auto func = this->get_function_view(iFn) ;

        Real coeff_iFn = coefficients[iFn] ;

        for (int jPt = 0 ; jPt < this->get_num_points() ; ++jPt)
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
    out << "ValueTable (num_functions=" << this->get_num_functions() << ",num_points=" << this->get_num_points() << ") :" << std::endl ;

    for (int iFunc = 0 ; iFunc < this->get_num_functions() ; iFunc++)
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



