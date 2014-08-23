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

//TODO: inline, check constructor
#include <igatools/utils/value_vector.h>
#include <igatools/base/tensor.h>
#include <igatools/base/exceptions.h>
//#include <boost/multi_array.hpp>

using std::vector ;

IGA_NAMESPACE_OPEN


template <class T>
ValueVector<T>::
ValueVector() : ValueVector<T>(0)
{}

template <class T>
ValueVector<T>::
ValueVector(const Index num_points)
    : std::vector<T>(num_points,T {})
{}

template <class T>
ValueVector<T>::
ValueVector(const vector<T> &vector_in)
    : std::vector<T> {vector_in}
{}

template <class T>
ValueVector<T>::
ValueVector(const std::initializer_list<T> &list)
    : std::vector<T>(list)
{}

template <class T>
ValueVector<T> &
ValueVector<T>::
operator=(const std::vector<T> &vector)
{
    std::vector<T>::operator=(vector);
    return (*this);
}



template <class T>
void
ValueVector<T>::
print_info(LogStream &out) const
{
    const int num_points = this->size() ;
    out << "ValueVector (num_points=" << num_points << ") :" << std::endl ;

    for (int iPt = 0 ; iPt < num_points ; iPt++)
        out << (*this)[ iPt ] << " " ;

    out << std::endl ;
}

template <class T>
void
ValueVector<T>::
zero()
{
    for (auto &value : (*this))
        value = T() ;
}


template< class T>
ValueVector<T>
operator*(const ValueVector<T> &a, const Real scalar)
{
    ValueVector<T> result = a ;

    for (auto &r : result)
        r *= scalar ;

    return result ;
}

/**
 * Performs the scalar-by-vector multiplication scalar * a
 */
template< class T>
ValueVector<T>
operator*(const Real scalar, const ValueVector<T> &a)
{
    ValueVector<T> result = a ;

    for (auto &r : result)
        r *= scalar ;

    return result ;
}

template< class T>
T&
ValueVector<T>::
operator[](const Index i)
{
	Assert(i >= 0 && i < this->get_num_points(),ExcIndexRange(i,0,this->get_num_points()));
	return vector<T>::operator[](i);
}

template< class T>
const T&
ValueVector<T>::
operator[](const Index i) const
{
	Assert(i >= 0 && i < this->get_num_points(),ExcIndexRange(i,0,this->get_num_points()));
	return vector<T>::operator[](i);
}

IGA_NAMESPACE_CLOSE



#include <igatools/utils/value_vector.inst>



