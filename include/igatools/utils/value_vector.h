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

#ifndef VALUE_VECTOR_H_
#define VALUE_VECTOR_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>


#include <vector>

IGA_NAMESPACE_OPEN

/**
 * @brief Container for objects of type T that refers to different evaluation points.
 *
 * This container class is derived from std::vector<T> with added a
 * function for printing its elements and a function for the reset of the element entries
 * using the proper default constructor.
 *
 * @tparam T Type of elements stored in the object.
 *
 * @author M.Martinelli
 * @date 2102,2013,2014
 */
template< class T >
class ValueVector :
    public std::vector<T>
{
public :
    /**
     * @name Constructors
     */
    ///@{
    explicit ValueVector() ;

    /**
     * Constructor. It builds a vector with num_points elements of type T, initialized using
     * the default constructor T().
     */
    explicit ValueVector(const Index num_points) ;

    // TODO (pauletti, Jul 11, 2014): it should be deleted after return by value in pf_accessor
    explicit ValueVector(const Size num_functions, const Size num_points)
        :ValueVector(num_points)
    {}
    /**
     * Constructor from a std::vector<T> object.
     * Performs a deep copy of the elements in @p vector_in.
     */
    explicit ValueVector(const std::vector<T> &vector_in) ;


    /**
     * Constructor from an initializer list.
     */
    ValueVector(const std::initializer_list<T> &list) ;

    /**
     * Copy constructor. Performs a deep copy of the content of the ValueVector object @p vector_in
     */
    ValueVector(const ValueVector<T> &vector_in) = default ;

    /**
     * Move constructor.
     */
    ValueVector(ValueVector<T> &&vector_in) = default ;

    /**
     * Destructor.
     */
    ~ValueVector() = default ;
    ///@}

    /**
     * @name Assignment operators
     */
    ///@{

    /**
     * Copy assignment operator. Performs a deep copy of the content of the ValueVector object.
     */
    ValueVector<T> &operator=(const ValueVector<T> &value_vector) = default ;

    /**
     * Copy assignment operator. Performs a deep copy of the content of the ValueVector object.
     */
    ValueVector<T> &operator=(const std::vector<T> &vector);

    /**
     * Move assignment operator.
     */
    ValueVector<T> &operator=(ValueVector<T> &&value_vector) = default ;

    ///@}

    /**
     * Returns the number of points.
     */
    Size get_num_points() const noexcept
    {
        return this->size();
    }
    /**
     * Returns the number of functions.
     */
    Size get_num_functions() const noexcept
    {
        return 1;
    }
    /**
     * @name Values initialization
     */
    ///@{

    /** Set all the values of the vector to zero. */
    void zero() ;

    ///@}


    /**
     * @name Printing info
     */
    ///@{
    /**
     * Prints the content of the ValueVector on the LogStream @p out.
     * Its use is intended mainly for testing and debugging purpose.
     */
    void print_info(LogStream &out) const ;
    ///@}

} ;


/**
 * Performs the scalar-by-vector multiplication <tt>scalar * a</tt>
 * @relates ValueVector
 */
template< class T>
ValueVector<T>
operator*(const ValueVector<T> &a, const Real scalar) ;

/**
 * Performs the scalar-by-vector multiplication <tt>scalar * a</tt>
 * @relates ValueVector
 */
template< class T>
ValueVector<T>
operator*(const Real scalar, const ValueVector<T> &a) ;

IGA_NAMESPACE_CLOSE


#endif /* VALUE_VECTOR_H_ */
