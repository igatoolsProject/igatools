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

#ifndef VALUE_TABLE_H_
#define __VALUE_TABLE_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/value_vector.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/utils/container_view.h>

#include <vector>

IGA_NAMESPACE_OPEN



/**
 * @class ValueTable
 * This class represents a 2-dimensional array for objects of type T.
 *
 * Each entry of the 2-dimensional array can be associated to a single function
 * at one single point.
 *
 * Internally the data are stored as single std::vector<T> of length num_functions * num_points
 * and the memory is ordered is made as num_functions chunks of length num_points objects of type T.
 * The array element with index i refers to the function index resulting by the integer division i/num_points
 * and to the point index resulting by the modulo operation i%num_points.
 *
 * In other words, the point-index runs faster than the function-index.
 *
 * @tparam T Type of the object to be stored in the table.
 * \todo Missing documentation.
 * @author M.Martinelli
 * @date 2013,2014
 */
template <class T>
class ValueTable :
    public DynamicMultiArray<T,2>
{
public :
    using iterator = typename DynamicMultiArray<T,2>::iterator ;
    using const_iterator = typename DynamicMultiArray<T,2>::const_iterator;
    using function_view = ContainerView<DynamicMultiArray<T,2>>;
    using const_function_view = ConstContainerView<DynamicMultiArray<T,2>>;


    /**
     * @name Constructors
     */
    ///@{
    /**
     * Default constructor.
     * Constructs an empty container, with no elements.
     */
    explicit ValueTable() ;

    /**
     * Constructor. Constructs a container for storing num_functions*num_points objects of type T.
     * @param[in] num_functions - Number of functions.
     * @param[in] num_points - Number of points.
     */
    explicit ValueTable(const Size num_functions, const Size num_points) ;

    /**
     * Copy constructor. Performs a deep copy of the ValueTable object @p table_in.
     * @warning Use this operator with caution because if the ValueTable to be copied is big,
     * this results in an CPU expensive operation.
     */
    ValueTable(const ValueTable<T> &table_in) = default ;

    /**
     * Move constructor.
     */
    ValueTable(ValueTable<T> &&table_in) = default ;


    /**
     * Destructor.
     */
    ~ValueTable() = default ;

    ///@}


    /**
     * @name Assignment operators
     */
    ///@{

    /**
     * Copy assignment operator. Performs a deep copy of the ValueTable object @p table_in.
     * @warning Use this operator with caution because if the ValueTable to be copied is big,
     * this results in an CPU expensive operation.
     */
    ValueTable<T> &operator=(const ValueTable<T> &table_in) = default ;


    /**
     * Move assignment operator.
     */
    ValueTable<T> &operator=(ValueTable<T> &&table_in) = default ;

    ///@}


    /**
     * @name Values initialization
     */
    ///@{

    /** Set all the values of the table to zero. */
    void zero() ;

    ///@}


    /**
     * @name Functions for resizing
     */
    ///@{

    /**
     * Resize the ValueTable in order to allocate space for @p num_functions functions and
     * @p num_points points.
     */
    void resize(const Size num_functions, const Size num_points) ;

    /**
     * Removes all elements from the ValueTable, leaving the container with a size of 0.
     */
    void clear() noexcept ;
    ///@}


    /**
     * @name Functions for getting size information
     */
    ///@{
    /**
     * Returns the number of elements in the ValueTable (= num_functions * num_points).
     */
    Size size() const;

    /**
     * Returns the number of points.
     */
    Size get_num_points() const noexcept ;

    /**
     * Returns the number of functions.
     */
    Size get_num_functions() const noexcept ;
    ///@}


    /**
     * @name Functions for getting values in the container
     */
    ///@{
    /**
     * Return a view of the elements relative to the i-th function. Non-const version.
     */
    function_view get_function_view(const int i);

    /**
     * Return a view of the elements relative to the i-th function. Const version.
     */
    const_function_view get_function_view(const int i) const;
    ///@}


    /**
     * @name Printing info
     */
    ///@{

    /**
     * Prints the content of the ValueTable on the LogStream @p out.
     */
    void print_info(LogStream &out) const ;
    ///@}


    /**
     * Returns the linear combination of the function values (at each evaluation points).
     * The size of vector of the @p coefficients must be equal to the number of functions
     * represented in the ValueTable.
     */
    ValueVector<T> evaluate_linear_combination(const std::vector<Real> &coefficients) const ;


private:

    /**
     * Number of functions for which the objects in the ValueTable refers to.
     */
    Size num_functions_ ;


    /**
     * Number of points for which the objects in the ValueTable refers to.
     */
    Size num_points_ ;
} ;


IGA_NAMESPACE_CLOSE


#endif /* VALUE_TABLE_H_ */
