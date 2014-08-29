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

#ifndef VALUE_CONTAINER_H_
#define VALUE_CONTAINER_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/utils/value_vector.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/utils/container_view.h>

#include <vector>

IGA_NAMESPACE_OPEN

template<class T>
class ValueContainer :
    public DynamicMultiArray<T,2>
{
public:
    /** Type for the iterator (non-const version). */
    using iterator = typename DynamicMultiArray<T,2>::iterator ;

    /** Type for the iterator (const version). */
    using const_iterator = typename DynamicMultiArray<T,2>::const_iterator;

    /** Type for the view (non-const version). */
    using view = ContainerView<DynamicMultiArray<T,2>>;


    /** Type for the view (const version). */
    using const_view = ConstContainerView<DynamicMultiArray<T,2>>;

    /** @name Constructors */
    ///@{
    /**
     * Constructor. Constructs a container for storing num_functions*num_points objects of type T.
     * @param[in] num_functions - Number of functions.
     * @param[in] num_points - Number of points.
     */
    explicit ValueContainer(const Size num_functions, const Size num_points)
        :
        DynamicMultiArray<T,2>(TensorSize<2>({num_points,num_functions})),
                      num_functions_ {num_functions},
    num_points_ {num_points}
    {
        Assert(num_functions >= 0, ExcLowerRange(num_functions,0));
        Assert(num_points >= 0, ExcLowerRange(num_points,0));
    }

    /**
     * Default constructor.
     * Initializes a container with no entries (num_points == num_functions == 0).
     */
    ValueContainer()
        :
        ValueContainer<T>(0,0)
    {}

    /** Copy constructor. */
    ValueContainer(const ValueContainer<T> &in) = default;

    /** Move constructor. */
    ValueContainer(ValueContainer<T> &&in) = default;

    ~ValueContainer() = default;
    ///@}


    /** @name Assignment operators*/
    ///@{
    /** Copy assignment operator.*/
    ValueContainer<T> &operator=(const ValueContainer<T> &in) = default;

    /** Move assignment operator.*/
    ValueContainer<T> &operator=(ValueContainer<T> &&in) = default;
    ///@}


    /**
     * @name Functions for getting size information
     */
    ///@{
    /**
     * Returns the number of elements in the ValueTable (= num_functions * num_points).
     */
    Size size() const
    {
        Assert(this->get_data().size() == this->flat_size(),
               ExcDimensionMismatch(this->get_data().size(),this->flat_size()));
        Assert(this->flat_size() == num_functions_ * num_points_,
               ExcDimensionMismatch(this->flat_size(), num_functions_ * num_points_)) ;

        return this->flat_size();
    }

    /**
     * Returns the number of points.
     */
    Size get_num_points() const noexcept
    {
        return num_points_;
    }

    /**
     * Returns the number of functions.
     */
    Size get_num_functions() const noexcept
    {
        return num_functions_;
    }
    ///@}

    /**
     * @name Values initialization
     */
    ///@{

    /**
     * Set all the values of the Container to zero.
     * @note The "zero" values means the default constructor T().
     */
    void zero()
    {
        for (auto &value : (*this))
            value = T() ;
    }

    ///@}

protected:

    /**
     * @name Functions for resizing
     */
    ///@{
    /**
     * Resize the ValueTable in order to allocate space for @p num_functions functions and
     * @p num_points points.
     */
    void resize(const Size num_functions, const Size num_points)
    {
        Assert(num_functions >= 0, ExcLowerRange(num_functions,0));
        Assert(num_points >= 0, ExcLowerRange(num_points,0));

        if (num_functions_ != num_functions ||
            num_points_ != num_points)
        {
            num_functions_ = num_functions;
            num_points_ = num_points;

            DynamicMultiArray<T,2>::resize(TensorSize<2>({num_points_,num_functions_}));
        }
    }

    /**
     * Removes all elements from the ValueTable, leaving the container with a size of 0.
     */
    void clear() noexcept
    {
        num_functions_ = 0;
        num_points_ = 0;
        DynamicMultiArray<T,2>::clear();
    }
    ///@}
    /**
     * Number of functions for which the objects in the ValueTable refers to.
     */
    Size num_functions_ ;


    /**
     * Number of points for which the objects in the ValueTable refers to.
     */
    Size num_points_ ;
};



IGA_NAMESPACE_CLOSE


#endif /* VALUE_CONTAINER_H_ */
