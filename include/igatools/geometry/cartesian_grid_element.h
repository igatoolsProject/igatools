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

#ifndef __CARTESIAN_GRID_ELEMENT_H_
#define __CARTESIAN_GRID_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/geometry/base_element.h>
#include <igatools/basis_functions/values_cache.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/geometry/cartesian_grid_iterator.h>
#include <igatools/utils/value_vector.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Element accessor for the CartesianGrid.
 *
 * The element can be queried for informations
 * that can be generated on-the-fly
 * (i.e. without the use of a cache) and for informations
 * that are obtained through a cache mechanism
 *
 *
 * See module (and the submodules) on \ref elements for a general overview.
 * @ingroup elements
 *
 * ### Quantities handled by the cache
 * - _Point i.e. evaluation points mapped in the parametric domain
 * - _W_Measure i.e. quadrature weights associated to each quadrature point,
 * multiplied by the element <tt>dim</tt>-dimensional measure.
 *
 * @author pauletti, 2012, 2013, 2014, 2015
 * @author martinelli, 2013, 2014, 2105
 *
 * @ingroup serializable
 */
template <int dim>
class CartesianGridElement
{
private:
    using self_t = CartesianGridElement<dim>;

public:
    /** Type required by the CartesianGridIterator templated iterator */
    using ContainerType = const CartesianGrid<dim>;


    using Point = Points<dim>;

    /** @name Constructors */
    ///@{
protected:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    CartesianGridElement() = default;

public:
    /**
     * Construct an accessor pointing to the element with
     * flat index @p elem_index of the CartesianGrid @p grid.
     */
    CartesianGridElement(const std::shared_ptr<ContainerType> grid,
                         const Index elem_index);


    /**
     * Copy constructor.
     * It can be used with different copy policies
     * (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a
     * classic copy constructor)
     * uses the <b>deep</b> copy.
     */
    CartesianGridElement(const self_t &elem,
                         const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    CartesianGridElement(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~CartesianGridElement() = default;
    ///@}

    /**
     * @name Functions for performing different kind of copy.
     */
    ///@{
    /**
     * Performs a deep copy of the input @p element,
     * i.e. a new local cache is built using the copy constructor on the local cache of @p element.
     *
     * @note In DEBUG mode, an assertion will be raised if the input local cache is not allocated.
     */
    void deep_copy_from(const self_t &element);

    /**
     * Performs a shallow copy of the input @p element. The current object will contain a pointer to the
     * local cache used by the input @p element.
     */
    void shallow_copy_from(const self_t &element);
    ///@}

    /** @name Assignment operators */
    ///@{
    /**
     * Copy assignment operator. Performs a <b>shallow copy</b> of the input @p element.
     *
     * @note Internally it uses the function shallow_copy_from().
     */
    self_t &operator=(const self_t &element);

    /**
     * Move assignment operator.
     */
    self_t &operator=(self_t &&elem) = default;
    ///@}


    /** Return the cartesian grid from which the element belongs.*/
    const std::shared_ptr<const CartesianGrid<dim>> get_grid() const;


    /** @name Functions related to the indices of the element in the cartesian grid. */
    ///@{
    /** Returns the index of the element in its flatten representation. */
    Index get_flat_index() const;

    /** Returns the index of the element in its tensor representation. */
    TensorIndex<dim> get_tensor_index() const;
    ///@}


    /**
     * @name Functions for managing/querying the element properties.
     */
    ///@{
    /**
     * Tests if a certain element @p property is TRUE.
     */
    bool is_property_true(const std::string &property) const;
    ///@}

    /**
     * @name Functions/operators for moving the element in the CartesianGrid.
     *
     * @note They should be called only by the CartesianGridIterator.
     */
    ///@{
    /**
     * Sets the index of the element using the flatten representation.
     * @note This function also updates the index for the tensor representation.
     * @warning This may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    void move_to(const Index flat_index);

    ///@}

    /**
     * @name Comparison operators
     * @note In order to be meaningful, the comparison must be performed on elements defined on
     * the <b>same grid</b>
     * (in the sense that the pointer to the grid held by the element must point to the same
     * grid object).
     */
    ///@{
    /**
     * True if the elements have the same index.
     *  @note In debug mode, it is also check they both refer to
     *  the same cartesian grid. No check is done on the cache.
     */
    bool operator==(const CartesianGridElement<dim> &elem) const;

    /**
     * True if the elements have different index.
     *  @note In debug mode, it is also check they both refer to
     *  the same cartesian grid. No check is done on the cache.
     */
    bool operator!=(const CartesianGridElement<dim> &elem) const;

    /**
     * True if the flat-index of the element on the left is smaller than
     * the flat-index of the element on the right.
     *  @note In debug mode, it is also check they both refer to
     *  the same cartesian grid. No check is done on the cache.
     */
    bool operator<(const CartesianGridElement<dim> &elem) const;

    /**
     * True if the flat-index of the element on the left is bigger than
     * the flat-index of the element on the right.
     *  @note In debug mode, it is also check they both refer to
     *  the same cartesian grid. No check is done on the cache.
     */
    bool operator>(const CartesianGridElement<dim> &elem) const;
    ///@}

    ///@name Query information that requires the use of the cache
    ///@{

    /**
     * Returns the lengths of the coordinate sides of the cartesian element.
     * For example in 2 dimensions
     * \code{.cpp}
       auto length = elem.coordinate_lenths();
       // length[0] is the length of the x-side of the element and
       // length[1] the length of the y-side of the element.
       \endcode
     */

    /**
       * Test if the element has a boundary face.
       */
    template<int k = (dim > 0) ? (dim-1) : 0 >
    bool is_boundary() const;

    /**
     * Test if the face identified by @p face_id on the current element is on the
     * boundary of the cartesian grid.
     */
    template<int k = (dim > 0) ? (dim-1) : 0>
    bool is_boundary(const Index sub_elem_id) const;
    ///@}

    /**
     * Return the @p i-th vertex
     */
    Point vertex(const int i) const;


    /**
     * Return the properties defined for the element.
     */
    SafeSTLVector<std::string> get_defined_properties() const;


    //TODO (martinelli, Apr 16, 2015): the returned value should be Points<topology_dim>, i.e.
    // the lengths should refer to the sub-element of dimension topology_dim
    //TODO (martinelli, Aug 13, 2015): maybe it is wothy to declare this function private
    template<int topology_dim>
    const Point get_side_lengths(const int topology_id) const;

    template <int topology_dim>
    Real get_measure(const int topology_id) const;

    /**
     * Returns the <tt>topology_dim</tt> dimensional topology_id-th sub-element measure
     * multiplied by the weights of the quadrature.
     */
    template <int topology_dim>
    ValueVector<Real> get_w_measures(const int topology_id) const;


    template <int topology_dim = dim>
    ValueVector<Point> get_points(const int topology_id = 0) const;

    ValueVector<Point> get_element_points() const;

    /**
     * Prints internal information about the CartesianGridElementAccessor.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out) const;

    void print_cache_info(LogStream &out) const;




private:
    template <class Accessor> friend class CartesianGridIteratorBase;
    friend class GridElementHandler<dim>;

    /** Cartesian grid from which the element belongs.*/
    std::shared_ptr<ContainerType> grid_;

    /** Tensor product indices of the current struct index @p flat_index_. */
    TensorIndex<dim> tensor_index_;

    /**
     * @name Types, data and methods for the cache.
     */
    ///@{

    /**
     * Alias used to define the container for the values in the cache.
     */
    using CType = boost::fusion::map<
                  boost::fusion::pair<    _Point,DataWithFlagStatus<ValueVector<Points<dim>>>>,
                  boost::fusion::pair<_W_Measure,DataWithFlagStatus<ValueVector<Real>>>
                  >;
    /**
     * Returns the flags that are valid to be used with this class.
     *
     * @note The valid flags are defined to be the ones that can be inferred from the ValueType(s)
     * used as key of the boost::fusion::map in CType.
     */
    static ValueFlags get_valid_flags();

public:
    using ValuesCache = FuncValuesCache<dim,CType>;

    using CacheType = AllSubElementsCache<ValuesCache>;

private:
    /** The local cache. */
    std::shared_ptr<CacheType> all_sub_elems_cache_;


    template <class ValueType, int topology_dim>
    const auto &
    get_values_from_cache(const int topology_id) const
    {
        Assert(all_sub_elems_cache_ != nullptr, ExcNullPtr());
        const auto &cache = all_sub_elems_cache_->template get_sub_elem_cache<topology_dim>(topology_id);
        return cache.template get_data<ValueType>();
    }
    ///@}

protected:
    /**
     * Performs a copy of the input @p element.
     * The type of copy (deep or shallow) is specified by the input parameter @p copy_policy.
     */
    void copy_from(const self_t &element,
                   const CopyPolicy &copy_policy);

    /**
     * ExceptionUnsupported Value Flag.
     */
    DeclException2(ExcFillFlagNotSupported, ValueFlags, ValueFlags,
                   << "The passed ValueFlag " << arg2
                   << " contains a non admissible flag " << (arg1 ^arg2));

    DeclException1(ExcCacheInUse, int,
                   << "The global cache is being used by " << arg1
                   << " iterator. Changing its value not allowed.");


private:
    /**
     * Creates a new object performing a deep copy of the current object using the CartesianGridElement
     * copy constructor.
     *
     * @note For information on the local cache treatment, see the documentation
     * of deep_copy_from().
     */
    std::shared_ptr<CartesianGridElement<dim> > clone() const;


#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif // SERIALIZATION
};

IGA_NAMESPACE_CLOSE

#endif
