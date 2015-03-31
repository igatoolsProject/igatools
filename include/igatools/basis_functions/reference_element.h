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


#ifndef REFERENCE_ELEMENT_H_
#define REFERENCE_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/reference_element_handler.h>

IGA_NAMESPACE_OPEN


template <int, int, int> class ReferenceSpace;

/**
 *
 * @ingroup elements
 */
template <int dim, int range, int rank>
class ReferenceElement : public SpaceElement<dim,0,range,rank>
{
public:
    /** Type for the grid accessor. */
    using GridAccessor = CartesianGridElement<dim>;

    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const ReferenceSpace<dim,range,rank> ;

    using Space = ReferenceSpace<dim,range,rank>;
    using ConstSpace = const ReferenceSpace<dim,range,rank>;

    using parent_t = SpaceElement<dim,0,range,rank>;

    using RefPoint = typename Space::RefPoint;
    using Point = typename Space::Point;
    using Value = typename Space::Value;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;

    using Div = typename Space::Div;

    ReferenceElement() = delete;

    ReferenceElement(const ReferenceElement<dim,range,rank> &elem,
                     const iga::CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Constructs an accessor to element number index of a
     * ReferenceSpace space.
     */
    ReferenceElement(const std::shared_ptr<ConstSpace> space,
                     const Index elem_index);

    /**
     * Constructs an accessor to element number index of a
     * Reference space.
     */
    ReferenceElement(const std::shared_ptr<ConstSpace> space,
                     const TensorIndex<dim> &elem_index);


    virtual ~ReferenceElement() = default;


    /** @name Functions/operators for moving the element in the ReferenceSpace.*/
    ///@{
    /**
     * Sets the index of the element using the flatten representation.
     * @note This function also updates the index for the tensor representation.
     * @warning This may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    virtual void move_to(const Index flat_index) ;
    ///@}


    /**
     * Creates a new object performing a deep copy of the current object using
     * the
     * copy constructor of the derived class.
     *
     * @note This function should be not called directly, but it should be
     * called its
     * specialization on a derived class. It would be better to define this
     * function
     * <em>pure virtual</em> but this will not allow to dereference an iterator
     * containing
     * a pointer to an object of kind ReferenceElement.
     */
    virtual std::shared_ptr<ReferenceElement<dim,range,rank> > clone() const
    {
        Assert(false,ExcMessage("This function must not be called. "
                                "You should call the clone() funtion of a derived base class."));
        return nullptr;
    }


    /**
     * This boost::mpl::map is used to select the return type of the template function evaluate_basis_at_points()
     * from its template parameter @p ValueType.
     *
     * @see evaluate_basis_at_points()
     */
    using map_ValueTypeId_ContainerType =
        boost::mpl::map<
        boost::mpl::pair<boost::mpl::int_<     _Value::id>,ValueTable<Value> >,
        boost::mpl::pair<boost::mpl::int_<  _Gradient::id>,ValueTable<Derivative<1>> >,
        boost::mpl::pair<boost::mpl::int_<   _Hessian::id>,ValueTable<Derivative<2>> >,
        boost::mpl::pair<boost::mpl::int_<_Divergence::id>,ValueTable<Div>>
        >;

    /**
     * This is a type-trait that convert a @p ValueType to the correspondent container returned
     * by the function evaluate_basis_at_points().
     *
     * @see evaluate_basis_at_points()
     */
    template <class ValueType>
    using ContType_from_ValueType = typename boost::mpl::at<
                                    map_ValueTypeId_ContainerType,boost::mpl::int_<ValueType::id>>::type;



    /**
     * @name Functions for the basis and field evaluations without the use of the cache.
     */
    ///@{
    /**
     * Returns a ValueTable with the quantity specified by the template parameter @p ValueType,
     * computed for all local basis function,
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_cache()/fill_cache().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    template <class ValueType>
    ContType_from_ValueType<ValueType>
    evaluate_basis_at_points(
        const Quadrature<dim> &points,
        const std::string &dofs_property)
    {
        auto elem_handler = ReferenceElementHandler<dim,range,rank>::create(this->space_);

        ValueFlags flags;
        if (ValueType::id == _Value::id)
            flags = ValueFlags::value;
        else if (ValueType::id == _Gradient::id)
            flags = ValueFlags::gradient;
        else if (ValueType::id == _Hessian::id)
            flags = ValueFlags::hessian;
        else if (ValueType::id == _Divergence::id)
            flags = ValueFlags::divergence;
        else
        {
            Assert(false,ExcNotImplemented());
        }

        elem_handler->reset_one_element(flags,points,this->get_flat_index());
        elem_handler->template init_cache<dim>(*this);
        elem_handler->template fill_cache<dim>(*this,0);

        //    Assert(false,ExcNotImplemented());

        return this->template get_basis<ValueType,dim>(0,dofs_property);

    }

    ///@}


    /**
     * Returns the <tt>k</tt> dimensional j-th sub-element measure
     * multiplied by the weights of the quadrature.
     */
    template <int k>
    ValueVector<Real> get_w_measures(const int j) const
    {
        return this->as_cartesian_grid_element_accessor().template get_w_measures<k>(j);
    }

    /**
     * Returns the gradient determinant of the identity map at the dilated quadrature points.
     */
    ValueVector<Real> get_element_w_measures() const
    {
        return this->template get_w_measures<dim>(0);
    }


    using OffsetTable = typename Space::template ComponentContainer<int>;
    using TensorSizeTable = typename Space::TensorSizeTable;

protected:

    /** Number of scalar basis functions along each direction, for all space components. */
    TensorSizeTable n_basis_direction_;

    /** Basis function ID offset between the different components. */
    OffsetTable comp_offset_;

    using Indexer = CartesianProductIndexer<dim>;
    using IndexerPtr = std::shared_ptr<Indexer>;
    using IndexerPtrTable = typename Space::template ComponentContainer<IndexerPtr>;

    /** Hash table for fast conversion between flat-to-tensor basis function ids. */
    IndexerPtrTable basis_functions_indexer_;

public:
    using parent_t::get_num_basis;

    /**
     * Returns the max. number of basis function that can have support on this element.
     */
    int get_num_basis() const override final;

    /**
     * Returns the basis function ID offset between the different components.
     */
    OffsetTable get_basis_offset() const;

    /**
     * Number of non-zero scalar basis functions associated
     * with the i-th space component on the element.
     * This makes sense as a reference B-spline space
     * is only allowed to be of the cartesian product type
     * V = V1 x V2 x ... X Vn.
     */
    int get_num_basis_comp(const int i) const;

    void print_info(LogStream &out) const;

protected:
    std::shared_ptr<const Space> space_;


public:
    std::shared_ptr<const Space> get_space() const
    {
        return space_;
    }
};



IGA_NAMESPACE_CLOSE


#endif // #ifndef REFERENCE_ELEMENT_H_

