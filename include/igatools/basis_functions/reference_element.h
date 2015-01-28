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

IGA_NAMESPACE_OPEN


template <int, int, int> class ReferenceSpace;

template <int dim, int range, int rank>
class ReferenceElement : public SpaceElement<ReferenceSpace<dim,range,rank>>
{
public:
    /** Type for the grid accessor. */
    using GridAccessor = CartesianGridElement<dim>;

    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const ReferenceSpace<dim,range,rank> ;

    using Space = ReferenceSpace<dim,range,rank>;
    using ConstSpace = const ReferenceSpace<dim,range,rank>;

    using parent_t = SpaceElement<ReferenceSpace<dim,range,rank>>;

    using RefPoint = typename Space::RefPoint;
    using Point = typename Space::Point;
    using Value = typename Space::Value;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;

    ReferenceElement()
    {
        Assert(false,ExcNotImplemented());
    }

    ReferenceElement(const ReferenceElement<dim,range,rank> &elem,
                     const iga::CopyPolicy &copy_policy = CopyPolicy::deep)
        :
        parent_t(elem,copy_policy)
    {};

    /**
     * Constructs an accessor to element number index of a
     * ReferenceSpace space.
     */
    ReferenceElement(const std::shared_ptr<ConstSpace> space,
                     const Index elem_index)
        :
        parent_t(space,elem_index)
    {
        Assert(this->get_space() != nullptr,ExcNullPtr());
    };

    /**
     * Constructs an accessor to element number index of a
     * Reference space.
     */
    ReferenceElement(const std::shared_ptr<ConstSpace> space,
                     const TensorIndex<dim> &elem_index)
        :
        parent_t(space,elem_index)
    {
        Assert(this->get_space() != nullptr,ExcNullPtr());
    };

    virtual ~ReferenceElement() = default;


    /** @name Functions/operators for moving the element in the ReferenceSpace.*/
    ///@{
    /**
     * Moves the element to the position that differs from the current one
     * for the quantity given by @p increment.
     *
     * If the resulting position after the movement is valid (i.e. within the grid), then the function
     * returns true, otherwise it returns false.
     */
    virtual bool jump(const TensorIndex<dim> &increment) ;

    /**
     * Sets the index of the element using the flatten representation.
     * @note This function also updates the index for the tensor representation.
     * @warning This may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    virtual void move_to(const Index flat_index) ;


    /**
     * Sets the index of the element using the tensor representation.
     * @note This function also updates the index for the flatten representation.
     * @warning this may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    virtual void move_to(const TensorIndex<dim> &tensor_index) ;

    /** Moves the element to the next active element in the CartesianGrid. */
    virtual void operator++() ;
    ///@}


    /**
     * Creates a new object performing a deep copy of the current object using the
     * copy constructor of the derived class.
     *
     * @note This function should be not called directly, but it should be called its
     * specialization on a derived class. It would be better to define this function
     * <em>pure virtual</em> but this will not allow to dereference an iterator containing
     * a pointer to an object of kind ReferenceElement.
     */
    virtual std::shared_ptr<ReferenceElement<dim,range,rank> > clone() const
    {
        Assert(false,ExcMessage("This function must not be called. "
                                "You should call the clone() funtion of a derived base class."));
        return nullptr;
    }



    /**
     * @name Functions for the basis and field evaluations without the use of
     * the cache.
     */
    ///@{
    /**
     * Returns a ValueTable with the <tt>deriv_order</tt>-th derivatives of all local basis function
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_cache()/fill_cache().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    template <int deriv_order>
    ValueTable<
    Conditional< deriv_order==0,
                 Value,
                 Derivative<deriv_order> > >
                 evaluate_basis_derivatives_at_points(const EvaluationPoints<dim> &points);

    ValueTable<Value>
    evaluate_basis_values_at_points(const EvaluationPoints<dim> &points)
    {
        return this->template evaluate_basis_derivatives_at_points<0>(points);
    }

    ValueTable<Derivative<1> >
    evaluate_basis_gradients_at_points(const EvaluationPoints<dim> &points)
    {
        return this->template evaluate_basis_derivatives_at_points<1>(points);
    }

    ValueTable<Derivative<2> >
    evaluate_basis_hessians_at_points(const EvaluationPoints<dim> &points)
    {
        return this->template evaluate_basis_derivatives_at_points<2>(points);
    }
    ///@}

};



IGA_NAMESPACE_CLOSE


#endif // #ifndef REFERENCE_ELEMENT_H_

