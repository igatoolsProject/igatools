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


#ifndef PHYSICAL_SPACE_ELEMENT_H
#define PHYSICAL_SPACE_ELEMENT_H

#include <igatools/base/config.h>

#include <igatools/base/evaluation_points.h>
#include <igatools/basis_functions/reference_element.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/nurbs_element.h>
#include <igatools/geometry/push_forward_element.h>

IGA_NAMESPACE_OPEN

template <class Accessor> class CartesianGridIterator;

/**
 *
 * @ingroup elements
 */
template<class PhysSpace>
class PhysicalSpaceElement
    :
    public SpaceElement<PhysSpace>
//,
    // todo: private PhysSpace::RefSpace::ElementAccessor,
//    private PhysSpace::PushForwardType::ElementAccessor
{
public :
    using parent_t = SpaceElement<PhysSpace>;
    // using parent_t::get_elem_cache;
    using parent_t::LocalCache;
    using parent_t::get_values;
//    using parent_t::is_boundary;

    /** Type required by the CartesianGridIterator templated iterator */
    using ContainerType = const PhysSpace;

    using Space = PhysSpace;
    using RefSpace = typename PhysSpace::RefSpace;
    using PushForwardType = typename PhysSpace::PushForwardType;
    using PfElemAccessor = typename PushForwardType::ElementAccessor;
    using RefElemAccessor = typename RefSpace::ElementAccessor;
    /*
        using PfElemAccessor::dim;
        using PfElemAccessor::space_dim;
        using PfElemAccessor::codim;
        using PfElemAccessor::type;
    //*/
    static const auto dim = PfElemAccessor::dim;
    static const auto space_dim = PfElemAccessor::space_dim;
    static const auto codim = PfElemAccessor::codim;
    static const auto type = PfElemAccessor::type;

#if 0
    using Value = typename PhysSpace::Value;
    using PhysPoint = typename PhysSpace::Point;
    using RefPoint = typename RefSpace::Point;


    using ValueMap = typename PfElemAccessor::MappingElementAccessor::ValueMap;

    /**
     * Typedef for specifying the derivatives of the basis function in the physical domain.
     * \tparam order - order of the derivative.
     */
    template <int order>
    using Derivative = typename PfElemAccessor::template PhysDerivative<RefSpace::range, RefSpace::rank, order>;

    using Div = typename PhysSpace::Div;
#endif
    /**
     * @name Constructors
     */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    PhysicalSpaceElement() = delete;

    PhysicalSpaceElement(const std::shared_ptr<ContainerType> space,
                         const Index index);


    PhysicalSpaceElement(const std::shared_ptr<ContainerType> space,
                         const TensorIndex<dim> &index);

    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a classic copy constructor)
     * uses the deep copy.
     */
    PhysicalSpaceElement(const PhysicalSpaceElement<PhysSpace> &in,
                         const CopyPolicy &copy_policy = CopyPolicy::deep);


    /**
     * Move constructor.
     */
    PhysicalSpaceElement(PhysicalSpaceElement<PhysSpace> &&in) = default;

    /**
     * Destructor.
     */
    ~PhysicalSpaceElement() = default;

    ///@}

    /**
     * @name Assignment operators
     */
    ///@{
    /**
     * Copy assignment operator. Performs a <b>shallow copy</b> of the input @p element.
     *
     * @note Internally it uses the function shallow_copy_from().
     */
    PhysicalSpaceElement<PhysSpace> &
    operator=(const PhysicalSpaceElement<PhysSpace> &in) = default;

    /**
     * Move assignment operator.
     */
    PhysicalSpaceElement<PhysSpace> &
    operator=(PhysicalSpaceElement<PhysSpace> &&in) = default;

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
    void deep_copy_from(const PhysicalSpaceElement<PhysSpace> &element);


    /**
     * Performs a shallow copy of the input @p element. The current object will contain a pointer to the
     * local cache used by the input @p element.
     */
    void shallow_copy_from(const PhysicalSpaceElement<PhysSpace> &element);
    ///@}

#if 0
    /** @name Functions for the basis and field evaluations without the use of the cache */
    ///@{

    /**
     * Returns a ValueTable with the <tt>deriv_order</tt>-th derivatives of all local basis function
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    template <int deriv_order>
    ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >
    evaluate_basis_derivatives_at_points(const ValueVector<RefPoint> &points) const;

    ///@}

#endif
    /**
     * @name Getting quantities that are geometry-related
     */
    ///@{

    /**
     * Returns the gradient determinant of the map at the dilated quadrature points.
     */
    /*
    using PhysSpace::PushForwardType::ElementAccessor::get_points;
    using PhysSpace::PushForwardType::ElementAccessor::get_measures;
    using PhysSpace::PushForwardType::ElementAccessor::get_w_measures;
    //*/

    /**
     * Returns the <tt>k</tt> dimensional j-th sub-element measure
     * multiplied by the weights of the quadrature.
     */
    template <int k>
    ValueVector<Real> get_w_measures(const int j) const
    {
        return push_fwd_element_->template get_w_measures<k>(j);
    }


    /**
     * Prints internal information about the BSplineElementAccessor.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out) const;

    void print_cache_info(LogStream &out) const;

private:

    /**
     * Return a const reference of this object as would be viewed as reference space element accessor.
     * This means that the returned object can be queried (but not modified) as the reference space
     * element accessor that is used as partial inheritance of the physical space element accessor.
     */
    const RefElemAccessor &get_ref_space_accessor() const;
    RefElemAccessor &get_ref_space_accessor();


    /**
     * Return a const reference of this object as would be viewed as push-forward element accessor.
     * This means that the returned object can be queried (but not modified) as the push-forward
     * element accessor that is used as partial inheritance of the physical space element accessor.
     */
    const PfElemAccessor &get_push_forward_accessor() const;
    PfElemAccessor &get_push_forward_accessor();

public:

    /**
     * @name Functions related to get the indices of the element.
     */
    ///@{
    /** Returns the index of the element in its flatten representation. */
    Index get_flat_index() const;

    /** Returns the index of the element in its tensor representation. */
    TensorIndex<dim> get_tensor_index() const;
    ///@}


    /** Return the cartesian grid from which the element belongs.*/
    const std::shared_ptr<const CartesianGrid<dim>> get_grid() const;


protected:

#if 0
    /**
     * For a given flags input argument identifies the face quantities and
     * returns a new ValueFlags variable containing only face quantities.
     * The output flags does not contain the word face.
     */
    ValueFlags get_face_flags(const ValueFlags fill_flag) const ;

#endif
    /** @name Functions/operators for moving the element in the CartesianGrid.*/
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


    /*
        bool operator==(const PhysicalSpaceElement<PhysSpace> &a) const;
        bool operator!=(const PhysicalSpaceElement<PhysSpace> &a) const;
        bool operator<(const PhysicalSpaceElement<PhysSpace> &a) const;
        bool operator>(const PhysicalSpaceElement<PhysSpace> &a) const;
    //*/

#if 0
    /**
     * This function returns the ValueFlags needed to be passed to the ReferenceSpacePhysicalAccessor
     * in order to compute the quantities specified by the input argument
     * @p fill_flag (i.e. the ValueFlags that refers to the PhysicalSpaceElement).
     */
    ValueFlags get_reference_space_accessor_fill_flags(const ValueFlags fill_flag) const;

    /**
     * This function returns the ValueFlags needed to be passed to the PushForwardAccessor
     * in order to compute the quantities specified by the input argument
     * @p fill_flag (i.e. the ValueFlags that refers to the PhysicalSpaceElement).
     */
    ValueFlags get_push_forward_accessor_fill_flags(const ValueFlags fill_flag) const;

#endif
    /**
     * Performs a copy of the input @p element.
     * The type of copy (deep or shallow) is specified by the input parameter @p copy_policy.
     */
    void copy_from(const PhysicalSpaceElement<PhysSpace> &element,
                   const CopyPolicy &copy_policy);


private:
    template <class Accessor> friend class CartesianGridIteratorBase;
    template <typename PSpace> friend class PhysSpaceElementHandler;

    std::shared_ptr<RefElemAccessor> ref_space_element_accessor_;


    std::shared_ptr<PfElemAccessor> push_fwd_element_;

    /**
     * Creates a new object performing a deep copy of the current object using the PhysicalSpaceElement
     * copy constructor.
     */
    std::shared_ptr<PhysicalSpaceElement<PhysSpace> > clone() const
    {
        auto elem = std::shared_ptr<PhysicalSpaceElement<PhysSpace> >(
                        new PhysicalSpaceElement(*this,CopyPolicy::deep));
        Assert(elem != nullptr, ExcNullPtr());
        return elem;
    }
};


IGA_NAMESPACE_CLOSE

#endif
