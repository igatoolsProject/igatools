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


#ifndef PHYSICAL_SPACE_ELEMENT_H
#define PHYSICAL_SPACE_ELEMENT_H

#include <igatools/base/config.h>

#include <igatools/base/quadrature.h>
#include <igatools/basis_functions/bspline_element.h>
//#include <igatools/basis_functions/nurbs_element_accessor.h>
#include <igatools/geometry/push_forward_element.h>

IGA_NAMESPACE_OPEN

template <typename Accessor> class GridForwardIterator;

template<class PhysSpace>
class PhysicalSpaceElement
    :
    public SpaceElement<PhysSpace>,
    // todo: private PhysSpace::RefSpace::ElementAccessor,
    private PhysSpace::PushForwardType::ElementAccessor
{
public :
    using parent_t = SpaceElement<PhysSpace>;
    // using parent_t::get_elem_cache;
    using parent_t::LocalCache;

    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const PhysSpace;

    using Space = PhysSpace;
    using RefSpace = typename PhysSpace::RefSpace;
    using PushForwardType = typename PhysSpace::PushForwardType;
    using PfElemAccessor = typename PushForwardType::ElementAccessor;
    using RefElemAccessor = typename RefSpace::ElementAccessor;

    using PfElemAccessor::dim;
    using PfElemAccessor::space_dim;
    using PfElemAccessor::codim;
    using PfElemAccessor::type;


    /** Type for the quadrature scheme. */
    using QuadratureType = Quadrature<dim>;

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
    using PhysSpace::PushForwardType::ElementAccessor::get_measures;

#if 0
    /**
     * Returns the gradient determinant of the map at the dilated quadrature points
     * on the face specified by @p face_id.
     */
    using PhysSpace::PushForwardType::ElementAccessor::get_face_measures;
#endif

    /**
     * Returns the quadrature weights multiplied by the
     * gradient determinant of map at the dilated quadrature points.
     */
    using PhysSpace::PushForwardType::ElementAccessor::get_w_measures;

#if 0
    /**
     * Returns the quadrature weights multiplied by the
     * gradient determinant of map at the dilated quadrature points
     * on the face specified by @p face_id.
     */
    using PhysSpace::PushForwardType::ElementAccessor::get_face_w_measures;


    /**
     * @todo Document this function
     */
    const PhysPoint &
    get_point(const Index qp,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns a const reference to the one-dimensional container with the values
     * of the map at the evaluation points.
     */
    const ValueVector< typename Mapping<dim,codim>::Value > &
    get_points() const;

    /**
     * Returns a const reference to the one-dimensional container with the values
     * of the map at the evaluation points on the face specified by @p face_id.
     */
    const ValueVector< typename Mapping<dim,codim>::Value > &
    get_face_points(const Index face_id) const;

    /**
     * \brief Return a const reference to the one-dimensional container with the gradients of the map (i.e. the Jacobian) at the evaluation points.
     * \return The const reference to the one-dimensional container with the values of the map at the evaluation points.
     * \author M.Martinelli
     * \date 03 Jun 2013
     */
    using PhysSpace::PushForwardType::ElementAccessor::get_map_gradients;


    /**
     * Return a const reference to the one-dimensional container with the normals at the face evaluation points.
     */
    using PhysSpace::PushForwardType::ElementAccessor::get_face_normals;

    /**
     * Test if the element has a boundary face.
     */
    bool is_boundary() const;

    /**
     * Test if the face @p face on the current element is on the boundary.
     */
    bool is_boundary(const Index face) const;

    ///@}


    /**
     * @name Functions for field evaluation (values, gradient and hessian)
     */
    ///@{
    //Fields related
    ValueVector< Value >
    evaluate_field(const vector<Real> &local_coefs,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    ValueVector< Derivative<1> >
    evaluate_field_gradients(const vector<Real> &local_coefs,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    ValueVector< Derivative<2> >
    evaluate_field_hessians(const vector<Real> &local_coefs,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;
    ///@}


    /**
     * Return a pointer to the physical space on which the element is defined.
     */
    std::shared_ptr<const PhysSpace> get_physical_space() const;


    using  push_forward_element_accessor = PushForwardElementAccessor< typename PhysSpace::PushForwardType>;


#endif

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
    RefElemAccessor &get_ref_space_accessor()
    {
        return ref_space_element_accessor_;
    }
#if 0

    /**
     * Return a const reference of this object as would be viewed as push-forward element accessor.
     * This means that the returned object can be queried (but not modified) as the push-forward
     * element accessor that is used as partial inheritance of the physical space element accessor.
     */
    const PfElemAccessor &get_push_forward_accessor() const;
#endif
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
     * Moves the element to the position that differs from the current one
     * for the quantity given by @p increment.
     *
     * If the resulting position after the movement is valid (i.e. within the grid), then the function
     * returns true, otherwise it returns false.
     */
    bool jump(const TensorIndex<dim> &increment);

    /**
     * Sets the index of the element using the flatten representation.
     * @note This function also updates the index for the tensor representation.
     * @warning This may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    void move_to(const Index flat_index);


    /**
     * Sets the index of the element using the tensor representation.
     * @note This function also updates the index for the flatten representation.
     * @warning this may be a dangerous function, be careful when using it
     * as it is easy to use incorrectly. Only use it if you know what you
     * are doing.
     */
    void move_to(const TensorIndex<dim> &tensor_index);

    /** Moves the element to the next valid element in the CartesianGrid. */
    void operator++();
    ///@}



    bool operator==(const PhysicalSpaceElement<PhysSpace> &a) const;
    bool operator!=(const PhysicalSpaceElement<PhysSpace> &a) const;
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
    template <typename Accessor> friend class GridForwardIterator;
    template <typename PSpace> friend class SpaceElementHandler;
    RefElemAccessor ref_space_element_accessor_;
};


IGA_NAMESPACE_CLOSE

#endif
