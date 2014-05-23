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


#ifndef PHYSICAL_SPACE_ELEMENT_ACCESSOR_H
#define PHYSICAL_SPACE_ELEMENT_ACCESSOR_H

#include <igatools/base/config.h>

#include <igatools/geometry/mapping.h>
#include <igatools/base/quadrature.h>
#include <igatools/base/function.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/basis_functions/nurbs_element_accessor.h>
#include <igatools/geometry/push_forward_element_accessor.h>

IGA_NAMESPACE_OPEN

template < typename Accessor > class GridForwardIterator;

//TODO: trilinos_vector.h should be called vector.h
//TODO: document this class

/**
 *
 * This class represents a single element in the physical space.
 *
 * Its main use is for the evaluation of quantities (at some evaluation points specified as input argument)
 * defined over the element. Namely, it provides the methods to evaluate:
 * - non-zero basis functions. The quantities that can be computed are:
 *   - values;
 *   - gradients;
 *   - hessians;
 *   - divergence.
 * - fields (i.e. a linear combination of the non-zero basis function). The quantities that can be computed are:
 *   - values;
 *   - gradients;
 *   - hessians;
 * - weights of the quadrature points;
 * - evaluation points in the physical space (i.e. the image of the evaluation points from the \f$ [0,1]^{\text{dim\_domain}} \f$ domain to
 * the physical element using the mapping defined in the PushForward.
 *
 * Many of the above quantities can be retrieved with two different interfaces:
 * - the first set of interfaces need as input argument the evaluation point index.
 *   The retrieved quantities relates to the evaluation point specified by the input index.
 * - the second set of interfaces retrieve the quantity of interests for all the evaluation point.
 *   This kind of interfaces return the reference of a container (no copy is made) and should be used when
 *   good performances are requested.
 *
 * See module on @ref accessors_iterators for a general overview.
 * @ingroup accessors
 *
 * @tparam PhysSpace - type for the space in the physical domain
 */
template<class PhysSpace>
class PhysicalSpaceElementAccessor
    :
public SpaceElementAccessor<
    PhysicalSpaceElementAccessor<PhysSpace>,PhysSpace,
    PhysSpace::RefSpace::dim,
    PhysSpace::PushForwardType::codim,
    PhysSpace::RefSpace::range,
    PhysSpace::RefSpace::rank>,
private PhysSpace::PushForwardType::ElementAccessor
{
public :
    using parent_t = SpaceElementAccessor<
                     PhysicalSpaceElementAccessor<PhysSpace>,PhysSpace,
                     PhysSpace::RefSpace::dim,
                     PhysSpace::PushForwardType::codim,
                     PhysSpace::RefSpace::range,
                     PhysSpace::RefSpace::rank>;


    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const PhysSpace;

    using RefSpace = typename PhysSpace::RefSpace;
    using PushForwardType = typename PhysSpace::PushForwardType;
    using PfElemAccessor = typename PushForwardType::ElementAccessor;
    using RefElemAccessor = typename RefSpace::ElementAccessor;

    using PfElemAccessor::dim;
    using PfElemAccessor::space_dim;
    using PfElemAccessor::codim;
    using PfElemAccessor::transformation_type;


    /** Type for the quadrature scheme. */
    using QuadratureType = Quadrature<dim>;
    using QuadratureFaceType = Quadrature<dim-1>;

    static const Size n_faces = UnitElement<dim>::faces_per_element;


    using Value = typename PfElemAccessor::template PhysValue<RefSpace::range, RefSpace::rank>;
    using ValueMap = typename PfElemAccessor::MappingElementAccessor::ValueMap;

    /**
     * Typedef for specifying the derivatives of the basis function in the physical domain.
     * \tparam order - order of the derivative.
     */
    template <int order>
    using Derivative = typename PfElemAccessor::template PhysDerivative<RefSpace::range, RefSpace::rank, order>;

    using Div = Values<dim, 1, 1>;

    /**
     * @name Constructors
     */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    PhysicalSpaceElementAccessor() = delete;

    PhysicalSpaceElementAccessor(const std::shared_ptr<ContainerType> space,
                                 const Index index);

    /**
     * Copy constructor.
     * @note Performs a deep copy of the PhysicalSpaceElementAccessor<PhysSpace> @p in,
     * except for the pointer to the PhysicalSpace.
     */
    PhysicalSpaceElementAccessor(const PhysicalSpaceElementAccessor<PhysSpace> &in) = default;

    /**
     * Move constructor.
     */
    PhysicalSpaceElementAccessor(PhysicalSpaceElementAccessor<PhysSpace> &&in) = default;

    /**
     * Destructor.
     */
    ~PhysicalSpaceElementAccessor() = default;

    ///@}

    /**
     * @name Assignment operators
     */
    ///@{
    /**
     * Copy assignment operator.
     * @note Performs a deep copy of the PhysicalSpaceElementAccessor<PhysSpace> @p in,
     * except for the pointer to the PhysicalSpace.
     */
    PhysicalSpaceElementAccessor<PhysSpace> &
    operator=(const PhysicalSpaceElementAccessor<PhysSpace> &in) = default;

    /**
     * Move assignment operator.
     */
    PhysicalSpaceElementAccessor<PhysSpace> &
    operator=(PhysicalSpaceElementAccessor<PhysSpace> &&in) = default;

    ///@}



    /**
     * @name Management of the cache used in PhysicalSpaceElementAccessor
     */
    ///@{
    void init_values(const ValueFlags fill_flag,
                     const QuadratureType &quad);

    void init_face_values(const Index face_id,
                          const ValueFlags fill_flag,
                          const QuadratureFaceType &quad);

    void fill_values(const TopologyId<dim> &topology_id = ElemTopology<dim>());

    void fill_face_values(const Index face_id);
    ///@}




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
    evaluate_basis_derivatives_at_points(const std::vector<Point<dim>> &points) const;

    ///@}


    /**
     * @name Getting quantities that are geometry-related
     */
    ///@{

    /**
     * Returns the gradient determinant of the map at the dilated quadrature points.
     */
    const ValueVector< Real > &
    get_measures(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the gradient determinant of the map at the dilated quadrature points
     * on the face specified by @p face_id.
     */
    const ValueVector< Real > &
    get_face_measures(const Index face_id) const;


    /**
     * Returns the quadrature weights multiplied by the
     * gradient determinant of map at the dilated quadrature points.
     */
    const ValueVector<Real> &
    get_w_measures(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the quadrature weights multiplied by the
     * gradient determinant of map at the dilated quadrature points
     * on the face specified by @p face_id.
     */
    const ValueVector<Real> &
    get_face_w_measures(const Index face_id) const;


    /**
     * @todo Document this function
     */
    const Point<space_dim> &
    get_point(const Index qp,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns a const reference to the one-dimensional container with the values
     * of the map at the evaluation points.
     */
    const ValueVector< typename Mapping<dim,codim>::ValueType > &
    get_points(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns a const reference to the one-dimensional container with the values
     * of the map at the evaluation points on the face specified by @p face_id.
     */
    const ValueVector< typename Mapping<dim,codim>::ValueType > &
    get_face_points(const Index face_id) const;

    /**
     * \brief Return a const reference to the one-dimensional container with the gradients of the map (i.e. the Jacobian) at the evaluation points.
     * \return The const reference to the one-dimensional container with the values of the map at the evaluation points.
     * \author M.Martinelli
     * \date 03 Jun 2013
     */
    const ValueVector< typename Mapping<dim,codim>::GradientType > &
    get_map_gradient_at_points(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;



    /**
     * Return a const reference to the one-dimensional container with the normals at the face evaluation points.
     */
    const ValueVector< typename Mapping<dim,codim>::ValueType > &
    get_face_normals(const Index face_id) const;

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
    evaluate_field(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    ValueVector< Derivative<1> >
    evaluate_field_gradients(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    ValueVector< Derivative<2> >
    evaluate_field_hessians(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;
    ///@}


    /**
     * Return a pointer to the physical space on which the element is defined.
     */
    std::shared_ptr<const PhysSpace>get_physical_space() const;


    using  push_forward_element_accessor = PushForwardElementAccessor< typename PhysSpace::PushForwardType>;


    /**
     * Prints internal information about the BSplineElementAccessor.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out, const VerbosityLevel verbosity_level = VerbosityLevel::normal) const;

    /**
     * Return a const reference of this object as would be viewed as reference space element accessor.
     * This means that the returned object can be queried (but not modified) as the reference space
     * element accessor that is used as partial inheritance of the physical space element accessor.
     */
    const RefElemAccessor &get_ref_space_accessor() const;

    /**
     * Return a const reference of this object as would be viewed as push-forward element accessor.
     * This means that the returned object can be queried (but not modified) as the push-forward
     * element accessor that is used as partial inheritance of the physical space element accessor.
     */
    const PfElemAccessor &get_push_forward_accessor() const;


    /** Returns the index of the element in its flatten representation. */
    Index get_flat_index() const ;



protected:


    /**
     * For a given flags input argument identifies the face quantities and
     * returns a new ValueFlags variable containing only face quantities.
     * The output flags does not contain the word face.
     */
    ValueFlags get_face_flags(const ValueFlags fill_flag) const ;

    void operator++();

    bool operator==(const PhysicalSpaceElementAccessor <PhysSpace> &a) const;

    bool operator!=(const PhysicalSpaceElementAccessor <PhysSpace> &a) const;

    /**
     * This function returns the ValueFlags needed to be passed to the ReferenceSpacePhysicalAccessor
     * in order to compute the quantities specified by the input argument
     * @p fill_flag (i.e. the ValueFlags that refers to the PhysicalSpaceElementAccessor).
     */
    ValueFlags get_reference_space_accessor_fill_flags(const ValueFlags fill_flag) const;

    /**
     * This function returns the ValueFlags needed to be passed to the PushForwardAccessor
     * in order to compute the quantities specified by the input argument
     * @p fill_flag (i.e. the ValueFlags that refers to the PhysicalSpaceElementAccessor).
     */
    ValueFlags get_push_forward_accessor_fill_flags(const ValueFlags fill_flag) const;

    // TODO (pauletti, Apr 23, 2014): why not private?
public :


private:
    template <typename Accessor> friend class GridForwardIterator;

    RefElemAccessor ref_space_element_accessor_;
};


IGA_NAMESPACE_CLOSE

#endif // PHYSICAL_SPACE_ELEMENT_ACCESSOR_H
