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


#ifndef __PHYSICALSPACEELEMENTACCESSOR_H
#define __PHYSICALSPACEELEMENTACCESSOR_H

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
//TODO: inline at least all the getters
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
 * @ingroup accessors_iterators
 *
 * @tparam PhysSpace - type for the space in the physical domain
 */
template<class PhysSpace>
class PhysicalSpaceElementAccessor
    :
public PhysSpace::RefSpace::ElementAccessor,
private PhysSpace::PushForwardType::ElementAccessor
{
public :
    using ContainerType = PhysSpace;

    using RefSpace = typename PhysSpace::RefSpace;
    using PushForwardType = typename PhysSpace::PushForwardType;
    using PfElemAccessor = typename PushForwardType::ElementAccessor;
    using RefElemAccessor = typename RefSpace::ElementAccessor;

    using RefElemAccessor::dim;
    using PfElemAccessor::space_dim;
    using PfElemAccessor::codim;
    using PfElemAccessor::transformation_type;



    /** Type for the quadrature scheme. */
    using QuadratureType = Quadrature<dim>;
    using QuadratureFaceType = Quadrature<dim-1>;

    static const Size n_faces = UnitElement<dim>::faces_per_element;


    using Value = typename PfElemAccessor::template PhysValue<RefSpace::dim_range, RefSpace::rank>;
    using ValueMap = typename PfElemAccessor::MappingElementAccessor::ValueMap;

    /**
     * Typedef for specifying the derivatives of the basis function in the physical domain.
     * \tparam order - order of the derivative.
     */
    template <int order>
    using Derivative = typename PfElemAccessor::template PhysDerivative<RefSpace::dim_range, RefSpace::rank, order>;

    using Div = Values<dim, 1, 1>;

    /**
     * @name Constructors
     */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    PhysicalSpaceElementAccessor() = delete;

    PhysicalSpaceElementAccessor(const PhysSpace &phys_space,
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

    int get_flat_index() const;


    /**
     * @name Management of the cache used in PhysicalSpaceElementAccessor
     */
    ///@{
    void init_values(const ValueFlags fill_flag,
                     const QuadratureType &quad);

    void init_face_values(const Index face_id,
                          const ValueFlags fill_flag,
                          const QuadratureFaceType &quad);

    void fill_values(const TopologyId &topology_id = ElemTopology());

    void fill_face_values(const Index face_id);
    ///@}


    /**
     * @name Getting the basis functions values, derivatives and hessians
     */
    ///@{
    ValueTable<Value> const &
    get_basis_values(const TopologyId &topology_id = ElemTopology()) const;

    /**
     * @brief Return the one-dimensional container with the
     * values of the i-th basis function at the evaluation points.
     */
    typename ValueTable<Value>::const_view
    get_basis_values(const Index func,const TopologyId &topology_id = ElemTopology()) const;

    const Value &
    get_basis_value(const Index func, const Index qp,const TopologyId &topology_id = ElemTopology()) const;


    ValueTable<Derivative<1> > const &
    get_basis_gradients(const TopologyId &topology_id = ElemTopology()) const;

    /**
     * \brief Return the one-dimensional container with the
     * gradients of the i-th basis function at the evaluation points.
     */
    typename ValueTable< Derivative<1> >::const_view
    get_basis_gradients(const Index func,const TopologyId &topology_id = ElemTopology()) const;

    const Derivative<1> &
    get_basis_gradient(const Index func, const Index qp,const TopologyId &topology_id = ElemTopology()) const;


    ValueTable<Derivative<2> > const &
    get_basis_hessians(const TopologyId &topology_id = ElemTopology()) const;

    /**
     * \brief Return the one-dimensional container with the
     * hessians of the i-th basis function at the evaluation points.
     */
    typename ValueTable< Derivative<2> >::const_view
    get_basis_hessians(const Index func,const TopologyId &topology_id = ElemTopology()) const;

    const Derivative<2> &
    get_basis_hessian(const Index func, const Index qp,const TopologyId &topology_id = ElemTopology()) const;

    Real get_basis_divergence(const Index func, const Index qp,const TopologyId &topology_id = ElemTopology()) const;
    ///@}


    /**
     * @name Getting quantities that are geometry-related
     */
    ///@{
    /**
     * Returns the quadrature weights multiplied by the
     * gradient determinant of map at the dilated quadrature points.
     */
    const ValueVector<Real> &
    get_w_measures(const TopologyId &topology_id = ElemTopology()) const;


    /**
     * @todo Document this function
     */
    const Point<space_dim> &
    get_point(const Index qp,const TopologyId &topology_id = ElemTopology()) const;

    /**
     * \brief Return a const reference to the one-dimensional container with the values of the map at the evaluation points.
     * \return The const reference to the one-dimensional container with the values of the map at the evaluation points.
     * \author M.Martinelli
     * \date 29 Jan 2013
     */
    const ValueVector< typename Mapping<dim,codim>::ValueType > &
    get_points(const TopologyId &topology_id = ElemTopology()) const;

    /**
     * \brief Return a const reference to the one-dimensional container with the gradients of the map (i.e. the Jacobian) at the evaluation points.
     * \return The const reference to the one-dimensional container with the values of the map at the evaluation points.
     * \author M.Martinelli
     * \date 03 Jun 2013
     */
    const ValueVector< typename Mapping<dim,codim>::GradientType > &
    get_map_gradient_at_points(const TopologyId &topology_id = ElemTopology()) const;



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
    evaluate_field(const std::vector<Real> &local_coefs,const TopologyId &topology_id = ElemTopology()) const;

    ValueVector< Derivative<1> >
    evaluate_field_gradients(const std::vector<Real> &local_coefs,const TopologyId &topology_id = ElemTopology()) const;

    ValueVector< Derivative<2> >
    evaluate_field_hessians(const std::vector<Real> &local_coefs,const TopologyId &topology_id = ElemTopology()) const;
    ///@}


    /**
     * Return a pointer to the physical space on which the element is defined.
     */
    const PhysSpace *get_physical_space() const;


    typedef PushForwardElementAccessor< typename PhysSpace::PushForwardType> push_forward_element_accessor;

private :

    /**
     * Typedef for specifying the derivatives of the basis function in the reference domain.
     * \tparam order - order of the derivative.
     */
    template< int order >
    using DerivativeRef_t = Derivatives<dim,RefSpace::dim_range,RefSpace::rank,order>;



protected:
    const PhysSpace *phys_space_;


    struct ValuesCache : CacheStatus
    {
        void reset(const int n_basis_per_element,
                   const QuadratureType &quad,
                   const BasisElemValueFlagsHandler &flags_handler);


        BasisElemValueFlagsHandler flags_handler_;

        Size n_points_ = 0;

        ValueTable<Value>         D0phi_;
        ValueTable<Derivative<1>> D1phi_;
        ValueTable<Derivative<2>> D2phi_;
    };


    struct ElementValuesCache : ValuesCache
    {
        void reset(const int n_basis_per_element,
                   const QuadratureType &quad,
                   const BasisElemValueFlagsHandler &flags_handler);
    };


    struct FaceValuesCache : ValuesCache
    {
        void reset(const Index face_id,
                   const int n_basis_per_element,
                   const QuadratureType &quad,
                   const BasisFaceValueFlagsHandler &flags_handler);


        void reset(const Index face_id,
                   const int n_basis_per_element,
                   const QuadratureFaceType &quad,
                   const BasisFaceValueFlagsHandler &flags_handler);
    };

    /**
     * For a given flags input argument identifies the face quantities and
     * returns a new ValueFlags variable containing only face quantities.
     * The output flags does not contain the word face.
     */
    ValueFlags get_face_flags(const ValueFlags fill_flag) const ;


    const ValuesCache &get_values_cache(const TopologyId &topology_id) const;
    ValuesCache &get_values_cache(const TopologyId &topology_id);


    ElementValuesCache elem_values_;

    std::array<FaceValuesCache, n_faces> face_values_;

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


public :
    template <typename Accessor> friend class GridForwardIterator;

};


IGA_NAMESPACE_CLOSE

#endif // __PHYSICALSPACEELEMENTACCESSOR_H
