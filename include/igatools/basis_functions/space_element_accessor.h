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


#ifndef SPACE_ELEMENT_ACCESSOR_H_
#define SPACE_ELEMENT_ACCESSOR_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/value_flags_handler.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/topology.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/utils/value_vector.h>
#include <igatools/utils/value_table.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/cartesian_product_indexer.h>

IGA_NAMESPACE_OPEN

template <typename Accessor> class GridForwardIterator;


/**
 * @brief Base class for element accessors on function spaces.
 *
 * This class is the base class for all element accessors on function spaces.
 * It collects all common functions for getting basis values (and derivatives)
 * and field values (and derivatives).
 *
 * Its design fulfills the
 * <a href="http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern"><em>Curiously Recurring Template Pattern (CRTP)</em></a>
 *  in order to achieve <em>static polymorphism</em>, which is an imitation of polymorphism
 * in programming code but which is resolved at compile time and thus does away with run-time
 * virtual-table lookups.
 *
 *
 * @ingroup accessors
 * @todo Complete the documentation
 * @author M. Martinelli
 * @date 13 May 2014
 */
// TODO (pauletti, Jun 5, 2014): from space one can get dim, codim, range and rank
template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
class SpaceElementAccessor : public CartesianGridElementAccessor<dim>
{
public:
    /** @name Types and aliases used and/or returned by the SpaceElementAccessor's methods. */
    ///@{
    /**
     * Typedef for specifying the value of the basis function.
     */
    using Value = typename Space::Value;

    using Point = typename Space::Point;

    /**
     * Typedef for specifying the divergence of the basis function.
     */
    using Div = typename Space::Div;

    /**
     * Typedef for specifying the derivatives of the basis function.
     */
    template <int order>
    using Derivative = typename Space::template Derivative<order>;


    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentContainer = typename Space::template ComponentContainer<T>;

    ///@}


    /** Number of faces per element. */
    static const Size n_faces = UnitElement<dim>::faces_per_element;



    /** Fill flags supported by this iterator */
    static const ValueFlags admisible_flag =
        ValueFlags::point|
        ValueFlags::measure |
        ValueFlags::w_measure |
        ValueFlags::face_point |
        ValueFlags::face_w_measure |
        ValueFlags::value |
        ValueFlags::gradient |
        ValueFlags::hessian |
        ValueFlags::divergence |
        ValueFlags::face_value |
        ValueFlags::face_gradient |
        ValueFlags::face_hessian |
        ValueFlags::face_divergence;


    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    SpaceElementAccessor() = delete;

    /**
     * Constructs an accessor to element number index of a
     * BsplineSpace space.
     */
    SpaceElementAccessor(const std::shared_ptr<const Space> space,
                         const int elem_index);


    /**
     * Copy constructor.
     * @note For the constructed object it
     * creates a new element cache, but it shares
     * the one dimensional cache with the copied element.
     */
    SpaceElementAccessor(const SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank> &elem)
        = default;

    /**
     * Move constructor.
     */
    SpaceElementAccessor(SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank> &&elem)
        = default;

    /**
     * Destructor.
     */
    ~SpaceElementAccessor() = default;
    ///@}



    /** @name Functions for the basis and field evaluations without the use of the cache */
    ///@{

    /**
     * Returns a ValueTable with the values of all local basis function
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueTable<Value>
    evaluate_basis_values_at_points(const std::vector<Point> &points) const;


    /**
     * Returns a ValueTable with the gradients of all local basis function
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueTable< Derivative<1> >
    evaluate_basis_gradients_at_points(const std::vector<Point> &points) const;


    /**
     * Returns a ValueTable with the hessians of all local basis function
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueTable< Derivative<2> >
    evaluate_basis_hessians_at_points(const std::vector<Point> &points) const;


    /**
     * Returns a ValueVector with the <tt>deriv_order</tt>-th derivatives of the field
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    template <int deriv_order>
    ValueVector< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >
    evaluate_field_derivatives_at_points(
        const std::vector<Real> &local_coefs,
        const std::vector<Point> &points) const;


    /**
     * Returns a ValueVector with the values of the field
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueVector<Value>
    evaluate_field_values_at_points(
        const std::vector<Real> &local_coefs,
        const std::vector<Point> &points) const;


    /**
     * Returns a ValueVector with the gradients of the field
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueVector< Derivative<1> >
    evaluate_field_gradients_at_points(
        const std::vector<Real> &local_coefs,
        const std::vector<Point> &points) const;


    /**
     * Returns a ValueVector with the hessians of the field
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueVector< Derivative<2> >
    evaluate_field_hessians_at_points(
        const std::vector<Real> &local_coefs,
        const std::vector<Point> &points) const;
    ///@}



    /** @name Functions returning the value of the basis functions. */
    ///@{
    /**
     * Returns the const reference to a ValueTable with the values of all local basis function
     * at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    ValueTable<Value> const &
    get_basis_values(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns a const view to the values of the <tt>i</tt>-th basis function at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    typename ValueTable<Value>::const_view
    get_basis_values(const Index i,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the const reference to the value of a local basis function
     * at one evaluation point.
     * @param[in] basis Local basis id.
     * @param[in] qp Point id.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    Value const &
    get_basis_value(const Index basis, const Index qp,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the const reference to a ValueTable with the values of all local basis function
     * at each evaluation point on the face specified by @p face_id.
     */
    ValueTable<Value> const &
    get_face_basis_values(const Index face_id) const;
    ///@}



    /** @name Functions returning the gradient of the basis functions. */
    ///@{
    /**
     * Returns the const reference to a ValueTable with the gradients of all local basis function
     * evaluated at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    ValueTable<Derivative<1>> const &get_basis_gradients(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns a const view to the gradients of the <tt>i</tt>-th basis function at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    typename ValueTable<Derivative<1> >::const_view
    get_basis_gradients(const Index i,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the const reference to the gradient of a local basis function
     * at one evaluation point.
     * @param[in] basis Local basis id.
     * @param[in] qp Point id.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    Derivative<1> const &get_basis_gradient(const Index basis, const Index qp,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the const reference to a ValueTable with the gradients of all local basis function
     * at each evaluation point on the face specified by @p face_id.
     */
    ValueTable<Derivative<1>> const &get_face_basis_gradients(const Index face_id) const;
    ///@}

    /** @name Functions returning the hessian of the basis functions. */
    ///@{
    /**
     * Returns the const reference to a ValueTable with hessians of all local basis function
     * at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    ValueTable<Derivative<2>> const &get_basis_hessians(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns a const view to the hessians of the <tt>i</tt>-th basis function at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    typename ValueTable<Derivative<2> >::const_view
    get_basis_hessians(const Index i,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the const reference to the hessian of a local basis function
     * at one evaluation point.
     * @param[in] basis Local basis id.
     * @param[in] qp Point id.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    Derivative<2> const &get_basis_hessian(const Index basis, const Index qp,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the const reference to a ValueTable with the hessians of all local basis function
     * at each evaluation point on the face specified by @p face_id.
     */
    ValueTable<Derivative<2>> const &get_face_basis_hessians(const Index face_id) const;
    ///@}

    /** @name Functions returning the divergence of the basis functions. */
    ///@{
    /**
     * Returns the const reference to a ValueTable with the values of all local basis function
     * at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    ValueTable<Div> const &get_basis_divergences(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns a const view to the divergences of the <tt>i</tt>-th basis function at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    typename ValueTable<Div>::const_view
    get_basis_divergences(const Index i,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the const reference to the divergence of a local basis function
     * at one evaluation point.
     * @param[in] basis Local basis id.
     * @param[in] qp Point id.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    Div const &get_basis_divergence(const Index basis, const Index qp,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the const reference to a ValueTable with the divergences of all local basis function
     * at each evaluation point on the face specified by @p face_id.
     */
    ValueTable<Div> const &get_face_basis_divergences(const Index face_id) const;
    ///@}




    /** @name Fields related */
    ///@{
    /**
     * Returns the ValueVector with the evaluation of the field @p local_coefs at the evaluation
     * points.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     * @see get_local_coefs
     */
    ValueVector<Value>
    evaluate_field(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;


    /**
     * Returns the ValueVector with the evaluation of the field @p local_coefs at the evaluation
     * points on the face specified by @p face_id.
     */
    ValueVector<Value>
    evaluate_face_field(const Index face_id, const std::vector<Real> &local_coefs) const;

    /**
     * Returns the ValueVector with the evaluation of the gradient of the field @p local_coefs
     * at the evaluation points.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     * @see get_local_coefs
     */
    ValueVector<Derivative<1> >
    evaluate_field_gradients(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the ValueVector with the evaluation of the gradient of the field @p local_coefs at the evaluation
     * points on the face specified by @p face_id.
     */
    ValueVector<Derivative<1> >
    evaluate_face_field_gradients(const Index face_id, const std::vector<Real> &local_coefs) const;

    /**
     * Returns the ValueVector with the evaluation of the hessians of the field @p local_coefs
     * at the evaluation points.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     * @see get_local_coefs
     */
    ValueVector<Derivative<2> >
    evaluate_field_hessians(const std::vector<Real> &local_coefs, const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the ValueVector with the evaluation of the hessian of the field @p local_coefs at the evaluation
     * points on the face specified by @p face_id.
     */
    ValueVector<Derivative<2> >
    evaluate_face_field_hessians(const Index face_id, const std::vector<Real> &local_coefs) const;

    /**
     * Returns the ValueVector with the evaluation of the divergences of the field @p local_coefs
     * at the evaluation points.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     * @see get_local_coefs
     */
    ValueVector<Div>
    evaluate_field_divergences(const std::vector<Real> &local_coefs, const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns the ValueVector with the evaluation of the divergence of the field @p local_coefs at the evaluation
     * points on the face specified by @p face_id.
     */
    ValueVector<Div>
    evaluate_face_field_divergences(const Index face_id, const std::vector<Real> &local_coefs) const;
    ///@}



    /** @name Query information without use of cache */
    ///@{
    /**
     *  Number of non zero basis functions over the current element.
     */
    Size get_num_basis() const;

    /**
     * Number of non-zero scalar basis functions associated
     * with the i-th space component on the element.
     * This makes sense as a reference B-spline space
     * is only allowed to be of the cartesian product type
     * V = V1 x V2 x ... X Vn.
     */
    int get_num_basis(const int i) const;

    /**
     * Returns the global dofs of the local (non zero) basis functions
     * on the element.
     * For example:
     * \code
       auto loc_to_glob = elem->get_local_to_global();
       // loc_to_glob[0] is the global id of the first element basis function
       // loc_to_glob[1] is the global id of the second element basis function
       // ...
      \endcode
     *
     */
    std::vector<Index> get_local_to_global() const;

    /**
     * Pointer to the BsplineSpace the accessor is iterating on.
     */
    std::shared_ptr<const Space> get_space() const;
    ///@}



    /**
     * Get the quadrature points used to initialize the element or a given element-face.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     * @see get_local_coefs
     */
    const Quadrature<dim> &
    get_quad_points(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;



protected:
    /**
     * Space for which the SpaceElementAccessor refers to.
     */
    std::shared_ptr<const Space> space_ = nullptr;


    /** Number of scalar basis functions along each direction, for all space components. */
    // const typename Space::SpaceDimensionTable &n_basis_direction_;

    /** Hash table for fast conversion between flat-to-tensor basis function ids. */
    ComponentContainer<std::shared_ptr<CartesianProductIndexer<dim> > > basis_functions_indexer_;

    /** Basis function ID offset between the different components. */
    ComponentContainer<int> comp_offset_;



    /**
     * Base class for the cache of the element values and for the cache of the face values.
     */
    class ValuesCache : public CacheStatus
    {
    public:
        /**
         * Allocate space for the values and derivatives
         * at quadrature points
         */
        void reset(const BasisElemValueFlagsHandler &flags_handler,
                   const ComponentContainer<TensorSize<dim> > &n_basis_direction,
                   const Quadrature<dim> &quad);

        /** Returns the values. */
        const ValueTable<Value> &get_values() const;

        /** Returns the gradients. */
        const ValueTable<Derivative<1>> &get_gradients() const;

        /** Returns the hessians. */
        const ValueTable<Derivative<2>> &get_hessians() const;

        /** Returns the divergences. */
        const ValueTable<Div> &get_divergences() const;


        //TODO: the member variables should be private
    public:

        BasisElemValueFlagsHandler flags_handler_;


        ValueTable<Value> phi_;
        ValueTable<Derivative<1>> D1phi_;
        ValueTable<Derivative<2>> D2phi_;

        ValueTable<Div> div_phi_;

        Quadrature<dim> quad_;
    };


    /**
     * Cache for the element values at quadrature points
     */
    class ElementValuesCache : public ValuesCache
    {
    public:
        /**
         * Allocate space for the values and derivatives
         * at quadrature points
         */
        void reset(const BasisElemValueFlagsHandler &flags_handler,
                   const ComponentContainer<TensorSize<dim> > &n_basis_direction,
                   const Quadrature<dim> &quad);

    };


    /**
     * Cache for the face values at quadrature points
     */
    class FaceValuesCache : public ValuesCache
    {
    public:
        /**
         * Allocate space for the values and derivatives
         * at quadrature points
         */
        void reset(const Index face_id,
                   const BasisFaceValueFlagsHandler &flags_handler,
                   const ComponentContainer<TensorSize<dim> > &n_basis_direction,
                   const Quadrature<dim> &quad);

        /**
         * Allocate space for the values and derivatives
         * at quadrature points for a specified face.
         */
        void reset(const Index face_id,
                   const BasisFaceValueFlagsHandler &flags_handler,
                   const ComponentContainer<TensorSize<dim> > &n_basis_direction,
                   const Quadrature<dim-1> &quad);

    };


    /**
     * Element cache to store the values and derivatives
     * of the basis functions on the element
     */
    ElementValuesCache elem_values_;

    /**
     * Face cache to store the values and derivatives
     * of the basis functions on the faces of the element
     */
    std::array<FaceValuesCache, n_faces> face_values_;


    /**
     * Initializes the element and faces cache according to
     * the quadrature number of point and the fill_flag.
     */
    void reset_element_and_faces_cache(
        const ValueFlags fill_flag,
        const Quadrature<dim> &quad);


public:

    /**
     * Fills the values cache of the <tt>face_id</tt>-th face, according to the evaluation points
     * and fill flags specifies in init_values.
     */
    void fill_face_values(const Index face_id);


    /**
     * @todo Document this function
     */
    const ValuesCache &get_values_cache(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * @todo Document this function
     */
    ValuesCache &get_values_cache(const TopologyId<dim> &topology_id = ElemTopology<dim>());


    /** Return a reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    CartesianGridElementAccessor<dim> &as_cartesian_grid_element_accessor();


    /** Return a const-reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    const CartesianGridElementAccessor<dim> &as_cartesian_grid_element_accessor() const;

private:

    /** Return a reference to "*this" as being an object of type DerivedElementAccessor.*/
    DerivedElementAccessor &as_derived_element_accessor();

    /** Return a const-reference to "*this" as being an object of type DerivedElementAccessor.*/
    const DerivedElementAccessor &as_derived_element_accessor() const;

};


IGA_NAMESPACE_CLOSE


#include <igatools/basis_functions/space_element_accessor-inline.h>

#endif // #ifndef SPACE_ELEMENT_ACCESSOR_
