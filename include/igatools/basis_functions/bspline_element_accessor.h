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


#ifndef BSPLINE_ELEMENT_ACCESSOR_H_
#define BSPLINE_ELEMENT_ACCESSOR_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/value_flags_handler.h>
#include <igatools/base/quadrature.h>
#include <igatools/base/function.h>
#include <igatools/geometry/topology.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/utils/value_vector.h>
#include <igatools/utils/value_table.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/cartesian_product_indexer.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/basis_functions/bernstein_basis.h>
#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>

IGA_NAMESPACE_OPEN

template <int dim, int range, int rank> class BSplineSpace;
template <typename Accessor> class GridForwardIterator;


/**
 * @todo Missing documentation
 */
template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
class SpaceElementAccessor : public CartesianGridElementAccessor<dim>
{
public:
    /** @name Types and aliases used and/or returned by the SpaceElementAccessor's methods. */
    ///@{
    /**
     * Typedef for specifying the value of the basis function.
     */
    using Value = Values<dim+codim, range, rank>;

    /**
     * Typedef for specifying the derivatives of the basis function.
     */
    template <int deriv_order>
    using Derivative = Derivatives<dim+codim, range, rank, deriv_order>;

    /**
     * Typedef for specifying the divergence of the basis function.
     */
    using Div = Values<dim+codim, 1, 1>;


    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentTable = StaticMultiArray<T,range,rank>;
    ///@}


    static const Size n_faces = UnitElement<dim>::faces_per_element;



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


    /** Return a reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    CartesianGridElementAccessor<dim> &as_cartesian_grid_element_accessor();


    /** Return a const-reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    const CartesianGridElementAccessor<dim> &as_cartesian_grid_element_accessor() const;


    /** Return a reference to "*this" as being an object of type DerivedElementAccessor.*/
    DerivedElementAccessor &as_derived_element_accessor();

    /** Return a const-reference to "*this" as being an object of type DerivedElementAccessor.*/
    const DerivedElementAccessor &as_derived_element_accessor() const;


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
    evaluate_basis_values_at_points(const std::vector<Point<dim>> &points) const;


    /**
     * Returns a ValueTable with the gradients of all local basis function
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueTable< Derivative<1> >
    evaluate_basis_gradients_at_points(const std::vector<Point<dim>> &points) const;


    /**
     * Returns a ValueTable with the hessians of all local basis function
     * at each point (in the unit domain) specified by the input argument <tt>points</tt>.
     * @note This function does not use the cache and therefore can be called any time without
     * needing to pre-call init_values()/fill_values().
     * @warning The evaluation <tt>points</tt> must belong to the unit hypercube
     * \f$ [0,1]^{\text{dim}} \f$ otherwise, in Debug mode, an assertion will be raised.
     */
    ValueTable< Derivative<2> >
    evaluate_basis_hessians_at_points(const std::vector<Point<dim>> &points) const;


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
        const std::vector<Point<dim>> &points) const;


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
        const std::vector<Point<dim>> &points) const;


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
        const std::vector<Point<dim>> &points) const;


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
        const std::vector<Point<dim>> &points) const;
    ///*}



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
    std::vector<Index> const &get_local_to_global() const;

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
    ComponentTable< TensorSize<dim> > n_basis_direction_;

    /** Hash table for fast conversion between flat-to-tensor basis function ids. */
    ComponentTable<std::shared_ptr<CartesianProductIndexer<dim> > > basis_functions_indexer_;

    /** Basis function ID offset between the different components. */
    ComponentTable<int> comp_offset_;



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
                   const StaticMultiArray<TensorSize<dim>,range,rank> &n_basis_direction,
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


        ValueTable<Value> phi_hat_;
        ValueTable<Derivative<1>> D1phi_hat_;
        ValueTable<Derivative<2>> D2phi_hat_;

        ValueTable<Div> div_phi_hat_;

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
                   const ComponentTable<TensorSize<dim> > &n_basis_direction,
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
                   const ComponentTable<TensorSize<dim> > &n_basis_direction,
                   const Quadrature<dim> &quad);

        /**
         * Allocate space for the values and derivatives
         * at quadrature points for a specified face.
         */
        void reset(const Index face_id,
                   const BasisFaceValueFlagsHandler &flags_handler,
                   const ComponentTable<TensorSize<dim> > &n_basis_direction,
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

public:
    /**
     * @todo Document this function
     */
    const ValuesCache &get_values_cache(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;


    /*
    private:
        const typename DerivedElementAccessor::ValuesCache &
        get_values_cache(const TopologyId<dim> &topology_id) const;
    //*/
};

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
SpaceElementAccessor(const std::shared_ptr<const Space> space,
                     const int elem_index)
    :
    CartesianGridElementAccessor<dim>(space->get_grid(), elem_index),
    space_(space)
{
    //--------------------------------------------------------------------------
    using Indexer = CartesianProductIndexer<dim>;

    for (int comp_id = 0; comp_id < Space::n_components; ++comp_id)
    {
        for (int j=0; j<dim; ++j)
            n_basis_direction_(comp_id)[j] = this->space_->get_degree()(comp_id)[j]+1;

        // creating the objects for fast conversion from flat-to-tensor indexing
        // (in practice it is an hash-table from flat to tensor indices)
        basis_functions_indexer_(comp_id) = std::shared_ptr<Indexer>(new Indexer(n_basis_direction_(comp_id)));
    }
    //--------------------------------------------------------------------------



    //--------------------------------------------------------------------------
    comp_offset_(0) = 0;
    for (int comp_id = 1; comp_id < Space::n_components; ++comp_id)
        comp_offset_(comp_id)= comp_offset_(comp_id-1) + n_basis_direction_(comp_id).flat_size();
    //--------------------------------------------------------------------------
};


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
CartesianGridElementAccessor<dim> &
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
as_cartesian_grid_element_accessor()
{
    return static_cast<CartesianGridElementAccessor<dim> &>(*this);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
const CartesianGridElementAccessor<dim> &
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
as_cartesian_grid_element_accessor() const
{
    return static_cast<const CartesianGridElementAccessor<dim> &>(*this);
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
DerivedElementAccessor &
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
as_derived_element_accessor()
{
    return static_cast<DerivedElementAccessor &>(*this);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
const DerivedElementAccessor &
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
as_derived_element_accessor() const
{
    return static_cast<const DerivedElementAccessor &>(*this);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_basis_values_at_points(const std::vector<Point<dim>> &points) const -> ValueTable<Value>
{
    return this->as_derived_element_accessor().template evaluate_basis_derivatives_at_points<0>(points);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_basis_gradients_at_points(const std::vector<Point<dim>> &points) const -> ValueTable<Derivative<1> >
{
    return this->as_derived_element_accessor().template evaluate_basis_derivatives_at_points<1>(points);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_basis_hessians_at_points(const std::vector<Point<dim>> &points) const -> ValueTable<Derivative<2> >
{
    return this->as_derived_element_accessor().template evaluate_basis_derivatives_at_points<2>(points);
}



template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
template <int deriv_order>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_field_derivatives_at_points(
    const std::vector<Real> &local_coefs,
    const std::vector<Point<dim>> &points) const ->
ValueVector< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >
{
    const auto &derived_element_accessor = this->as_derived_element_accessor();
    Assert(derived_element_accessor.get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(derived_element_accessor.get_num_basis(),local_coefs.size()));

    const auto derivatives_phi_hat =
    derived_element_accessor.template evaluate_basis_derivatives_at_points<deriv_order>(points);
    Assert(derivatives_phi_hat.get_num_functions() == derived_element_accessor.get_num_basis(),
    ExcDimensionMismatch(derivatives_phi_hat.get_num_functions(), derived_element_accessor.get_num_basis())) ;

    return derivatives_phi_hat.evaluate_linear_combination(local_coefs) ;
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_field_values_at_points(
    const std::vector<Real> &local_coefs,
    const std::vector<Point<dim>> &points) const -> ValueVector<Value>
{
    return this->evaluate_field_derivatives_at_points<0>(local_coefs,points);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_field_gradients_at_points(
    const std::vector<Real> &local_coefs,
    const std::vector<Point<dim>> &points) const -> ValueVector<Derivative<1> >
{
    return this->evaluate_field_derivatives_at_points<1>(local_coefs,points);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_field_hessians_at_points(
    const std::vector<Real> &local_coefs,
    const std::vector<Point<dim>> &points) const -> ValueVector<Derivative<2> >
{
    return this->evaluate_field_derivatives_at_points<2>(local_coefs,points);
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_values(const TopologyId<dim> &topology_id) const -> ValueTable<Value> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.values_filled(), ExcCacheNotFilled());

    return cache.phi_hat_;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_face_basis_values(const Index face_id) const -> ValueTable<Value> const &
{
    return this->get_basis_values(FaceTopology<dim>(face_id));
}



template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_values(const Index i,const TopologyId<dim> &topology_id) const -> typename ValueTable<Value>::const_view
{
    return this->get_basis_values(topology_id).get_function_view(i);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_value(const Index basis, const Index qp,const TopologyId<dim> &topology_id) const -> Value const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_values(basis,topology_id)[qp];
}




template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_gradients(const TopologyId<dim> &topology_id) const -> ValueTable<Derivative<1>> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.gradients_filled(), ExcCacheNotFilled());

    return cache.D1phi_hat_;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_face_basis_gradients(const Index face_id) const -> ValueTable<Derivative<1>> const &
{
    return this->get_basis_gradients(FaceTopology<dim>(face_id));
}



template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_gradients(const Index i,const TopologyId<dim> &topology_id) const -> typename ValueTable<Derivative<1>>::const_view
{
    return this->get_basis_gradients(topology_id).get_function_view(i);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_gradient(const Index basis, const Index qp,const TopologyId<dim> &topology_id) const -> Derivative<1> const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_gradients(basis,topology_id)[qp];
}





template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_hessians(const TopologyId<dim> &topology_id) const -> ValueTable<Derivative<2>> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.hessians_filled(), ExcCacheNotFilled());

    return cache.D2phi_hat_;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_face_basis_hessians(const Index face_id) const -> ValueTable<Derivative<2>> const &
{
    return this->get_basis_hessians(FaceTopology<dim>(face_id));
}



template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_hessians(const Index i,const TopologyId<dim> &topology_id) const -> typename ValueTable<Derivative<2>>::const_view
{
    return this->get_basis_hessians(topology_id).get_function_view(i);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_hessian(const Index basis, const Index qp,const TopologyId<dim> &topology_id) const -> Derivative<2> const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_hessians(basis,topology_id)[qp];
}




template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_divergences(const TopologyId<dim> &topology_id) const -> ValueTable<Div> const &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_filled() == true, ExcCacheNotFilled());
    Assert(cache.flags_handler_.divergences_filled(), ExcCacheNotFilled());

    return cache.div_phi_hat_;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_face_basis_divergences(const Index face_id) const -> ValueTable<Div> const &
{
    return this->get_basis_divergences(FaceTopology<dim>(face_id));
}



template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_divergences(const Index i,const TopologyId<dim> &topology_id) const -> typename ValueTable<Div>::const_view
{
    return this->get_basis_divergences(topology_id).get_function_view(i);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_basis_divergence(const Index basis, const Index qp,const TopologyId<dim> &topology_id) const -> Div const &
{
    Assert(qp >= 0 && qp < this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size(),
           ExcIndexRange(qp,0,this->get_values_cache(topology_id).quad_.get_num_points_direction().flat_size()));
    return this->get_basis_divergences(basis,topology_id)[qp];
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_values_cache(const TopologyId<dim> &topology_id) const -> const ValuesCache &
{
    Assert(topology_id.is_element() || topology_id.is_face(),
           ExcMessage("Only element or face topology is allowed."));
    if (topology_id.is_element())
    {
        return elem_values_;
    }
    else
    {
        Assert(topology_id.get_id()>=0 && topology_id.get_id() < n_faces,
               ExcIndexRange(topology_id.get_id(),0,n_faces));
        return face_values_[topology_id.get_id()];
    }
}





template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
void
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
ValuesCache::
reset(const BasisElemValueFlagsHandler &flags_handler,
      const ComponentTable<TensorSize<dim> > &n_basis_direction,
      const Quadrature<dim> &quad)
{
    quad_ = quad;
    const TensorSize<dim> n_points_direction = quad_.get_num_points_direction();

    flags_handler_ = flags_handler;

    //--------------------------------------------------------------------------
    // computing the total number of basis functions and the number of evaluation points
    const int total_n_points = n_points_direction.flat_size();

    int total_n_basis = 0;
    for (int i = 0; i < Space::n_components; ++i)
        total_n_basis += n_basis_direction(i).flat_size();

    Assert(total_n_points > 0, ExcLowerRange(total_n_points,1));
    Assert(total_n_basis > 0, ExcLowerRange(total_n_basis,1));
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // resizing the containers for the basis functions
    if (flags_handler_.fill_values())
    {
        if (phi_hat_.get_num_points() != total_n_points ||
            phi_hat_.get_num_functions() != total_n_basis)
        {
            phi_hat_.resize(total_n_basis,total_n_points);
            phi_hat_.zero();
        }
    }
    else
    {
        phi_hat_.clear();
    }

    if (flags_handler_.fill_gradients())
    {
        if (D1phi_hat_.get_num_points() != total_n_points ||
            D1phi_hat_.get_num_functions() != total_n_basis)
        {
            D1phi_hat_.resize(total_n_basis,total_n_points);
            D1phi_hat_.zero();
        }
    }
    else
    {
        D1phi_hat_.clear();
    }


    if (flags_handler_.fill_divergences())
    {
        Assert(flags_handler_.fill_gradients(),
               ExcMessage("Divergence requires gradient to be filled."));

        if (div_phi_hat_.get_num_points() != total_n_points ||
            div_phi_hat_.get_num_functions() != total_n_basis)
        {
            div_phi_hat_.resize(total_n_basis,total_n_points);
            div_phi_hat_.zero();
        }
    }
    else
    {
        div_phi_hat_.clear();
    }



    if (flags_handler_.fill_hessians())
    {
        if (D2phi_hat_.get_num_points() != total_n_points ||
            D2phi_hat_.get_num_functions() != total_n_basis)
        {
            D2phi_hat_.resize(total_n_basis,total_n_points);
            D2phi_hat_.zero();
        }
    }
    else
    {
        D2phi_hat_.clear();
    }
    //--------------------------------------------------------------------------


    this->set_initialized(true);
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
void
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
ElementValuesCache::
reset(const BasisElemValueFlagsHandler &flags_handler,
      const ComponentTable<TensorSize<dim> > &n_basis_direction,
      const Quadrature<dim> &quad)
{
    ValuesCache::reset(flags_handler, n_basis_direction,quad);
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
void
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
FaceValuesCache::
reset(const Index face_id,
      const BasisFaceValueFlagsHandler &flags_handler,
      const ComponentTable<TensorSize<dim> > &n_basis_direction,
      const Quadrature<dim> &quad1)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    const auto quad = quad1.collapse_to_face(face_id);

    ValuesCache::reset(flags_handler, n_basis_direction,quad);
}



template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
void
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
FaceValuesCache::
reset(const Index face_id,
      const BasisFaceValueFlagsHandler &flags_handler,
      const ComponentTable<TensorSize<dim> > &n_basis_direction,
      const Quadrature<dim-1> &quad)
{
    Assert(false,ExcNotImplemented()) ;
    AssertThrow(false,ExcNotImplemented()) ;
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
ValuesCache::
get_values() const -> const ValueTable<Value> &
{
    return phi_hat_;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
ValuesCache::
get_gradients() const -> const ValueTable<Derivative<1>> &
{
    return D1phi_hat_;
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
ValuesCache::
get_hessians() const -> const ValueTable<Derivative<2>> &
{
    return D2phi_hat_;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
ValuesCache::
get_divergences() const -> const ValueTable<Div> &
{
    return div_phi_hat_;
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_field(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const
-> ValueVector<Value>
{
    Assert(this->get_values_cache(topology_id).is_filled() == true, ExcCacheNotFilled());
    Assert(this->get_values_cache(topology_id).flags_handler_.fill_values() == true, ExcCacheNotFilled());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D0phi_hat = this->get_basis_values(topology_id) ;
    Assert(D0phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D0phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D0phi_hat.evaluate_linear_combination(local_coefs) ;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_face_field(const Index face_id, const std::vector<Real> &local_coefs) const
-> ValueVector<Value>
{
    return this->evaluate_field(local_coefs,FaceTopology<dim>(face_id));
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_field_gradients(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const
-> ValueVector< Derivative<1> >
{
    Assert(this->get_values_cache(topology_id).is_filled() == true, ExcCacheNotFilled());
    Assert(this->get_values_cache(topology_id).flags_handler_.fill_gradients() == true, ExcCacheNotFilled());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D1phi_hat = this->get_basis_gradients(topology_id) ;
    Assert(D1phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D1phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D1phi_hat.evaluate_linear_combination(local_coefs) ;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_face_field_gradients(const Index face_id, const std::vector<Real> &local_coefs) const
-> ValueVector< Derivative<1> >
{
    return this->evaluate_field_gradients(local_coefs,FaceTopology<dim>(face_id));
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_field_divergences(
    const std::vector<Real> &local_coefs,
    const TopologyId<dim> &topology_id) const -> ValueVector<Div>
{
    Assert(this->get_values_cache(topology_id).is_filled() == true, ExcCacheNotFilled());
    Assert(this->get_values_cache(topology_id).flags_handler_.fill_divergences() == true, ExcCacheNotFilled());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &div_phi_hat = this->get_basis_divergences(topology_id) ;
    Assert(div_phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(div_phi_hat.get_num_functions(), this->get_num_basis())) ;

    return div_phi_hat.evaluate_linear_combination(local_coefs) ;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_face_field_divergences(const Index face_id, const std::vector<Real> &local_coefs) const
-> ValueVector<Div>
{
    return this->evaluate_field_divergences(local_coefs,FaceTopology<dim>(face_id));
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_field_hessians(const std::vector<Real> &local_coefs,const TopologyId<dim> &topology_id) const -> ValueVector< Derivative<2> >
{
    Assert(this->get_values_cache(topology_id).is_filled() == true, ExcCacheNotFilled());
    Assert(this->get_values_cache(topology_id).flags_handler_.fill_hessians() == true, ExcCacheNotFilled());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D2phi_hat = this->get_basis_hessians(topology_id) ;
    Assert(D2phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D2phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D2phi_hat.evaluate_linear_combination(local_coefs) ;
}

template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
evaluate_face_field_hessians(const Index face_id, const std::vector<Real> &local_coefs) const
-> ValueVector< Derivative<2> >
{
    return this->evaluate_field_hessians(local_coefs,FaceTopology<dim>(face_id));
}





template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
Size
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_num_basis() const
{
    return this->space_->get_num_basis_per_element();
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
int
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_num_basis(const int i) const
{
    const auto &degree_comp = this->space_->get_degree()(i);
    int component_num_basis = 1;
    for (int j = 0; j < dim; ++j)
        component_num_basis *= degree_comp[j] + 1;

    return component_num_basis;
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_local_to_global() const -> std::vector<Index> const &
{
    return space_->get_element_global_dofs()[this->get_flat_index()];
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_space() const -> std::shared_ptr<const Space>
{
    return space_;
}


template<class DerivedElementAccessor,class Space,int dim,int codim,int range,int rank>
inline
auto
SpaceElementAccessor<DerivedElementAccessor,Space,dim,codim,range,rank>::
get_quad_points(const TopologyId<dim> &topology_id) const -> const Quadrature<dim> &
{
    const auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_initialized(), ExcNotInitialized());

    return cache.quad_;
}



/**
 * See module on \ref accessors_iterators for a general overview.
 * @ingroup accessors_iterators
 */
template <int dim, int range, int rank>
class BSplineElementAccessor :
    public SpaceElementAccessor<
    BSplineElementAccessor<dim,range,rank>,BSplineSpace<dim,range,rank>,dim,0,range,rank>
{
public:
    using parent_t = SpaceElementAccessor<
                     BSplineElementAccessor<dim,range,rank>,BSplineSpace<dim, range, rank>,dim,0,range,rank>;

    /** Type for the grid accessor. */
    using GridAccessor = CartesianGridElementAccessor<dim>;

    /** Type required by the GridForwardIterator templated iterator */
    using ContainerType = const BSplineSpace<dim, range, rank> ;

    /** Type required for the generic algorithm on the spaces (plots??) */
    using Space = BSplineSpace<dim, range, rank> ;

    /** Number of faces of the element. */
    using parent_t::n_faces;


    using ValuesCache = typename parent_t::ValuesCache;


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


public:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    BSplineElementAccessor() = delete;

    /**
     * Constructs an accessor to element number index of a
     * BsplineSpace space.
     */
    BSplineElementAccessor(const std::shared_ptr<ContainerType> space,
                           const int elem_index);

    /**
     * Copy constructor.
     * @note For the constructed object it
     * creates a new element cache, but it shares
     * the one dimensional cache with the copied element.
     */
    BSplineElementAccessor(const BSplineElementAccessor<dim, range, rank> &elem)
        = default;

    /**
     * Move constructor.
     */
    BSplineElementAccessor(BSplineElementAccessor<dim, range, rank> &&elem)
        = default;

    /**
     * Destructor.
     */
    ~BSplineElementAccessor() = default;
    ///@}


    /** @name Assignment operators */
    ///@{

    /**
     * Copy assignment operator.
     * @note Creates a new element cache, but it shares
     * the one dimensional cache with the copied element.
     */
    BSplineElementAccessor<dim, range, rank> &
    operator=(const BSplineElementAccessor<dim, range, rank> &elem)
        = default;

    /**
     * Move assignment operator.
     */
    BSplineElementAccessor<dim, range, rank> &
    operator=(BSplineElementAccessor<dim, range, rank> &&elem)
        = default;
    ///@}



    /** @name Query information that requires the use of the cache */
    ///@{

    /**
     * Initializes the internal cache for the efficient
     * computation of the values requested in
     * the fill_flag on the given quadrature points.
     * This implies a uniform quadrature scheme
     * (i.e. the same for all elements).
     * @note This function should be called before fill_values()
     */
    void init_values(const ValueFlags fill_flag,
                     const Quadrature<dim> &quad);

    /**
     * For a given face quadrature.
     */
    void init_face_values(const Index face_id,
                          const ValueFlags fill_flag,
                          const Quadrature<dim-1> &quad);

    /**
     * Fills the element values cache according to the evaluation points
     * and fill flags specifies in init_values.
     */
    void fill_values();

    void fill_face_values(const Index face_id);

    /** Reset the global cache */
    void reset_global_cache();

    /**
     * Typedef for specifying the value of the basis function in the
     * reference domain.
     */
    using Value = Values<dim, range, rank>;

    /**
     * Typedef for specifying the derivatives of the basis function in the
     * reference domain.
     */
    template <int deriv_order>
    using Derivative = Derivatives<dim, range, rank, deriv_order>;

    /**
     * Typedef for specifying the divergence of the basis function in the
     * reference domain.
     */
    using Div = Values<dim, 1, 1>;

protected:

public:

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
     * Prints internal information about the BSplineElementAccessor.
     * Its main use is for testing and debugging.
     */
    void print_info(LogStream &out, const VerbosityLevel verbosity_level = VerbosityLevel::normal) const;

private:
    /**
     * Initilizes (reserves memory) the
     * univariate basis values
     * and derivatives at quadrature points cache
     * for efficient use of computations.
     * This function implies a uniform quadrature schemes
     * (the same for each element).
     * The fill_flag provides what information to compute.
     */
    void reset_univariate_cache(const Quadrature<dim> &quad,
                                const int max_der);

    /**
     * Initializes the element cache according to
     * the quadrature number of point and the fill_flag.
     */
    void reset_element_cache(const ValueFlags fill_flag,
                             const Quadrature<dim> &quad);



    /**
     * @name Containers for the cache of the element values and for the
     * cache of the face values
     */
    ///@{

    /**
     * This type store the values, first second derivatives
     * of a 1D Bspline functions, i.e BasisValues1d[k]
     * stores the values of the k-th derivative of the (p+1) basis function on a given interval
     * at the quadrature points.
     * BasisValues1d[k] is a (p+1) x n_qp matrix
     */
    using BasisValues1d = std::vector<DenseMatrix>;

protected:

    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentTable = StaticMultiArray<T,range,rank>;

    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentDirectionTable =
        StaticMultiArray<CartesianProductArray<T,dim>, range, rank>;




private:

    ComponentTable<
    DynamicMultiArray<std::shared_ptr<BSplineElementScalarEvaluator<dim>>,dim>> scalar_evaluators_;


    using univariate_values_t = ComponentTable<std::array<const BasisValues1d *,dim>>;

    /**
     * Fills the cache (accordingly with the flags_handler status)
     * from the univariate values (and derivatives up to the order
     * specified by @p max_deriv_order).
     *
     *
     * @note The BSplineElementAccessor @p elem is needed in order to call the function
     * elem.evaluate_bspline_derivatives<p>()
     * that computes the @p p-th order derivatives of a BSpline from the univariate values.
     */
    void fill_values_cache_from_univariate(const int max_deriv_order,
                                           const univariate_values_t &values_1D,
                                           ValuesCache &cache);




    ///@}

public:
    /**
     * For a given flags input argument identifies the face quantities and
     * returns a new ValueFlags variable containing only face quantities.
     * The output flags does not contain the word face.
     */
    ValueFlags get_face_flags(const ValueFlags fill_flag) const ;


private:

    /**
     * Computes the k-th order derivative of the non-zero B-spline basis
     * functions over the current element,
     *   at the evaluation points pre-allocated in the cache.
     *
     * \warning If the output result @p derivatives_phi_hat is not correctly pre-allocated,
     * an exception will be raised.
     */
    template <int deriv_order>
    void evaluate_bspline_derivatives(const ComponentTable<std::array<const BasisValues1d *, dim> > &elem_values,
                                      const ValuesCache &cache,
                                      ValueTable<
                                      Conditional<(deriv_order==0),Value,Derivative<deriv_order> >
                                      > &derivatives_phi_hat) const;





    class GlobalCache : public CacheStatus
    {
    public:
        int max_deriv_order_ = 0;
    };

    /**
     * Cache for the efficient use of Bspline basis on a uniform
     * quadrature scheme on all elements.
     */
    class GlobalElemCache : public GlobalCache
    {
    public:
        /**
         * Allocates space for and compute the values and
         * derivatives at quadrature points for one
         * dimensional B-splines for each component
         * and each direction of a component and each interval
         * of a direction.
         */
        void reset(const Space &space,
                   const Quadrature<dim> &quad,
                   const int max_der);

        /**
         * univariate B-splines values and derivatives at
         * quadrature points
         * splines1d_cache_data_[comp][dir][interval][order][function][point]
         */
        ComponentDirectionTable<BasisValues1d> splines1d_cache_data_;

        ComponentDirectionTable<const BasisValues1d *> splines1d_cache_;

    };

    class GlobalFaceCache : public GlobalCache
    {
    public:
        /**
         * Allocates space for and compute the values and
         * derivatives at quadrature points for one
         * dimensional B-splines for each component
         * and each direction of a component and each interval
         * of a direction.
         */
        void reset(const Space &space,
                   const Quadrature<dim> &quad1,
                   const Index face_id,
                   const int max_der);

        /**
         * univariate B-splines values and derivatives at
         * quadrature points
         * splines1d_cache_data_[comp][interval][order][function][point]
         */
        ComponentTable<BasisValues1d> splines1d_cache_data_;

        ComponentTable<const BasisValues1d *> splines1d_cache_;

    };

    /**
     * Tensor product style space sized cache for
     * storing the values and derivatives of the
     * one dimensional values of the basis
     * function at the quadrature points
     */
    std::shared_ptr< GlobalElemCache > values_1d_elem_ = nullptr;

    std::array<std::shared_ptr<GlobalFaceCache>, n_faces> values_1d_faces_;



protected:


    /** Returns the Bezier extraction operator relative to the current element. */
    ComponentTable< std::array< const DenseMatrix *,dim> > get_bezier_extraction_operator() const;


private:


    template <typename Accessor> friend class GridForwardIterator;




public:
    const ComponentTable<
    DynamicMultiArray<
    std::shared_ptr<
    BSplineElementScalarEvaluator<dim>>,dim> > &get_scalar_evaluators() const;

};




IGA_NAMESPACE_CLOSE

#endif
