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


#ifndef SPACE_ELEMENT_H_
#define SPACE_ELEMENT_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/flags_handler.h>

#include <igatools/base/quadrature.h>
#include <igatools/geometry/cartesian_grid_element.h>

#include <igatools/utils/value_vector.h>
#include <igatools/utils/value_table.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/cartesian_product_indexer.h>

IGA_NAMESPACE_OPEN

template <class Accessor, class Allocator> class CartesianGridIterator;




template<class Space>
class SpaceElement : public CartesianGridElement<Space::dim>
{
protected:
    using base_t =  CartesianGridElement<Space::dim>;
private:
    using self_t = SpaceElement<Space>;

public:
    using DerivedElementAccessor = typename Space::ElementAccessor;

    using RefPoint = typename Space::RefPoint;
    using Point = typename Space::Point;
    using Value = typename Space::Value;
    template <int order>
    using Derivative = typename Space::template Derivative<order>;
    using Div = typename Space::Div;

    static const int dim       = Space::dim;
    static const int codim     = Space::codim;
    static const int space_dim = Space::space_dim;
    static const int range     = Space::range;
    static const int rank      = Space::rank;

    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentContainer = typename Space::template ComponentContainer<T>;
    using SpaceDimensionTable = typename Space::SpaceDimensionTable;
    ///@}


    /** @name Constructors */
    ///@{
    /**
     * Default constructor.
     */
    SpaceElement()
    {
        Assert(false,ExcNotImplemented());
    }

    /**
     * Constructs an accessor to element number index of a
     * function space.
     */
    SpaceElement(const std::shared_ptr<const Space> space,
                 const Index elem_index);

    /**
     * Constructs an accessor to element number index of a
     * function space.
     */
    SpaceElement(const std::shared_ptr<const Space> space,
                 const TensorIndex<dim> &elem_index);

    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a classic copy constructor)
     * uses the deep copy.
     */
    SpaceElement(const self_t &elem,
                 const CopyPolicy &copy_policy = CopyPolicy::deep);

    /**
     * Move constructor.
     */
    SpaceElement(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~SpaceElement() = default;
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

    template<int order = 0, int k = dim>
    auto
    get_values(const int j = 0) const
    {
        Assert(local_cache_ != nullptr, ExcNullPtr());
        const auto &cache = local_cache_->template get_value_cache<k>(j);
        Assert(cache.is_filled() == true, ExcCacheNotFilled());
        return cache.template get_der<order>();
    }

    template<int order, int k>
    auto
    linear_combination(const vector<Real> &loc_coefs, const int id) const
    {
        const auto &basis_values =
            this->template get_values<order, k>(id);
        return basis_values.evaluate_linear_combination(loc_coefs) ;
    }

#if 0
    /** @name Functions returning the value of the basis functions. */
    ///@{
    /**
     * Returns the const reference to a ValueTable with the values of all local basis function
     * at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    ValueTable<Value> const &get_basis_values() const;

    /**
     * Returns a const view to the values of the <tt>i</tt>-th basis function at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    typename ValueTable<Value>::const_view
    get_basis_values(const Index i) const;

    /**
     * Returns the const reference to the value of a local basis function
     * at one evaluation point.
     * @param[in] basis Local basis id.
     * @param[in] qp Point id.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    Value const &get_basis_value(const Index basis, const Index qp);

    /**
     * Returns the const reference to a ValueTable with the values of all local basis function
     * at each evaluation point on the face specified by @p face_id.
     */
    ValueTable<Value> const &get_face_basis_values(const Index face_id) const;
    ///@}



    /** @name Functions returning the gradient of the basis functions. */
    ///@{
    /**
     * Returns the const reference to a ValueTable with the gradients of all local basis function
     * evaluated at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    ValueTable<Derivative<1>> const &get_basis_gradients() const;

    /**
     * Returns a const view to the gradients of the <tt>i</tt>-th basis function at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    typename ValueTable<Derivative<1> >::const_view
    get_basis_gradients(const Index i,) const;

    /**
     * Returns the const reference to the gradient of a local basis function
     * at one evaluation point.
     * @param[in] basis Local basis id.
     * @param[in] qp Point id.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    Derivative<1> const &get_basis_gradient(const Index basis, const Index qp)

    /**
     * Returns the const reference to a ValueTable with the gradients of all local basis function
     * at each evaluation point on the face specified by @p face_id.
     */
    ValueTable<Derivative<1>> const
                           &get_face_basis_gradients(const Index face_id) const;
    ///@}

    /** @name Functions returning the hessian of the basis functions. */
    ///@{
    /**
     * Returns the const reference to a ValueTable with hessians of all local basis function
     * at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    ValueTable<Derivative<2>> const &get_basis_hessians() const;

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
    Derivative<2> const &get_basis_hessian(const Index basis, const Index qp) const;

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
    ValueTable<Div> const &get_basis_divergences() const;

    /**
     * Returns a const view to the divergences of the <tt>i</tt>-th basis function at each evaluation point.
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     */
    typename ValueTable<Div>::const_view
    get_basis_divergences(const Index i) const;

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
    evaluate_field(const vector<Real> &local_coefs) const;


    /**
     * Returns the ValueVector with the evaluation of the field @p local_coefs at the evaluation
     * points on the face specified by @p face_id.
     */
    ValueVector<Value>
    evaluate_face_field(const Index face_id, const vector<Real> &local_coefs) const;

    /**
     * Returns the ValueVector with the evaluation of the gradient of the field @p local_coefs
     * at the evaluation points.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     * @see get_local_coefs
     */
    ValueVector<Derivative<1> >
    evaluate_field_gradients(const vector<Real> &local_coefs) const;

    /**
     * Returns the ValueVector with the evaluation of the gradient of the field @p local_coefs at the evaluation
     * points on the face specified by @p face_id.
     */
    ValueVector<Derivative<1> >
    evaluate_face_field_gradients(const Index face_id,
                                  const vector<Real> &local_coefs) const;

    /**
     * Returns the ValueVector with the evaluation of the hessians of the field @p local_coefs
     * at the evaluation points.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     * @see get_local_coefs
     */
    ValueVector<Derivative<2> >
    evaluate_field_hessians(const vector<Real> &local_coefs) const;

    /**
     * Returns the ValueVector with the evaluation of the hessian of the field @p local_coefs at the evaluation
     * points on the face specified by @p face_id.
     */
    ValueVector<Derivative<2> >
    evaluate_face_field_hessians(const Index face_id,
                                 const vector<Real> &local_coefs) const;

    /**
     * Returns the ValueVector with the evaluation of the divergences of the field @p local_coefs
     * at the evaluation points.
     *
     * @note The @p topology_id parameter can be used to select values on the element
     * (it's the default behaviour if @p topology_id is not specified) or on a element-face. See the TopologyId documentation).
     * @see get_local_coefs
     */
    ValueVector<Div>
    evaluate_field_divergences(const vector<Real> &local_coefs) const;

    /**
     * Returns the ValueVector with the evaluation of the divergence of the field @p local_coefs at the evaluation
     * points on the face specified by @p face_id.
     */
    ValueVector<Div>
    evaluate_face_field_divergences(const Index face_id,
                                    const vector<Real> &local_coefs) const;
    ///@}

#endif

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
     * Returns the basis function ID offset between the different components.
     */
    ComponentContainer<int> get_basis_offset() const;

    /**
     * Returns the global dofs of the local (non zero) basis functions
     * on the element.
     * For example:
     * \code
       auto loc_to_glob = elem->get_local_to_global();
       // loc_to_glob[0] is the global id of the first basis function on the element
       // loc_to_glob[1] is the global id of the second basis function on the element
       // ...
      \endcode
     *
     */
    vector<Index> get_local_to_global() const;

    /**
     * Returns the patch dofs of the local (non zero) basis functions
     * on the element.
     */
    vector<Index> get_local_to_patch() const;

    /**
     * Pointer to the @p Space upon which the accessor is iterating on.
     */
    std::shared_ptr<const Space> get_space() const;
    ///@}


    void print_info(LogStream &out) const;

    void print_cache_info(LogStream &out) const;


private:
    /**
     * Space for which the SpaceElement refers to.
     */
    std::shared_ptr<const Space> space_ = nullptr;


protected:
    /** Number of scalar basis functions along each direction, for all space components. */
    typename Space::SpaceDimensionTable n_basis_direction_;


    /** Hash table for fast conversion between flat-to-tensor basis function ids. */
    ComponentContainer<std::shared_ptr<CartesianProductIndexer<dim> > > basis_functions_indexer_;

    /** Basis function ID offset between the different components. */
    ComponentContainer<int> comp_offset_;

    /**
     * Base class for the cache of the element values and
     * for the cache of the face values.
     */
    class ValuesCache : public CacheStatus
    {
    public:
        /**
         * Allocate space for the values and derivatives
         * of the element basis functions at quadrature points
         * as specify by the flag
         */
        void resize(const FunctionFlags &flags_handler,
                    const Size total_n_points,
                    const Size n_basis);


        /** Returns the divergences. */
        //const ValueTable<Div> &get_divergences() const;


    public:
        void print_info(LogStream &out) const;

        FunctionFlags flags_handler_;

        std::tuple<ValueTable<Value>,
            ValueTable<Derivative<1>>,
            ValueTable<Derivative<2>>> values_;

        template<int k>
        auto &get_der()
        {
            return std::get<k>(values_);
        }

        template<int k>
        const auto &get_der() const
        {
            return std::get<k>(values_);
        }

        template<int k>
        void resize_der(const int n_basis, const int n_points)
        {
            auto &value = std::get<k>(values_);
            if (value.get_num_points() != n_points ||
                value.get_num_functions() != n_basis)
            {
                value.resize(n_basis, n_points);
                value.zero();
            }
        }

        template<int k>
        void clear_der()
        {
            auto &value = std::get<k>(values_);
            value.clear();
        }

    };


    class LocalCache
    {
    public:
        LocalCache() = default;

        LocalCache(const LocalCache &in) = default;

        LocalCache(LocalCache &&in) = default;

        ~LocalCache() = default;


        LocalCache &operator=(const LocalCache &in) = delete;

        LocalCache &operator=(LocalCache &&in) = delete;

        void print_info(LogStream &out) const;

        template <int k>
        ValuesCache &
        get_value_cache(const int j)
        {
            return std::get<k>(values_)[j];
        }

        template <int k>
        const ValuesCache &
        get_value_cache(const int j) const
        {
            return std::get<k>(values_)[j];
        }

        CacheList<ValuesCache, dim> values_;
    };

    /** The local (element and face) cache. */
    std::shared_ptr<LocalCache> local_cache_;

public:
    std::shared_ptr<LocalCache> &get_local_cache()
    {
        return local_cache_;
    }

protected:
    /**
     * Performs a copy of the input @p element.
     * The type of copy (deep or shallow) is specified by the input parameter @p copy_policy.
     */
    void copy_from(const self_t &element,
                   const CopyPolicy &copy_policy);



public:

    /** Return a reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    CartesianGridElement<dim> &as_cartesian_grid_element_accessor();


    /** Return a const-reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    const CartesianGridElement<dim> &as_cartesian_grid_element_accessor() const;

private:
#if 0
    /** Return a reference to "*this" as being an object of type DerivedElementAccessor.*/
    DerivedElementAccessor &as_derived_element_accessor();

    /** Return a const-reference to "*this" as being an object of type DerivedElementAccessor.*/
    const DerivedElementAccessor &as_derived_element_accessor() const;
#endif
};


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space_element-inline.h>


#endif // #ifndef SPACE_ELEMENT_ACCESSOR_

