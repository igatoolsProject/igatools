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

template <class Grid> class FunctionSpaceOnGrid;

template <class Accessor> class CartesianGridIterator;

template <int dim>
class SpaceElementBase : private CartesianGridElement<dim>
{
protected:
    using base_t = CartesianGridElement<dim>;
private:
    using self_t = SpaceElementBase<dim>;

    using Space = FunctionSpaceOnGrid<CartesianGrid<dim>>;
public:
    using base_t::get_flat_index;
    using base_t::get_tensor_index;
    using base_t::get_grid;
    using base_t::is_boundary;

    using Topology = typename base_t::Topology;

    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    SpaceElementBase() = delete;

    /**
     * Constructs an accessor to element number index of a
     * function space.
     */
    SpaceElementBase(const std::shared_ptr<const Space> space,
                 const Index elem_index)
    :
    	base_t(space->get_grid(), elem_index),
		space0_(space)
    {
    	Assert(space0_ != nullptr, ExcNullPtr());
    }



    /**
     * Copy constructor.
     * It can be used with different copy policies (i.e. deep copy or shallow copy).
     * The default behaviour (i.e. using the proper interface of a classic copy constructor)
     * uses the deep copy.
     */
    SpaceElementBase(const self_t &elem,
                 const CopyPolicy &copy_policy = CopyPolicy::deep)
    :
    	base_t(elem,copy_policy),
		space0_(elem.space0_)
    {}

    /**
     * Move constructor.
     */
    SpaceElementBase(self_t &&elem) = default;

    /**
     * Destructor.
     */
    ~SpaceElementBase() = default;
    ///@}


    /** Return a reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    CartesianGridElement<dim> &as_cartesian_grid_element_accessor()
	{
	    return static_cast<CartesianGridElement<dim> &>(*this);
	}



    /** Return a const-reference to "*this" as being an object of type CartesianGridElementAccessor.*/
    const CartesianGridElement<dim> &as_cartesian_grid_element_accessor() const
    {
        return static_cast<const CartesianGridElement<dim> &>(*this);
    }

    void print_info(LogStream &out) const
    {
    	base_t::print_info(out);
    }

    void print_cache_info(LogStream &out) const
    {
        base_t::print_cache_info(out);
    }

    void
    copy_from(const SpaceElementBase<dim> &elem,const CopyPolicy &copy_policy)
    {
        if (this != &elem)
        {
        	CartesianGridElement<dim>::copy_from(elem,copy_policy);
        	space0_ = elem.space0_;
        }
    }


private:
    std::shared_ptr<const Space> space0_;


public:

    /**
     * @name Comparison operators
     *
     * @warning To be comparable, two SpaceElement objects must be defined on the same space,
     * otherwise an assertion will be raised (in Debug mode).
     */
    ///@{
    bool operator==(const self_t &a) const
    {
        Assert(space0_ == a.space0_,
               ExcMessage("Comparison between elements defined on different spaces"));
        return this->as_cartesian_grid_element_accessor() == a.as_cartesian_grid_element_accessor();
    }


    bool operator!=(const self_t &a) const
    {
        Assert(space0_ == a.space0_,
               ExcMessage("Comparison between elements defined on different spaces"));
        return this->as_cartesian_grid_element_accessor() != a.as_cartesian_grid_element_accessor();
    }

    bool operator<(const self_t &a) const
    {
        Assert(space0_ == a.space0_,
               ExcMessage("Comparison between elements defined on different spaces"));
        return this->as_cartesian_grid_element_accessor() < a.as_cartesian_grid_element_accessor();
    }

    bool operator>(const self_t &a) const
    {
        Assert(space0_ == a.space0_,
               ExcMessage("Comparison between elements defined on different spaces"));
        return this->as_cartesian_grid_element_accessor() > a.as_cartesian_grid_element_accessor();
    }
    ///@}
};


template<class Space>
class SpaceElement : public SpaceElementBase<Space::dim>
{
protected:
    using base_t =  SpaceElementBase<Space::dim>;
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

    /*
    using base_t::get_flat_index;
    using base_t::get_tensor_index;
    using base_t::get_grid;
    using base_t::is_boundary;
//*/
    using Topology = typename base_t::Topology;

    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentContainer = typename Space::template ComponentContainer<T>;
    using TensorSizeTable = typename Space::TensorSizeTable;
    ///@}


    /** @name Constructors */
    ///@{
    /**
     * Default constructor. Not allowed to be used.
     */
    SpaceElement() = delete;

    /**
     * Constructs an accessor to element number index of a
     * function space.
     */
    SpaceElement(const std::shared_ptr<const Space> space,
                 const Index elem_index);
#if 0
    /**
     * Constructs an accessor to element number index of a
     * function space.
     */
    SpaceElement(const std::shared_ptr<const Space> space,
                 const TensorIndex<dim> &elem_index);
#endif
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
    get_values(const int j, const std::string &dofs_property = DofProperties::active) const
    {
        Assert(local_cache_ != nullptr, ExcNullPtr());
        const auto &cache = local_cache_->template get_value_cache<k>(j);
        Assert(cache.is_filled() == true, ExcCacheNotFilled());
        const auto values_all_elem_dofs = cache.template get_der<order>();

        //--------------------------------------------------------------------------------------
        // filtering the values that correspond to the dofs with the given property --- begin
        vector<Index> dofs_global;
        vector<Index> dofs_local_to_patch;
        vector<Index> dofs_local_to_elem;

        this->space_->get_element_dofs(
            this->as_cartesian_grid_element_accessor(),
            dofs_global,
            dofs_local_to_patch,
            dofs_local_to_elem,
            dofs_property);

        const auto n_active_dofs = dofs_local_to_elem.size();
        const auto n_pts = values_all_elem_dofs.get_num_points();

        decltype(values_all_elem_dofs) values_active_elem_dofs(n_active_dofs,n_pts);

        int fn = 0;
        for (const auto loc_dof : dofs_local_to_elem)
        {
            const auto values_all_elem_dofs_fn = values_all_elem_dofs.get_function_view(loc_dof);

            const auto values_active_elem_dofs_fn = values_active_elem_dofs.get_function_view(fn);

            std::copy(values_all_elem_dofs_fn.begin(),
                      values_all_elem_dofs_fn.end(),
                      values_active_elem_dofs_fn.begin());

            ++fn;
        }
        // filtering the values that correspond to the dofs with the given property --- end
        //--------------------------------------------------------------------------------------

        return values_active_elem_dofs;
    }

    auto
    get_element_values(const std::string &dofs_property = DofProperties::active) const
    {
        return this->template get_values<0,dim>(0,dofs_property);
    }

    template<int order, int k>
    auto
    linear_combination(const vector<Real> &loc_coefs,
                       const int id,
                       const std::string &dofs_property) const
    {
        const auto &basis_values =
            this->template get_values<order, k>(id,dofs_property);
        return basis_values.evaluate_linear_combination(loc_coefs) ;
    }


    template<int k = dim>
    ValueTable<Div> get_divergences(const int id,
                                    const std::string &dofs_property) const
    {
        /*
        Assert(local_cache_ != nullptr, ExcNullPtr());
        const auto &cache = local_cache_->template get_value_cache<k>(id);
        Assert(cache.is_filled() == true, ExcCacheNotFilled());

        Assert(cache.flags_handler_.gradients_filled() == true, ExcCacheNotFilled());
        //*/
        const auto &basis_gradients =
            this->template get_values<1,k>(id,dofs_property);

        const int n_basis = basis_gradients.get_num_functions();
        const int n_pts   = basis_gradients.get_num_points();

        ValueTable<Div> div(n_basis,n_pts);

        auto div_it = div.begin();
        for (const auto &grad : basis_gradients)
        {
            *div_it = trace(grad);
            ++div_it;
        }

        return div;
    }



    ValueTable<Div> get_element_divergences(const std::string &dofs_property) const
    {
        return get_divergences<dim>(0,dofs_property);
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
     *  Number of non zero basis functions with the given @p dofs_property,
     *  over the current element.
     */
    Size get_num_basis(const std::string &dofs_property) const;

    virtual Size get_num_basis() const = 0;


    /**
     * Returns the global dofs of the local (non zero) basis functions
     * on the element.
     *
     * @note The dofs can be filtered invoking the function with the argument @p dof_property.
     * If @p dof_property is equal to DofProperties::active, then no filter is applied.
     *
     * For example:
     * \code
       auto loc_to_glob_all = elem->get_local_to_global(DofProperties::active);
       // loc_to_glob_all[0] is the global id of the first basis function on the element
       // loc_to_glob_all[1] is the global id of the second basis function on the element
       // ...
       auto loc_to_glob_active = elem->get_local_to_global(DofProperties::active);
       // loc_to_glob_active[0] is the global id of the first active basis function on the element
       // loc_to_glob_active[1] is the global id of the second active basis function on the element
       // ...
      \endcode
     *
     */
    vector<Index> get_local_to_global(const std::string &dofs_property = DofProperties::active) const;

    /**
     * Returns the patch dofs of the local (non zero) basis functions
     * on the element.
     *
     * @note The dofs can be filtered invoking the function with the argument @p dof_property.
     * If @p dof_property is equal to DofProperties::active, then no filter is applied.
     *
     */
    vector<Index> get_local_to_patch(const std::string &dofs_property) const;

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
#ifndef NDEBUG
            // TODO (pauletti, Mar 17, 2015): bad checking, should be k independent
            if (k == 0)
            {
                Assert(flags_handler_.values_filled(),
                       ExcMessage("Values cache is not filled."));
            }
            else if (k == 1)
            {
                Assert(flags_handler_.gradients_filled(),
                       ExcMessage("Gradients cache is not filled."));
            }
            else if (k == 2)
            {
                Assert(flags_handler_.hessians_filled(),
                       ExcMessage("Hessians cache is not filled."));
            }
            else
            {
                Assert(false,ExcMessage("Derivative order >=3 is not supported."));
            }
#endif


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
    // TODO (pauletti, Mar 17, 2015): this cannot be public, if needed it measn wrong desing
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




};


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space_element-inline.h>


#endif // #ifndef SPACE_ELEMENT_ACCESSOR_

