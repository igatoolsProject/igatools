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
//#include <igatools/base/value_types.h>
//#include <igatools/base/cache_status.h>
//#include <igatools/base/flags_handler.h>

//#include <igatools/base/function.h>
#include <igatools/basis_functions/values_cache.h>

#include <igatools/base/quadrature.h>

//#include <igatools/utils/value_vector.h>
//#include <igatools/utils/value_table.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/cartesian_product_indexer.h>

#include <igatools/basis_functions/spline_space.h>

#include <igatools/basis_functions/space_element_base.h>



IGA_NAMESPACE_OPEN




template<int dim,int codim,int range,int rank>
class SpaceElement : public SpaceElementBase<dim>
{
protected:
    using base_t =  SpaceElementBase<dim>;
private:
    using self_t = SpaceElement<dim,codim,range,rank>;

public:

    using Func = Function<dim,codim,range,rank>;

    using RefPoint = typename Func::RefPoint;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    template <int order>
    using Derivative = typename Func::template Derivative<order>;
    using Div = typename Func::Div;

    static const int space_dim = Func::space_dim;

    using Topology = typename base_t::Topology;

    /**
     * For each component gives a product array of the dimension
     */
    template<class T>
    using ComponentContainer = typename SplineSpace<dim,range,rank>::template ComponentContainer<T>;
    using TensorSizeTable = typename SplineSpace<dim,range,rank>::TensorSizeTable;
    ///@}


    /** @name Constructors */
    ///@{
    using SpaceElementBase<dim>::SpaceElementBase;

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
    virtual ~SpaceElement() = default;
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


    template <class ValueType, int topology_dim = dim>
    auto
    get_basis(const int topology_id, const std::string &dofs_property = DofProperties::active) const
    {
        Assert(local_cache_ != nullptr, ExcNullPtr());
        const auto &cache = local_cache_->template get_value_cache<topology_dim>(topology_id);
        const auto values_all_elem_dofs = cache.template get_der<ValueType>();

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

        const auto n_filtered_dofs = dofs_local_to_elem.size();
        const auto n_pts = values_all_elem_dofs.get_num_points();

        decltype(values_all_elem_dofs) values_filtered_elem_dofs(n_filtered_dofs,n_pts);

        int fn = 0;
        for (const auto loc_dof : dofs_local_to_elem)
        {
            const auto values_all_elem_dofs_fn = values_all_elem_dofs.get_function_view(loc_dof);

            const auto values_filtered_elem_dofs_fn = values_filtered_elem_dofs.get_function_view(fn);

            std::copy(values_all_elem_dofs_fn.begin(),
                      values_all_elem_dofs_fn.end(),
                      values_filtered_elem_dofs_fn.begin());

            ++fn;
        }
        // filtering the values that correspond to the dofs with the given property --- end
        //--------------------------------------------------------------------------------------

        return values_filtered_elem_dofs;
    }


    template <class ValueType>
    auto
    get_basis_element(const std::string &dofs_property = DofProperties::active) const
    {
        return this->template get_basis<ValueType,dim>(0,dofs_property);
    }

    template <class ValueType, int k = dim>
    auto
    linear_combination(const vector<Real> &loc_coefs,
                       const int topology_id,
                       const std::string &dofs_property) const
    {
        const auto &basis_values =
            this->template get_basis<ValueType, k>(topology_id,dofs_property);
        return basis_values.evaluate_linear_combination(loc_coefs) ;
    }





    void print_info(LogStream &out) const;

    void print_cache_info(LogStream &out) const;

    using CType = boost::fusion::map<
            boost::fusion::pair<     _Value,ValueTable<Value>>,
            boost::fusion::pair<  _Gradient,ValueTable<Derivative<1>>>,
            boost::fusion::pair<   _Hessian,ValueTable<Derivative<2>>>,
            boost::fusion::pair<_Divergence,ValueTable<Div>>
            >;

    using Cache = BasisValuesCache<dim,CType,FunctionFlags>;

protected:



    /** The local (element and face) cache. */
    std::shared_ptr<LocalCache<Cache>> local_cache_;

public:
    // TODO (pauletti, Mar 17, 2015): this cannot be public, if needed it means wrong desing
    std::shared_ptr<LocalCache<Cache>> &
                                    get_local_cache()
    {
        return this->local_cache_;
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



#endif // #ifndef SPACE_ELEMENT_H_

