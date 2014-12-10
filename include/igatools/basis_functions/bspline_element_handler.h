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

#ifndef BSPLINE_ELEMENT_HANDLER_H_
#define BSPLINE_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/flags_handler.h>
#include <igatools/base/quadrature.h>

//TODO(pauletti, Sep 9, 2014): should we instantiate the cartesian product instead
#include <igatools/utils/cartesian_product_array-template.h>

#include <igatools/utils/value_table.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>


IGA_NAMESPACE_OPEN

/**
 * Global BSplineSpace uniform quadrature
 * computational optimization cache, storing the interval length
 * in each direction
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class BSplineElementHandler : public GridElementHandler<dim_>
{
    using base_t = GridElementHandler<dim_>;
    using Space = BSplineSpace<dim_,range_,rank_>;
    static const Size n_components =  Space::n_components;


    template<class T>
    using ComponentContainer = typename Space::template ComponentContainer<T>;

    template<class T>
    using ComponentDirectionTable = ComponentContainer<CartesianProductArray<T,dim_>>;

    template<class T>
    using ComponentDirectionContainer = ComponentContainer<std::array<T,dim_>>;

    using SpaceDimensionTable = typename Space::SpaceDimensionTable;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;

    using Value = typename Space::Value;

protected:
    using ElementIterator = typename Space::ElementIterator;
    using ElementAccessor = typename Space::ElementAccessor;



public:
    static const int dim = dim_;

    //Allocates and fill the (global) cache
    BSplineElementHandler(std::shared_ptr<const Space> space);

    template<int k>
    void reset(const ValueFlags flag, const Quadrature<k> &quad);

//protected:
    template <int k>
    void fill_cache(ElementAccessor &elem, const int j);

    template <int k>
    void init_cache(ElementAccessor &elem);

//    void init_all_caches(ElementAccessor &elem);
public:
    template <int k>
    void fill_cache(ElementIterator &elem, const int j)
    {
        fill_cache<k>(*elem, j);
    }

    template <int k>
    void init_cache(ElementIterator &elem)
    {
        init_cache<k>(*elem);
    }

//    void init_all_caches(ElementIterator &elem)
//    {
//        init_all_caches(*elem);
//    }


    //Allocates the ElementIterator element_cache
    void init_element_cache(ElementIterator &elem);

    //Fill the ElementIterator element_cache
    void fill_element_cache(ElementIterator &elem);

    void print_info(LogStream &out) const;

private:
    const Quadrature<dim> &get_quad() const;

public:

    void
    copy_to_inactive_components_values(const vector<Index> &inactive_comp,
                                       const std::array<Index, n_components> &active_map,
                                       ValueTable<Value> &D_phi) const;

    template <int order>
    void
    copy_to_inactive_components(const vector<Index> &inactive_comp,
                                const std::array<Index, n_components> &active_map,
                                ValueTable<Derivative<order>> &D_phi) const;


private:

    /**
     * Computes the k-th order derivative of the non-zero B-spline basis
     * functions over the current element,
     *   at the evaluation points pre-allocated in the cache.
     *
     * \warning If the output result @p derivatives_phi_hat is not correctly pre-allocated,
     * an exception will be raised.
     */
    template <int order>
    void evaluate_bspline_derivatives(
        const  ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
        ValueTable<Derivative<order>> &D_phi) const;

    void evaluate_bspline_values(
        const  ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
        ValueTable<Value> &D_phi) const;

private:
    std::shared_ptr<const Space> space_;

    SpaceDimensionTable n_basis_;

    ComponentContainer<Size> comp_offset_;

    std::array<FunctionFlags, dim + 1> flags_;

    template <class T>
    using DirectionTable = CartesianProductArray<T, dim_>;
    using BasisValues = ComponentContainer<BasisValues1d>;
    /**
     * B-splines values and derivatives at quadrature points.
     * The values are stored in the un tensor product way.
     *
     * splines1d_[dir][interval][comp][order][function][point]
     */
    class GlobalCache : public DirectionTable<BasisValues>
    {
    public:
        using DirectionTable<BasisValues>::DirectionTable;
        // TODO (pauletti, Sep 23, 2014): document and split definition
        auto get_element_values(const TensorIndex<dim> &id)
        {
            ComponentContainer<TensorProductFunctionEvaluator<dim> >
            result((this->entry(0,0)).get_comp_map());
            for (auto c : result.get_active_components_id())
            {
                for (int i = 0; i < dim; ++i)
                    result[c][i] =
                        BasisValues1dConstView((this->entry(i, id[i]))[c]);
                result[c].update_size();

            }
            return result;
        }
    };

    CacheList<GlobalCache, dim> splines1d_;
};

IGA_NAMESPACE_CLOSE


#endif
