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

#ifndef BSPLINE_UNIFORM_QUAD_CACHE_H_
#define BSPLINE_UNIFORM_QUAD_CACHE_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/value_flags_handler.h>
#include <igatools/base/quadrature.h>
#include <igatools/utils/cartesian_product_array-template.h>
#include <igatools/utils/value_table.h>
#include <igatools/geometry/grid_uniform_quad_cache.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_scalar_evaluator.h>

IGA_NAMESPACE_OPEN

/**
 * Global BSplineSpace uniform quadrature
 * computational optimization cache, storing the interval length
 * in each direction.
 *
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class BSplineUniformQuadCache : public GridUniformQuadCache<dim_>
{
    using base_t = GridUniformQuadCache<dim_>;
    using Space = BSplineSpace<dim_,range_,rank_>;
    using ElementIterator = typename Space::ElementIterator;

    template<class T>
    using ComponentContainer = typename Space::template ComponentContainer<T>;

    template<class T>
    using ComponentDirectionTable = ComponentContainer<CartesianProductArray<T,dim_>>;

    template<class T>
    using ComponentDirectionContainer = ComponentContainer<std::array<T,dim_>>;

    using SpaceDimensionTable = typename Space::SpaceDimensionTable;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;
    //using typename parent_t::Point;
    using Value = typename Space::Value;
protected:
    using ElementAccessor = typename Space::ElementAccessor;
    void init_element_cache(ElementAccessor &elem);
    void fill_element_cache(ElementAccessor &elem);
public:
    static const int dim = dim_;

    //Allocates and fill the (global) cache
    BSplineUniformQuadCache(std::shared_ptr<const Space> space,
                            const ValueFlags flag,
                            const Quadrature<dim> &quad);

    //Allocates the ElementIterator element_cache
    void init_element_cache(ElementIterator &elem);

    //Fill the ElementIterator element_cache
    void fill_element_cache(ElementIterator &elem);

    void print_info(LogStream &out) const;

    class BasisValues1d
    {
    public:
        BasisValues1d()
        {}
        BasisValues1d(const int max_der_order, const int n_func, const int n_points)
            :
            values_(max_der_order, DenseMatrix(n_func, n_points))
        {}

        void resize(const int max_der_order, const int n_func, const int n_points)
        {
            values_.resize(max_der_order);
            for (auto matrix: values_)
                matrix.resize(n_func, n_points);
        }

        void print_info(LogStream &out) const
        {
            values_.print_info(out);
        }

        auto &get_derivative(const int order)
        {
        	return values_[order];
        }
    private:
        vector<DenseMatrix> values_;
    };

private:
    template <int order>
    using Val = Conditional<(order==0), Value, Derivative<order> >;

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
    		const ComponentDirectionContainer<const BasisValues1d *> &elem_values,
    		ValueTable<Val<order>>& D_phi) const;

private:
    const int n_derivatives = 3;

    std::shared_ptr<const Space> space_;

    SpaceDimensionTable n_basis_;

    ComponentContainer<Size> comp_offset_;

    BasisElemValueFlagsHandler flags_;

    Quadrature<dim> quad_;

    /**
     * univariate B-splines values and derivatives at
     * quadrature points
     * splines1d_[comp][dir][interval][order][function][point]
     */
    ComponentDirectionTable<BasisValues1d> splines1d_;

    ComponentContainer<DynamicMultiArray<std::shared_ptr<BSplineElementScalarEvaluator<dim>>,dim>> scalar_evaluators_;


    using univariate_values_t = ComponentContainer<std::array<const BasisValues1d *,dim>>;

};

IGA_NAMESPACE_CLOSE

#endif /* BSPLINE_UNIFORM_QUAD_CACHE_H_ */
