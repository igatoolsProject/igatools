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


#ifndef NURBS_ELEMENT_HANDLER_H_
#define NURBS_ELEMENT_HANDLER_H_

#include <igatools/base/config.h>

#ifdef NURBS

#include <igatools/base/cache_status.h>
#include <igatools/base/new_flags_handler.h>
#include <igatools/base/quadrature.h>

//TODO(pauletti, Sep 9, 2014): should we instantiate the cartesian product instead
//#include <igatools/utils/cartesian_product_array-template.h>

#include <igatools/utils/value_table.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/basis_functions/nurbs_space.h>


IGA_NAMESPACE_OPEN

/**
 * Global NURBSSpace uniform quadrature
 * computational optimization cache.
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class NURBSElementHandler : public BSplineElementHandler<dim_>
{
    using base_t = GridElementHandler<dim_>;
    using Space = NURBSSpace<dim_,range_,rank_>;
    static const Size n_components =  Space::n_components;

    template<class T>
    using ComponentContainer = typename Space::template ComponentContainer<T>;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;

    using Value = typename Space::Value;

protected:
    using ElementIterator = typename Space::ElementIterator;
    using ElementAccessor = typename Space::ElementAccessor;


public:
    static const int dim = dim_;

    NURBSElementHandler(std::shared_ptr<const Space> space);

    template<int k>
    void reset(const NewValueFlags flag, const Quadrature<k> &quad);

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
        fill_cache<k>(elem.get_accessor(), j);
    }

    template <int k>
    void init_cache(ElementIterator &elem)
    {
        init_cache<k>(elem.get_accessor());
    }

//    void init_all_caches(ElementIterator &elem)
//    {
//        init_all_caches(elem.get_accessor());
//    }


    //Allocates the ElementIterator element_cache
    void init_element_cache(ElementIterator &elem);

    //Fill the ElementIterator element_cache
    void fill_element_cache(ElementIterator &elem);

    void print_info(LogStream &out) const;

private:
    std::shared_ptr<const Space> space_;

    std::array<FunctionFlags, dim + 1> flags_;

//
//
//    using ElementCache = typename BSplineElementAccessor<dim_,range_,rank_>::ValuesCache;
//
//    /**
//     * Computes the 0-th order derivative of the non-zero NURBS basis functions over the element
//     * at the evaluation points, from the BSpline values contained in <tt>bspline_cache</tt>
//     * and the NURBS weights local to the element @p element_weights.
//     * \warning If the output result @p D0_phi_hat is not correctly pre-allocated,
//     * an exception will be raised.
//     */
//    void
//    evaluate_nurbs_values(
//        const ElementCache &bspline_cache,
//        const vector<Real> &element_weights,
//        const ComponentContainer<int> &elem_basis_offset,
//        ValueTable<Value> &D0_phi_hat) const ;
//
//    /**
//     * Computes the 1-st order derivative of the non-zero NURBS basis functions over the element
//     * at the evaluation points, from the BSpline values contained in <tt>bspline_cache</tt>
//     * and the NURBS weights local to the element @p element_weights.
//     * \warning If the output result @p D1_phi_hat is not correctly pre-allocated,
//     * an exception will be raised.
//     */
//    void
//    evaluate_nurbs_gradients(
//        const ElementCache &bspline_cache,
//        const vector<Real> &element_weights,
//        const ComponentContainer<int> &elem_basis_offset,
//        ValueTable< Derivative<1> > &D1_phi_hat) const ;
//
//    /**
//     * Computes the 2-st order derivative of the non-zero NURBS basis functions over the element,
//     * at the evaluation points, from the BSpline values contained in <tt>bspline_cache</tt>
//     * and the NURBS weights local to the element @p element_weights.
//     * \warning If the output result @p D1_phi_hat is not correctly pre-allocated,
//     * an exception will be raised.
//     */
//    void
//    evaluate_nurbs_hessians(
//        const ElementCache &bspline_cache,
//        const vector<Real> &element_weights,
//        const ComponentContainer<int> &elem_basis_offset,
//        ValueTable< Derivative<2> > &D2_phi_hat) const ;


};




IGA_NAMESPACE_CLOSE


#endif // #ifndef NURBS_ELEMENT_HANDLER_H_

#endif
