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

#ifndef NURBS_UNIFORM_QUAD_CACHE_H_
#define NURBS_UNIFORM_QUAD_CACHE_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/value_flags_handler.h>
#include <igatools/base/quadrature.h>

//TODO(pauletti, Sep 9, 2014): should we instantiate the cartesian product instead
//#include <igatools/utils/cartesian_product_array-template.h>

#include <igatools/utils/value_table.h>
#include <igatools/geometry/grid_uniform_quad_cache.h>
#include <igatools/basis_functions/nurbs_space.h>

IGA_NAMESPACE_OPEN

/**
 * Global NURBSSpace uniform quadrature
 * computational optimization cache.
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class NURBSUniformQuadCache : public GridUniformQuadCache<dim_>
{
    using base_t = GridUniformQuadCache<dim_>;
    using Space = NURBSSpace<dim_,range_,rank_>;
    static const Size n_components =  Space::n_components;
    using ElementIterator = typename Space::ElementIterator;

    template <int order>
    using Derivative = typename Space::template Derivative<order>;

    using Value = typename Space::Value;


protected:
    using ElementAccessor = typename Space::ElementAccessor;
    void init_element_cache(ElementAccessor &elem);
    void fill_element_cache(ElementAccessor &elem);
    void fill_face_cache(ElementAccessor &elem, const int face);

public:
    static const int dim = dim_;

    //Allocates and fill the (global) cache
    NURBSUniformQuadCache(std::shared_ptr<const Space> space,
                          const ValueFlags flag,
                          const Quadrature<dim> &quad);

    //Allocates the ElementIterator element_cache
    void init_element_cache(ElementIterator &elem);

    //Fill the ElementIterator element_cache
    void fill_element_cache(ElementIterator &elem);

    /**
     * Fills the ElementIterator face_cache
     * element dependent part
     */
    void fill_face_cache(ElementIterator &elem, const int face);

    void print_info(LogStream &out) const;

private:
    // TODO (pauletti, Sep 8, 2014): to be passed some how
//    const int n_derivatives = 3;

    std::shared_ptr<const Space> space_;

//    SpaceDimensionTable n_basis_;

//    ComponentContainer<Size> comp_offset_;

    BasisElemValueFlagsHandler flags_;

    BasisFaceValueFlagsHandler face_flags_;

//    Quadrature<dim> quad_;

    BSplineUniformQuadCache<dim_,range_,rank_> bspline_uniform_quad_cache_;


    using ElementCache = typename BSplineElementAccessor<dim_,range_,rank_>::ValuesCache;

    /**
     * Computes the 0-th order derivative of the non-zero NURBS basis functions over the element
     * at the evaluation points, from the BSpline values contained in <tt>bspline_cache</tt>
     * and the NURBS weights local to the element @p element_weights.
     * \warning If the output result @p D0_phi_hat is not correctly pre-allocated,
     * an exception will be raised.
     */
    void
    evaluate_nurbs_values(
        const ElementCache &bspline_cache,
        const vector<Real> &element_weights,
        ValueTable<Value> &D0_phi_hat) const ;

    /**
     * Computes the 1-st order derivative of the non-zero NURBS basis functions over the element
     * at the evaluation points, from the BSpline values contained in <tt>bspline_cache</tt>
     * and the NURBS weights local to the element @p element_weights.
     * \warning If the output result @p D1_phi_hat is not correctly pre-allocated,
     * an exception will be raised.
     */
    void
    evaluate_nurbs_gradients(
        const ElementCache &bspline_cache,
        const vector<Real> &element_weights,
        ValueTable< Derivative<1> > &D1_phi_hat) const ;

    /**
     * Computes the 2-st order derivative of the non-zero NURBS basis functions over the element,
     * at the evaluation points, from the BSpline values contained in <tt>bspline_cache</tt>
     * and the NURBS weights local to the element @p element_weights.
     * \warning If the output result @p D1_phi_hat is not correctly pre-allocated,
     * an exception will be raised.
     */
    void
    evaluate_nurbs_hessians(
        const ElementCache &bspline_cache,
        const vector<Real> &element_weights,
        ValueTable< Derivative<2> > &D2_phi_hat) const ;


};




IGA_NAMESPACE_CLOSE


#endif // #ifdef NURBS_UNIFORM_QUAD_CACHE_H_
