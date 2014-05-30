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

#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/basis_functions/bernstein_basis.h>
#include <igatools/base/exceptions.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/base/function.h>
#include <igatools/utils/multi_array_utils.h>

#include <algorithm>
#include <numeric>
#include <memory>

#include <boost/numeric/ublas/io.hpp>

using std::reverse;
using std::accumulate;
using std::sort;

using std::shared_ptr;
using std::make_shared;

using std::array;
using std::vector;


IGA_NAMESPACE_OPEN


//#define NOT_OPTIMIZED

namespace
{

#ifndef NOT_OPTIMIZED

/**
 * @class DerivativeEvaluationSymmetryManager
 *
 * This class is used to manage the symmetries in the evaluation of a k-th order
 *  derivative of a scalar function \f$ f \mathbb{R}^d -to \mathbb{R} \f$
 *  according with the Schwarz theorem (equality of mixed partials).
 *
 * The total number of derivatives is given by \f$ d^k \f$ where the number of
 * different values
 * is \f$ \binorm{d-1+k}{d-1} = \binorm{d-1+k}{k} \f$.
 *
 */
// TODO (pauletti, Nov 1, 2013): Document how this class work and its internals
template<int dim, int order>
class DerivativeSymmetryManager
{
public:
    static const int num_entries_total = pow(dim,order);
    static const int num_entries_eval = constexpr_binomial_coefficient(dim-1+order,order);
    static const int num_entries_copy = num_entries_total - num_entries_eval;

    typedef Derivatives<dim,1,1,order> Derivative_t;

    /** @name Constructors */
    ///@{
    /** Constructor */
    DerivativeSymmetryManager();

    /** Copy constructor */
    DerivativeSymmetryManager(
        const DerivativeSymmetryManager<dim,order> &in) = default;

    /** Move constructor */
    DerivativeSymmetryManager(
        DerivativeSymmetryManager<dim,order> &&in) = default;

    ~DerivativeSymmetryManager() = default;
    ///@}


    /** @name Assignment operators */
    ///@{
    /** Copy assignment operator */
    DerivativeSymmetryManager<dim,order> &operator=(
        const DerivativeSymmetryManager<dim,order> &) = default;

    /** Move assignment operator */
    DerivativeSymmetryManager<dim,order> &operator=(
        DerivativeSymmetryManager<dim,order> &&) = default;
    ///@}


    static const int new_deriv_order = order>0?order:1;
    const array<int,num_entries_eval> &get_entries_flat_id_evaluate() const
    {
        return entries_flat_id_evaluate_;
    }

    const array<int,num_entries_copy> &get_entries_flat_id_copy_to() const
    {
        return entries_flat_id_copy_to_;
    }

    const array<int,num_entries_copy> &get_entries_flat_id_copy_from() const
    {
        return entries_flat_id_copy_from_;
    }

    const array<TensorIndex<new_deriv_order>,num_entries_total> &get_entries_tensor_id() const
    {
        return entries_tensor_id_;
    }
private:

    bool test_if_evaluate(const array<int,order> &tensor_index) const;

    /** Tensor ids of all the entries */

    array<TensorIndex<new_deriv_order> ,num_entries_total> entries_tensor_id_;

    /** Flat ids of the entries that need to be computed */
    array<int,num_entries_eval> entries_flat_id_evaluate_;

    /** Flat ids of the destination entries that need to be copied */
    array<int,num_entries_copy> entries_flat_id_copy_to_;

    /** Flat ids of the source entries that need to be copied */
    array<int,num_entries_copy> entries_flat_id_copy_from_;

    TensorSize<order> size_deriv_index_;
};




template<int dim, int order>
DerivativeSymmetryManager<dim,order>::
DerivativeSymmetryManager()
{
    size_deriv_index_.fill(dim);

    typedef MultiArrayUtils<order> MAUtils;


    Derivatives<dim,1,1,order> derivative;

    int eval_id = 0;
    int copy_id = 0;
    if (order == 0)
    {
        for (int flat_id = 0; flat_id < num_entries_total; ++flat_id)
        {
            entries_flat_id_evaluate_[eval_id] = flat_id;
            eval_id++;
        }
    }
    else
    {
        auto weights = MAUtils::compute_weight(size_deriv_index_);

        for (Index flat_id = 0; flat_id < num_entries_total; ++flat_id)
        {
            entries_tensor_id_[flat_id] = derivative.flat_to_tensor_index(flat_id);

            auto tensor_id = MAUtils::flat_to_tensor_index(flat_id,weights);

            if (test_if_evaluate(tensor_id))
            {
                entries_flat_id_evaluate_[eval_id] = flat_id;
                eval_id++;
            }
            else
            {
                entries_flat_id_copy_to_[copy_id] = flat_id;
                sort(tensor_id.begin(),tensor_id.end());
                reverse(tensor_id.begin(),tensor_id.end());
                entries_flat_id_copy_from_[copy_id] = MAUtils::tensor_to_flat_index(tensor_id,weights);
                copy_id++;
            }
        }
    }
    Assert(eval_id == num_entries_eval, ExcDimensionMismatch(eval_id,num_entries_eval));
    Assert(copy_id == num_entries_copy, ExcDimensionMismatch(copy_id,num_entries_copy));

}


template<int dim, int order>
inline
bool
DerivativeSymmetryManager<dim,order>::
test_if_evaluate(const array<int,order> &tensor_index) const
{
    bool test_result = true;
    for (int i = 0; i < order-1; ++i)
    {
        if (tensor_index[i+1] > tensor_index[i])
        {
            test_result = false;
            break;
        }
    }
    return test_result;
}


template<int order>
class DerivativeSymmetryManager<0,order>
{
public:
    static const int num_entries_total = 1;
    static const int num_entries_eval  = 1;
    static const int num_entries_copy  = num_entries_total - num_entries_eval;

    typedef Derivatives<0,1,1,order> Derivative_t;

    DerivativeSymmetryManager()
    {
        entries_tensor_id_[0] = {{}};
        entries_flat_id_evaluate_[0] = 0;
    }

    const array<int,num_entries_eval> &get_entries_flat_id_evaluate() const
    {
        return entries_flat_id_evaluate_;
    }

    const array<int,num_entries_copy> &get_entries_flat_id_copy_to() const
    {
        return entries_flat_id_copy_to_;
    }

    const array<int,num_entries_copy> &get_entries_flat_id_copy_from() const
    {
        return entries_flat_id_copy_from_;
    }

    const array<TensorIndex<order>,num_entries_total> &get_entries_tensor_id() const
    {
        return entries_tensor_id_;
    }

private:
    /** Tensor ids of all the entries */
    array<TensorIndex<order>,num_entries_total> entries_tensor_id_;

    /** Flat ids of the entries that need to be computed */
    array<int,num_entries_eval> entries_flat_id_evaluate_;

    /** Flat ids of the destination entries that need to be copied */
    array<int,num_entries_copy> entries_flat_id_copy_to_;

    /** Flat ids of the source entries that need to be copied */
    array<int,num_entries_copy> entries_flat_id_copy_from_;

};


};
#endif



template <int dim, int range, int rank>
BSplineElementAccessor<dim, range, rank>::
BSplineElementAccessor(const std::shared_ptr<ContainerType> space,
                       const int index)
    :
    SpaceElementAccessor<
    BSplineElementAccessor<dim,range,rank>,BSplineSpace<dim, range, rank>,dim,0,range,rank>(space,index)
{}



template <int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
init_values(const ValueFlags fill_flag,
            const Quadrature<dim> &quad)
{
    Assert((fill_flag|admisible_flag) == admisible_flag,
           typename CartesianGridElementAccessor<dim>::ExcFillFlagNotSupported(admisible_flag, fill_flag));

    // initalizing the cache of the CartesianGridElementAccessor
    {
        ValueFlags grid_flag = ValueFlags::none;
        if (contains(fill_flag , ValueFlags::point))
            grid_flag |= ValueFlags::point;
        if (contains(fill_flag , ValueFlags::measure))
            grid_flag |= ValueFlags::measure;
        if (contains(fill_flag , ValueFlags::w_measure))
            grid_flag |= ValueFlags::w_measure;
        if (contains(fill_flag , ValueFlags::face_point))
            grid_flag |= ValueFlags::face_point;
        if (contains(fill_flag , ValueFlags::face_w_measure))
            grid_flag |= ValueFlags::face_w_measure;
        CartesianGridElementAccessor<dim>::init_values(grid_flag,quad);
    }

    auto f_flag = fill_flag;

    if (contains(f_flag, ValueFlags::divergence))
        f_flag |= ValueFlags::gradient;

    if (contains(f_flag, ValueFlags::face_divergence))
        f_flag |= ValueFlags::face_gradient;

    int max_der_order = -1;

    if (contains(f_flag, ValueFlags::value))
        max_der_order=std::max(max_der_order,0);

    if (contains(f_flag, ValueFlags::face_value))
        max_der_order=std::max(max_der_order,0);

    if (contains(f_flag, ValueFlags::gradient))
        max_der_order=std::max(max_der_order,1);

    if (contains(f_flag, ValueFlags::face_gradient))
        max_der_order=std::max(max_der_order,1);

    if (contains(f_flag, ValueFlags::hessian))
        max_der_order=std::max(max_der_order,2);

    if (contains(f_flag, ValueFlags::face_hessian))
        max_der_order=std::max(max_der_order,2);


    Assert(max_der_order>=0, ExcMessage("Not a right ValueFlag"));

    reset_univariate_cache(quad, max_der_order);

    this->reset_element_and_faces_cache(f_flag, quad);
}



template <int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
init_face_values(const Index face_id,
                 const ValueFlags fill_flag,
                 const Quadrature<dim-1> &quad)
{
    AssertThrow(false,ExcNotImplemented()) ;
}



template <int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
reset_global_cache()
{
    values_1d_elem_.reset();

    for (int f = 0; f < n_faces; ++f)
        values_1d_faces_[f].reset();
}

template <int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
reset_univariate_cache(const Quadrature<dim> &quad, const int max_der)
{
    Assert(values_1d_elem_.use_count() < 2,
           ExcMessage("Resetting a shared cache, use force if this is what you Really want."));

    if (values_1d_elem_.use_count() == 0)
        values_1d_elem_= make_shared<GlobalElemCache>();

    values_1d_elem_->reset(*(this->space_), quad, max_der);

    for (int f = 0; f < n_faces; ++f)
    {
        Assert(values_1d_faces_[f].use_count() < 2,
               ExcMessage("Resetting a shared cache, use force if this is what you Really want."));

        if (values_1d_faces_[f].use_count() == 0)
            values_1d_faces_[f]= make_shared<GlobalFaceCache>();
        values_1d_faces_[f]->reset(*(this->space_), quad, f, max_der);
    }
}






// TODO (pauletti, May 30, 2014): we may not need this
#if 0
template <int dim, int range, int rank>
auto
BSplineElementAccessor<dim, range, rank>::
get_bezier_extraction_operator() const -> ComponentTable< std::array< const DenseMatrix *,dim> >
{
    ComponentTable< std::array< const DenseMatrix *,dim> > bezier_op;

    const int n_components = bezier_op.flat_size();

    TensorIndex<dim> element_tid = this->get_tensor_index();

    for (int comp_id = 0 ; comp_id < n_components ; ++comp_id)
    {
        for (int dir = 0 ; dir < dim ; ++dir)
            bezier_op(comp_id)[dir] = this->space_->bezier_op_(comp_id).get_data_direction(dir)[element_tid[dir]];
    }

    return bezier_op;

}
#endif


template < int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
GlobalCache::
reset(const Space &space,
      const Quadrature<dim> &quad,
      const int max_der,
      const std::array<std::vector<int>,dim> &intervals_id)
{
    //------------------------------------------------------------------------------------------
    // resizing the structures for the one dimensional splines
    this->max_deriv_order_ = max_der;

    for (int i = 0 ; i < dim ; ++i)
    {
        n_intervals_(i) = intervals_id[i].size();

        Assert(n_intervals_(i) == space.get_grid()->get_num_elements_dim()(i) ||
               n_intervals_(i) == 1,
               ExcMessage("Invalid number of intervals along the direction " + std::to_string(i)));
    }


    // resizing the structures for the one dimensional splines
    const int max_der_plus_one = this->max_deriv_order_ + 1;

    for (int iComp = 0; iComp < space.num_active_components_; ++iComp)
    {
        this->splines1d_cache_data_(iComp).resize(n_intervals_);
        for (int i = 0 ; i < dim ; ++i)
            for (int j = 0 ; j < n_intervals_(i) ; ++j)
                this->splines1d_cache_data_(iComp).entry(i,j).resize(max_der_plus_one);
    }

    for (int iComp = 0; iComp < Space::n_components; ++iComp)
    {
        Assert(n_intervals_ ==
               this->splines1d_cache_data_(space.map_component_to_active_data_(iComp)).tensor_size(),
               ExcMessage("Not same size"));
        this->splines1d_cache_(iComp).resize(this->n_intervals_);
    }

    this->set_initialized(true);
    //------------------------------------------------------------------------------------------





    /*
     * For each component and each direction we consider the number
     * of intervals in the space.
     * Then in each interval we compute the values and derivatives of
     * the one dimensional B-splines on each quadrature point.
     */

    const auto &bezier_op_ = space.bezier_op_;
    const auto &degree_    = space.degree_;
    const auto &eval_points = quad.get_points();

    const auto lengths = space.get_grid()->get_element_lengths();
    BasisValues1d bernstein_values(max_deriv_order_+1);
    for (int iComp = 0; iComp < space.num_active_components_; ++iComp)
    {
        for (int jDim = 0; jDim < dim; ++jDim)
        {
            const auto &intervals_id_dir = intervals_id[jDim];

            const int num_intervals = intervals_id_dir.size();
            Assert(num_intervals == n_intervals_[jDim],
                   ExcDimensionMismatch(num_intervals, n_intervals_[jDim]));
            const int degree = degree_(iComp)[jDim];
            const vector<Real> &pt_coords = eval_points.get_data_direction(jDim);


            // fill values and derivatives of the Bernstein's polynomials at
            // quad points in [0,1]

            for (int deriv_order = 0; deriv_order <= max_deriv_order_; ++deriv_order)
                bernstein_values[ deriv_order ] =
                    BernsteinBasis::derivative(deriv_order, degree, pt_coords);

            const auto &bez_iComp_jDim = bezier_op_(iComp).get_data_direction(jDim);
            const auto &lengths_jDim = lengths.get_data_direction(jDim);

            // compute the one dimensional B-splines at quad point on the reference interval
            for (int i = 0 ; i < num_intervals ; ++i)
            {
                const auto &M = *(bez_iComp_jDim[intervals_id_dir[i]]);
                const Real one_div_size = 1.0 / lengths_jDim[intervals_id_dir[i]];
                BasisValues1d &basis = splines1d_cache_data_(iComp).entry(jDim,i);

                for (int deriv_order = 0; deriv_order <= max_deriv_order_; ++deriv_order)
                {
                    const Real scaling_coef = std::pow(one_div_size, deriv_order);
                    basis[ deriv_order ] = scaling_coef * prec_prod(M, bernstein_values[ deriv_order ]);
                } //end loop deriv_order

            } // end loop interval
        } // end loop jDim
    } // end loop iComp


    //------------------------------------------------------------------------------------------
    // assign the basis1D data to the proper component/interval through their memory address
    for (int iComp = 0; iComp < Space::n_components; iComp++)
    {
        const int active_comp = space.map_component_to_active_data_(iComp);

        const auto &spline1d_active_data = splines1d_cache_data_(active_comp);

        const auto n_intervals_multi_d = spline1d_active_data.tensor_size();

        for (int jDim = 0; jDim < dim; jDim++)
        {
            const Size n_intervals_1d = n_intervals_multi_d[jDim];

            const auto &data = spline1d_active_data.get_data_direction(jDim);

            for (Size i = 0; i < n_intervals_1d; ++i)
                splines1d_cache_(iComp).entry(jDim,i) = &(data[i]);
        } // end loop jDim
    } // end loop iComp
    //-------------------------------------------------------------------------


    this->set_filled(true);
}

template < int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
GlobalElemCache::
reset(const Space &space,
      const Quadrature<dim> &quad,
      const int max_der)
{

    //------------------------------------------------------------------------------------------
    // resizing the structures for the one dimensional splines
    const auto grid = space.get_grid();
    array<vector<int>,dim> intervals_id;
    for (int i = 0 ; i < dim ; ++i)

    {
        auto &intervals_id_direction = intervals_id[i];

        const int n_intervals = grid->get_num_elements_dim()(i);

        intervals_id_direction.resize(n_intervals);

        for (int id = 0 ; id < n_intervals ; ++id)
            intervals_id_direction[id] = id;
    }

    GlobalCache::reset(space,quad,max_der,intervals_id);

#ifndef NDEBUG
    for (int i = 0 ; i < dim ; ++i)
        Assert(this->n_intervals_(i) == grid->get_num_elements_dim()(i),
               ExcDimensionMismatch(this->n_intervals_(i),grid->get_num_elements_dim()(i)));
#endif
    //------------------------------------------------------------------------------------------

}


template <int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
GlobalFaceCache::
reset(const Space &space,
      const Quadrature<dim> &quad1,
      const Index face_id,
      int max_der)
{
    const auto quad = quad1.collapse_to_face(face_id);

    const int const_dir = UnitElement<dim>::face_constant_direction[face_id];


    array<vector<int>,dim> intervals_id;
    for (int i = 0 ; i < dim ; ++i)

    {
        auto &intervals_id_direction = intervals_id[i];
        if (i != const_dir)
        {
            const int n_intervals = space.get_grid()->get_num_elements_dim()(i);

            intervals_id_direction.resize(n_intervals);
            for (int id = 0 ; id < n_intervals ; ++id)
                intervals_id_direction[id] = id;
        }
        else
        {
            intervals_id_direction.clear();
            if (face_id % 2 == 0)
                intervals_id_direction.push_back(0);
            else
                intervals_id_direction.push_back(space.get_grid()->get_num_elements_dim()(i)-1);

            Assert(intervals_id_direction.size() == 1,
                   ExcDimensionMismatch(intervals_id_direction.size(),1));
        }
    }

    GlobalCache::reset(space,quad,max_der,intervals_id);

#ifndef NDEBUG
    for (int i = 0 ; i < dim ; ++i)
    {
        if (i != const_dir)
        {
            Assert(this->n_intervals_(i) == space.get_grid()->get_num_elements_dim()(i),
                   ExcDimensionMismatch(this->n_intervals_(i),space.get_grid()->get_num_elements_dim()(i)));
        }
        else
        {
            Assert(this->n_intervals_(i) == 1,ExcDimensionMismatch(this->n_intervals_(i),1));
        }
    }
#endif
}


template <int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
fill_values_cache_from_univariate(const int max_deriv_order,
                                  const univariate_values_t &univariate_values,
                                  ValuesCache &cache)
{
    //--------------------------------------------------------------------------
    vector<std::array<Values1DConstView,dim>> values1D(max_deriv_order+1);

    TensorSize<dim> n_basis_direction;

    const auto &degree = this->get_space()->get_degree();

    for (int comp = 0; comp < Space::n_components; ++comp)
    {
        for (int i = 0; i < dim ; ++i)
            n_basis_direction(i) = degree(comp)[i]+1;


        auto &scalar_evaluator_comp = scalar_evaluators_(comp);

        scalar_evaluator_comp.resize(n_basis_direction);

        const auto &univariate_values_comp = univariate_values(comp);

        const Size n_basis = scalar_evaluator_comp.flat_size();

        for (Index flat_basis_id = 0 ; flat_basis_id < n_basis ; ++flat_basis_id)
        {
            const auto tensor_basis_id = scalar_evaluator_comp.flat_to_tensor(flat_basis_id);

            for (int dir = 0 ; dir < dim ; ++dir)
            {
                const auto &basis_with_ders = univariate_values_comp[dir];

                Assert(values1D.size() == basis_with_ders->size(),
                       ExcDimensionMismatch(values1D.size(),basis_with_ders->size()));

                for (int order = 0 ; order <= max_deriv_order ; ++order)
                {
                    const DenseMatrix &funcs = (*basis_with_ders)[order];
                    values1D[order][dir] = Values1DConstView(funcs,tensor_basis_id[dir]);
                } //end order loop

            } // end dir loop


            scalar_evaluator_comp(flat_basis_id) =
                shared_ptr<BSplineElementScalarEvaluator<dim>>(
                    new BSplineElementScalarEvaluator<dim>(values1D));

        } // end flat_basis_id loop

    } // end icomp loop
    //--------------------------------------------------------------------------






    //--------------------------------------------------------------------------
    if (cache.flags_handler_.fill_values())
    {
        evaluate_bspline_derivatives<0>(univariate_values,
                                        cache,
                                        cache.phi_);

        cache.flags_handler_.set_values_filled(true);
    }

    if (cache.flags_handler_.fill_gradients())
    {
        evaluate_bspline_derivatives<1>(univariate_values,
                                        cache,
                                        cache.D1phi_);

        cache.flags_handler_.set_gradients_filled(true);
    }

    if (cache.flags_handler_.fill_hessians())
    {
        evaluate_bspline_derivatives<2>(univariate_values,
                                        cache,
                                        cache.D2phi_);

        cache.flags_handler_.set_hessians_filled(true);
    }

    if (cache.flags_handler_.fill_divergences())
    {
        Assert(cache.flags_handler_.gradients_filled(),
               ExcMessage("Divergence requires gradient to be filled."));

        auto D1  = cache.D1phi_.begin();
        auto div = cache.div_phi_.begin();
        auto end = cache.D1phi_.end();
        for (; D1 != end; ++D1, ++div)
            *div = trace(*D1);

        cache.flags_handler_.set_divergences_filled(true);
    }
    //--------------------------------------------------------------------------

    cache.set_filled(true);
}



template <int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
fill_values(const TopologyId<dim> &topology_id)
{
    Assert(topology_id.is_element() || topology_id.is_face(),
           ExcMessage("Only element or face topology is allowed."));

    auto &cache = this->get_values_cache(topology_id);
    Assert(cache.is_initialized(), ExcNotInitialized());

    const auto &element_tensor_id = this->get_tensor_index();
    ComponentContainer<array<const BasisValues1d *, dim>> elem_univariate_values;


    CartesianGridElementAccessor<dim>::fill_values(topology_id);

    if (topology_id.is_element())
    {
        const auto global_elem_cache = values_1d_elem_;
        Assert(global_elem_cache->is_filled(), ExcCacheNotFilled());

        for (int iComp=0; iComp< Space::n_components; ++iComp)
        {
            const auto &values_1D_comp = global_elem_cache->splines1d_cache_(iComp);
            for (int i = 0; i < dim; ++i)
                elem_univariate_values(iComp)[i] = values_1D_comp.get_data_direction(i)[element_tensor_id[i]];
        }
    } // if (topology_id.is_element())
    else // if (topology_id.is_face())
    {
        const int face_id = topology_id.get_id();

        const auto global_face_cache = values_1d_faces_[face_id];
        Assert(global_face_cache->is_filled(), ExcCacheNotFilled());

        const int const_dir = UnitElement<dim>::face_constant_direction[face_id];

        for (int iComp=0; iComp < Space::n_components ; ++iComp)
        {
            const auto &values_1D_comp = global_face_cache->splines1d_cache_(iComp);
            for (int i = 0; i < dim; ++i)
            {
                if (i != const_dir)
                    elem_univariate_values(iComp)[i] =
                        values_1D_comp.get_data_direction(i)[element_tensor_id[i]];
                else
                    elem_univariate_values(iComp)[i] =
                        values_1D_comp.get_data_direction(i)[0];
            } // end loop i
        } // end loop iComp

        Assert(values_1d_elem_->max_deriv_order_ == global_face_cache->max_deriv_order_,
               ExcDimensionMismatch(values_1d_elem_->max_deriv_order_,global_face_cache->max_deriv_order_));
    } // if (topology_id.is_face())


    this->fill_values_cache_from_univariate(
        values_1d_elem_->max_deriv_order_,
        elem_univariate_values,
        cache);
}








template <int dim, int range, int rank>
template < int deriv_order >
void
BSplineElementAccessor<dim, range, rank>::
evaluate_bspline_derivatives(const ComponentContainer<std::array<const BasisValues1d *, dim> > &elem_values,
                             const ValuesCache &cache,
                             ValueTable< Conditional<(deriv_order==0),Value,Derivative<deriv_order> > >
                             &derivatives_phi_hat) const
{
    Assert(derivatives_phi_hat.size() > 0, ExcEmptyObject());
    Assert(derivatives_phi_hat.get_num_functions() == this->get_num_basis(),
           ExcDimensionMismatch(derivatives_phi_hat.get_num_functions(),this->get_num_basis()));

    const TensorSize<dim> n_points_direction = cache.quad_.get_num_points_direction();
    const Size num_points = n_points_direction.flat_size();
    Assert(derivatives_phi_hat.get_num_points() == num_points,
           ExcDimensionMismatch(derivatives_phi_hat.get_num_points(),num_points));

    /*
     * This code computes any order of derivatives for a multivariate
     * B-spline on the current element
     * We use the formula
     * \partial_(\alpha_1,...,\alpha_n) B(qp) = \Pi d^{\alpha_i} B_i(qp_i)
     */


    Assert(rank < 2, ExcMessage("For rank> 1 the basis function are not implemented/tested."));




    if (deriv_order == 0)
    {
        TensorIndex<dim> zero_tensor_id; // [0,0,..,0] tensor index
        for (int iComp = 0; iComp < this->space_->num_active_components_; ++iComp)
        {
            const int n_basis = this->get_num_basis(iComp);
            Assert(n_basis == this->n_basis_direction_(iComp).flat_size(), ExcMessage("different sizes"));

            const Size comp_offset_i = this->comp_offset_(iComp);

            DynamicMultiArray<Real,dim> derivative_scalar_component(n_points_direction);
            for (int func_flat_id = 0; func_flat_id < n_basis; ++func_flat_id)
            {
                auto derivatives_phi_hat_ifn = derivatives_phi_hat.get_function_view(comp_offset_i+func_flat_id);

                const auto &scalar_bspline = *scalar_evaluators_(iComp)(func_flat_id);

                //TODO: remove this if!!! (Maybe re-think about the BSplineSpace for dim==0)
                if (dim > 0)
                {
                    scalar_bspline.evaluate_derivative_at_points(zero_tensor_id,derivative_scalar_component);

                    for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
                        derivatives_phi_hat_ifn[point_flat_id](iComp) = derivative_scalar_component(point_flat_id);
                }
                else
                {
                    for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
                        derivatives_phi_hat_ifn[point_flat_id](iComp) = 1.0;
                }

            } // end func_flat_id loop

        } // end iComp loop


        if (this->space_->homogeneous_range_)
        {
            const auto n_basis = this->space_->get_num_basis_per_element(0);
            for (int comp = 1; comp < Space::n_components; ++comp)
            {
                const Size offset = this->comp_offset_(comp);
                for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
                {
                    const auto values_phi_hat_copy_from = derivatives_phi_hat.get_function_view(basis_i);
                    auto values_phi_hat_copy_to = derivatives_phi_hat.get_function_view(offset+basis_i);

                    for (int qp = 0; qp < num_points; ++qp)
                        values_phi_hat_copy_to[qp](comp) = values_phi_hat_copy_from[qp](0);

                } //end loop basis_i
            } // end loop comp
        } // end if (space_->homogeneous_range_)

    } // end if (deriv_order == 0)
    else // if (deriv_order != 0)
    {
        Assert(deriv_order >= 1, ExcLowerRange(deriv_order,1));

        using DerSymmMngr_t = DerivativeSymmetryManager<dim,deriv_order>;
        DerSymmMngr_t derivative_symmetry_manager;
        const auto &derivatives_flat_id_evaluate = derivative_symmetry_manager.get_entries_flat_id_evaluate();
        const auto &derivatives_flat_id_copy_to = derivative_symmetry_manager.get_entries_flat_id_copy_to();
        const auto &derivatives_flat_id_copy_from = derivative_symmetry_manager.get_entries_flat_id_copy_from();

        const auto &derivatives_tensor_id = derivative_symmetry_manager.get_entries_tensor_id();

        const Size n_derivatives_eval = DerSymmMngr_t::num_entries_eval;
        const Size n_derivatives_copy = DerSymmMngr_t::num_entries_copy;

        using der_t = Conditional<deriv_order==0,Values<dim,range,rank>,Derivatives<dim,range,rank,deriv_order>>;

        for (int iComp = 0; iComp < this->space_->num_active_components_; ++iComp)
        {
            const int n_basis = this->get_num_basis(iComp);
            Assert(n_basis == this->n_basis_direction_(iComp).flat_size(), ExcMessage("different sizes"));

            const Size comp_offset_i = this->comp_offset_(iComp);

            DynamicMultiArray<Real,dim> derivative_scalar_component(n_points_direction);
            for (int func_flat_id = 0; func_flat_id < n_basis; ++func_flat_id)
            {
                auto derivatives_phi_hat_ifn = derivatives_phi_hat.get_function_view(comp_offset_i+func_flat_id);

                const auto &scalar_bspline = *scalar_evaluators_(iComp)(func_flat_id);

                for (int entry_id = 0; entry_id < n_derivatives_eval; ++entry_id)
                {
                    const int entry_flat_id = derivatives_flat_id_evaluate[entry_id];
                    const auto entry_tensor_id = derivatives_tensor_id[entry_flat_id];

                    // from the entry_tensor_id we get the right derivative order
                    TensorIndex<dim> deriv_order_tensor_id; // [0,0,..,0] tensor index
                    for (int i = 0; i < deriv_order; ++i)
                        ++(deriv_order_tensor_id[entry_tensor_id[i]]);

                    //TODO: remove this if!!! (Maybe re-think about the BSplineSpace for dim==0)
                    if (dim > 0)
                    {
                        scalar_bspline.evaluate_derivative_at_points(deriv_order_tensor_id,derivative_scalar_component);

                        for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
                            derivatives_phi_hat_ifn[point_flat_id](entry_flat_id)(iComp) = derivative_scalar_component(point_flat_id);
                    }
                    else
                    {
                        for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
                            derivatives_phi_hat_ifn[point_flat_id](entry_flat_id)(iComp) = 1.0;
                    }

                } // end entry_id loop

            } // end func_flat_id loop

            for (int func_flat_id = 0; func_flat_id < n_basis; ++func_flat_id)
            {
                auto derivatives_phi_hat_ifn = derivatives_phi_hat.get_function_view(comp_offset_i+func_flat_id);
                for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
                {
                    der_t &derivative = derivatives_phi_hat_ifn[point_flat_id];

                    // here we copy the computed quantities to the symmetric part of the tensor
                    for (int entry_id = 0; entry_id < n_derivatives_copy; ++entry_id)
                        derivative(derivatives_flat_id_copy_to[entry_id])(iComp) =
                            derivative(derivatives_flat_id_copy_from[entry_id])(iComp);
                } // end point_flat_id loop
            } //end func_flat_id loop

        } // end iComp loop

        if (this->space_->homogeneous_range_)
        {
            const Size n_ders = Derivative<deriv_order>::size;
            const auto n_basis = this->space_->get_num_basis_per_element(0);
            for (int comp = 1; comp < Space::n_components; ++comp)
            {
                const Size offset = this->comp_offset_(comp);
                for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
                {
                    const auto derivatives_phi_hat_copy_from = derivatives_phi_hat.get_function_view(basis_i);
                    auto derivatives_phi_hat_copy_to = derivatives_phi_hat.get_function_view(offset+basis_i);
                    for (int qp = 0; qp < num_points; ++qp)
                    {
                        const der_t &values_0 = derivatives_phi_hat_copy_from[qp];
                        der_t &values = derivatives_phi_hat_copy_to[qp];

                        for (int der = 0; der < n_ders; ++der)
                            values(der)(comp) = values_0(der)(0);
                    }
                } //end loop basis_i
            } // end loop comp
        } // end if (space_->homogeneous_range_)

    } // end if (deriv_order > 0)

}








template <int dim, int range, int rank>
auto
BSplineElementAccessor<dim, range, rank>::
get_scalar_evaluators() const ->
const ComponentContainer<
DynamicMultiArray<
shared_ptr<BSplineElementScalarEvaluator<dim>>,dim>> &
{
    return scalar_evaluators_;
}


template <int dim, int range, int rank>
template<int deriv_order>
auto
BSplineElementAccessor<dim, range, rank>::
evaluate_basis_derivatives_at_points(const vector<Point<dim>> &points) const ->
ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >
{
    using return_t = ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >;

    const Size n_basis  = this->get_num_basis();
    const Size n_points = points.size();

    return_t derivatives_phi_hat(n_basis,n_points);

    const auto bezier_op = this->get_bezier_extraction_operator();

    if (deriv_order == 0)
    {
        for (Size pt_id = 0 ; pt_id < n_points ; ++pt_id)
        {
            const Point<dim> point = points[pt_id];

            auto derivatives_phi_hat_ipt = derivatives_phi_hat.get_point_view(pt_id);

#ifndef NDEBUG
            for (int dir = 0 ; dir < dim ; ++dir)
                Assert(point[dir] >= 0.0 && point[dir] <= 1.0,
                ExcMessage("Evaluation point " + std::to_string(pt_id) + " not in the unit-domain."));
#endif

            for (int iComp = 0; iComp < this->space_->num_active_components_; ++iComp)
            {
                //------------------------------------------------------------------------------
                // evaluation of the values/derivarives of the 1D Bernstein polynomials -- begin
                array<boost::numeric::ublas::vector<Real>,dim> bernstein_values;
                const TensorSize<dim> basis_component_t_size = this->n_basis_direction_(iComp);
                for (int dir = 0 ; dir < dim ; ++dir)
                {
                    const int n_basis_1D = basis_component_t_size(dir);
                    const int degree = n_basis_1D - 1 ;
                    bernstein_values[dir] = BernsteinBasis::derivative(0,degree,point[dir]);
                }
                // evaluation of the values/derivarives of the 1D Bernstein polynomials -- end
                //------------------------------------------------------------------------------


                //--------------------------------------------------------------------------------
                // apply the Bezier extraction operator for the functions on this element -- begin
                const auto &bezier_op_comp = bezier_op(iComp);

                array<boost::numeric::ublas::vector<Real>,dim> bspline_basis;
                for (int dir = 0 ; dir < dim ; ++dir)
                {
                    const auto &M = *(bezier_op_comp[dir]);

                    bspline_basis[dir] = prec_prod(M, bernstein_values[dir]);
                }
                // apply the Bezier extraction operator for the functions on this element -- end
                //--------------------------------------------------------------------------------


                //--------------------------------------------------------------------------------
                // multiply the spline 1D in order to have the multi-d value -- begin

                const Size comp_offset_i = this->comp_offset_(iComp);

                const auto &basis_flat_to_tensor = *(this->basis_functions_indexer_)(iComp);

                const int n_basis_component = basis_component_t_size.flat_size();

                for (Size basis_fid = 0 ; basis_fid < n_basis_component ; ++basis_fid)
                {
                    const TensorIndex<dim> basis_tid = basis_flat_to_tensor(basis_fid);

                    Real value = 1.0;
                    for (int dir = 0 ; dir < dim ; ++dir)
                        value *= bspline_basis[dir][basis_tid[dir]];

                    derivatives_phi_hat_ipt[basis_fid+comp_offset_i](iComp) = value;
                }
                // multiply the spline 1D in order to have the multi-d value -- end
                //--------------------------------------------------------------------------------

            } // end iComp loop
        } // end pt_id loop

        if (this->space_->homogeneous_range_)
        {
            const auto n_basis = this->space_->get_num_basis_per_element(0);
            for (int comp = 1; comp < Space::n_components; ++comp)
            {
                const Size offset = this->comp_offset_(comp);
                for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
                {
                    const auto values_phi_hat_copy_from = derivatives_phi_hat.get_function_view(basis_i);
                    auto values_phi_hat_copy_to = derivatives_phi_hat.get_function_view(offset+basis_i);

                    for (int qp = 0; qp < n_points; ++qp)
                        values_phi_hat_copy_to[qp](comp) = values_phi_hat_copy_from[qp](0);

                } //end loop basis_i
            } // end loop comp
        } // end if (space_->homogeneous_range_)

    } // end if (deriv_order == 0)
    else
    {

        Assert(deriv_order >= 1, ExcLowerRange(deriv_order,1));

        using DerSymmMngr_t = DerivativeSymmetryManager<dim,deriv_order>;
        DerSymmMngr_t derivative_symmetry_manager;
        const auto &derivatives_flat_id_evaluate = derivative_symmetry_manager.get_entries_flat_id_evaluate();
        const auto &derivatives_flat_id_copy_to = derivative_symmetry_manager.get_entries_flat_id_copy_to();
        const auto &derivatives_flat_id_copy_from = derivative_symmetry_manager.get_entries_flat_id_copy_from();

        const auto &derivatives_tensor_id = derivative_symmetry_manager.get_entries_tensor_id();

        const Size n_derivatives_eval = DerSymmMngr_t::num_entries_eval;
        const Size n_derivatives_copy = DerSymmMngr_t::num_entries_copy;

        const array<Real,dim> elem_lengths = CartesianGridElement<dim>::get_coordinate_lengths();

        for (Size pt_id = 0 ; pt_id < n_points ; ++pt_id)
        {
            const Point<dim> point = points[pt_id];

            auto derivatives_phi_hat_ipt = derivatives_phi_hat.get_point_view(pt_id);

#ifndef NDEBUG
            for (int dir = 0 ; dir < dim ; ++dir)
                Assert(point[dir] >= 0.0 && point[dir] <= 1.0,
                ExcMessage("Evaluation point " + std::to_string(pt_id) + " not in the unit-domain."));
#endif

            for (int iComp = 0; iComp < this->space_->num_active_components_; ++iComp)
            {
                //------------------------------------------------------------------------------
                // evaluation of the values/derivarives of the 1D Bernstein polynomials -- begin
                array<array<boost::numeric::ublas::vector<Real>,dim>,deriv_order+1> bernstein_values;
                const TensorSize<dim> basis_component_t_size = this->n_basis_direction_(iComp);

                for (int order = 0 ; order <= deriv_order ; ++order)
                {
                    for (int dir = 0 ; dir < dim ; ++dir)
                    {
                        const Real scaling_coef = pow(1.0/elem_lengths[dir],order);

                        const int n_basis_1D = basis_component_t_size(dir);
                        const int degree = n_basis_1D - 1 ;
                        bernstein_values[order][dir] =
                        scaling_coef * BernsteinBasis::derivative(order,degree,point[dir]);
                    }
                }
                // evaluation of the values/derivarives of the 1D Bernstein polynomials -- end
                //------------------------------------------------------------------------------



                //--------------------------------------------------------------------------------
                // apply the Bezier extraction operator for the functions on this element -- begin
                const auto &bezier_op_comp = bezier_op(iComp);

                array<array<boost::numeric::ublas::vector<Real>,dim>,deriv_order+1> bspline_basis;
                for (int order = 0 ; order <= deriv_order ; ++order)
                {
                    for (int dir = 0 ; dir < dim ; ++dir)
                    {
                        const auto &M = *(bezier_op_comp[dir]);

                        bspline_basis[order][dir] = prec_prod(M, bernstein_values[order][dir]);
                    }
                }
                // apply the Bezier extraction operator for the functions on this element -- end
                //--------------------------------------------------------------------------------



                //--------------------------------------------------------------------------------
                // multiply the spline 1D in order to have the multi-d value -- begin
                const Size comp_offset_i = this->comp_offset_(iComp);

                const auto &basis_flat_to_tensor = *(this->basis_functions_indexer_)(iComp);

                const int n_basis_component = basis_component_t_size.flat_size();

                for (Size basis_fid = 0 ; basis_fid < n_basis_component ; ++basis_fid)
                {

                    const TensorIndex<dim> basis_tid = basis_flat_to_tensor(basis_fid);

                    auto &deriv = derivatives_phi_hat_ipt[basis_fid+comp_offset_i];

                    //--------------------------------------------------------------------------------
                    for (int entry_id = 0; entry_id < n_derivatives_eval; ++entry_id)
                    {
                        const auto entry_eval_fid = derivatives_flat_id_evaluate[entry_id];
                        const auto entry_eval_tid = derivatives_tensor_id[entry_eval_fid];

                        // from the entry_tid we get the right derivative order
                        TensorIndex<dim> deriv_order_tid; // [0,0,..,0] tensor index
                        for (int i = 0; i < deriv_order; ++i)
                            ++(deriv_order_tid[entry_eval_tid[i]]);

                        //TODO: remove this if!!! (Maybe re-think about the BSplineSpace for dim==0)
                        if (dim > 0)
                        {
                            Real value = bspline_basis[deriv_order_tid[0]][0][basis_tid[0]];
                            for (int dir = 1 ; dir < dim ; ++dir)
                                value *= bspline_basis[deriv_order_tid[dir]][dir][basis_tid[dir]];

                            deriv(entry_eval_fid)(iComp) = value;
                        }
                        else
                        {
                            deriv(entry_eval_fid)(iComp) = 1.0;
                        }

                    } // end entry_id loop
                    //--------------------------------------------------------------------------------


                    //--------------------------------------------------------------------------------
                    // here we copy the computed quantities to the symmetric part of the tensor
                    for (int entry_id = 0; entry_id < n_derivatives_copy; ++entry_id)
                        deriv(derivatives_flat_id_copy_to[entry_id])(iComp) =
                            deriv(derivatives_flat_id_copy_from[entry_id])(iComp);
                    //--------------------------------------------------------------------------------

                } // end basis_fid loop
                // multiply the spline 1D in order to have the multi-d value -- end
                //--------------------------------------------------------------------------------

            } // end iComp loop

        } // end pt_id loop

        if (this->space_->homogeneous_range_)
        {
            const Size n_ders = Derivative<deriv_order>::size;
            const auto n_basis = this->space_->get_num_basis_per_element(0);
            for (int comp = 1; comp < Space::n_components; ++comp)
            {
                const Size offset = this->comp_offset_(comp);
                for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
                {
                    const auto derivatives_phi_hat_copy_from = derivatives_phi_hat.get_function_view(basis_i);
                    auto derivatives_phi_hat_copy_to = derivatives_phi_hat.get_function_view(offset+basis_i);
                    for (int qp = 0; qp < n_points; ++qp)
                    {
                        const auto &values_0 = derivatives_phi_hat_copy_from[qp];
                        auto &values = derivatives_phi_hat_copy_to[qp];

                        for (int der = 0; der < n_ders; ++der)
                            values(der)(comp) = values_0(der)(0);
                    } // end loop qp
                } //end loop basis_i
            } // end loop comp

        } // end if (space_->homogeneous_range_)

    } // end if (deriv_order != 0)

    return derivatives_phi_hat;
}

template <int dim, int range, int rank>
void
BSplineElementAccessor<dim, range, rank>::
print_info(LogStream &out,const VerbosityLevel verbosity_level) const
{
    using std::endl;

    const std::string tab = "   ";


    out << "BSplineElementAccessor<" << dim << "," << range << "," << rank << "> info:" << endl;
    out.push(tab);

    CartesianGridElementAccessor<dim>::print_info(out,verbosity_level);

    if (contains(verbosity_level,VerbosityLevel::debug))
    {
        out << "Element cache memory address = " << &this->elem_values_ << endl;
        this->elem_values_.flags_handler_.print_info(out);

        for (int i = 0 ; i < n_faces ; ++i)
        {
            out << "Face[" << i << "] cache memory address = " << &this->face_values_[i] << endl;
            this->face_values_[i].flags_handler_.print_info(out);
        }
    }

    out.pop();

}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_element_accessor.inst>


