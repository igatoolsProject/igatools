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
//TODO(pauletti, Mar 3, 2014): rename dim_domain to dim and deriv_order to order
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



template <int dim_domain, int dim_range, int rank>
BSplineElementAccessor<dim_domain, dim_range, rank>::
BSplineElementAccessor(const Space_t &space, const int index)
    :
    CartesianGridElementAccessor<dim_domain>(*(space.get_grid()), index),
    space_(&space)
{}



//TODO: inline this
template <int dim_domain, int dim_range, int rank>
int
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_num_basis() const
{
    return (this->space_->get_num_basis_per_element());
}


template <int dim_domain, int dim_range, int rank>
int
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_component_num_basis(const int i) const
{
    const auto &degree_comp = this->space_->get_degree()(i);
    int component_num_basis = 1;
    for (int j = 0; j < dim_domain; ++j)
        component_num_basis *= degree_comp[j] + 1;

    return component_num_basis;
}



template <int dim_domain, int dim_range, int rank>
vector<Index> const &
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_local_to_global() const
{
    return space_->element_global_dofs_[this->get_flat_index()];
}



template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
init_values(const ValueFlags fill_flag,
            const Quadrature<dim_domain> &quad)
{
    Assert((fill_flag|admisible_flag) == admisible_flag,
           typename CartesianGridElementAccessor<dim_domain>::ExcFillFlagNotSupported(admisible_flag, fill_flag));

    // initalizing the cache of the CartesianGridElementAccessor
    {
        ValueFlags grid_flag = ValueFlags::none;
        if (contains(fill_flag , ValueFlags::point))
            grid_flag |= ValueFlags::point;
        if (contains(fill_flag , ValueFlags::w_measure))
            grid_flag |= ValueFlags::w_measure;
        if (contains(fill_flag , ValueFlags::face_point))
            grid_flag |= ValueFlags::face_point;
        if (contains(fill_flag , ValueFlags::face_w_measure))
            grid_flag |= ValueFlags::face_w_measure;
        CartesianGridElementAccessor<dim_domain>::init_values(grid_flag,quad);
    }

    auto f_flag = fill_flag;
    if (contains(f_flag, ValueFlags::divergence))
        f_flag |= ValueFlags::gradient;

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
    reset_element_cache(f_flag, quad);
}



template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
init_face_values(const Index face_id,
                 const ValueFlags fill_flag,
                 const Quadrature<dim_domain-1> &quad)
{
    AssertThrow(false,ExcNotImplemented()) ;
}



template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
reset_global_cache()
{
    values_1d_data_.reset();

    for (int f = 0; f < n_faces; ++f)
        values_1d_faces_[f].reset();
}

template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
reset_univariate_cache(const Quadrature<dim_domain> &quad, const int max_der)
{
    Assert(values_1d_data_.use_count() < 2,
           ExcMessage("Resetting a shared cache, use force if this is what you Really want."));

    if (values_1d_data_.use_count() == 0)
        values_1d_data_= make_shared<UniformQuadCache>();

    values_1d_data_->reset(*space_, quad, max_der);

    for (int f = 0; f < n_faces; ++f)
    {
        Assert(values_1d_faces_[f].use_count() < 2,
               ExcMessage("Resetting a shared cache, use force if this is what you Really want."));

        if (values_1d_faces_[f].use_count() == 0)
            values_1d_faces_[f]= make_shared<GlobalFaceCache>();
        values_1d_faces_[f]->reset(*space_, quad, f, max_der);
    }
}



template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
reset_element_cache(const ValueFlags fill_flag,
                    const Quadrature<dim_domain> &quad)
{
    //--------------------------------------------------------------------------
    StaticMultiArray<TensorSize<dim_domain>, dim_range, rank> n_basis_direction;
    for (int i = 0; i < Space_t::n_components; ++i)
        for (int j=0; j<dim_domain; ++j)
            n_basis_direction(i)[j] = space_->degree_(i)[j]+1;
    //--------------------------------------------------------------------------


    elem_values_.reset(fill_flag, n_basis_direction, quad);

    Index face_id = 0 ;
    const auto face_fill_flag = get_face_flags(fill_flag) ;

//    if ( !contains(face_fill_flag, ValueFlags::none) )
    if (face_fill_flag != ValueFlags::none)
        for (auto& face_value : face_values_)
            face_value.reset(face_id++, face_fill_flag, n_basis_direction, quad);
}




template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
ElementValuesCache::
reset(const ValueFlags fill_flag,
      const StaticMultiArray<TensorSize<dim_domain>, dim_range, rank> &n_basis_direction,
      const Quadrature<dim_domain> &quad)
{
    ValuesCache::reset(fill_flag, n_basis_direction,quad.get_num_points_direction());
}


template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
ValuesCache::
reset(const ValueFlags fill_flag,
      const StaticMultiArray<TensorSize<dim_domain>, dim_range, rank> &n_basis_direction,
      const TensorSize<dim_domain> &n_points_direction)
{
    this->size_.reset(n_basis_direction,
                      n_points_direction);


    //--------------------------------------------------------------------------
    // computing the total number of basis functions and the number of evaluation points
    const int total_n_points = this->size_.n_points_direction_.flat_size();

    int total_n_basis = 0;
    for (int i = 0; i < Space_t::n_components; ++i)
        total_n_basis += this->size_.n_basis_direction_(i).flat_size();

    Assert(total_n_points > 0, ExcLowerRange(total_n_points,1));
    Assert(total_n_basis > 0, ExcLowerRange(total_n_basis,1));
    //--------------------------------------------------------------------------


//    using std::cout;
//    using std::endl;

    //--------------------------------------------------------------------------
    // resizing the containers for the basis functions

    Assert(contains(fill_flag, ValueFlags::none),
           ExcMessage("nothing to reset"));

    int max_der_order = -1;
    if (contains(fill_flag , ValueFlags::value))
    {
//      cout << "BSplineElementAccessor::ValuesData::reset() ---> value" <<endl;
        fill_values_ = true;

        if (D0phi_hat_.get_num_points() != total_n_points ||
            D0phi_hat_.get_num_functions() != total_n_basis)
        {
            D0phi_hat_.resize(total_n_basis,total_n_points);
            D0phi_hat_.zero();
        }

        if (phi_hat_.get_num_points() != total_n_points ||
            phi_hat_.get_num_functions() != total_n_basis)
        {
            phi_hat_.resize(total_n_basis,total_n_points);
            phi_hat_.zero();
        }

        max_der_order=std::max(max_der_order,0);

    }
    else
    {
        fill_values_ = false;
        D0phi_hat_.clear();
        phi_hat_.clear();
    }


    if (contains(fill_flag , ValueFlags::divergence))
    {
        fill_divs_ = true;
        if (div_phi_hat_.get_num_points() != total_n_points ||
            div_phi_hat_.get_num_functions() != total_n_basis)
        {
            div_phi_hat_.resize(total_n_basis,total_n_points);
            div_phi_hat_.zero();
        }

        max_der_order=std::max(max_der_order,1);
    }
    else
    {
        fill_divs_ = false;
        D1phi_hat_.clear();
    }


    if (contains(fill_flag , ValueFlags::gradient))
    {
//      cout << "BSplineElementAccessor::ValuesData::reset() ---> gradient" <<endl;
        fill_gradients_ = true;
        if (D1phi_hat_.get_num_points() != total_n_points ||
            D1phi_hat_.get_num_functions() != total_n_basis)
        {
            D1phi_hat_.resize(total_n_basis,total_n_points);
            D1phi_hat_.zero();
        }

        max_der_order=std::max(max_der_order,1);
    }
    else
    {
        fill_gradients_ = false;
        D1phi_hat_.clear();
    }

    if (contains(fill_flag , ValueFlags::hessian))
    {
        fill_hessians_ = true;
        if (D2phi_hat_.get_num_points() != total_n_points ||
            D2phi_hat_.get_num_functions() != total_n_basis)
        {
            D2phi_hat_.resize(total_n_basis,total_n_points);
            D2phi_hat_.zero();
        }

        max_der_order=std::max(max_der_order,2);
    }
    else
    {
        fill_hessians_ = false;
        D2phi_hat_.clear();
    }

//    Assert((max_der_order>=0), ExcMessage("Not a right ValueFlag"));
    //--------------------------------------------------------------------------


    this->set_initialized(true);
//    this->set_filled(false);
}

template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
ValuesCache::
get_values() const -> const ValueTable<Value> &
{
    return phi_hat_;
}

template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
ValuesCache::
get_gradients() const -> const ValueTable<Derivative<1>> &
{
    return D1phi_hat_;
}


template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
ValuesCache::
get_hessians() const -> const ValueTable<Derivative<2>> &
{
    return D2phi_hat_;
}

template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
ValuesCache::
get_divergences() const -> const ValueTable<Div> &
{
    return div_phi_hat_;
}


template <int dim_domain, int dim_range, int rank>
ValueFlags
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_flags(const ValueFlags fill_flag) const
{

    ValueFlags face_fill_flag = ValueFlags::none ;

    if (contains(fill_flag , ValueFlags::face_value))
        face_fill_flag |= ValueFlags::value ;

    if (contains(fill_flag , ValueFlags::face_divergence))
        face_fill_flag |= ValueFlags::divergence ;

    if (contains(fill_flag , ValueFlags::face_gradient))
        face_fill_flag |= ValueFlags::gradient ;

    if (contains(fill_flag , ValueFlags::face_hessian))
        face_fill_flag |= ValueFlags::hessian ;

    return face_fill_flag ;
}


template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
FuncPointSize::
reset(StaticMultiArray<TensorSize<dim_domain>, dim_range, rank> n_basis_direction,
      TensorSize<dim_domain> n_points_direction)
{
    n_basis_direction_ = n_basis_direction;
    n_points_direction_ = n_points_direction;

    //--------------------------------------------------------------------------
    // creating the objects for fast conversion from flat-to-tensor indexing
    // in practice is an hash-table from flat to tensor indices
    using Indexer = CartesianProductIndexer<dim_domain>;

    const int n_components = StaticMultiArray<TensorIndex<dim_domain>, dim_range, rank>::n_entries;
    for (int iComp = 0; iComp < n_components; ++iComp)
        basis_functions_indexer_(iComp) = shared_ptr<Indexer>(new Indexer(n_basis_direction_(iComp)));

    points_indexer_ = shared_ptr<Indexer>(new Indexer(n_points_direction_));
    //--------------------------------------------------------------------------
}



template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
FaceValuesCache::
reset(const Index face_id,
      const ValueFlags fill_flag,
      const StaticMultiArray<TensorSize<dim_domain>, dim_range, rank> &n_basis_direction,
      const Quadrature<dim_domain> &quad1)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    const auto quad = quad1.collapse_to_face(face_id);

    ValuesCache::reset(fill_flag, n_basis_direction,quad.get_num_points_direction());
}



template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
FaceValuesCache::
reset(const Index face_id,
      const ValueFlags fill_flag,
      const StaticMultiArray<TensorSize<dim_domain>, dim_range, rank> &n_basis_direction,
      const Quadrature<dim_domain-1> &quad)
{
    AssertThrow(false,ExcNotImplemented()) ;
}



template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
GlobalFaceCache::
reset(const Space_t &space,
      const Quadrature<dim_domain> &quad1,
      const Index face_id,
      int max_der)
{
    const auto quad = quad1.collapse_to_face(face_id);
    max_deriv_order = max_der;
    // resizing the structures for the one dimensional splines
    const int n_active_components = space.num_active_components_;
    const int max_der_plus_one = max_deriv_order + 1;

    for (int iComp = 0; iComp < n_active_components; ++iComp)
    {
        splines1d_cache_data_(iComp).resize(max_der_plus_one);
    }
    this->set_initialized(true);

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
    BasisValues1d bernstein_values(max_der_plus_one);
    for (int iComp = 0; iComp < n_active_components; ++iComp)
    {
        BasisValues1d &basis = splines1d_cache_data_(iComp);

        const int jDim = UnitElement<dim_domain>::face_constant_direction[face_id];//const_dir;
        {
            const int degree = degree_(iComp)[ jDim ];
            const std::vector<Real> &pt_coords = eval_points.get_data_direction(jDim);

            // fill values and derivatives of the Bernstein's polynomials at
            // quad points in [0,1]

            for (int deriv_order = 0; deriv_order < max_der_plus_one; ++deriv_order)
                bernstein_values[ deriv_order ] =
                    BernsteinBasis::derivative(deriv_order, degree, pt_coords);


            const auto &bez_iComp_jDim = bezier_op_(iComp).get_data_direction(jDim);
            const auto &lengths_jDim = lengths.get_data_direction(jDim);

            // compute the one dimensional B-splines at quad point on the reference interval
            const int interval = 0;
            {
                const boost::numeric::ublas::matrix<Real> &M = *(bez_iComp_jDim[interval]);
                const Real one_div_size = iga::Real(1.0) / lengths_jDim[interval];
                for (int deriv_order = 0; deriv_order < max_der_plus_one; ++deriv_order)
                {
                    const Real scaling_coef =
                        std::pow(one_div_size, deriv_order);
                    basis[ deriv_order ] = scaling_coef *
                                           prod(M, bernstein_values[ deriv_order ]);
                }
            }
        } // end loop jDim
    } // end loop iComp


    // assign the basis1D data to the proper component/interval through their memory address
    for (int iComp = 0; iComp < Space_t::n_components; iComp++)
    {
        const int active_comp = space.map_component_to_active_data_(iComp);
        splines1d_cache_(iComp) = &(splines1d_cache_data_(active_comp));
    }

    this->set_filled(true);
}



template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
fill_values()
{
    Assert(values_1d_data_->is_filled(), ExcNotInitialized());
    Assert(elem_values_.is_initialized(), ExcNotInitialized());

    CartesianGridElementAccessor<dim_domain>::fill_values();

    const auto  &element_tensor_id = this->get_tensor_index();
    StaticMultiArray<array<const BasisValues1d *, dim_domain>, dim_range, rank>
    elem_univariate_values;
    for (int iComp=0; iComp<space_->num_active_components_; ++iComp)
    {
        auto &univariate_values = values_1d_data_->splines1d_cache_(iComp);
        for (int i = 0; i < dim_domain; ++i)
            elem_univariate_values(iComp)[i] = univariate_values.get_data_direction(i)[element_tensor_id[i]];
    }

    if (elem_values_.fill_values_)
    {
        evaluate_bspline_derivatives<0>(elem_values_.size_,
                                        elem_univariate_values,
                                        elem_values_.D0phi_hat_);
        auto phi_hat = elem_values_.phi_hat_.begin();
        for (auto &D0phi_hat : elem_values_.D0phi_hat_)
        {
            *phi_hat = (D0phi_hat)(0);
            ++phi_hat;
        }
    }

    if (elem_values_.fill_gradients_)
        evaluate_bspline_derivatives<1>(elem_values_.size_,
                                        elem_univariate_values,
                                        elem_values_.D1phi_hat_);

    if (elem_values_.fill_hessians_)
        evaluate_bspline_derivatives<2>(elem_values_.size_,
                                        elem_univariate_values,
                                        elem_values_.D2phi_hat_);

    if (elem_values_.fill_divs_)
    {
        auto D1  = elem_values_.D1phi_hat_.begin();
        auto div = elem_values_.div_phi_hat_.begin();
        auto end = elem_values_.D1phi_hat_.end();
        for (; D1 != end; ++D1, ++div)
            *div = trace(*D1);
    }

    elem_values_.set_filled(true);
}



template <int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
fill_face_values(const Index face_id)
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));

    auto &face_value = face_values_[face_id] ;

    Assert(face_value.is_initialized(), ExcNotInitialized());

    CartesianGridElementAccessor<dim_domain>::fill_face_values(face_id);
    const int const_dir = UnitElement<dim_domain>::face_constant_direction[face_id];

    const auto &element_tensor_id = this->get_tensor_index();
    StaticMultiArray<array<const BasisValues1d *,dim_domain>,dim_range,rank>
    elem_univariate_values;
    for (int iComp=0; iComp<space_->num_active_components_; ++iComp)
    {
        for (int i = 0; i < int (dim_domain); ++i)
        {
            if (i==const_dir)
                elem_univariate_values(iComp)[i] = values_1d_faces_[face_id]->splines1d_cache_(iComp);
            else
                elem_univariate_values(iComp)[i] =
                    values_1d_data_->splines1d_cache_(iComp).get_data_direction(i)[element_tensor_id[i]];
        }
    }
    if (face_value.fill_values_)
    {
        evaluate_bspline_derivatives<0>(face_value.size_,
                                        elem_univariate_values,
                                        face_value.D0phi_hat_);
        auto phi_hat = face_values_[face_id].phi_hat_.begin();
        for (const auto &D0phi_hat : face_value.D0phi_hat_)
        {
            *phi_hat = (D0phi_hat)(0);
            ++phi_hat;
        }
    }

    if (face_value.fill_gradients_)
        evaluate_bspline_derivatives<1>(face_value.size_,
                                        elem_univariate_values,
                                        face_value.D1phi_hat_);

    if (face_value.fill_hessians_)
        evaluate_bspline_derivatives<2>(face_value.size_,
                                        elem_univariate_values,
                                        face_value.D2phi_hat_);

    face_value.set_filled(true);
}



template < int dim_domain, int dim_range, int rank>
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
UniformQuadCache::
reset(const Space_t &space,
      const Quadrature<dim_domain> &quad,
      const int max_der)
{
    //------------------------------------------------------------------------------------------
    max_deriv_order = max_der;
    // resizing the structures for the one dimensional splines
    const int n_active_components = space.num_active_components_;
    const auto &n_elem = space.get_grid()->get_num_elements_dim();
    const int max_der_plus_one = max_deriv_order + 1;

    for (int iComp = 0; iComp < n_active_components; ++iComp)
    {
        Assert(n_elem == space.bezier_op_(iComp).tensor_size(),
               ExcMessage("Not same size"));
        splines1d_cache_data_(iComp).resize(n_elem);
        for (int i=0; i<dim_domain; ++i)
            for (int j=0; j<n_elem[i]; ++j)
                splines1d_cache_data_(iComp).entry(i,j).resize(max_der_plus_one);
    }


    for (int iComp = 0; iComp < Space_t::n_components; ++iComp)
    {
        Assert(n_elem == splines1d_cache_data_(space.map_component_to_active_data_(iComp)).tensor_size(),
               ExcMessage("Not same size"));
        splines1d_cache_(iComp).resize(n_elem);
    }
    this->set_initialized(true);
    //------------------------------------------------------------------------------------------


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
    BasisValues1d bernstein_values(max_der_plus_one);
    for (int iComp = 0; iComp < n_active_components; ++iComp)
    {
        for (int jDim = 0; jDim < dim_domain; ++jDim)
        {
            const int num_intevals = n_elem[jDim];
            const int degree = degree_(iComp)[ jDim ];
            const vector<Real> &pt_coords = eval_points.get_data_direction(jDim);

            // fill values and derivatives of the Bernstein's polynomials at
            // quad points in [0,1]

            for (int deriv_order = 0; deriv_order < max_der_plus_one; ++deriv_order)
                bernstein_values[ deriv_order ] =
                    BernsteinBasis::derivative(deriv_order, degree, pt_coords);

            const auto &bez_iComp_jDim = bezier_op_(iComp).get_data_direction(jDim);
            const auto &lengths_jDim = lengths.get_data_direction(jDim);

            // compute the one dimensional B-splines at quad point on the reference interval
            for (int interval = 0; interval < num_intevals; interval++)
            {
                const auto &M = *(bez_iComp_jDim[interval]);
                const Real one_div_size = 1.0 / lengths_jDim[interval];
                BasisValues1d &basis = splines1d_cache_data_(iComp).entry(jDim,interval);

                for (int deriv_order = 0; deriv_order < max_der_plus_one; ++deriv_order)
                {
                    const Real scaling_coef = std::pow(one_div_size, deriv_order);
                    basis[ deriv_order ] = scaling_coef * prod(M, bernstein_values[ deriv_order ]);
                }
            }
        } // end loop jDim
    } // end loop iComp
    //------------------------------------------------------------------------------------------


    //------------------------------------------------------------------------------------------
    // assign the basis1D data to the proper component/interval through their memory address
    for (int iComp = 0; iComp < Space_t::n_components; iComp++)
    {
        const int active_comp = space.map_component_to_active_data_(iComp);

        const auto &spline1d_active_data = splines1d_cache_data_(active_comp);

        const auto n_intervals_multi_d = spline1d_active_data.tensor_size();

        for (int jDim = 0; jDim < dim_domain; jDim++)
        {
            const Size n_intervals_1d = n_intervals_multi_d[jDim];

            const auto &data = spline1d_active_data.get_data_direction(jDim);

            for (Size i = 0; i < n_intervals_1d; ++i)
                splines1d_cache_(iComp).entry(jDim,i) = &(data[i]);
        }
    }

    this->set_filled(true);
    //-------------------------------------------------------------------------

}



template <int dim_domain, int dim_range, int rank>
template < int deriv_order >
void
BSplineElementAccessor<dim_domain, dim_range, rank>::
evaluate_bspline_derivatives(const FuncPointSize &size,
                             StaticMultiArray<std::array<const BasisValues1d *, dim_domain>, dim_range, rank> &elem_values,
                             ValueTable< Derivative<deriv_order> > &derivatives_phi_hat) const
{
    Assert(derivatives_phi_hat.size() > 0, ExcEmptyObject());
    Assert(derivatives_phi_hat.get_num_functions() == this->get_num_basis(),
           ExcDimensionMismatch(derivatives_phi_hat.get_num_functions(),this->get_num_basis()));
    Assert(derivatives_phi_hat.get_num_points() ==
           size.n_points_direction_.flat_size(),
           ExcDimensionMismatch(derivatives_phi_hat.get_num_points(),
                                size.n_points_direction_.flat_size()));

    /*
     * This code computes any order of derivatives for a multivariate
     * B-spline on the current element
     * We use the formula
     * \partial_(\alpha_1,...,\alpha_n) B(qp) = \Pi d^{\alpha_i} B_i(qp_i)
     */


    Assert(rank < 2, ExcMessage("For rank> 1 the basis function are not implemented/tested."));

    array<int, Space_t::n_components> comp_offset;
    comp_offset[0] = 0;
    for (int i = 1; i < Space_t::n_components; ++i)
        comp_offset[i]= comp_offset[i-1] + size.n_basis_direction_(i).flat_size();

    typedef Derivatives<dim_domain,dim_range,rank,deriv_order> der_t;

    TensorIndex<dim_domain> deriv_order_tensor_id;


#ifndef NOT_OPTIMIZED
    typedef DerivativeSymmetryManager<dim_domain,deriv_order> DerSymmMngr_t;
    DerSymmMngr_t derivative_symmetry_manager;
    const auto &derivatives_flat_id_evaluate = derivative_symmetry_manager.get_entries_flat_id_evaluate();
    const auto &derivatives_flat_id_copy_to = derivative_symmetry_manager.get_entries_flat_id_copy_to();
    const auto &derivatives_flat_id_copy_from = derivative_symmetry_manager.get_entries_flat_id_copy_from();

    const auto &derivatives_tensor_id = derivative_symmetry_manager.get_entries_tensor_id();

    const int n_derivatives_eval = DerSymmMngr_t::num_entries_eval;
    const int n_derivatives_copy = DerSymmMngr_t::num_entries_copy;
#endif


    const int num_points = size.n_points_direction_.flat_size();

    const auto &points_indexer = *(size.points_indexer_);

//    LogStream out;
//    out << &derivatives_phi_hat << std::endl ;
//    derivatives_phi_hat.print_info(out) ;
    for (int iComp=0; iComp < space_->num_active_components_; ++iComp)
    {
        const auto &splines1d_direction = elem_values(iComp);


        const int n_basis = get_component_num_basis(iComp);
        Assert(n_basis == size.n_basis_direction_(iComp).flat_size(), ExcMessage("different sizes"));

        const auto &functions_indexer = *(size.basis_functions_indexer_(iComp));

        const int comp_offset_i = comp_offset[iComp];

#ifdef NOT_OPTIMIZED
        // NOT_OPTIMIZED branch

        const int n_partial_ders = der_t::size;

        for (int func_flat_id = 0; func_flat_id < n_basis; ++func_flat_id)
        {
            const auto func_tensor_id = functions_indexer.get_tensor_index(func_flat_id);

            auto derivatives_phi_hat_ifn = derivatives_phi_hat.get_function_view(comp_offset_i+func_flat_id);

            for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
            {

                const auto &point_tensor_id = points_indexer.get_tensor_index(point_flat_id);

                der_t &derivative = derivatives_phi_hat_ifn[point_flat_id];

                for (int entry_flat_id = 0; entry_flat_id < n_partial_ders; ++entry_flat_id)
                {
                    const auto entry_tensor_id = derivative.flat_to_tensor_index(entry_flat_id);

                    // from the entry id we compute the right derivative order
                    for (int i = 0; i < dim_domain; ++i)
                        deriv_order_tensor_id[i] = 0;

                    for (int i = 0; i < deriv_order; ++i)
                        ++(deriv_order_tensor_id[entry_tensor_id[i]]);

                    // Main computation
                    Real partial_der = Real(1.0);
                    for (int i = 0; i < dim_domain; ++i)
                    {
                        const auto &splines1d = *(splines1d_direction[i]);
                        partial_der *= splines1d[deriv_order_tensor_id[i]](func_tensor_id[i],point_tensor_id[i]);
                    }
                    derivative(entry_flat_id)(iComp) = partial_der;

                } // end entry_flat_id loop
            } // end point_flat_id loop
        } //end func_flat_id loop

        // end NOT_OPTIMIZED branch
#else
        // OPTIMIZED branch

        for (int func_flat_id = 0; func_flat_id < n_basis; ++func_flat_id)
        {
            const auto func_tensor_id = functions_indexer.get_tensor_index(func_flat_id);

//            out << "comp_offset_i=" << comp_offset_i
//                << "  func_flat_id=" << func_flat_id
//                << "  func_id=" << comp_offset_i+func_flat_id<<std::endl;
            auto derivatives_phi_hat_ifn = derivatives_phi_hat.get_function_view(comp_offset_i+func_flat_id);

            for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
            {
                const auto &point_tensor_id = points_indexer.get_tensor_index(point_flat_id);

                der_t &derivative = derivatives_phi_hat_ifn[point_flat_id];

                for (int entry_id = 0; entry_id < n_derivatives_eval; ++entry_id)
                {
                    const int entry_flat_id = derivatives_flat_id_evaluate[entry_id];
                    const auto entry_tensor_id = derivatives_tensor_id[entry_flat_id];

                    // from the entry id we compute the right derivative order
                    for (int i = 0; i < dim_domain; ++i)
                        deriv_order_tensor_id[i] = 0;

                    for (int i = 0; i < deriv_order; ++i)
                        ++(deriv_order_tensor_id[entry_tensor_id[i]]);

                    if (dim_domain != 0)
                    {
                        // Main computation
                        const auto &splines1d = *(splines1d_direction[0]);
                        Real partial_der = splines1d[deriv_order_tensor_id[0]](func_tensor_id[0], point_tensor_id[0]);
                        for (int i = 1; i < dim_domain; ++i)
                        {
                            const auto &splines1d = *(splines1d_direction[i]);
                            partial_der *= splines1d[deriv_order_tensor_id[i]](func_tensor_id[i], point_tensor_id[i]);
                        }
                        derivative(entry_flat_id)(iComp) = partial_der;
                    }
                    else
                    {
                        derivative(entry_flat_id)(iComp) = 1.0;
                    }

                } // end entry_id loop

                // here we copy the computed quantities to the symmetric part of the tensor
                for (int entry_id = 0; entry_id < n_derivatives_copy; ++entry_id)
                    derivative(derivatives_flat_id_copy_to[entry_id])(iComp) = derivative(derivatives_flat_id_copy_from[entry_id])(iComp);

            } // end point_flat_id loop
        } //end func_flat_id loop

        //end OPTIMIZED branch
#endif
    } // end iComp loop

    if (space_->homogeneous_range_)
    {
        const int n_ders = Derivative<deriv_order>::size;
        const auto n_basis = space_->get_component_num_basis_per_element(0);
        for (int iComp = 1; iComp < Space_t::n_components; ++iComp)
        {
            const int offset = comp_offset[iComp];
            for (int basis_i = 0; basis_i < n_basis;  ++basis_i)
            {
                const auto derivatives_phi_hat_copy_from = derivatives_phi_hat.get_function_view(basis_i);
                auto derivatives_phi_hat_copy_to = derivatives_phi_hat.get_function_view(offset+basis_i);
                for (int qp = 0; qp < num_points; ++qp)
                {
                    const der_t &values_0 = derivatives_phi_hat_copy_from[qp];
                    der_t &values = derivatives_phi_hat_copy_to[qp];
                    for (int j = 0; j < n_ders; ++j)
                    {
                        values(j)(iComp) = values_0(j)(0);
                    }
                }
            }
        } // end loop iComp
    } // end if (space_->homogeneous_range_)
}






template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_space() const -> const Space_t *
{
    return (space_);
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_values() const -> ValueTable<Value> const &
{
    Assert(elem_values_.is_filled() == true, ExcCacheNotFilled());
    Assert(elem_values_.fill_values_, ExcCacheNotFilled());

    return elem_values_.phi_hat_;
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_values(const Index i) const -> typename ValueTable<Value>::const_view
{
    return this->get_basis_values().get_function_view(i);
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_divergences() const -> ValueTable<Div> const &
{
    Assert(elem_values_.is_filled() == true, ExcCacheNotFilled());
    Assert(elem_values_.fill_divs_, ExcCacheNotFilled());

    return elem_values_.div_phi_hat_;
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_divergences(const Index i) const -> typename ValueTable<Div>::const_view
{
    return this->get_basis_divergences().get_function_view(i);
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_gradients() const -> ValueTable<Derivative<1>> const &
{
    Assert(elem_values_.is_filled() == true, ExcCacheNotFilled());
    Assert(elem_values_.fill_gradients_, ExcCacheNotFilled());

    return elem_values_.D1phi_hat_;
}

template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_gradients(const Index i) const -> typename ValueTable<Derivative<1>>::const_view
{
    return this->get_basis_gradients().get_function_view(i);
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_hessians() const -> ValueTable<Derivative<2>> const &
{
    Assert(elem_values_.is_filled() == true, ExcCacheNotFilled());
    Assert(elem_values_.fill_hessians_, ExcCacheNotFilled());

    return elem_values_.D2phi_hat_;
}


template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_hessians(const Index i) const -> typename ValueTable<Derivative<2>>::const_view
{
    return this->get_basis_hessians().get_function_view(i);
}


template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_value(const Index basis, const Index qp) const -> Value const &
{
    /*
        const auto &data = elem_values_.phi_hat_;

        Assert(basis >= 0 && basis < int(data.get_num_functions()),
               ExcIndexRange(basis,0,int(data.get_num_functions())));
        Assert(qp >= 0 && qp < int(data.get_num_points()),
               ExcIndexRange(qp,0,int(data.get_num_points())));

        return data.get_function_view(basis)[qp];
     //*/

    Assert(qp >= 0 && qp < elem_values_.size_.n_points_direction_.flat_size(),
           ExcIndexRange(qp,0,elem_values_.size_.n_points_direction_.flat_size()));
    return this->get_basis_values(basis)[qp];
}


template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_divergence(const Index basis, const Index qp) const -> Div const &
{
    /*
        const auto &data = elem_values_.phi_hat_;

        Assert(basis >= 0 && basis < int(data.get_num_functions()),
               ExcIndexRange(basis,0,int(data.get_num_functions())));
        Assert(qp >= 0 && qp < int(data.get_num_points()),
               ExcIndexRange(qp,0,int(data.get_num_points())));

        return data.get_function_view(basis)[qp];
     //*/

    Assert(qp >= 0 && qp < elem_values_.size_.n_points_direction_.flat_size(),
           ExcIndexRange(qp,0,elem_values_.size_.n_points_direction_.flat_size()));
    return this->get_basis_divergences(basis)[qp];
}

template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_gradient(const Index basis, const Index qp) const -> Derivative<1> const &
{
    /*
        const auto &data = elem_values_.D1phi_hat_;

        Assert(basis >= 0 && basis < int(data.get_num_functions()),
               ExcIndexRange(basis,0,int(data.get_num_functions())));
        Assert(qp >= 0 && qp < int(data.get_num_points()),
               ExcIndexRange(qp,0,int(data.get_num_points())));

        return data.get_function_view(basis)[qp];
    //*/
    Assert(qp >= 0 && qp < elem_values_.size_.n_points_direction_.flat_size(),
           ExcIndexRange(qp,0,elem_values_.size_.n_points_direction_.flat_size()));
    return this->get_basis_gradients(basis)[qp];
}

template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_basis_hessian(const Index basis, const Index qp) const -> Derivative<2> const &
{
    /*
        const auto &data = elem_values_.D2phi_hat_;

        Assert(basis >= 0 && basis < int(data.get_num_functions()),
               ExcIndexRange(basis,0,int(data.get_num_functions())));
        Assert(qp >= 0 && qp < int(data.get_num_points()),
               ExcIndexRange(qp,0,int(data.get_num_points())));

        return data.get_function_view(basis)[qp];
    //*/
    Assert(qp >= 0 && qp < elem_values_.size_.n_points_direction_.flat_size(),
           ExcIndexRange(qp,0,elem_values_.size_.n_points_direction_.flat_size()));
    return this->get_basis_hessians(basis)[qp];
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_values(const Index face_id) const -> ValueTable<Value> const &
{
    Assert(face_id >= 0 && face_id < n_faces, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled() == true, ExcCacheNotFilled());
    Assert(face_values_[face_id].fill_values_, ExcCacheNotFilled());

    return face_values_[face_id].phi_hat_;
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_values(const Index face_id, const Index i) const -> typename ValueTable<Value>::const_view
{
    return this->get_face_basis_values(face_id).get_function_view(i);
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_divergences(const Index face_id) const -> ValueTable<Div> const &
{
    Assert(face_id >= 0 && face_id < n_faces, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled() == true, ExcCacheNotFilled());
    Assert(face_values_[face_id].fill_divs_, ExcCacheNotFilled());

    return face_values_[face_id].div_phi_hat_;
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_divergences(const Index face_id, const Index i) const -> typename ValueTable<Div>::const_view
{
    return this->get_face_basis_divergences(face_id).get_function_view(i);
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_gradients(const Index face_id) const -> ValueTable<Derivative<1>> const &
{
    Assert(face_id >= 0 && face_id < n_faces, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled() == true, ExcCacheNotFilled());
    Assert(face_values_[face_id].fill_gradients_, ExcCacheNotFilled());

    return face_values_[face_id].D1phi_hat_;
}

template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_gradients(const Index face_id, const Index i) const -> typename ValueTable<Derivative<1>>::const_view
{
    return this->get_face_basis_gradients(face_id).get_function_view(i);
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_hessians(const Index face_id) const -> ValueTable<Derivative<2>> const &
{
    Assert(face_id >= 0 && face_id < n_faces, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled() == true, ExcCacheNotFilled());
    Assert(face_values_[face_id].fill_hessians_, ExcCacheNotFilled());

    return face_values_[face_id].D2phi_hat_;
}


template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_hessians(const Index face_id, const Index i) const -> typename ValueTable<Derivative<2>>::const_view
{
    return this->get_face_basis_hessians(face_id).get_function_view(i);
}


template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_value(const Index face_id, const Index basis, const Index qp) const -> Value const &
{

    Assert(qp >= 0 && qp < elem_values_.size_.n_points_direction_.flat_size(),
           ExcIndexRange(qp,0,elem_values_.size_.n_points_direction_.flat_size()));
    return this->get_face_basis_values(face_id, basis)[qp];
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_divergence(const Index face_id, const Index basis, const Index qp) const -> Div const &
{
    Assert(qp >= 0 && qp < elem_values_.size_.n_points_direction_.flat_size(),
           ExcIndexRange(qp,0,elem_values_.size_.n_points_direction_.flat_size()));
    return this->get_face_basis_divergences(face_id, basis)[qp];
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_gradient(const Index face_id, const Index basis, const Index qp) const -> Derivative<1> const &
{
    Assert(qp >= 0 && qp < elem_values_.size_.n_points_direction_.flat_size(),
           ExcIndexRange(qp,0,elem_values_.size_.n_points_direction_.flat_size()));
    return this->get_face_basis_gradients(face_id, basis)[qp];
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
get_face_basis_hessian(const Index face_id, const Index basis, const Index qp) const -> Derivative<2> const &
{
    Assert(qp >= 0 && qp < elem_values_.size_.n_points_direction_.flat_size(),
           ExcIndexRange(qp,0,elem_values_.size_.n_points_direction_.flat_size()));
    return this->get_face_basis_hessians(face_id, basis)[qp];
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
evaluate_field(const std::vector<Real> &local_coefs) const
-> ValueVector<Value>
{
    Assert(elem_values_.is_filled() == true, ExcCacheNotFilled());
    Assert(elem_values_.fill_values_ == true, ExcInvalidState());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D0phi_hat = this->get_basis_values() ;
    Assert(D0phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D0phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D0phi_hat.evaluate_linear_combination(local_coefs) ;
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
evaluate_field_gradients(const std::vector<Real> &local_coefs) const -> ValueVector< Derivative<1> >
{
    Assert(elem_values_.is_filled() == true, ExcCacheNotFilled());
    Assert(elem_values_.fill_gradients_ == true, ExcInvalidState());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D1phi_hat = this->get_basis_gradients() ;
    Assert(D1phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D1phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D1phi_hat.evaluate_linear_combination(local_coefs) ;
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
evaluate_field_hessians(const std::vector<Real> &local_coefs) const -> ValueVector< Derivative<2> >
{
    Assert(elem_values_.is_filled() == true, ExcCacheNotFilled());
    Assert(elem_values_.fill_hessians_ == true, ExcInvalidState());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D2phi_hat = this->get_basis_hessians() ;
    Assert(D2phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D2phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D2phi_hat.evaluate_linear_combination(local_coefs) ;
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
evaluate_face_field(const Index face_id, const std::vector<Real> &local_coefs) const
-> ValueVector<Value>
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled() == true, ExcCacheNotFilled());
    Assert(face_values_[face_id].fill_values_ == true, ExcInvalidState());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D0phi_hat = this->get_face_basis_values(face_id) ;
    Assert(D0phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D0phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D0phi_hat.evaluate_linear_combination(local_coefs) ;
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
evaluate_face_field_gradients(const Index face_id, const std::vector<Real> &local_coefs) const -> ValueVector< Derivative<1> >
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled() == true, ExcCacheNotFilled());
    Assert(face_values_[face_id].fill_gradients_ == true, ExcInvalidState());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D1phi_hat = this->get_face_basis_gradients(face_id) ;
    Assert(D1phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D1phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D1phi_hat.evaluate_linear_combination(local_coefs) ;
}



template <int dim_domain, int dim_range, int rank>
auto
BSplineElementAccessor<dim_domain, dim_range, rank>::
evaluate_face_field_hessians(const Index face_id, const std::vector<Real> &local_coefs) const -> ValueVector< Derivative<2> >
{
    Assert(face_id < n_faces && face_id >= 0, ExcIndexRange(face_id,0,n_faces));
    Assert(face_values_[face_id].is_filled() == true, ExcCacheNotFilled());
    Assert(face_values_[face_id].fill_hessians_ == true, ExcInvalidState());
    Assert(this->get_num_basis() == local_coefs.size(),
    ExcDimensionMismatch(this->get_num_basis(),local_coefs.size()));

    const auto &D2phi_hat = this->get_face_basis_hessians(face_id) ;
    Assert(D2phi_hat.get_num_functions() == this->get_num_basis(),
    ExcDimensionMismatch(D2phi_hat.get_num_functions(), this->get_num_basis())) ;

    return D2phi_hat.evaluate_linear_combination(local_coefs) ;
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_element_accessor.inst>


