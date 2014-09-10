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

#include <igatools/basis_functions/bspline_uniform_quad_cache.h>
#include <igatools/basis_functions/bernstein_basis.h>
#include <igatools/utils/multi_array_utils.h>

#include <algorithm>
using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace
{
//TODO(pauletti, Sep 9, 2014): this class seems to have a more general use put in own
// file under group tensor product values utilities
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
                std::sort(tensor_id.begin(),tensor_id.end());
                std::reverse(tensor_id.begin(),tensor_id.end());
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

template<int dim_, int range_ , int rank_>
BSplineUniformQuadCache<dim_, range_, rank_>::
BSplineUniformQuadCache(shared_ptr<const Space> space,
                        const ValueFlags flag,
                        const Quadrature<dim> &quad)
    :
    GridUniformQuadCache<dim_>(space->get_grid(), flag, quad),
    space_(space),
    n_basis_(space_->get_num_all_element_basis()),
    flags_(flag),
    face_flags_(flag),
    quad_(quad),
    splines1d_(space->get_grid()->get_num_intervals(),
    		BasisValues(space->get_components_map()))
    {


	// Compute the component offsets
    comp_offset_[0] = 0;
    for (int j = 1; j < Space::n_components; ++j)
        comp_offset_[j] = comp_offset_[j-1] + n_basis_.comp_dimension[j-1];


    const auto &n_inter = space->get_grid()->get_num_intervals();
    const auto &n_points = quad.get_num_points_direction();

    // Allocate space for the BasisValues1D
    for (int dir = 0 ; dir < dim ; ++dir)
    {
        const auto &n_pts = n_points[dir];
        for (int j = 0 ; j < n_inter[dir] ; ++j)
        {
            auto &splines1d = splines1d_.entry(dir, j);
            for (auto comp : splines1d.get_active_components_id())
            	splines1d[comp].resize(n_derivatives, n_basis_[comp][dir], n_pts);
        }
    }

    /*
     * For each direction, interval and component we compute the 1D bspline
     * basis evaluate at the 1D component of the tensor product quadrature
     */
    const auto &degree      = space->get_degree();
    const auto &bezier_op   = space_->operators_;
    const auto &points      = quad_.get_points();
    const auto &lengths = this->lengths_;

    BasisValues bernstein_values(n_basis_.get_comp_map());

    for (int dir = 0 ; dir < dim ; ++dir)
    {
    	// fill values and derivatives of the Bernstein's polynomials at
        // quad points in [0,1]
    	for (auto comp : bernstein_values.get_active_components_id())
        {
            const int deg = degree[comp][dir];
            bernstein_values[comp].resize(n_derivatives, deg+1, n_points[dir]);
            const auto &pt_coords = points.get_data_direction(dir);
            for (int order = 0; order < n_derivatives; ++order)
            	bernstein_values[comp].get_derivative(order) =
            			BernsteinBasis::derivative(order, deg, pt_coords);
        }

    	const auto& inter_lengths = lengths.get_data_direction(dir);
    	for (int j = 0 ; j < n_inter[dir] ; ++j)
    	{
            auto &splines1d = splines1d_.entry(dir, j);
            for (auto comp : splines1d.get_active_components_id())
            {
                const auto &berns_values = bernstein_values[comp];
                auto &basis = splines1d[comp];
                const auto &oper = bezier_op.get_operator(comp,dir)[j];
                const Real one_div_size = 1.0 / inter_lengths[j];
                for (int order = 0; order < n_derivatives; ++order)
                {
                    const Real scale = std::pow(one_div_size, order);
                    const auto &b_values = berns_values.get_derivative(order);
                    basis.get_derivative(order) =
                            scale * prec_prod(oper, b_values);
                }
            }
        }

    }

//        for (auto comp : splines1d_.get_active_components_id())
//    {
//        auto &splines1d = splines1d_[comp];
//
//        {
//            const int num_intervals = n_intervals[jDim];
//            const int deg = degree[comp][jDim];
//            BasisValues1d bernstein_values(n_derivatives, deg+1, n_points[jDim]);
//
//            const auto &pt_coords = points.get_data_direction(jDim);
//
//            // fill values and derivatives of the Bernstein's polynomials at
//            // quad points in [0,1]
//            for (int order = 0; order < n_derivatives; ++order)
//                bernstein_values.get_derivative(order) =
//                    BernsteinBasis::derivative(order, deg, pt_coords);
//
//            const auto &bez_iComp_jDim = bezier_op.get_operator(comp,jDim);
//            const auto &lengths_jDim = lengths.get_data_direction(jDim);
//
//            // compute the one dimensional B-splines at quad point on the reference interval
//            for (int i = 0 ; i < num_intervals ; ++i)
//            {
//                const auto &M = bez_iComp_jDim[i];
//                const Real one_div_size = 1.0 / lengths_jDim[i];
//                BasisValues1d &basis = splines1d.entry(jDim,i);
//
//                for (int order = 0; order < n_derivatives; ++order)
//                {
//                    const Real scaling_coef = std::pow(one_div_size, order);
//                    basis.get_derivative(order) = scaling_coef * prec_prod(M, bernstein_values.get_derivative(order));
//                } //end loop order
//
//            } // end loop interval
//        } // end loop jDim
//    } // end loop comp
}



template<int dim_, int range_ , int rank_>
void
BSplineUniformQuadCache<dim_, range_, rank_>::
init_element_cache(ElementAccessor &elem)
{
    base_t::init_element_cache(elem);
    auto n_basis = space_->get_num_all_element_basis();
    auto &cache = elem.get_elem_cache();
    cache.resize(flags_, quad_, n_basis);

    auto &face_cache = elem.face_values_;
    for (auto f: base_t::faces)
    {
        auto &f_cache = face_cache[f];
        f_cache.resize(face_flags_, quad_, n_basis, f);
    }
}



template<int dim_, int range_ , int rank_>
void
BSplineUniformQuadCache<dim_, range_, rank_>::
init_element_cache(ElementIterator &elem)
{
    init_element_cache(elem.get_accessor());
}



template <int dim, int range, int rank>
template <int order>
void
BSplineUniformQuadCache<dim, range, rank>::
evaluate_bspline_derivatives(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Val<order>> &D_phi) const
{

    Assert(D_phi.size() > 0, ExcEmptyObject());
//  Assert(D_phi.get_num_functions() == this->get_num_basis(),
//           ExcDimensionMismatch(D_phi.get_num_functions(),this->get_num_basis()));

    const auto n_points_direction = quad_.get_num_points_direction();
    const Size num_points = n_points_direction.flat_size();
    Assert(D_phi.get_num_points() == num_points,
           ExcDimensionMismatch(D_phi.get_num_points(),num_points));

    /*
     * This code computes any order of derivatives for a multivariate
     * B-spline on the current element
     * We use the formula
     * \partial_(\alpha_1,...,\alpha_n) B(qp) = \Pi d^{\alpha_i} B_i(qp_i)
     */

    if (order == 0)
    {
        const TensorIndex<dim> zero; // [0,0,..,0] tensor index
        for (int comp : elem_values.get_active_components_id())
        {
        	auto &values = elem_values[comp];
        	//const auto &n_basis_dir = n_basis_[comp];
        	const int total_n_basis = n_basis_.comp_dimension[comp];
            const Size offset = comp_offset_[comp];

            for (int func_id = 0; func_id < total_n_basis; ++func_id)
            {
            	auto D_phi_i = D_phi.get_function_view(offset + func_id);
            	auto const &func = values.func_flat_to_tensor(func_id);

            	for (int point_id = 0; point_id < num_points; ++point_id)
            	{
            	    auto const &pts  = values.points_flat_to_tensor(point_id);

                    D_phi_i[point_id](comp) = values.evaluate(zero, func, pts);
            	}
            } // end func_id loop
        } // end comp loop

        for (int comp : elem_values.get_inactive_components_id())
        {
            const auto n_basis = n_basis_.comp_dimension[comp];
            const auto act_comp = elem_values.active(comp);
            const Size act_offset = comp_offset_[act_comp];
            const Size offset = comp_offset_[comp];
            for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
            {
                const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
                auto inact_D_phi = D_phi.get_function_view(offset+basis_i);

                for (int qp = 0; qp < num_points; ++qp)
                    inact_D_phi[qp](comp) = act_D_phi[qp](act_comp);
            }
        }

    } // end if (order == 0)
#if 0
    else // if (order != 0)
    {
        Assert(order >= 1, ExcLowerRange(order,1));

        using DerSymmMngr_t = DerivativeSymmetryManager<dim,order>;
        DerSymmMngr_t derivative_symmetry_manager;
        const auto &derivatives_flat_id_evaluate = derivative_symmetry_manager.get_entries_flat_id_evaluate();
        const auto &derivatives_flat_id_copy_to = derivative_symmetry_manager.get_entries_flat_id_copy_to();
        const auto &derivatives_flat_id_copy_from = derivative_symmetry_manager.get_entries_flat_id_copy_from();

        const auto &derivatives_tensor_id = derivative_symmetry_manager.get_entries_tensor_id();

        const Size n_derivatives_eval = DerSymmMngr_t::num_entries_eval;
        const Size n_derivatives_copy = DerSymmMngr_t::num_entries_copy;

        using der_t = Conditional<order==0,Value,Derivative<order>>;

        for (int comp : scalar_evaluators_.get_active_components_id())
        {
            const auto n_basis = n_basis_.comp_dimension[comp];
            const Size comp_offset_i = comp_offset_[comp];

            DynamicMultiArray<Real,dim> derivative_scalar_component(n_points_direction);
            for (int func_flat_id = 0; func_flat_id < n_basis; ++func_flat_id)
            {
                auto D_phi_i = D_phi.get_function_view(comp_offset_i+func_flat_id);

                const auto &scalar_bspline = *scalar_evaluators_[comp][func_flat_id];

                for (int entry_id = 0; entry_id < n_derivatives_eval; ++entry_id)
                {
                    const int entry_flat_id = derivatives_flat_id_evaluate[entry_id];
                    const auto entry_tensor_id = derivatives_tensor_id[entry_flat_id];

                    // from the entry_tensor_id we get the right derivative order
                    TensorIndex<dim> deriv_order_tensor_id; // [0,0,..,0] tensor index
                    for (int i = 0; i < order; ++i)
                        ++(deriv_order_tensor_id[entry_tensor_id[i]]);

                    //TODO: remove this if!!! (Maybe re-think about the BSplineSpace for dim==0)
                    if (dim > 0)
                    {
                        scalar_bspline.evaluate_derivative_at_points(deriv_order_tensor_id,derivative_scalar_component);

                        for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
                            D_phi_i[point_flat_id](entry_flat_id)(comp) = derivative_scalar_component[point_flat_id];
                    }
                    else
                    {
                        for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
                            D_phi_i[point_flat_id](entry_flat_id)(comp) = 1.0;
                    }

                } // end entry_id loop

            } // end func_flat_id loop

            for (int func_flat_id = 0; func_flat_id < n_basis; ++func_flat_id)
            {
                auto D_phi_i = D_phi.get_function_view(comp_offset_i+func_flat_id);
                for (int point_flat_id = 0; point_flat_id < num_points; ++point_flat_id)
                {
                    der_t &derivative = D_phi_i[point_flat_id];

                    // here we copy the computed quantities to the symmetric part of the tensor
                    for (int entry_id = 0; entry_id < n_derivatives_copy; ++entry_id)
                        derivative(derivatives_flat_id_copy_to[entry_id])(comp) =
                            derivative(derivatives_flat_id_copy_from[entry_id])(comp);
                } // end point_flat_id loop
            } //end func_flat_id loop

        } // end comp loop

        for (int comp : scalar_evaluators_.get_inactive_components_id())
        {
            const Size n_ders = Derivative<order>::size;


            const auto n_basis = n_basis_.comp_dimension[comp];
            const auto scalar_eval_act_comp = scalar_evaluators_.active(comp);
            const Size act_offset = comp_offset_[scalar_eval_act_comp];
            const Size offset = comp_offset_[comp];
            for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
            {
                const auto derivatives_phi_hat_copy_from = D_phi.get_function_view(act_offset+basis_i);
                auto derivatives_phi_hat_copy_to = D_phi.get_function_view(offset+basis_i);
                for (int qp = 0; qp < num_points; ++qp)
                {
                    const der_t &values_0 = derivatives_phi_hat_copy_from[qp];
                    der_t &values = derivatives_phi_hat_copy_to[qp];

                    for (int der = 0; der < n_ders; ++der)
                        values(der)(comp) = values_0(der)(scalar_eval_act_comp);
                }
            } //end loop basis_i
        } // end loop comp
    } // end if (order > 0)
#endif
}



template<int dim_, int range_ , int rank_>
void
BSplineUniformQuadCache<dim_, range_, rank_>::
fill_element_cache(ElementIterator &elem)
{
    fill_element_cache(elem.get_accessor());
}



template<int dim_, int range_ , int rank_>
void
BSplineUniformQuadCache<dim_, range_, rank_>::
fill_element_cache(ElementAccessor &elem)
{
    base_t::fill_element_cache(elem);
    auto &cache = elem.get_elem_cache();

    const auto &elem_id = elem.get_tensor_index();
    auto values = splines1d_.get_element_values(elem_id);


    //--------------------------------------------------------------------------
    if (cache.flags_handler_.fill_values())
    {
        evaluate_bspline_derivatives<0>(values, cache.phi_);
        cache.flags_handler_.set_values_filled(true);
    }
    if (cache.flags_handler_.fill_gradients())
    {
        evaluate_bspline_derivatives<1>(values, cache.D1phi_);
        cache.flags_handler_.set_gradients_filled(true);
    }

    if (cache.flags_handler_.fill_hessians())
    {
        evaluate_bspline_derivatives<2>(values, cache.D2phi_);
        cache.flags_handler_.set_hessians_filled(true);
    }

    if (cache.flags_handler_.fill_divergences())
    {
        //TODO(pauletti, Sep 7, 2014): create a specialize exception
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



template<int dim_, int range_ , int rank_>
void
BSplineUniformQuadCache<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    out.begin_item("Grid Cache:");
    base_t::print_info(out);
    out.end_item();


    out.begin_item("One dimensional splines cache:");
    splines1d_.print_info(out);

    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_uniform_quad_cache.inst>
