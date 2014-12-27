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


#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bernstein_basis.h>

#include <igatools/utils/multi_array_utils.h>

#include <algorithm>
#include <numeric>
#include <memory>

using std::reverse;
using std::accumulate;
using std::sort;

using std::shared_ptr;
using std::make_shared;
using std::array;


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
BSplineElement<dim, range, rank>::
BSplineElement(const std::shared_ptr<ContainerType> space,
               const Index index)
    :
    parent_t(space,index)
{}



template <int dim, int range, int rank>
BSplineElement<dim, range, rank>::
BSplineElement(const std::shared_ptr<ContainerType> space,
               const TensorIndex<dim> &index)
    :
    parent_t(space,index)
{}


template <int dim, int range, int rank>
BSplineElement<dim, range, rank>::
BSplineElement(const self_t &elem,
               const CopyPolicy &copy_policy)
    :
    parent_t(elem,copy_policy)
{}




template <int dim, int range, int rank>
auto
BSplineElement<dim, range, rank>::
evaluate_univariate_derivatives_at_points(
    const int deriv_order,
    const std::array<vector<Real>,dim> &points) const
-> ComponentContainer<std::array<ValueTable<Real>,dim> >
{
    TensorSize<dim> n_points_direction;
    for (int i = 0 ; i < dim ; ++i)
    {
        Assert(points[i].empty() == false,ExcEmptyObject());

        n_points_direction[i] = points[i].size();
    }

    const auto &element_tensor_id = this->get_tensor_index();

    ComponentContainer< array<ValueTable<Real>,dim> > funcs1D_table(this->space_->get_components_map());

    const auto degree_table = this->space_->get_degree();
    const auto &bezier_op_ = this->space_->operators_;

    const auto element_lengths = CartesianGridElement<dim>::template get_coordinate_lengths<dim>(0);

    for (int comp : funcs1D_table.get_active_components_id())
    {
        auto &funcs1D_comp = funcs1D_table[comp];

        const auto &degree_comp = degree_table[comp];

        auto n_basis_direction = TensorSize<dim>(degree_comp+1);

        for (int i = 0 ; i < dim ; ++i)
        {
            const auto &M = bezier_op_.get_operator(comp,i)[element_tensor_id[i]];

            const auto lengths_dir = element_lengths[i];
            const Real one_div_size = Real(1.0) / lengths_dir;
            const Real scaling_coef = std::pow(one_div_size, deriv_order);

            // compute the one dimensional Bernstein at quad point on the unit interval
            const auto B = BernsteinBasis::derivative(deriv_order,degree_comp[i],points[i]);

            // compute the one dimensional B-splines at quad point on the reference interval
            const auto basis = scaling_coef * prec_prod(M,B);

            auto &funcs1D_comp_dir = funcs1D_comp[i];
            funcs1D_comp_dir.resize(n_basis_direction[i],n_points_direction[i]);

            for (int fn = 0 ; fn < n_basis_direction[i] ; ++fn)
            {
                auto fn_view = funcs1D_comp_dir.get_function_view(fn);

                for (int pt = 0 ; pt < n_points_direction[i] ; ++pt)
                    fn_view[pt] = basis(fn,pt);
            } // end fn loop
        } // end dir loop
    } // end comp loop


    return funcs1D_table;
}


template <int dim, int range, int rank>
auto
BSplineElement<dim, range, rank>::
evaluate_univariate_derivatives_at_points(const int deriv_order, const Quadrature<dim> &quad) const
-> ComponentContainer<std::array<ValueTable<Real>,dim> >
{
    std::array<vector<Real>,dim> points_coords;
    for (int i = 0 ; i < dim ; ++i)
        points_coords[i] = quad.get_points().get_data_direction(i);

    return this->evaluate_univariate_derivatives_at_points(deriv_order,points_coords);
}

template <int dim, int range, int rank>
auto
BSplineElement<dim, range, rank>::
evaluate_univariate_derivatives_at_points(
    const int deriv_order,
    const ValueVector<Point> &points) const -> ComponentContainer<std::array<ValueTable<Real>,dim> >
{
    std::array<vector<Real>,dim> points_coords;
    for (const auto &pt : points)
        for (int i = 0 ; i < dim ; ++i)
            points_coords[i].push_back(pt[i]);

    return this->evaluate_univariate_derivatives_at_points(deriv_order,points_coords);
}







// TODO (pauletti, Jun 3, 2014): the name and arguments should be more consistent
// this function and evaluate_bspline_derivatives

template <int dim, int range, int rank>
template<int deriv_order>
auto
BSplineElement<dim, range, rank>::
evaluate_basis_derivatives_at_points(const ValueVector<Point> &points) const ->
ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >
{
    using return_t = ValueTable< Conditional< deriv_order==0,Value,Derivative<deriv_order> > >;

    const Size n_basis  = this->get_num_basis();
    const Size n_points = points.size();

    return_t D_phi(n_basis,n_points);

    const auto bezier_op = this->space_->operators_.get_element_operators(this->get_tensor_index());

    if (deriv_order == 0)
    {
        for (Size pt_id = 0 ; pt_id < n_points ; ++pt_id)
        {
            const Point point = points[pt_id];

            auto derivatives_phi_hat_ipt = D_phi.get_point_view(pt_id);

#ifndef NDEBUG
            for (int dir = 0 ; dir < dim ; ++dir)
                Assert(point[dir] >= 0.0 && point[dir] <= 1.0,
                ExcMessage("Evaluation point " + std::to_string(pt_id) + " not in the unit-domain."));
#endif
//            auto n_basis = this->n_basis_direction_;
            auto degree = this->space_->get_degree();
            for (int iComp : bezier_op.get_active_components_id())
            {
                //------------------------------------------------------------------------------
                // evaluation of the values/derivarives of the 1D Bernstein polynomials -- begin
                array<boost::numeric::ublas::vector<Real>,dim> bernstein_values;
                for (int dir = 0 ; dir < dim ; ++dir)
                    bernstein_values[dir] = BernsteinBasis::derivative(0,degree[iComp][dir],point[dir]);
                // evaluation of the values/derivarives of the 1D Bernstein polynomials -- end
                //------------------------------------------------------------------------------


                //--------------------------------------------------------------------------------
                // apply the Bezier extraction operator for the functions on this element -- begin
                const auto &bezier_op_comp = bezier_op[iComp];

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

                const Size comp_offset_i = this->comp_offset_[iComp];

                const auto &basis_flat_to_tensor = *(this->basis_functions_indexer_)[iComp];

                const int n_basis_component = this->get_num_basis(iComp);

                for (Size basis_fid = 0 ; basis_fid < n_basis_component ; ++basis_fid)
                {
                    const auto &basis_tid = basis_flat_to_tensor[basis_fid];

                    Real value = 1.0;
                    for (int dir = 0 ; dir < dim ; ++dir)
                        value *= bspline_basis[dir][basis_tid[dir]];

                    derivatives_phi_hat_ipt[basis_fid+comp_offset_i](iComp) = value;
                }
                // multiply the spline 1D in order to have the multi-d value -- end
                //--------------------------------------------------------------------------------

            } // end iComp loop
        } // end pt_id loop

        for (int comp : bezier_op.get_inactive_components_id())
        {
            const int n_basis = this->get_num_basis(comp);
            const Size offset = this->comp_offset_[comp];
            const Size act_offset = this->comp_offset_[bezier_op.active(comp)];
            for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
            {
                const auto values_phi_hat_copy_from = D_phi.get_function_view(act_offset+basis_i);
                auto values_phi_hat_copy_to = D_phi.get_function_view(offset+basis_i);

                for (int qp = 0; qp < n_points; ++qp)
                    values_phi_hat_copy_to[qp](comp) = values_phi_hat_copy_from[qp](0);

            } //end loop basis_i
        } // end loop comp

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

        const auto elem_lengths = CartesianGridElement<dim>::get_coordinate_lengths();

        for (Size pt_id = 0 ; pt_id < n_points ; ++pt_id)
        {
            const Point point = points[pt_id];

            auto derivatives_phi_hat_ipt = D_phi.get_point_view(pt_id);

#ifndef NDEBUG
            for (int dir = 0 ; dir < dim ; ++dir)
                Assert(point[dir] >= 0.0 && point[dir] <= 1.0,
                ExcMessage("Evaluation point " + std::to_string(pt_id) + " not in the unit-domain."));
#endif
//            auto n_basis = this->n_basis_direction_;
            auto degree = this->space_->get_degree();
            for (int iComp : bezier_op.get_active_components_id())
            {
                //------------------------------------------------------------------------------
                // evaluation of the values/derivarives of the 1D Bernstein polynomials -- begin
                array<array<boost::numeric::ublas::vector<Real>,dim>,deriv_order+1> bernstein_values;
                for (int order = 0 ; order <= deriv_order ; ++order)
                {
                    for (int dir = 0 ; dir < dim ; ++dir)
                    {
                        const Real scaling_coef = pow(1.0/elem_lengths[dir],order);

                        bernstein_values[order][dir] =
                        scaling_coef * BernsteinBasis::derivative(order,degree[iComp][dir],point[dir]);
                    }
                }
                // evaluation of the values/derivarives of the 1D Bernstein polynomials -- end
                //------------------------------------------------------------------------------



                //--------------------------------------------------------------------------------
                // apply the Bezier extraction operator for the functions on this element -- begin
                const auto &bezier_op_comp = bezier_op[iComp];

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
                const Size comp_offset_i = this->comp_offset_[iComp];

                const auto &basis_flat_to_tensor = *(this->basis_functions_indexer_)[iComp];

                const int n_basis_component = this->get_num_basis(iComp);

                for (Size basis_fid = 0 ; basis_fid < n_basis_component ; ++basis_fid)
                {

                    const auto &basis_tid = basis_flat_to_tensor[basis_fid];

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

        for (int comp : bezier_op.get_inactive_components_id())
        {
            const Size n_ders = Derivative<deriv_order>::size;
            const auto n_basis =  this->get_num_basis(comp);
            const Size act_offset = this->comp_offset_[bezier_op.active(comp)];
            const Size offset = this->comp_offset_[comp];
            for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
            {
                const auto derivatives_phi_hat_copy_from = D_phi.get_function_view(act_offset+basis_i);
                auto derivatives_phi_hat_copy_to = D_phi.get_function_view(offset+basis_i);
                for (int qp = 0; qp < n_points; ++qp)
                {
                    const auto &values_0 = derivatives_phi_hat_copy_from[qp];
                    auto &values = derivatives_phi_hat_copy_to[qp];

                    for (int der = 0; der < n_ders; ++der)
                        values(der)(comp) = values_0(der)(0);
                } // end loop qp
            } //end loop basis_i
        } // end loop comp

    } // end if (deriv_order != 0)

    return D_phi;
}



IGA_NAMESPACE_CLOSE

//#include <igatools/basis_functions/bspline_element.inst>


