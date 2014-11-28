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

#include <igatools/basis_functions/nurbs_element_handler.h>
#include <igatools/basis_functions/nurbs_element.h>

using std::shared_ptr;

#ifdef NURBS

IGA_NAMESPACE_OPEN

template<int dim_, int range_ , int rank_>
NURBSElementHandler<dim_, range_, rank_>::
NURBSElementHandler(shared_ptr<const Space> space)
    :
    base_t(space->get_spline_space()),
    space_(space)
{}


#if 0
template<int dim_, int range_ , int rank_>
NURBSElementHandler<dim_, range_, rank_>::
NURBSElementHandler(shared_ptr<const Space> space)
    :
    base_t(space->get_grid()),
    space_(space),
    n_basis_(space_->get_num_all_element_basis())
{

    // Compute the component offsets
    comp_offset_[0] = 0;
    for (int j = 1; j < Space::n_components; ++j)
        comp_offset_[j] = comp_offset_[j-1] + n_basis_.comp_dimension[j-1];
}
#endif


template<int dim_, int range_ , int rank_>
template<int k>
void
NURBSElementHandler<dim_, range_, rank_>::
reset(const NewValueFlags flag,
      const Quadrature<k> &quad1)
{
    base_t::template reset<k>(flag, quad1);
    flags_[k] = flag;
}



template<int dim_, int range_ , int rank_>
template<int k>
void
NURBSElementHandler<dim_, range_, rank_>::
init_cache(ElementAccessor &elem)
{
    base_t::template init_cache<k>(elem.bspline_elem_);

    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    const auto n_basis = elem.get_num_basis();

    for (auto &s_id: UnitElement<dim>::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_value_cache<k>(s_id);
        const auto n_points = this->template get_num_points<k>();
        s_cache.resize(flags_[k], n_points, n_basis);
    }
}





template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
init_element_cache(ElementIterator &elem)
{
    init_cache<dim>(elem.get_accessor());
}


#if 0
template <int dim, int range, int rank>
void
NURBSElementHandler<dim, range, rank>::
copy_to_inactive_components_values(const vector<Index> &inactive_comp,
                                   const std::array<Index, n_components> &active_map,
                                   ValueTable<Value> &D_phi) const
{
    const Size n_points = D_phi.get_num_points();
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis = n_basis_.comp_dimension[comp];
        const Size act_offset = comp_offset_[act_comp];
        const Size offset     = comp_offset_[comp];
        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
        {
            const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
            auto inact_D_phi = D_phi.get_function_view(offset+basis_i);
            for (int qp = 0; qp < n_points; ++qp)
                inact_D_phi[qp](comp) = act_D_phi[qp](act_comp);
        }
    }
}



template <int dim, int range, int rank>
template <int order>
void
NURBSElementHandler<dim, range, rank>::
copy_to_inactive_components(const vector<Index> &inactive_comp,
                            const std::array<Index, n_components> &active_map,
                            ValueTable<Derivative<order>> &D_phi) const
{
    const Size n_points = D_phi.get_num_points();
    const Size n_ders = Derivative<order>::size;
    for (int comp : inactive_comp)
    {
        const auto act_comp = active_map[comp];
        const auto n_basis = n_basis_.comp_dimension[comp];
        const Size act_offset = comp_offset_[act_comp];
        const Size offset     = comp_offset_[comp];
        for (Size basis_i = 0; basis_i < n_basis;  ++basis_i)
        {
            const auto act_D_phi = D_phi.get_function_view(act_offset+basis_i);
            auto     inact_D_phi = D_phi.get_function_view(offset+basis_i);
            for (int qp = 0; qp < n_points; ++qp)
                for (int der = 0; der < n_ders; ++der)
                    inact_D_phi[qp](der)(comp) = act_D_phi[qp](der)(act_comp);
        }
    }
}


template <int dim, int range, int rank>
void
NURBSElementHandler<dim, range, rank>::
evaluate_bspline_values(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Value> &D_phi) const
{
    const Size n_points = D_phi.get_num_points();
    const TensorIndex<dim> der_tensor_id; // [0,0,..,0] tensor index
    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int total_n_basis = n_basis_.comp_dimension[comp];
        const Size offset = comp_offset_[comp];

        for (int func_id = 0; func_id < total_n_basis; ++func_id)
        {
            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func = values.func_flat_to_tensor(func_id);
            for (int point_id = 0; point_id < n_points; ++point_id)
            {
                auto const &pts  = values.points_flat_to_tensor(point_id);
                D_phi_i[point_id](comp) = values.evaluate(der_tensor_id, func, pts);
            }
        } // end func_id loop
    } // end comp loop

    copy_to_inactive_components_values(elem_values.get_inactive_components_id(),
                                       elem_values.get_comp_map(), D_phi);
}



template <int dim, int range, int rank>
template <int order>
void
NURBSElementHandler<dim, range, rank>::
evaluate_bspline_derivatives(
    const ComponentContainer<TensorProductFunctionEvaluator<dim>> &elem_values,
    ValueTable<Derivative<order>> &D_phi) const
{
    /*
     * This code computes any order of derivatives for a multivariate
     * B-spline on the current element
     * We use the formula
     * \partial_(\alpha_1,...,\alpha_n) B(qp) = \Pi d^{\alpha_i} B_i(qp_i)
     */

    Assert(D_phi.size() > 0, ExcEmptyObject());
    //  Assert(D_phi.get_num_functions() == this->get_num_basis(),
    //           ExcDimensionMismatch(D_phi.get_num_functions(),this->get_num_basis()));
    const Size n_points = D_phi.get_num_points();


    TensorFunctionDerivativesSymmetry<dim,order> sym;
    const auto n_der =  TensorFunctionDerivativesSymmetry<dim,order>::num_entries_eval;

    const auto &univariate_order = sym.univariate_order ;
    const auto &copy_indices = sym.copy_indices;

    for (int comp : elem_values.get_active_components_id())
    {
        auto &values = elem_values[comp];
        const int total_n_basis = n_basis_.comp_dimension[comp];
        const Size offset = comp_offset_[comp];

        for (int func_id = 0; func_id < total_n_basis; ++func_id)
        {
            auto D_phi_i = D_phi.get_function_view(offset + func_id);
            auto const &func = values.func_flat_to_tensor(func_id);
            for (int der_id = 0; der_id<n_der; ++der_id)
            {
                const auto &copy_indices_der = copy_indices[der_id];
                const auto copy_indices_der_size = copy_indices_der.size();

                auto const &der_tensor_id = univariate_order[der_id];
                for (int point_id = 0; point_id < n_points; ++point_id)
                {
                    auto const &pts  = values.points_flat_to_tensor(point_id);
                    auto &der = D_phi_i[point_id];
                    der(copy_indices_der[0])(comp) = values.evaluate(der_tensor_id, func, pts);
                    for (int k=1; k<copy_indices_der_size; ++k)
                        der(copy_indices_der[k])(comp) = der(copy_indices_der[0])(comp);
                }
            }

        }
    } // end comp loop

    copy_to_inactive_components<order>(elem_values.get_inactive_components_id(),
                                       elem_values.get_comp_map(), D_phi);
}
#endif



template<int dim_, int range_ , int rank_>
template <int k>
void
NURBSElementHandler<dim_, range_, rank_>::
fill_cache(ElementAccessor &elem, const int j)
{
    base_t::template fill_cache<k>(elem.bspline_elem_, j);


    Assert(elem.local_cache_ != nullptr, ExcNullPtr());
    auto &cache = elem.local_cache_->template get_value_cache<k>(j);

    auto &flags = cache.flags_handler_;
    if (flags.fill_values())
    {
        auto &values = cache.template get_der<0>();
        evaluate_nurbs_from_bspline(elem.bspline_elem_,elem.weight_elem_,values);
        flags.set_values_filled(true);
    }
    if (flags.fill_gradients())
    {
        auto &values = cache.template get_der<1>();
        evaluate_nurbs_from_bspline(elem.bspline_elem_,elem.weight_elem_,values);
        flags.set_gradients_filled(true);
    }
    if (flags.fill_hessians())
    {
        auto &values = cache.template get_der<2>();
        evaluate_nurbs_from_bspline(elem.bspline_elem_,elem.weight_elem_,values);
        flags.set_hessians_filled(true);
    }

    cache.set_filled(true);
}



template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
fill_element_cache(ElementIterator &elem)
{
    fill_cache<dim>(elem.get_accessor(), 0);
}




template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
#if 0

    out.begin_item("Grid Cache:");
    base_t::print_info(out);
    out.end_item();

    cacheutils::print_caches(splines1d_, out);
    //out.begin_item("One dimensional splines cache:");
    //splines1d_.print_info(out);
    //out.end_item();
#endif
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
evaluate_nurbs_from_bspline(
    const typename Space::SpSpace::ElementAccessor &bspline_elem,
    const typename Space::WeightFunction::ElementAccessor &weight_elem,
    ValueTable<Value> &phi) const
{
    Assert(false,ExcNotImplemented());
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
evaluate_nurbs_from_bspline(
    const typename Space::SpSpace::ElementAccessor &bspline_elem,
    const typename Space::WeightFunction::ElementAccessor &weight_elem,
    ValueTable<Derivative<1>> &D1_phi) const
{
    Assert(false,ExcNotImplemented());
}

template<int dim_, int range_ , int rank_>
void
NURBSElementHandler<dim_, range_, rank_>::
evaluate_nurbs_from_bspline(
    const typename Space::SpSpace::ElementAccessor &bspline_elem,
    const typename Space::WeightFunction::ElementAccessor &weight_elem,
    ValueTable<Derivative<2>> &D2_phi) const
{
    Assert(false,ExcNotImplemented());
}


IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/nurbs_element_handler.inst>

#endif // #ifdef NURBS
