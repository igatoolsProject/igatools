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


using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim_, int range_ , int rank_>
BSplineUniformQuadCache<dim_, range_, rank_>::
BSplineUniformQuadCache(shared_ptr<const Space> space,
                        const ValueFlags flag,
                        const Quadrature<dim> &quad)
    :
    GridUniformQuadCache<dim_>(space->get_grid(), flag, quad),
    space_(space),
    flags_(flag),
    quad_(quad),
    splines1d_(space->get_components_map(),
               CartesianProductArray<BasisValues1d, dim_>(space->get_grid()->get_num_intervals()))
{
    const int n_derivatives = 3;
    const auto grid = space->get_grid();
    //const auto n_intervals = grid->get_num_intervals();
    const auto degree = space->get_degree();
//    const auto n_funcs = space->get_num_basis_per_element_table();
    const auto n_points = quad.get_num_points_direction();
    for (auto comp : splines1d_.get_active_components_id())
    {
        auto &splines1d = splines1d_[comp];
        const auto size = splines1d.tensor_size();
        for (int dir = 0 ; dir < dim ; ++dir)
        {
            auto n_fun = degree[comp][dir]+1;
            auto n_pts = n_points[dir];
            for (int j = 0 ; j < size[dir] ; ++j)
                splines1d.entry(dir,j).resize(n_derivatives, n_fun, n_pts);
        }
    }


}



template<int dim_, int range_ , int rank_>
void
BSplineUniformQuadCache<dim_, range_, rank_>::
init_element_cache(ElementIterator &elem)
{
//    // TODO (pauletti, Aug 14, 2014): create get_cache in accessor
//    auto &cache = elem.get_accessor().elem_values_;
//    cache.resize(flags_, quad_);
}



template<int dim_, int range_ , int rank_>
void
BSplineUniformQuadCache<dim_, range_, rank_>::
fill_element_cache(ElementIterator &elem)
{
//    auto &cache = elem.get_accessor().elem_values_;
//    auto meas = lengths_.tensor_product(elem->get_tensor_index());
//    cache.fill(meas);
//    cache.set_filled(true);
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
    // TODO (pauletti, Aug 21, 2014): This should just be splines1d_.print_info
    for (auto spline : splines1d_)
    {
        for (int dir = 0 ; dir < dim ; ++dir)
            for (auto basis : spline.get_data_direction(dir))
                basis.print_info(out);
    }
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_uniform_quad_cache.inst>
