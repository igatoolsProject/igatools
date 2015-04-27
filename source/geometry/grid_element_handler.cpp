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

#include <igatools/geometry/grid_element_handler.h>
//#include <igatools/geometry/unit_element.h>


using std::shared_ptr;
using std::array;

IGA_NAMESPACE_OPEN


template <int dim>
GridElementHandler<dim>::
GridElementHandler(shared_ptr<GridType> grid)
    :
    grid_(grid)
{}

template <int dim>
std::shared_ptr<GridElementHandler<dim> >
GridElementHandler<dim>::
create(std::shared_ptr<GridType> grid)
{
    using ElemHandler = GridElementHandler<dim>;
    auto elem_handler = std::make_shared<ElemHandler>(grid);
    Assert(elem_handler != nullptr,ExcNullPtr());
    return elem_handler;
}


template <int dim>
template<int k>
void
GridElementHandler<dim>::
reset(const ValueFlags flag,
      const Quadrature<k> &quad)
{
    const auto valid_flags = ElementAccessor::get_valid_flags();
    auto grid_flag = flag & valid_flags;

    if (contains(flag, ValueFlags::value))
        grid_flag |= ValueFlags::point;

    flags_[k] = grid_flag;

    cacheutils::extract_sub_elements_data<k>(quad_all_sub_elems_) = quad;
}

#if 0
template <int dim>
void
GridElementHandler<dim>::
init_all_caches(ElementAccessor &elem)
{
    auto &cache = elem.all_sub_elems_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::CacheType;
        cache = shared_ptr<Cache>(new Cache);
    }

    const auto &quad = cacheutils::extract_sub_elements_data<dim>(quad_all_sub_elems_);

    boost::fusion::for_each(cache->cache_all_sub_elems_,
                            [&](auto & value_dim) -> void
    {
        using PairType = typename std::remove_reference<decltype(value_dim)>::type;
        const int topology_dim = PairType::first_type::value;
        auto &cache_same_topology_dim = value_dim.second;
        int topology_id = 0;
        for (auto &cache_same_topology_id : cache_same_topology_dim)
        {
            cache_same_topology_id.resize(
                flags_[dim],
                quad.template collapse_to_sub_element<topology_dim>(topology_id));
            ++topology_id;
        }
    }
                           );
}
#endif

template <int dim>
template <int k>
void
GridElementHandler<dim>::
init_cache(ElementAccessor &elem)
{
    auto &cache = elem.all_sub_elems_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::CacheType;
        cache = std::make_shared<Cache>();
    }

    for (auto &s_id: UnitElement<dim>::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_sub_elem_cache<k>(s_id);
//        s_cache.resize(flags_[k], extend_sub_elem_quad<k, dim>(
//                           cacheutils::extract_sub_elements_data<k>(quad_all_sub_elems_), s_id));
        s_cache.resize(flags_[k], this->template get_num_points<k>());

    }
}




template <int dim>
template <int k>
void
GridElementHandler<dim>::
fill_cache(ElementAccessor &elem, const int j)
{
    Assert(elem.all_sub_elems_cache_ != nullptr, ExcNullPtr());
    auto &cache = elem.all_sub_elems_cache_->template get_sub_elem_cache<k>(j);

    const auto &quadrature = this->template get_quadrature<k>();

    if (cache.template status_fill<_Point>())
    {
        const auto unit_points = quadrature.get_points();

        const auto translate = elem.vertex(0);
        const auto dilate    = elem.template get_coordinate_lengths<k>(j);

        const auto n_pts = unit_points.get_num_points();

        const auto &sub_unit_elem = UnitElement<dim>::template get_elem<k>(j);
        auto &ref_pts = cache.template get_data<_Point>();
        for (int pt = 0 ; pt < n_pts ; ++pt)
        {
            const auto &unit_pt = unit_points[pt];
            auto &ref_pt = ref_pts[pt];

            int sub_elem_dir = 0;
            for (const auto active_dir : sub_unit_elem.active_directions)
            {
                ref_pt[active_dir] = translate[active_dir] +
                                     dilate[active_dir] * unit_pt[sub_elem_dir] ;
                ++sub_elem_dir;
            }

            sub_elem_dir = 0;
            for (const auto constant_dir : sub_unit_elem.constant_directions)
            {
                ref_pt[constant_dir] = translate[constant_dir] +
                                       dilate[constant_dir] * sub_unit_elem.constant_values[sub_elem_dir];
                ++sub_elem_dir;
            }
        }

        cache.template set_status_filled<_Point>(true);
    }

    if (cache.template status_fill<_W_Measure>())
    {
        cache.template get_data<_W_Measure>() = elem.template get_measure<k>(j) * quadrature.get_weights();
        cache.template set_status_filled<_W_Measure>(true);
    }


    cache.set_filled(true);
}



template <int dim>
auto
GridElementHandler<dim>::
get_grid() const -> std::shared_ptr<const GridType>
{
    return grid_;
}


template <int dim>
void
GridElementHandler<dim>::
print_info(LogStream &out) const
{
    out.begin_item("Quadrature cache for all dimensions:");

    boost::fusion::for_each(quad_all_sub_elems_,
                            [&](const auto & data_same_topology_dim)
    {
        using PairType = typename std::remove_reference<decltype(data_same_topology_dim)>::type;
        using SubDimType = typename PairType::first_type;

        const auto &quad_same_subdim = data_same_topology_dim.second;

        out.begin_item("Quadrature cache for dimension: " + std::to_string(SubDimType::value));
        quad_same_subdim.print_info(out);
        out.end_item();
    }
                           );
    out.end_item();

}


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_element_handler.inst>
