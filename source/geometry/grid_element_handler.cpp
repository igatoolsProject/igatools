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
    auto elem_handler = std::shared_ptr<ElemHandler>(new ElemHandler(grid));
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
    flags_[k] = flag & valid_flags;

    cacheutils::extract_sub_elements_data<k>(quad_) = quad;
}


template <int dim>
void
GridElementHandler<dim>::
init_all_caches(ElementAccessor &elem)
{
    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::CacheType;
        cache = shared_ptr<Cache>(new Cache);
    }

    const auto &quad = cacheutils::extract_sub_elements_data<dim>(quad_);

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
//#endif
//    Assert(false,ExcNotImplemented());
}


template <int dim>
template <int k>
void
GridElementHandler<dim>::
init_cache(ElementAccessor &elem)
{
    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::CacheType;
        cache = shared_ptr<Cache>(new Cache);
    }

    for (auto &s_id: Topology::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_sub_elem_cache<k>(s_id);
//        auto &quad = std::get<k>(quad_);
//        s_cache.resize(flags_[k], extend_sub_elem_quad<k, dim>(quad, s_id));
        s_cache.resize(flags_[k], extend_sub_elem_quad<k, dim>(
                           cacheutils::extract_sub_elements_data<k>(quad_), s_id));
    }
}




template <int dim>
template <int k>
void
GridElementHandler<dim>::
fill_cache(ElementAccessor &elem, const int j)
{
    Assert(elem.local_cache_ != nullptr, ExcNullPtr());
    auto &cache = elem.local_cache_->template get_sub_elem_cache<k>(j);

//    auto &flags = cache.flags_handler_;

    if (cache.template status_fill<_Point>())
    {
        auto translate = elem.vertex(0);
        auto dilate    = elem.template get_coordinate_lengths<k>(j);

        const int n_pts = cache.unit_points_.get_num_points();

        const auto &unit_pts = cache.unit_points_;
        auto &ref_pts = cache.template get_data<_Point>();
        for (int pt = 0 ; pt < n_pts ; ++pt)
        {
            const auto &unit_pt = unit_pts[pt];
            auto &ref_pt = ref_pts[pt];

            for (const auto dir : Topology::active_directions)
                ref_pt[dir] = unit_pt[dir] * dilate[dir] + translate[dir];
        }

        cache.template set_status_filled<_Point>(true);
    }

    if (cache.template status_fill<_W_Measure>())
    {
        cache.template get_data<_W_Measure>() = elem.template get_measure<k>(j) * cache.unit_weights_;
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


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_element_handler.inst>
