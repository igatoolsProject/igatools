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

#if 0
namespace
{
struct UniformQuadFunc
{
    void func(auto &val_cache, const GridFlags &flag, const auto &quad)
    {
        val_cache.resize(flag, quad);
    }
};

template<class Quad, class... Args>
void
init_unif_caches(const GridFlags &flag, const Quad &quad, std::tuple<Args...> &t)
//init_unif_caches(const GridFlags &flag, const Quad &quad, boost::fusion::vector<Args...> &t)
{
    const int dim = Quad::dim;
    const int low = dim==0? 0 : dim-num_sub_elem;
    UniformQuadFunc f;
    TupleFunc1< UniformQuadFunc, GridFlags, Quad, decltype(t), sizeof...(Args), low>::apply_func(f, flag, quad, t);
}
};
#endif


template <int dim>
GridElementHandler<dim>::
GridElementHandler(shared_ptr<GridType> grid)
    :
    grid_(grid)
//  ,
//    lengths_(grid->get_element_lengths())
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
    flags_[k] = flag;
//    auto &quad_k = std::get<k>(quad_);
//    quad_k = quad;
    boost::fusion::at_c<k>(quad_) = quad;
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
//    init_unif_caches(flags_[dim], std::get<dim>(quad_), cache->values_);
//    init_unif_caches(flags_[dim], boost::fusion::at_c<dim>(quad_), cache->values_);
    const auto &quad = boost::fusion::at_c<dim>(quad_);

    boost::fusion::for_each(cache->values_,
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
        auto &s_cache = cache->template get_value_cache<k>(s_id);
//        auto &quad = std::get<k>(quad_);
//        s_cache.resize(flags_[k], extend_sub_elem_quad<k, dim>(quad, s_id));
        s_cache.resize(flags_[k], extend_sub_elem_quad<k, dim>(boost::fusion::at_c<k>(quad_), s_id));
    }
}




template <int dim>
template <int k>
void
GridElementHandler<dim>::
fill_cache(ElementAccessor &elem, const int j)
{
    Assert(elem.local_cache_ != nullptr, ExcNullPtr());
    auto &cache = elem.local_cache_->template get_value_cache<k>(j);

    auto &flags = cache.flags_handler_;

    if (flags.template fill<_Point>())
    {
        auto translate = elem.vertex(0);
        auto dilate    = elem.template get_coordinate_lengths<k>(j);

        const int n_pts = cache.unit_points_.get_num_points();

        const auto &unit_pts = cache.unit_points_;
        auto &ref_pts = cache.template get_der<_Point>();
        for (int pt = 0 ; pt < n_pts ; ++pt)
        {
            const auto &unit_pt = unit_pts[pt];
            auto &ref_pt = ref_pts[pt];

            for (const auto dir : Topology::active_directions)
                ref_pt[dir] = unit_pt[dir] * dilate[dir] + translate[dir];
        }

        flags.template set_filled<_Point>(true);
    }

    if (flags.template fill<_W_Measure>())
    {
        cache.template get_der<_W_Measure>() = elem.template get_measure<k>(j) * cache.unit_weights_;
        flags.template set_filled<_W_Measure>(true);
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
