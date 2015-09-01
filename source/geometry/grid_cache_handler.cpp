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

#include <igatools/geometry/grid_cache_handler.h>

using std::shared_ptr;

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
set_flags(const Flags &flag)
{
  flags_[k] = flag;

#if 0
  cacheutils::extract_sub_elements_data<k>(quad_all_sub_elems_) = quad;
#endif
}

template <int dim>
void
GridElementHandler<dim>::
set_flags(const topology_variant &sdim,
          const Flags &flag)
{
  auto disp = SetFlagsDispatcher(flag, flags_);
  boost::apply_visitor(disp, sdim);
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
template <int sdim>
void
GridElementHandler<dim>::
init_cache(ElementAccessor &elem,
           std::shared_ptr<const Quadrature<sdim>> quad) const
{
    Assert(quad != nullptr,ExcNullPtr());

    auto &q = elem.quad_list_.template get_quad<sdim>();
    q = quad;

    auto &cache = elem.all_sub_elems_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::CacheType;
        cache = std::make_shared<Cache>();
    }

    for (auto &s_id: UnitElement<dim>::template elems_ids<sdim>())
    {
        auto &s_cache = cache->template get_sub_elem_cache<sdim>(s_id);
        s_cache.resize(flags_[sdim], quad->get_num_points());

    }
}




template <int dim>
template <int sdim>
void
GridElementHandler<dim>::
fill_cache(ElementAccessor &elem, const int j) const
{
  using _Point = typename ElementAccessor::_Point;
  using _Weight = typename ElementAccessor::_Weight;
  Assert(elem.all_sub_elems_cache_ != nullptr, ExcNullPtr());
  auto &cache = elem.all_sub_elems_cache_->template get_sub_elem_cache<sdim>(j);

  const auto &s_quad = elem.quad_list_.template get_quad<sdim>();

  if (cache.template status_fill<_Point>())
  {
    auto quad = extend_sub_elem_quad<sdim,dim>(*s_quad, j);

    const auto translate = elem.vertex(0);
    const auto dilate    = elem.template get_side_lengths<dim>(0);
    quad.dilate(dilate);
    quad.translate(translate);
    auto &ref_pts = cache.template get_data<_Point>();
    ref_pts = quad.get_points();
    cache.template set_status_filled<_Point>(true);
  }

  if (cache.template status_fill<_Weight>())
  {
    cache.template get_data<_Weight>() =
      elem.template get_measure<sdim>(j) * s_quad->get_weights();
    cache.template set_status_filled<_Weight>(true);
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
  out.begin_item("Flags for each subdimension");
  // flags_.print_info(out);
  out.end_item();

}


#ifdef SERIALIZATION
template <int dim>
template<class Archive>
void
GridElementHandler<dim>::
serialize(Archive &ar, const unsigned int version)
{
  using namespace boost::serialization;
  auto non_const_grid = std::const_pointer_cast<CartesianGrid<dim>>(grid_);
  ar &boost::serialization::make_nvp("grid_",non_const_grid);
  grid_ = non_const_grid;
  Assert(grid_ != nullptr,ExcNullPtr());

  ar &make_nvp("flags_",flags_);
}
#endif // SERIALIZATION


IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_cache_handler.inst>
