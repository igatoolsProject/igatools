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
#include <igatools/geometry/unit_element.h>


using std::shared_ptr;
using std::array;

IGA_NAMESPACE_OPEN

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
void init_unif_caches(const GridFlags &flag, const Quad &quad, std::tuple<Args...> &t)
{
    const int dim = Quad::dim;
    const int low = dim==0? 0 : dim-num_sub_elem;
    UniformQuadFunc f;
    TupleFunc1<UniformQuadFunc, GridFlags, Quad, decltype(t), sizeof...(Args), low>::apply_func(f, flag, quad, t);
}
};



template <int dim_>
GridElementHandler<dim_>::
GridElementHandler(shared_ptr<GridType> grid)
    :
    grid_(grid),
    lengths_(grid->get_element_lengths())
{}



template <int dim_>
template<int k>
void
GridElementHandler<dim_>::
reset(const ValueFlags flag,
      const Quadrature<k> &quad)
{
    flags_[k] = flag;
    auto &quad_k = std::get<k>(quad_);
    quad_k = quad;
}



template <int dim_>
void
GridElementHandler<dim_>::
init_all_caches(ElementAccessor &elem)
{
    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }
    init_unif_caches(flags_[dim], std::get<dim>(quad_), cache->values_);
}



template <int dim_>
template <int k>
void
GridElementHandler<dim_>::
init_cache(ElementAccessor &elem)
{
    auto &cache = elem.local_cache_;
    if (cache == nullptr)
    {
        using Cache = typename ElementAccessor::LocalCache;
        cache = shared_ptr<Cache>(new Cache);
    }

    for (auto &s_id: UnitElement<dim>::template elems_ids<k>())
    {
        auto &s_cache = cache->template get_value_cache<k>(s_id);
        auto &quad = std::get<k>(quad_);
        s_cache.resize(flags_[k], extend_sub_elem_quad<k, dim>(quad, s_id));
    }
}



template <int dim_>
void
GridElementHandler<dim_>::
init_element_cache(ElementAccessor &elem)
{
    init_cache<dim>(elem);
}



template <int dim_>
void
GridElementHandler<dim_>::
init_element_cache(ElementIterator &elem)
{
    init_element_cache(*elem);
}



template <int dim_>
template <int k>
void
GridElementHandler<dim_>::
fill_cache(ElementAccessor &elem, const int j)
{
    Assert(elem.local_cache_ != nullptr, ExcNullPtr());
    auto &cache = elem.local_cache_->template get_value_cache<k>(j);

    const auto &index = elem.get_tensor_index();
    const TensorIndex<k> active(UnitElement<dim>::template get_elem<k>(j).active_directions);

    auto &flags = cache.flags_handler_;

    auto meas = lengths_.template sub_tensor_product<k>(index, active);

    if (flags.fill_measures())
    {
        cache.measure_ = meas;
        flags.set_measures_filled(true);
    }
    if (flags.fill_lengths())
    {
        cache.lengths_ = lengths_.cartesian_product(index);
        flags.set_lengths_filled(true);
    }

    cache.set_filled(true);
}



template <int dim_>
void
GridElementHandler<dim_>::
fill_element_cache(ElementAccessor &elem)
{
    fill_cache<dim>(elem, 0);
}



template <int dim_>
void
GridElementHandler<dim_>::
fill_element_cache(ElementIterator &elem)
{
    fill_element_cache(*elem);
}



template <int dim_>
void
GridElementHandler<dim_>::
print_info(LogStream &out) const
{
    out.begin_item("Lengths:");
    lengths_.print_info(out);
    out.end_item();
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/grid_element_handler.inst>
