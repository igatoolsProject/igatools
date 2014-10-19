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

#ifndef MAPPING_UNIFORM_QUAD_CACHE_H_
#define MAPPING_UNIFORM_QUAD_CACHE_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/value_flags_handler.h>
#include <igatools/base/quadrature.h>
#include <igatools/utils/value_table.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/geometry/mapping.h>

IGA_NAMESPACE_OPEN

/**
 * Global BSplineSpace uniform quadrature
 * computational optimization cache, storing the interval length
 * in each direction
 */
template<int dim_, int codim_ = 0>
class MappingUniformQuadCache : public GridElementHandler<dim_>
{
private:
    using base_t = GridElementHandler<dim_>;
    using Map = Mapping<dim_, codim_>;
    using ElementIterator = typename Map::ElementIterator;
public:
    static const int dim = dim_;

    MappingUniformQuadCache(std::shared_ptr<const Map> map,
                            const ValueFlags flag,
                            const Quadrature<dim> &quad);


protected:
    using ElementAccessor = typename Map::ElementAccessor;

    void init_element_cache(ElementAccessor &elem);

    void fill_element_cache(ElementAccessor &elem);

    //  void fill_face_cache(ElementAccessor &elem, const int face);

public:
    void init_element_cache(ElementIterator &elem);
    void fill_element_cache(ElementIterator &elem);
    /**
     * Fills the ElementIterator face_cache
     * element dependent part
     */
    //void fill_face_cache(ElementIterator &elem, const int face);

    //void print_info(LogStream &out) const;

private:
    ValueFlags flags_;

    Quadrature<dim> quad_;
};

IGA_NAMESPACE_CLOSE

#endif
