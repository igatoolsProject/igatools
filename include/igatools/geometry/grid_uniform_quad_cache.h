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

#ifndef GRID_UNIFORM_QUAD_CACHE_H_
#define GRID_UNIFORM_QUAD_CACHE_H_

#include <igatools/base/config.h>
#include <igatools/base/cache_status.h>
#include <igatools/base/value_flags_handler.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/cartesian_grid.h>

template <int dim_>
class GridUniformQuadCache : public CacheStatus
{
    using GridType = CartesianGrid<dim_>;
    using ElementIterator = typename GridType::ElementIterator;

    /**
     * @brief Global CartesianGrid cache, storing the interval length
     * in each direction.
     *
     * For now only a uniform quad is taken care of.
     */
public:
    static const int dim = dim_;

    //Allocates and fill the (global) cache
    GridUniformQuadCache(std::shared_ptr<const GridType> grid,
                         const ValueFlags flag,
                         const Quadrature<dim> &quad)
:
    grid_(grid),
    flags_(flag),
    length_(grid->get_element_lengths())
{}

    //Allocates the ElementIterator element_cache
    void init_element_cache(ElementIterator &elem)
    {
        const auto n_points_direction = quad.get_num_points_direction();
        const Size n_points = n_points_direction.flat_size();

        flags_handler_ = flags_handler;

        if (flags_handler_.fill_points())
        {
            this->unit_points_ = quad.get_points();
            flags_handler_.set_points_filled(true);
        }

        if (flags_handler_.fill_w_measures())
        {
            if (this->w_measure_.size() != n_points)
                this->w_measure_.resize(n_points);

            this->unit_weights_ = quad.get_weights().get_flat_tensor_product();
        }
        else
        {
            Assert(True, ExcMessage("Should not get here"));
            this->w_measure_.clear() ;
            this->unit_weights_.clear() ;
        }
    }

    //Fill the ElementIterator element_cache
    void fill_element_cache(ElementIterator &elem)
    {

    }

private:
    std::shared_ptr<const GridType> grid_;

    GridElemValueFlagsHandler flags_;

    CartesianProductArray<Real, dim> lengths_;

    Quadrature<dim> quad_;
};



#endif /* GRID_UNIFORM_QUAD_CACHE_H_ */
