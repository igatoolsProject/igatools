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

#ifndef MULTIPLICITY_H_
#define MULTIPLICITY_H_

#include <igatools/base/config.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/geometry/cartesian_grid.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Multiplicity of knot vectors without repetition.
 *
  *
 * @author pauletti, 2013-2014
 *
 */
template<int dim, int range, int rank>
class Multiplicity :
        public StaticMultiArray<CartesianProductArray<Size, dim>,range,rank>
{
private:
    using Grid     = CartesianGrid<dim>;

    using T =  CartesianProductArray<Size, dim>;


public:
    using Degrees  = TensorIndex<dim>;
    using DegreeTable = StaticMultiArray<Degrees, range, rank>;


    using parent_t = StaticMultiArray<CartesianProductArray<Size, dim>,range,rank>;


public:
    using parent_t::StaticMultiArray;

    /**
     * Maximun regularity multiplicity vectors associated with
     * the grid and degrees
     */
    explicit Multiplicity(parent_t &mult, const DegreeTable &deg)
    :
      parent_t::StaticMultiArray(mult),
      deg_(deg)
    {}

    /**
     * Maximun regularity multiplicity vectors associated with
     * the grid and degrees
     */
    explicit Multiplicity(std::shared_ptr<const Grid> knots,
                          const DegreeTable &deg,
                          const bool max_reg);

private:
    /**
     * Fill the multiplicy for the maximum possible regularity
     *  of the given number of knots
     */
    void fill_max_regularity();
//
//    /**
//     * Fill the multiplicy for the maximum possible regularity
//     *  of the given number of knots
//     */
//    void fill_max_regularity(const int_array<dim> degree);
//
//    /**
//     * Computes the cumulative multiplicity
//     */
//    parent_t accumulate() ;
public:
    /**
      * @todo (Nov 9 2013, antolin): document.
     */
    parent_t compute_index_space_offset() ;

private:
    DegreeTable deg_;
};



IGA_NAMESPACE_CLOSE





#endif // #ifndef MULTIPLICITY_H_
