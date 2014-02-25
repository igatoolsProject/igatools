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

IGA_NAMESPACE_OPEN

/**
 * @brief Multiplicity of knot vectors.
 *
 * It is used in conjunction with the knot vectors
 * to initiate a BSpline space.
 *
 * @author pauletti, 2013
 *
 */
template<int dim>
class Multiplicity : public CartesianProductArray<Size,dim>
{
public:
    /** We use the father's constructors. */
    using CartesianProductArray<Size,dim>::CartesianProductArray;

    /**
     * Fill the multiplicy for the maximum possible regularity
     *  of the given number of knots
     */
    void fill_max_regularity(const int degree);

    /**
     * Fill the multiplicy for the maximum possible regularity
     *  of the given number of knots
     */
    void fill_max_regularity(const int_array<dim> degree);

    /**
     * Computes the cumulative multiplicity
     */
    Multiplicity<dim> accumulate() ;

    /**
      * @todo (Nov 9 2013, antolin): document.
     */
    Multiplicity<dim> compute_index_space_offset(const std::array<int,dim> &degree) ;

};

IGA_NAMESPACE_CLOSE


#endif // #ifndef MULTIPLICITY_H_
