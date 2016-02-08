//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#ifndef BBOX_H_
#define BBOX_H_

#include <igatools/base/config.h>
#include <igatools/base/exceptions.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/safe_stl_array.h>

IGA_NAMESPACE_OPEN

/**
 * @brief Bounding Box, a dim-dimensional rectangular
 * box described by the intervals.
 *
 * eg. BBox<2> box {{0,0.5},{1,2}}.
 */
template<int dim>
class BBox : public SafeSTLArray<SafeSTLArray<Real, 2>, dim>
{
public:
  using SafeSTLArray<SafeSTLArray<Real, 2>, dim>::SafeSTLArray;

  /**
   * @brief Default constructor. It builds a box representing the <tt>dim</tt>-dimensional
   * unit hypercube \f$[0,1]^{dim}\f$.
   */
  BBox();

  /**
   * @brief Translates the BBox by the amount specified by <tt>translation_amount</tt>.
   * @param translation_amount
   */
  void
  translate(const Points<dim> &translation_amount);

  /**
   * @brief Dilates the intervals defining the bounding box.
   *
   * @note The argument <tt>dilation_factor</tt> must be positive in each coordinate direction,
   * otherwise (in Debug mode) an assertion will be raised.
   */
  void
  dilate(const Points<dim> &dilation_factor);


  /**
   * @brief Returns true if the bounding box represents the unit <tt>dim</tt>-dimensional hypercube
   * \f$ [0,1]^{\text{dim}}\f$.
   */
  bool is_unit() const;


  /**
   * @brief Returns the lengths of the box side.
   */
  SafeSTLArray<Real,dim> get_side_lengths() const;

  /**
   * @brief Returns TRUE if the <tt>point</tt> is inside or on the boundary.
   *
   * @note The optional non-negative parameter <tt>eps</tt> is to perform the test on an box that is enlarged
   * by an amount <tt>eps</tt> along all the coordinate direction,
   * in order to fix numerical problems when the <tt>point</tt> is very close to the box boundary.
   */
  bool is_point_inside(const Points<dim> &point,const Real eps = 0.0) const;
};


IGA_NAMESPACE_CLOSE


#endif // BBOX_H
