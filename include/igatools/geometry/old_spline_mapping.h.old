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

#ifndef SPLINE_MAPPING_H_
#define SPLINE_MAPPING_H_

#include <igatools/base/config.h>
#include <igatools/geometry/mapping.h>

IGA_NAMESPACE_OPEN

template <class RefSpace>
class SplineMapping
    : public Mapping<RefSpace::dim, RefSpace::range - RefSpace::dim>
{
protected:
    using base_t = Mapping<RefSpace::dim, RefSpace::range - RefSpace::dim>;

    using base_t::dim;
    using base_t::codim;
    using base_t::space_dim;

    using typename base_t::Point;
    using typename base_t::Value;
    using typename base_t::Gradient;
    using typename base_t::Hessian;
    using typename base_t::GridType;
    using typename base_t::ElementIterator;

    using self_t = SplineMapping<RefSpace>;

public:
    /** Inheriting the constructor from the Mapping class. */
    using base_t::Mapping;


    /** @name Function used to modify the position of the control points */
    ///@{
    /**
     * Sets the control points defining the map.
     * @param[in] control_points - Coordinates of the control points in the Euclidean space.
     */
    virtual void set_control_points(const vector<Real> &control_points) = 0;
    ///@}


    /**
     * Returns a vector containing a copy of the control point values.
     */
    virtual vector<Real> get_control_points() const = 0;



    /** Return the space used to define (or defined by) the SplineMapping.*/
    virtual std::shared_ptr<RefSpace> get_iga_space() = 0;
};

IGA_NAMESPACE_CLOSE

#endif // #ifndef SPLINE_MAPPING_H_



