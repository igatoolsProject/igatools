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


#ifndef __SPACE_TOOLS_H_
#define __SPACE_TOOLS_H_

#include <igatools/base/config.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_slice.h>
#include <igatools/base/function.h>

#include <igatools/basis_functions/physical_space.h>
#include <boost/optional.hpp>


#include <memory>
#include <map>

IGA_NAMESPACE_OPEN


/**
 * This namespace collect functions that work on the
 * physical spaces and fields or functions.
 * Such as:
 * - projections
 * - interpolations
 * - error computations
 */
namespace space_tools
{

/**
 * Determine the knot span index.
 *
 * @return The knot span index of the value @p u in the knot vector @p U.
 * @param[in] p Degree.
 * @param[in] u Knot values for which the span is requested.
 * @param[in] U Knot vector with repeated values.
 *
 * @note The implementation of this function is based on "The NURBS Book" Algorithm A2.1
 */
Index find_span(
    const int p,
    const Real u,
    const std::vector<Real> &U);


//TODO the order of parameters should be consistent
/**
 * Computes an integral norm of the difference between two functions.
 * In this case one function is a Function and the other one an IG field.
 * @note mostly use to compute the convergence rates when the exact solution is known.
 * @todo document a little more
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
Real integrate_difference(std::shared_ptr<const typename Space::Func> exact_solution,
                          std::shared_ptr<const Space> space,
                          const Quadrature< Space::dim > &quad,
                          const Norm &norm_flag,
                          const Vector<la_pack> &solution_coefs,
                          std::vector< Real > &element_error);




// TODO (pauletti, Jun 18, 2014):use space::Function
// TODO (pauletti, Jun 18, 2014): use a quadrature table, and use a default quad if none provided
/**
 * Perform an (L2)-Projection the function @p func
 * onto the space @p space using the quadrature rule @p quad.
 *  The projection is a numerical vector (the coefficients of
 *  the projected function)
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
Vector<la_pack>
projection_l2(
    const typename Space::Func &func,
    std::shared_ptr<const Space> space,
    const Quadrature<Space::dim> &quad
);

/**
 * Projects (using the L2 scalar product) a function to the whole or part
 * of the boundary of the domain.
 * The piece of the domain is indicated by the boundary ids and the
 * projection is computed using the provided quadrature rule.
 *
 * The projected function is returned in boundary_values, a map containing all
 * indices of degrees of freedom at the boundary and the computed coefficient value
 * for this degree of freedom.
 *
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
void project_boundary_values(
    const typename Space::Func &func,
    std::shared_ptr<const Space> space,
    const Quadrature<Space::dim-1> &quad,
    const std::set<boundary_id>  &boundary_ids,
    std::map<Index, Real>  &boundary_values);

/**
 * See documentation above.
 */
template<class Space, LAPack la_pack = LAPack::trilinos>
void project_boundary_values(
    const typename Space::Func &func,
    std::shared_ptr<const Space> space,
    const Quadrature<Space::dim-1> &quad,
    const boundary_id bdry_id,
    std::map<Index,Real>  &boundary_values) ;

//TODO: who uses the next function? delete
/**
 * Transform a set (CartesianProductArray) of points from the unit reference cube [0,1]^{dim} to the
 * element-based reference domain. It returns also the interval id of the elements that contains the points.
 * @param[in] reference_patch - Reference patch.
 * @param[in] points_ref - Set of points in the unit cube [0,1]^{dim}
 * @param[out] points_element - Set of point in the element-based reference domain.
 * @param knot_interval_id - ID of the intervals owning the points.
 */
template < int dim >
void reference_to_element(
    const CartesianGrid< dim > &reference_patch,
    const CartesianProductArray< Real, dim> &points_ref,
    CartesianProductArray< Real, dim> &points_element,
    CartesianProductArray< int, dim> &knot_interval_id) ;

} ;



IGA_NAMESPACE_CLOSE


#endif /* __ERROR_EVALUATION_H_ */
