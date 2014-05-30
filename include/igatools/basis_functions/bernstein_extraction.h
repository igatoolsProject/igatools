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

#ifndef __BERNSTEIN_EXTRACTION_H_
#define __BERNSTEIN_EXTRACTION_H_

#include <igatools/base/config.h>

#include <igatools/basis_functions/space_spec.h>
#include <boost/numeric/ublas/matrix.hpp>

IGA_NAMESPACE_OPEN

/**
 * A spline function restricted to each interval determined by
 * the knots is a polynomial of order m.
 *
 * In particular each B-spline can be expressed a linear combination
 * of the Berstein polynomial.
 *
 * This class computes and stores theses coefficients.
 */
template<int dim, int range = 1, int rank = 1>
class BernsteinExtraction
{
public:
    using matrix = typename boost::numeric::ublas::matrix<Real>;
    using Space = SplineSpace<dim, range, rank>;
    using DegreeTable = typename Space::DegreeTable;
    using KnotsTable = typename Space::KnotsTable;
    using MultiplicityTable = typename Space::MultiplicityTable;

private:
    using Operators = CartesianProductArray<matrix, dim>;
    using OperatorsTable = typename Space::template ComponentContainer<Operators>;

public:
    /**
     * Construct the extraction operators.
     */
    BernsteinExtraction(std::shared_ptr<CartesianGrid<dim> > grid,
                       const KnotsTable &rep_knots,
                       const MultiplicityTable &acum_mult,
                       const DegreeTable &deg);

    /**
     * Print the class content
     */
    void print_info(LogStream &out) const;

    const std::vector<matrix> &get_operator(const int comp, const int  dir) const
    {
        return ext_operators(comp).get_data_direction(dir);
    }
private:
    std::vector<matrix>
    fill_extraction(const int m,
                    const std::vector<Real>    &knots,
                    const std::vector<Real>    &rep_knots,
                    const std::vector<Index>   &acum_mult);

    matrix compute(const matrix &M_j_1,
                   typename std::vector<Real>::const_iterator  y,
                   const Real a,
                   const Real b);

private:
    OperatorsTable ext_operators;

};


IGA_NAMESPACE_CLOSE

#endif
