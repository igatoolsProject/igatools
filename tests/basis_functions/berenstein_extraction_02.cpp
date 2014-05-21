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
/*
 *  Test for bezier extraction
 *  author: pauletti
 *  date:
 *
 */

#include "../tests.h"
#include <igatools/basis_functions/space_spec.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using boost::numeric::ublas::matrix;
using boost::numeric::ublas::matrix_row;


template<int dim, int range = 1, int rank = 1>
class BersteinExtraction
{
public:
    using Space = SpaceSpec<dim, range, rank>;
    using DegreeTable = typename Space::DegreeTable;
    using KnotsTable = typename Space::KnotsTable;
    using MultiplicityTable = typename Space::MultiplicityTable;

private:
    using Operators = CartesianProductArray<matrix<Real>, dim>;
    using OperatorsTable = typename Space::template ComponentContainer<Operators>;

public:

    BersteinExtraction(std::shared_ptr<CartesianGrid<dim> > grid,
                       const KnotsTable &rep_knots,
                       const MultiplicityTable &acum_mult,
                       const DegreeTable &deg);


    vector<matrix<Real>>
    fill_extraction(const int m,
                    const vector<Real>    &knots,
                    const vector<Real>    &rep_knots,
                    const vector<Index>   &acum_mult);

    void print_info(LogStream &out) const;

private:
    matrix<Real> compute(const matrix<Real> &M_j_1,
                         typename vector<Real>::const_iterator  y,
                         const Real a,
                         const Real b);

private:
    OperatorsTable ext_operators;

};



template<int dim, int range, int rank>
auto
BersteinExtraction<dim, range, rank>::
compute(const matrix<Real> &M_j_1,
        typename vector<Real>::const_iterator  y,
        const Real a,
        const Real b) -> matrix<Real>
{
    const int j = M_j_1.size1() + 1;
    matrix<Real> M_j(j,j);

    vector<Real> alpha(j);
    vector<Real> one_alpha(j,1);
    vector<Real> beta(j, b-a);

    for (int k = 0; k < j; ++k)
    {
        alpha[k] = (y[k+j] - a)/(y[k+j]-y[k]);
        one_alpha[k] -= alpha[k];

        beta[k] /= (y[k+j]-y[k]);
    }

    for (int l = 0; l < j-1; ++l)
    {
        //k = 0
        M_j(0, l) = alpha[0] * M_j_1(0, l);
        //k = 1,...,j-2
        for (int k = 1; k < j-1; ++k)
        {
            M_j(k, l) = alpha[k] * M_j_1(k, l) + one_alpha[k] * M_j_1(k-1, l);
        }
        //k = j-1
        M_j(j-1, l) = one_alpha[j-1] * M_j_1(j-2, l);
    }


    const int l = j-1;

    //k = 0
    M_j(0, l) = M_j(0, l-1) - beta[0] * M_j_1(0, l-1);
    //k = 1,...,j-2
    for (int k = 1; k < j-1; ++k)
    {
        M_j(k, l) = M_j(k, l-1) + beta[k] * (M_j_1(k-1, l-1) - M_j_1(k, l-1) );
    }
    //k = j-1
    M_j(j-1, l) = M_j(j-1, l-1) + beta[j-1] * M_j_1(j-2, j-2);

    return M_j;
}



template<int dim, int range, int rank>
void
BersteinExtraction<dim, range, rank>::
print_info(LogStream &out) const
{
    int c=0;
    for (const auto &comp : ext_operators)
    {
        out << "Component[" << c++ << "]: " << endl;
        for (int j = 0; j < dim; ++j)
        {
            out << "Direction[" << j << "]:" << endl;
            for (const auto &M : comp.get_data_direction(j))
                out << M << endl;
        }
    }
}



template<int dim, int range, int rank>
auto
BersteinExtraction<dim, range, rank>::
fill_extraction(const int m,
                const vector<Real>    &knots,
                const vector<Real>    &rep_knots,
                const vector<Index>   &acum_mult) -> vector<matrix<Real>>
{
    const int n_elem = knots.size()-1;

    vector<matrix<Real>>  operators(n_elem, matrix<Real>(m,m));
    const auto &x = knots;
    const auto &y = rep_knots;

    for (int n=0; n < n_elem; ++n)
    {
        const auto a = x[n];
        const auto b = x[n+1];

        matrix<Real> M(1,1);
        M(0,0) = 1/(b-a);
        for(int j = 2; j<=m; ++j)
        {
            const int s = acum_mult[n+1] - j;

            auto M1 = compute(M, y.begin()+s, a, b);
            M.assign_temporary(M1);
        }

        //Normalized
        auto M2(M);
        const int s = acum_mult[n+1] - m;
        for (int k = 0; k < m; ++k)
        {
            matrix_row<matrix<double> > mr(M2, k);
            mr *= (y[s+k+m]-y[s+k]);
        }
        operators[n] = M2;
    }
    return operators;
}



template<int dim, int range, int rank>
BersteinExtraction<dim, range, rank>::
BersteinExtraction(std::shared_ptr<CartesianGrid<dim> > grid,
                   const KnotsTable &rep_knots,
                   const MultiplicityTable &acum_mult,
                   const DegreeTable &deg)
{
    for (int i = 0; i < Space::n_components; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            const int m = deg(i)[j] + 1;
            auto opers =
                    fill_extraction(m,
                                    grid->get_knot_coordinates(j),
                                    rep_knots(i).get_data_direction(j),
                                    acum_mult(i).get_data_direction(j));
            ext_operators(i).copy_data_direction(j,opers);
        }
    }
}



int main()
{
    out.depth_console(10);

    {
        const int dim = 1;
        int degree = 1;

        CartesianProductArray<Real, dim> knots({{0,1,2,3}});
        auto grid = CartesianGrid<dim>::create(knots);

        typename BersteinExtraction<dim>::KnotsTable rep_knots({{0,0,1,2,3,3}});
        typename BersteinExtraction<dim>::MultiplicityTable acum_mult ({{{0,2,3,4,6}}});
        typename BersteinExtraction<dim>::DegreeTable deg{{degree}};
        BersteinExtraction<dim> operators(grid, rep_knots, acum_mult, deg);
        operators.print_info(out);
    }


    {
        const int dim = 2;
        int degree = 1;

        CartesianProductArray<Real, dim> knots({{0,1,2,3}, {3,4,5}});
        auto grid = CartesianGrid<dim>::create(knots);

        typename BersteinExtraction<dim>::KnotsTable
        rep_knots({{0,0,1,2,3,3},{3,3,4,5,5}});
        typename BersteinExtraction<dim>::MultiplicityTable
        acum_mult ({{{0,2,3,4,6}, {0,2,3,5}}});
        typename BersteinExtraction<dim>::DegreeTable deg{{degree, degree}};
        BersteinExtraction<dim> operators(grid, rep_knots, acum_mult, deg);
        operators.print_info(out);
    }




    //
//    {
//        int degree = 2;
//        vector<Real>    knots = {0,1};
//        vector<Real>    rep_knots = {0,0,0,1,1,1};
//        vector<Index>   acum_mult = {0,3,6};
//
//        fill_extraction( degree,knots,rep_knots, acum_mult);
//    }
//
//
//    {
//        int degree = 3;
//        vector<Real>    knots = {0,1};
//        vector<Real>    rep_knots = {0,0,0,0,1,1,1,1};
//        vector<Index>   acum_mult = {0,4,8};
//
//        fill_extraction( degree,knots,rep_knots, acum_mult);
//    }


    return 0;
}
