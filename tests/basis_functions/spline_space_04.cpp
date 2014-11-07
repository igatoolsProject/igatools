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
#include <igatools/basis_functions/spline_space.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using boost::numeric::ublas::matrix;
using boost::numeric::ublas::matrix_row;

matrix<Real> compute(const matrix<Real> &M_j_1,
                     typename vector<Real>::const_iterator  y,
                     const Real a,
                     const Real b)
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
        M_j(k, l) = M_j(k, l-1) + beta[k] * (M_j_1(k-1, l-1) - M_j_1(k, l-1));
    }
    //k = j-1
    M_j(j-1, l) = M_j(j-1, l-1) + beta[j-1] * M_j_1(j-2, j-2);

    return M_j;
}

void fill_extraction(const int degree,
                     const vector<Real>    &knots,
                     const vector<Real>    &rep_knots,
                     const vector<Index>   &acum_mult)
//,                    vector<matrix<Real>>  &extraction_operators)
{
    // interval n
    const int n=0;
    const int m = degree+1;

    const auto &x = knots;
    const auto &y = rep_knots;

    const auto a = x[n];
    const auto b = x[n+1];

    matrix<Real> M(1,1);
    M(0,0) = 1/(b-a);
    for (int j = 2; j<=m; ++j)
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
    out << M << endl;

}




int main()
{
    out.depth_console(10);

    {
        int degree = 1;
        vector<Real>    knots = {0,1};
        vector<Real>    rep_knots = {0,0,1,1};
        vector<Index>   acum_mult = {0,2,4};

        fill_extraction(degree,knots,rep_knots, acum_mult);
    }


    {
        int degree = 2;
        vector<Real>    knots = {0,1};
        vector<Real>    rep_knots = {0,0,0,1,1,1};
        vector<Index>   acum_mult = {0,3,6};

        fill_extraction(degree,knots,rep_knots, acum_mult);
    }


    {
        int degree = 3;
        vector<Real>    knots = {0,1};
        vector<Real>    rep_knots = {0,0,0,0,1,1,1,1};
        vector<Index>   acum_mult = {0,4,8};

        fill_extraction(degree,knots,rep_knots, acum_mult);
    }


    return 0;
}
