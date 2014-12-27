//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

#include <igatools/basis_functions/bernstein_extraction.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using std::endl;
using std::array;

using std::shared_ptr;
using std::make_shared;
using boost::numeric::ublas::matrix_row;

IGA_NAMESPACE_OPEN

template<int dim, int range, int rank>
auto
BernsteinExtraction<dim, range, rank>::
compute(const matrix &M_j_1,
        typename vector<Real>::const_iterator  y,
        const Real a,
        const Real b) -> matrix
{
    const int j = M_j_1.size1() + 1;
    matrix M_j(j,j);

    vector<Real> alpha(j);
    vector<Real> one_alpha(j,1);
    vector<Real> beta(j, b-a);

    for (int k = 0; k < j; ++k)
    {
        alpha[k] = (y[k+j] - a) / (y[k+j]-y[k]);
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



template<int dim, int range, int rank>
void
BernsteinExtraction<dim, range, rank>::
print_info(LogStream &out) const
{
    int c=0;
    for (const auto &comp : ext_operators_)
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
BernsteinExtraction<dim, range, rank>::
fill_extraction(const int m,
                const vector<Real>    &knots,
                const vector<Real>    &rep_knots,
                const vector<Index>   &acum_mult) -> vector<matrix>
{
    const int n_elem = knots.size()-1;

    vector<matrix>  operators(n_elem, matrix(m,m));
    // const auto &x = knots;
    const auto &y = rep_knots;

    auto x = knots;
    x[0] = y[m-1];
    x[n_elem] = *(y.end()-m);
    for (int n=0; n < n_elem; ++n)
    {
        const auto a = x[n];
        const auto b = x[n+1];

        matrix M(1,1);
        M(0,0) = 1/(b-a);
        for (int k = m-2; k>=0; --k)
        {
            const int s = acum_mult[n] + k;

            auto M1 = compute(M, y.begin()+s, a, b);
            M.assign_temporary(M1);
        }

        //Normalized
        auto M2(M);
        const int s = acum_mult[n];
        for (int k = 0; k < m; ++k)
        {
            matrix_row<matrix> mr(M2, k);
            mr *= (y[s+k+m]-y[s+k]);
        }
        operators[n] = M2;
    }
    return operators;
}



template<int dim, int range, int rank>
BernsteinExtraction<dim, range, rank>::
BernsteinExtraction(std::shared_ptr<CartesianGrid<dim> > grid,
                    const KnotsTable &rep_knots,
                    const MultiplicityTable &acum_mult,
                    const DegreeTable &deg)
{
    for (const int &i : ext_operators_.get_active_components_id())
    {
        for (int j = 0; j < dim; ++j)
        {
            const int m = deg[i][j] + 1;
            auto opers =
                fill_extraction(m,
                                grid->get_knot_coordinates(j),
                                rep_knots[i].get_data_direction(j),
                                acum_mult[i].get_data_direction(j));
            ext_operators_[i].copy_data_direction(j,opers);
        }
    }
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bernstein_extraction.inst>

