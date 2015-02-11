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


#include <igatools/basis_functions/bernstein_basis.h>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using boost::numeric::ublas::row;
using boost::numeric::ublas::matrix;
using boost::math::binomial_coefficient;

IGA_NAMESPACE_OPEN



boost::numeric::ublas::vector<Real>
BernsteinBasis::evaluate(const int p, const Real x)
{
    Assert(x >= 0.0 && x <= 1.0,
           ExcMessage("Point not in the unit interval [0,1]"));
    Assert(p >= 0, ExcLowerRange(p,0));

    const int n_basis = p + 1 ;

    boost::numeric::ublas::vector<Real> B(n_basis);

    boost::numeric::ublas::scalar_vector<Real> ones(n_basis,1.0);
    boost::numeric::ublas::vector<Real> t(ones);
    boost::numeric::ublas::vector<Real> one_t(ones);

    for (int k = 1 ; k < n_basis ; ++k)
        for (int i = k ; i < n_basis ; ++i)
        {
            t(i)     *= x;
            one_t(i) *= 1.-x;
        }

    for (int i = 0 ; i < n_basis ; ++i)
    {
        Real C = binomial_coefficient<Real>(p, i);
        B(i) = C * t(i) * one_t(p-i) ;
    }

    return B;
}

boost::numeric::ublas::vector<Real>
BernsteinBasis::derivative(
    const int order,
    const int p,
    const Real x)
{
    Assert(x >= 0.0 && x <= 1.0,
           ExcMessage("Point not in the unit interval [0,1]"));

    Assert(p >= 0, ExcLowerRange(p,0));

    Assert(order >= 0, ExcLowerRange(order,0));

    const int n_basis = p + 1 ;


    if (order > 0)
    {
        /*
         * To compute derivatives we use the recusion formula
         * dB^k = p* ( dB^{k-1}_{i-1} - dB^{k-1}_{i}).
         * To stop the recusion we specialize to the function evaluation in
         * derivative<0>.
         */
        if (p==0)
            return boost::numeric::ublas::zero_vector<Real>(n_basis);

        boost::numeric::ublas::vector<Real> dB(n_basis);
        boost::numeric::ublas::vector<Real> B = BernsteinBasis::derivative(order-1,p-1,x);

        dB(0) = - B(0);
        dB(p) =   B(p-1);
        for (int i = 1 ; i < p ; ++i)
            dB(i) = B(i-1) - B(i);

        dB *= p;

        return dB;
    } // end if (order > 0)
    else
    {
        return BernsteinBasis::evaluate(p,x);
    } // end if (order == 0)
}



matrix<Real>
BernsteinBasis::evaluate(const int p,  const vector<Real> &points)
{
    /*
     * First we compute 2 table where each colum we store
     * 1 1 1, t t t, t^2 t^2 t^2, ...
     * 1 1 1, 1-t 1-t 1-t, (1-t)^2 (1-t)^2 (1-t)^2, ...
     *
     */
    const int n_points = points.size() ;
    const int n_basis  = p + 1 ;

#ifndef NDEBUG
    for (int i = 0 ; i < n_points ; ++i)
        Assert(points[i] >= 0.0 && points[i] <= 1.0,
               ExcMessage("Point " + std::to_string(i) + "not in the unit interval [0,1]"));
#endif

    matrix<Real> B(n_basis, n_points);

    boost::numeric::ublas::scalar_matrix<Real> ones(n_basis,n_points,1.);
    matrix<Real> t(ones);
    matrix<Real> one_t(ones);


    for (int k = 1 ; k < n_basis ; ++k)
        for (int i = k ; i < n_basis ; ++i)
            for (int j = 0 ; j < n_points ; ++j)
            {
                t(i,j)   *= points[j];
                one_t(i,j) *= 1.-points[j];
            }


    for (int i = 0 ; i < n_basis ; ++i)
    {
        Real C = binomial_coefficient<Real>(p, i);
        for (int j = 0 ; j < n_points ; ++j)
            B(i,j) = C * t(i,j) * one_t(p-i,j) ;
    }

    return (B);
}

matrix<Real>
BernsteinBasis::derivative(
    const int order,
    const int p,
    const vector< Real > &points)
{
    if (order > 0)
    {
        /*
         * To compute derivatives we use the recusion formula
         * dB^k = p* ( dB^{k-1}_{i-1} - dB^{k-1}_{i}).
         * To stop the recusion we specialize to the function evaluation in
         * derivative<0>.
         */
        const int n_points = points.size() ;

#ifndef NDEBUG
        for (int i = 0 ; i < n_points ; ++i)
            Assert(points[i] >= 0.0 && points[i] <= 1.0,
                   ExcMessage("Point " + std::to_string(i) + "not in the unit interval [0,1]"));
#endif

        if (p==0)
            return (boost::numeric::ublas::zero_matrix<Real>(p+1,n_points));

        matrix<Real> dB(p+1, n_points);
        matrix<Real> B(p, n_points);
        B = derivative(order-1,p-1,points);

        row(dB,0) = - row(B,0);
        row(dB,p) =   row(B,p-1);
        for (int i = 1 ; i < p ; ++i)
            row(dB,i) = row(B,i-1) - row(B,i);

        dB *= p;
        return (dB);
    }
    else
    {
        return evaluate(p, points) ;
    }
}


IGA_NAMESPACE_CLOSE
