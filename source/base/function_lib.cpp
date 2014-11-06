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

#include <igatools/base/function_lib.h>
#include <igatools/base/function_element.h>

IGA_NAMESPACE_OPEN

namespace functions
{

template<int dim, int codim, int range, int rank>
ConstantFunction<dim, codim, range, rank>::
ConstantFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
                 const Value &b,
                 const NewValueFlags &flag,
                 const Quadrature<dim> &quad)
    :
    parent_t::FormulaFunction(grid, flag, quad),
    b_(b)
{}



template<int dim, int codim, int range, int rank>
auto
ConstantFunction<dim, codim, range, rank>::
create(std::shared_ptr<const CartesianGrid<dim>> grid,
       const Value &b,
       const NewValueFlags &flag,
       const Quadrature<dim> &quad) ->  std::shared_ptr<base_t>
{
    return std::shared_ptr<base_t>(new self_t(grid, b, flag, quad));
}



template<int dim, int codim, int range, int rank>
auto
ConstantFunction<dim, codim, range, rank>::
evaluate_0(const ValueVector<Point> &points,
           ValueVector<Value> &values) const -> void
{
    for (auto &val : values)
        val =  b_;
}



template<int dim, int codim, int range, int rank>
auto
ConstantFunction<dim, codim, range, rank>::
evaluate_1(const ValueVector<Point> &points,
           ValueVector<Derivative<1>> &values) const -> void
{
    for (auto &val : values)
        val = 0.;
}



template<int dim, int codim, int range, int rank>
auto
ConstantFunction<dim, codim, range, rank>::
evaluate_2(const ValueVector<Point> &points,
           ValueVector<Derivative<2>> &values) const -> void
{
    for (auto &val : values)
        val = 0.;
}



//------------------------------------------------------------------------------
template<int dim, int codim, int range>
LinearFunction<dim, codim, range>::
LinearFunction(std::shared_ptr<const CartesianGrid<dim>> grid,
               const Gradient &A, const Value &b)
    :
    parent_t::FormulaFunction(grid),
    A_(A),
    b_(b)
{}



template<int dim, int codim, int range>
auto
LinearFunction<dim, codim, range>::
create(std::shared_ptr<const CartesianGrid<dim>> grid,
       const Gradient &A, const Value &b) ->  std::shared_ptr<base_t>
{
    return std::shared_ptr<base_t>(new self_t(grid, A, b));
}



template<int dim, int codim, int range>
auto
LinearFunction<dim, codim, range>::
evaluate_0(const ValueVector<Point> &points,
           ValueVector<Value> &values) const -> void
{
    auto point = points.begin();
    for (auto &val : values)
    {
        val = action(A_, *point) + b_;
        ++point;
    }
}



template<int dim, int codim, int range>
auto
LinearFunction<dim, codim, range>::
evaluate_1(const ValueVector<Point> &points,
           ValueVector<Derivative<1>> &values) const -> void
{
    for (auto &val : values)
        val = A_;
}

template<int dim, int codim, int range>
auto
LinearFunction<dim, codim, range>::
evaluate_2(const ValueVector<Point> &points,
           ValueVector<Derivative<2>> &values) const -> void
{
    for (auto &val : values)
        val = 0.;
}


//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

template<int dim>
BallFunction<dim>::
BallFunction(std::shared_ptr<const CartesianGrid<dim>> grid)
    :
    parent_t::FormulaFunction(grid)
{}



template<int dim>
auto
BallFunction<dim>::
create(std::shared_ptr<const CartesianGrid<dim>> grid) ->  std::shared_ptr<base_t>
{
    return std::shared_ptr<base_t>(new self_t(grid));
}



template<int dim>
template<int order>
auto
BallFunction<dim>::get_aux_vals(const ValueVector<Point> &points) const
{
    std::array<std::array<vector<std::array<double, dim> >, order>, 2> val_table;
    auto &cos_val = val_table[0];
    auto &sin_val = val_table[1];

    const int n_points = points.size();

    for (int der = 0; der < order; ++der)
    {
        cos_val[der].resize(n_points);
        sin_val[der].resize(n_points);
    }

    for (int qp = 0; qp < n_points; ++qp)
    {
        sin_val[0][qp][0] = points[qp][0];
        for (int i = 1; i < dim; ++i)
        {
            sin_val[0][qp][i]   = sin(points[qp][i]);
            cos_val[0][qp][i-1] = cos(points[qp][i]);
        }
        cos_val[0][qp][dim-1] = 1;

        for (int der = 1; der < order; ++der)
        {
            auto res = std::div(der,2);
            sin_val[der][qp][0] = der>1? 0. : 1.;
            for (int i = 1; i < dim; ++i)
            {
                sin_val[der][qp][i] =
                    std::pow(-1, res.quot) *
                    (res.rem == 0? sin_val[0][qp][i]: cos_val[0][qp][i-1]);
                cos_val[der][qp][i-1] = -sin_val[der-1][qp][i];
            }
            cos_val[der][qp][dim-1] = 1.;
        }
    }
    return val_table;
}


template<int dim>
auto
BallFunction<dim>::
evaluate_0(const ValueVector<Point> &points,
           ValueVector<Value> &values) const -> void
{
    const auto val_table = get_aux_vals<1>(points);
    auto &cos_val = val_table[0];
    auto &sin_val = val_table[1];

    const int der = 0;
    const auto &s = sin_val[der];
    const auto &c = cos_val[der];
    const int n_points = points.size();

    for (int qp = 0; qp < n_points; ++qp)
    {
        auto &x = values[qp];
        double y = 1.;
        for (int i = 0; i < dim; ++i)
        {
            y *= s[qp][i];
            x[i] = y * c[qp][i];
        }
    }
}



template<int dim>
auto
BallFunction<dim>::
evaluate_1(const ValueVector<Point> &points,
           ValueVector<Derivative<1>> &values) const -> void
{
    const auto val_table = get_aux_vals<2>(points);
    auto &cos_val = val_table[0];
    auto &sin_val = val_table[1];

    const auto &s = sin_val[0];
    const auto &c = cos_val[0];
    const auto &s_p = sin_val[1];
    const auto &c_p = cos_val[1];
    const int n_points = points.size();

    for (int qp = 0; qp < n_points; ++qp)
    {
        auto &grad = values[qp];
        grad = 0.;

        for (int i = 0; i < dim-1; ++i)
        {
            for (int j = 0; j < i+2; ++j)
            {
                double djy = 1.;
                for (int k = 0; k < i+1; ++k)
                    djy *= k!=j ? s[qp][k] : s_p[qp][k];
                grad[j][i] = djy * (i+1!=j ? c[qp][i] : c_p[qp][i]);
            }
        }

        const int i = dim-1;
        for (int j = 0; j < dim; ++j)
        {
            double djy = 1.;
            for (int k = 0; k < i+1; ++k)
                djy *= k!=j ? s[qp][k] : s_p[qp][k];
            grad[j][i] = djy;
        }
    }
}



template<int dim>
auto
BallFunction<dim>::
evaluate_2(const ValueVector<Point> &points,
           ValueVector<Derivative<2>> &values) const -> void
{
    const auto val_table = get_aux_vals<3>(points);
    auto &cos_val = val_table[0];
    auto &sin_val = val_table[1];

    const auto &s = sin_val[0];
    const auto &c = cos_val[0];
    const auto &s_p = sin_val[1];
    const auto &c_p = cos_val[1];
    const auto &s_2p = sin_val[2];
    const auto &c_2p = cos_val[2];
    const int n_points = points.size();

    for (int qp = 0; qp < n_points; ++qp)
    {
        auto &hessian = values[qp];
        hessian = 0.;
        for (int i = 0; i < dim-1; ++i)
        {
            for (int j = 0; j < i+2; ++j)
            {
                for (int k = 0; k < j+1; ++k)
                {
                    double d2jy = 1.;
                    for (int l = 0; l < i+1; ++l)
                    {
                        double factor;
                        if (j==k)
                            factor = l==j ? s_2p[qp][l] : s[qp][l];
                        else
                            factor = (l==j || l==k) ?  s_p[qp][l] : s[qp][l];

                        d2jy *= factor;
                    }
                    double factor;
                    if (j==k)
                        factor = (i+1)==j ? c_2p[qp][i] : c[qp][i];
                    else
                        factor = ((i+1)==j || (i+1)==k) ?
                                 c_p[qp][i] : c[qp][i];

                    hessian[j][k][i] = d2jy * factor;
                }
            }
        }

        const int i = dim-1;
        for (int j = 0; j < dim; ++j)
            for (int k = 0; k < j+1; ++k)
            {
                double d2jy = 1.;
                for (int l = 0; l < dim; ++l)
                {
                    double factor;
                    if (j==k)
                        factor = l==j ? s_2p[qp][l] : s[qp][l];
                    else
                        factor = (l==j || l==k) ?  s_p[qp][l] : s[qp][l];

                    d2jy *= factor;
                }
                hessian[j][k][i] = d2jy;
            }


        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                for (int k = 0; k< j; ++k)
                {
                    hessian[k][j][i] = hessian[j][k][i];
                }
    }
}



//------------------------------------------------------------------------------
} // of namespace functions.

IGA_NAMESPACE_CLOSE

#include <igatools/base/function_lib.inst>
