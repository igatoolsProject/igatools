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

#include <igatools/geometry/mapping_lib.h>

#include <igatools/base/exceptions.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>

using std::vector;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
LinearMapping<dim_, codim_>::
LinearMapping(const GradientType &A, const ValueType &b)
    :
    base_t(CartesianGrid<dim>::create()),
    A_(A),
    b_(b)
{}



template<int dim_, int codim_>
LinearMapping<dim_, codim_>::
LinearMapping(const shared_ptr<GridType> grid,
              const GradientType &A, const ValueType &b)
    :
    base_t(grid),
    A_(A),
    b_(b)
{}



template<int dim_, int codim_>
auto
LinearMapping<dim_, codim_>::
create(const GradientType &A, const ValueType &b) -> shared_ptr<base_t>
{
    return (shared_ptr<base_t>(new self_t(A,b)));
}


template<int dim_, int codim_>
auto
LinearMapping<dim_, codim_>::
create(const shared_ptr<GridType> grid,
       const GradientType &A, const ValueType &b) -> shared_ptr<base_t>
{
    return (shared_ptr<base_t>(new self_t(grid, A,b)));
}

template<int dim_, int codim_>
shared_ptr< Mapping<dim_,codim_> >
LinearMapping<dim_, codim_>::
clone() const
{
    return (shared_ptr<Mapping<dim_,codim_>>(new self_t(*this)));
}

template<int dim_, int codim_>
ValueFlags
LinearMapping<dim_, codim_>::
required_flags() const
{
    return ValueFlags::point;
}



template<int dim_, int codim_>
void
LinearMapping<dim_, codim_>::
set_element(const CartesianGridElementAccessor<dim> &elem)
{
    points_ = elem.get_points();
}



template<int dim_, int codim_>
void
LinearMapping<dim_, codim_>::
set_face_element(const Index face_id,
                 const CartesianGridElementAccessor<dim> &elem)
{
    Assert(face_id < UnitElement<dim_>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim_>::faces_per_element));
    face_points_[face_id] = elem.get_face_points(face_id);
}



template<int dim_, int codim_>
void
LinearMapping<dim_, codim_>::
evaluate(vector<ValueType> &values) const
{
    const int num_points = points_.size();
    for (int i = 0; i<num_points; i++)
    {
        const auto &x = points_[i];
        values[i] = action(A_,x);
        values[i] += b_;
    }
}



template<int dim_, int codim_>
void
LinearMapping<dim_, codim_>::
evaluate_gradients(vector<GradientType> &gradients) const
{
    const int num_points = points_.size();
    for (int i = 0; i<num_points; i++)
        gradients[i] = A_;
}



template<int dim_, int codim_>
void
LinearMapping<dim_, codim_>::
evaluate_hessians(vector<HessianType> &hessians) const
{
    const int num_points = points_.size();
    for (int i = 0; i<num_points; i++)
        hessians[i] = 0.;
}



template<int dim_, int codim_>
void
LinearMapping<dim_, codim_>::
evaluate_face(const Index face_id, vector<ValueType> &values) const
{
    AssertThrow(false,ExcNotImplemented());
}



template<int dim_, int codim_>
void
LinearMapping<dim_, codim_>::
evaluate_face_gradients(const Index face_id, vector<GradientType> &gradients) const
{
    AssertThrow(false,ExcNotImplemented());
}



template<int dim_, int codim_>
void
LinearMapping<dim_, codim_>::
evaluate_face_hessians(const Index face_id, vector<HessianType> &hessians) const
{
    AssertThrow(false,ExcNotImplemented());
}



/******************************************************************************/
template<int dim_>
BallMapping<dim_>::
BallMapping(const shared_ptr<GridType> grid)
    : base_t(grid)
{}



template<int dim_>
auto
BallMapping<dim_>::
create(const shared_ptr<GridType> grid) -> shared_ptr<base_t>
{
    return (shared_ptr<base_t> (new self_t(grid)));
}

template<int dim_>
shared_ptr< Mapping<dim_,0> >
BallMapping<dim_>::
clone() const
{
    return (shared_ptr<Mapping<dim_,0>>(new self_t(*this)));
}



template<int dim_>
ValueFlags
BallMapping<dim_>::required_flags() const
{
    return ValueFlags::point;
}



template<int dim_>
void
BallMapping<dim_>::set_element(const CartesianGridElementAccessor<dim> &elem)
{
    points_ = elem.get_points();
    const int n_points = points_.size();

    for (int der = 0; der < order; ++der)
    {
        cos_val[der].resize(n_points);
        sin_val[der].resize(n_points);
    }

    for (int qp = 0; qp < n_points; ++qp)
    {
        sin_val[0][qp][0] = points_[qp][0];
        for (int i = 1; i < dim; ++i)
        {
            sin_val[0][qp][i]   = sin(points_[qp][i]);
            cos_val[0][qp][i-1] = cos(points_[qp][i]);
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
}



template<int dim_>
void
BallMapping<dim_>::set_face_element(const Index face_id,
                                    const CartesianGridElementAccessor<dim> &elem)
{
    Assert(face_id < UnitElement<dim_>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim_>::faces_per_element));
    face_points_[face_id] = elem.get_face_points(face_id);
    auto &points = face_points_[face_id] ;
    const int n_points = points.size();

    auto &f_cos_val = face_cos_val[face_id] ;
    auto &f_sin_val = face_sin_val[face_id] ;

    for (int der = 0; der < order; ++der)
    {
        f_cos_val[der].resize(n_points);
        f_sin_val[der].resize(n_points);
    }

    for (int qp = 0; qp < n_points; ++qp)
    {
        f_sin_val[0][qp][0] = points[qp][0];
        for (int i = 1; i < dim; ++i)
        {
            f_sin_val[0][qp][i]   = sin(points[qp][i]);
            f_cos_val[0][qp][i-1] = cos(points[qp][i]);
        }
        f_cos_val[0][qp][dim-1] = 1;

        for (int der = 1; der < order; ++der)
        {
            auto res = std::div(der,2);
            f_sin_val[der][qp][0] = der>1? 0. : 1.;
            for (int i = 1; i < dim; ++i)
            {
                f_sin_val[der][qp][i] =
                    std::pow(-1, res.quot) *
                    (res.rem == 0? f_sin_val[0][qp][i]: f_cos_val[0][qp][i-1]);
                f_cos_val[der][qp][i-1] = -f_sin_val[der-1][qp][i];
            }
            f_cos_val[der][qp][dim-1] = 1.;
        }
    }
}



template<int dim_>
void
BallMapping<dim_>::evaluate(vector<ValueType> &values) const
{
    const int der = 0;
    const auto &s = sin_val[der];
    const auto &c = cos_val[der];
    const int n_points = points_.size();

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



template<int dim_>
void
BallMapping<dim_>::
evaluate_gradients(vector<GradientType> &gradients) const
{
    const auto &s = sin_val[0];
    const auto &c = cos_val[0];
    const auto &s_p = sin_val[1];
    const auto &c_p = cos_val[1];
    const int n_points = points_.size();

    for (int qp = 0; qp < n_points; ++qp)
    {
        auto &grad = gradients[qp];
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



template<int dim_>
void
BallMapping<dim_>::
evaluate_hessians(vector<HessianType> &hessians) const
{
    const auto &s = sin_val[0];
    const auto &c = cos_val[0];
    const auto &s_p = sin_val[1];
    const auto &c_p = cos_val[1];
    const auto &s_2p = sin_val[2];
    const auto &c_2p = cos_val[2];
    const int n_points = points_.size();

    for (int qp = 0; qp < n_points; ++qp)
    {
        auto &hessian = hessians[qp];
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



template<int dim_>
void
BallMapping<dim_>::
evaluate_face(const Index face_id, vector<ValueType> &values) const
{
    AssertThrow(false,ExcNotImplemented());
}



template<int dim_>
void
BallMapping<dim_>::
evaluate_face_gradients(const Index face_id, vector<GradientType> &gradients) const
{
    AssertThrow(false,ExcNotImplemented());
}



template<int dim_>
void
BallMapping<dim_>::
evaluate_face_hessians(const Index face_id, vector<HessianType> &hessians) const
{
    AssertThrow(false,ExcNotImplemented());
}



/******************************************************************************/
template<int dim_>
SphereMapping<dim_>::
SphereMapping(const shared_ptr<GridType> grid, const Real R_)
    : base_t(grid),
      R(R_)
{}



template<int dim_>
auto
SphereMapping<dim_>::
create(const shared_ptr<GridType> grid) -> shared_ptr<base_t>
{
    return (shared_ptr<base_t>(new self_t(grid)));
}



template<int dim_>
shared_ptr< Mapping<dim_,1> >
SphereMapping<dim_>::
clone() const
{
    return (shared_ptr<Mapping<dim_,1>>(new self_t(*this)));
}


template<int dim_>
ValueFlags
SphereMapping<dim_>::required_flags() const
{
    return ValueFlags::point;
}



template<int dim_>
void
SphereMapping<dim_>::set_element(const CartesianGridElementAccessor<dim> &elem)
{
    points_ = elem.get_points();
    const int n_points = points_.size();

    for (int der = 0; der < order; ++der)
    {
        cos_val[der].resize(n_points);
        sin_val[der].resize(n_points);
    }

    for (int qp = 0; qp < n_points; ++qp)
    {
        sin_val[0][qp][0] = R;
        for (int i = 1; i < base_t::space_dim; ++i)
        {
            sin_val[0][qp][i]   = sin(points_[qp][i-1]);
            cos_val[0][qp][i-1] = cos(points_[qp][i-1]);
        }
        cos_val[0][qp][space_dim-1] = 1.;

        for (int der = 1; der < order; ++der)
        {
            auto res = std::div(der,2);
            for (int i = 1; i < base_t::space_dim; ++i)
            {
                sin_val[der][qp][i] =
                    std::pow(-1, res.quot) *
                    (res.rem == 0? sin_val[0][qp][i]: cos_val[0][qp][i-1]);
                cos_val[der][qp][i-1] = -sin_val[der-1][qp][i];
            }
            cos_val[der][qp][base_t::space_dim-1] = 1.;
        }
    }
}



template<int dim_>
void
SphereMapping<dim_>::set_face_element(const Index face_id,
                                      const CartesianGridElementAccessor<dim> &elem)
{
    Assert(face_id < UnitElement<dim_>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<dim_>::faces_per_element));
    face_points_[face_id] = elem.get_face_points(face_id);
    auto &points = face_points_[face_id] ;
    const int n_points = points_.size();

    auto &f_cos_val = face_cos_val[face_id] ;
    auto &f_sin_val = face_sin_val[face_id] ;

    for (int der = 0; der < order; ++der)
    {
        f_cos_val[der].resize(n_points);
        f_sin_val[der].resize(n_points);
    }

    for (int qp = 0; qp < n_points; ++qp)
    {
        f_sin_val[0][qp][0] = R;
        for (int i = 1; i < base_t::space_dim; ++i)
        {
            f_sin_val[0][qp][i]   = sin(points[qp][i-1]);
            f_cos_val[0][qp][i-1] = cos(points[qp][i-1]);
        }
        f_cos_val[0][qp][space_dim-1] = 1.;

        for (int der = 1; der < order; ++der)
        {
            auto res = std::div(der,2);
            for (int i = 1; i < base_t::space_dim; ++i)
            {
                f_sin_val[der][qp][i] =
                    std::pow(-1, res.quot) *
                    (res.rem == 0? f_sin_val[0][qp][i]: f_cos_val[0][qp][i-1]);
                f_cos_val[der][qp][i-1] = -f_sin_val[der-1][qp][i];
            }
            f_cos_val[der][qp][base_t::space_dim-1] = 1.;
        }
    }
}



template<int dim_>
void
SphereMapping<dim_>::evaluate(vector<ValueType> &values) const
{
    const int der = 0;
    const auto &s = sin_val[der];
    const auto &c = cos_val[der];
    const int n_points = points_.size();

    for (int qp = 0; qp < n_points; ++qp)
    {
        auto &x = values[qp];
        double y = 1.;
        for (int i = 0; i < base_t::space_dim; ++i)
        {
            y *= s[qp][i];
            x[i] = y * c[qp][i];
        }
    }
}


template<int dim_>
void
SphereMapping<dim_>::
evaluate_gradients(vector<GradientType> &gradients) const
{
    const auto &s = sin_val[0];
    const auto &c = cos_val[0];
    const auto &s_p = sin_val[1];
    const auto &c_p = cos_val[1];
    const int n_points = points_.size();

    for (int qp = 0; qp < n_points; ++qp)
    {
        auto &grad = gradients[qp];
        grad = 0.;

        for (int i = 0; i < base_t::space_dim-1; ++i)
        {
            for (int j = 1; j < i+2; ++j)
            {
                double djy = 1.;
                for (int k = 0; k < i+1; ++k)
                    djy *= k!=j ? s[qp][k] : s_p[qp][k];
                grad[j-1][i] = djy * (i+1!=j ? c[qp][i] : c_p[qp][i]);
            }
        }

        const int i = base_t::space_dim-1;
        for (int j = 1; j < base_t::space_dim; ++j)
        {
            double djy = 1.;
            for (int k = 0; k < i+1; ++k)
                djy *= k!=j ? s[qp][k] : s_p[qp][k];
            grad[j-1][i] = djy;
        }
    }
}


template<int dim_>
void
SphereMapping<dim_>::
evaluate_hessians(vector <HessianType> &hessians) const
{
    const auto &s = sin_val[0];
    const auto &c = cos_val[0];
    const auto &s_p = sin_val[1];
    const auto &c_p = cos_val[1];
    const auto &s_2p = sin_val[2];
    const auto &c_2p = cos_val[2];
    const int n_points = points_.size();

    for (int qp = 0; qp < n_points; ++qp)
    {
        auto &hessian = hessians[qp];
        hessian = 0.;
        for (int i = 0; i < base_t::space_dim-1; ++i)
        {
            for (int j = 1; j < i+2; ++j)
            {
                for (int k = 1; k < j+1; ++k)
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

                    hessian[j-1][k-1][i] = d2jy * factor;
                }
            }
        }

        const int i = base_t::space_dim-1;
        for (int j = 1; j < base_t::space_dim; ++j)
            for (int k = 1; k < j+1; ++k)
            {
                double d2jy = 1.;
                for (int l = 1; l < base_t::space_dim; ++l)
                {
                    double factor;
                    if (j==k)
                        factor = l==j ? s_2p[qp][l] : s[qp][l];
                    else
                        factor = (l==j || l==k) ?  s_p[qp][l] : s[qp][l];

                    d2jy *= factor;
                }
                hessian[j-1][k-1][i] = d2jy;
            }


        for (int i = 0; i < base_t::space_dim; ++i)
            for (int j = 0; j < dim; ++j)
                for (int k = 0; k< j; ++k)
                {
                    hessian[k][j][i] = hessian[j][k][i];
                }
    }
}



template<int dim_>
void
SphereMapping<dim_>::
evaluate_face(const Index face_id, vector<ValueType> &values) const
{
    AssertThrow(false,ExcNotImplemented());
}



template<int dim_>
void
SphereMapping<dim_>::
evaluate_face_gradients(const Index face_id, vector<GradientType> &gradients) const
{
    AssertThrow(false,ExcNotImplemented());
}



template<int dim_>
void
SphereMapping<dim_>::
evaluate_face_hessians(const Index face_id, vector<HessianType> &hessians) const
{
    AssertThrow(false,ExcNotImplemented());
}


/******************************************************************************/

CylindricalAnnulus::CylindricalAnnulus(
    const Real r0,
    const Real r1,
    const Real h0,
    const Real h1,
    const Real theta0,
    const Real theta1)
    :
    AnalyticalMapping<3,0>(CartesianGrid<dim>::create()),
    r0_(r0),
    r1_(r1),
    h0_(h0),
    h1_(h1),
    theta0_(theta0),
    theta1_(theta1),
    dR_(r1_-r0_),
    dT_(theta1_-theta0_),
    dH_(h1_-h0_)
{}



CylindricalAnnulus::CylindricalAnnulus(
    const Real r0,
    const Real r1,
    const Real h1,
    const Real theta1)
    :
    CylindricalAnnulus(r0, r1, 0.0, h1, 0.0, theta1)
{}



ValueFlags
CylindricalAnnulus::
required_flags() const
{
    return ValueFlags::point;
}



void
CylindricalAnnulus::
set_element(const CartesianGridElementAccessor<3> &elem)
{
    points_ = elem.get_points();
}



void
CylindricalAnnulus::
set_face_element(const Index face_id,
                 const CartesianGridElementAccessor<3> &elem)
{
    Assert(face_id < UnitElement<3>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<3>::faces_per_element));
    face_points_[face_id] = elem.get_face_points(face_id);
}



void CylindricalAnnulus::evaluate(vector<ValueType> &values) const
{
    const Size num_points = points_.size();

    Assert(Size(values.size()) == num_points,
           ExcDimensionMismatch(values.size(),num_points));
    for (int iPt = 0; iPt < num_points; iPt++)
    {
        auto &F = values[iPt];

        const auto &pt = points_[iPt];

        const Real theta = pt[0];
        const Real r     = pt[1];
        const Real z     = pt[2];

        F[0] = (dR_ * r + r0_) * cos(dT_ * theta);
        F[1] = (dR_ * r + r0_) * sin(dT_ * theta);
        F[2] = h0_ + z * dH_;
    }
}



void CylindricalAnnulus::
evaluate_gradients(vector<GradientType> &gradients) const
{
    const Size num_points = points_.size();

    Assert(Size(gradients.size()) == num_points,
           ExcDimensionMismatch(gradients.size(),num_points));
    for (int iPt = 0; iPt < num_points; iPt++)
    {
        auto &dF = gradients[iPt];

        const auto &pt = points_[iPt];

        const Real theta = pt[0];
        const Real r     = pt[1];

        dF[0][0] = - dT_ * (dR_ * r + r0_) * sin(dT_ * theta);
        dF[0][1] =   dT_ * (dR_ * r + r0_) * cos(dT_ * theta);
        dF[0][2] = 0.0;

        dF[1][0] = dR_ * cos(dT_ * theta);
        dF[1][1] = dR_ * sin(dT_ * theta);
        dF[1][2] = 0.0;

        dF[2][0] = 0.0;
        dF[2][1] = 0.0;
        dF[2][2] = dH_;

    }
}


void CylindricalAnnulus::
evaluate_hessians(vector<HessianType> &hessians) const
{
    const Size num_points = points_.size();

    Assert(Size(hessians.size()) == num_points,
           ExcDimensionMismatch(hessians.size(),num_points));
    for (int iPt = 0; iPt < num_points; iPt++)
    {
        auto &d2F = hessians[iPt];

        const auto &pt = points_[iPt];

        const Real theta = pt[0];
        const Real r     = pt[1];

        d2F[0][0][0] = - dT_ * dT_ * (dR_ * r + r0_) * cos(dT_ * theta);
        d2F[0][0][1] = - dT_ * dT_ * (dR_ * r + r0_) * sin(dT_ * theta);
        d2F[0][0][2] = 0.0;

        d2F[1][0][0] = -dT_ * dR_ * sin(dT_ * theta);
        d2F[1][0][1] =  dT_ * dR_ * cos(dT_ * theta);
        d2F[1][0][2] = 0.0;

        d2F[2][0][0] = 0.0;
        d2F[2][0][1] = 0.0;
        d2F[2][0][2] = 0.0;


        d2F[0][1][0] = - dT_ * dR_ * sin(dT_ * theta);
        d2F[0][1][1] =   dT_ * dR_ * cos(dT_ * theta);
        d2F[0][1][2] = 0.0;

        d2F[1][1][0] = 0.0;
        d2F[1][1][1] = 0.0;
        d2F[1][1][2] = 0.0;

        d2F[2][1][0] = 0.0;
        d2F[2][1][1] = 0.0;
        d2F[2][1][2] = 0.0;


        d2F[0][2][0] = 0.0;
        d2F[0][2][1] = 0.0;
        d2F[0][2][2] = 0.0;

        d2F[1][2][0] = 0.0;
        d2F[1][2][1] = 0.0;
        d2F[1][2][2] = 0.0;

        d2F[2][2][0] = 0.0;
        d2F[2][2][1] = 0.0;
        d2F[2][2][2] = 0.0;
    }
}



void CylindricalAnnulus::
evaluate_face(const Index face_id, vector<ValueType> &values) const
{
    AssertThrow(false,ExcNotImplemented());
}



void CylindricalAnnulus::
evaluate_face_gradients(const Index face_id, vector<GradientType> &gradients) const
{
    AssertThrow(false,ExcNotImplemented());
}



void CylindricalAnnulus::
evaluate_face_hessians(const Index face_id, vector<HessianType> &hessians) const
{
    AssertThrow(false,ExcNotImplemented());
}




auto
CylindricalAnnulus::
create(
    const Real r0,
    const Real r1,
    const Real h0,
    const Real h1,
    const Real theta0,
    const Real theta1) -> shared_ptr<base_t>
{
    return (shared_ptr<base_t>(new self_t(r0, r1, h0, h1, theta0, theta1)));
}


shared_ptr< Mapping<3,0> >
CylindricalAnnulus::
clone() const
{
    return shared_ptr<Mapping<3,0>>(new self_t(*this));
}

IGA_NAMESPACE_CLOSE


#include <igatools/geometry/mapping_lib.inst>
