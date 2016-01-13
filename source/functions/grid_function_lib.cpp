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

#include <igatools/functions/grid_function_lib.h>

IGA_NAMESPACE_OPEN

namespace grid_functions
{
//------------------------------------------------------------------------------
template<int dim, int range>
ConstantGridFunction<dim,range>::
ConstantGridFunction(const SharedPtrConstnessHandler<GridType> &domain,
                     const Value &b)
  :
  parent_t(domain),
  b_(b)
{}



template<int dim, int range>
auto
ConstantGridFunction<dim,range>::
create(const std::shared_ptr<GridType> &domain,
       const Value &b) ->  std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(domain),b));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;

}

template<int dim, int range>
auto
ConstantGridFunction<dim,range>::
const_create(const std::shared_ptr<const GridType> &domain,
             const Value &b) ->  std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(domain), b));
}

#ifdef MESH_REFINEMENT
template<int dim, int range>
void
ConstantGridFunction<dim,range>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  this->grid_function_previous_refinement_ =
    self_t::const_create(
      this->get_grid()->get_grid_pre_refinement(),b_);
}
#endif // MESH_REFINEMENT


template<int dim, int range>
auto
ConstantGridFunction<dim,range>::
evaluate_0(const ValueVector<GridPoint> &points,
           ValueVector<Value> &values) const -> void
{
  auto point = points.begin();
  for (auto &val : values)
  {
    val = b_;
    ++point;
  }
}



template<int dim, int range>
auto
ConstantGridFunction<dim,range>::
evaluate_1(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<1>> &values) const -> void
{
  for (auto &val : values)
    val = 0.0;
}

template<int dim, int range>
auto
ConstantGridFunction<dim,range>::
evaluate_2(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<2>> &values) const -> void
{
  for (auto &val : values)
    val = 0.0;
}

template<int dim, int range>
void
ConstantGridFunction<dim,range>::
print_info(LogStream &out) const
{
  out.begin_item("ConstantGridFunction<"
                 + std::to_string(dim) + ","
                 + std::to_string(range) + ">");

  out.begin_item("b:");
  out << b_ ;
  out.end_item();

  out << "Name: " << this->name_ << std::endl;

  out.end_item();
}

template<int dim, int range>
auto
ConstantGridFunction<dim,range>::
get_constant_value() const -> const Value &
{
    return b_;
}


//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
template<int dim, int range>
LinearGridFunction<dim,range>::
LinearGridFunction(const SharedPtrConstnessHandler<GridType> &domain,
                   const Derivative<1> &A,
                   const Value &b)
  :
  parent_t(domain),
  A_(A),
  b_(b)
{}



template<int dim, int range>
auto
LinearGridFunction<dim,range>::
create(const std::shared_ptr<GridType> &domain,
       const Derivative<1> &A,
       const Value &b) ->  std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(domain),A,b));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;

}

template<int dim, int range>
auto
LinearGridFunction<dim,range>::
const_create(const std::shared_ptr<const GridType> &domain,
             const Derivative<1> &A,
             const Value &b) ->  std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(domain), A, b));
}

#ifdef MESH_REFINEMENT
template<int dim, int range>
void
LinearGridFunction<dim,range>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  this->grid_function_previous_refinement_ =
    self_t::const_create(
      this->get_grid()->get_grid_pre_refinement(),A_,b_);
}
#endif // MESH_REFINEMENT


template<int dim, int range>
auto
LinearGridFunction<dim,range>::
evaluate_0(const ValueVector<GridPoint> &points,
           ValueVector<Value> &values) const -> void
{
  auto point = points.begin();
  for (auto &val : values)
  {
    val = action(A_, *point) + b_;
    ++point;
  }
}



template<int dim, int range>
auto
LinearGridFunction<dim,range>::
evaluate_1(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<1>> &values) const -> void
{
  for (auto &val : values)
    val = A_;
}

template<int dim, int range>
auto
LinearGridFunction<dim,range>::
evaluate_2(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<2>> &values) const -> void
{
  for (auto &val : values)
    val = 0.;
}

template<int dim, int range>
void
LinearGridFunction<dim,range>::
print_info(LogStream &out) const
{
  out.begin_item("LinearGridFunction<"
                 + std::to_string(dim) + ","
                 + std::to_string(range) + ">");

  out.begin_item("A:");
  out << A_ ;
  out.end_item();

  out.begin_item("b:");
  out << b_ ;
  out.end_item();

  out << "Name: " << this->name_ << std::endl;

  out.end_item();
}



template<int dim, int range>
auto
LinearGridFunction<dim,range>::
get_A () const -> const Derivative<1> &
{
    return A_;
}



template<int dim, int range>
auto
LinearGridFunction<dim,range>::
get_b () const -> const Value &
{
    return b_;
}


//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
template<int dim>
IdentityGridFunction<dim>::
IdentityGridFunction(const SharedPtrConstnessHandler<GridType> &domain)
  :
  parent_t(domain)
{}



template<int dim>
auto
IdentityGridFunction<dim>::
create(const std::shared_ptr<GridType> &domain)
->  std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(domain)));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;
}

template<int dim>
auto
IdentityGridFunction<dim>::
const_create(const std::shared_ptr<const GridType> &domain)
->  std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(domain)));
}


template<int dim>
auto
IdentityGridFunction<dim>::
evaluate_0(const ValueVector<GridPoint> &points,
           ValueVector<Value> &values) const -> void
{
  auto point = points.begin();
  for (auto &val : values)
  {
    val = *point;
    ++point;
  }
}

#ifdef MESH_REFINEMENT
template<int dim>
void
IdentityGridFunction<dim>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  this->grid_function_previous_refinement_ =
    self_t::const_create(
      this->get_grid()->get_grid_pre_refinement());
}
#endif // MESH_REFINEMENT


template<int dim>
auto
IdentityGridFunction<dim>::
evaluate_1(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<1>> &values) const -> void
{
  for (auto &val : values)
  {
    val = 0.0;
    for (int i = 0 ; i < dim ; ++i)
      val[i][i] = 1.0;
  }
}

template<int dim>
auto
IdentityGridFunction<dim>::
evaluate_2(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<2>> &values) const -> void
{
  for (auto &val : values)
    val = 0.0;
}

template<int dim>
void
IdentityGridFunction<dim>::
print_info(LogStream &out) const
{
  out.begin_item("IdentityGridFunction<"
                 + std::to_string(dim) + ">");
  out.end_item();

  out << "Name: " << this->name_ << std::endl;
}


//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template<int dim>
BallGridFunction<dim>::
BallGridFunction(const SharedPtrConstnessHandler<GridType> &grid)
  :
  parent_t::FormulaGridFunction(grid)
{}



template<int dim>
auto
BallGridFunction<dim>::
create(std::shared_ptr<GridType> grid) ->  std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid)));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;
}

template<int dim>
auto
BallGridFunction<dim>::
const_create(std::shared_ptr<const GridType> grid) ->  std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<GridType>(grid)));
}

#ifdef MESH_REFINEMENT
template<int dim>
void
BallGridFunction<dim>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  this->grid_function_previous_refinement_ =
    self_t::const_create(
      this->get_grid()->get_grid_pre_refinement());
}
#endif // MESH_REFINEMENT


template<int dim>
template<int order>
auto
BallGridFunction<dim>::get_aux_vals(const ValueVector<GridPoint> &points) const
{
  SafeSTLArray<SafeSTLArray<SafeSTLVector<SafeSTLArray<double, dim> >, order>, 2> val_table;
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
    const auto &point = points[qp];

    sin_val[0][qp][0] = point[0];
    for (int i = 1; i < dim; ++i)
    {
      sin_val[0][qp][i]   = sin(point[i]);
      cos_val[0][qp][i-1] = cos(point[i]);
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
BallGridFunction<dim>::
evaluate_0(const ValueVector<GridPoint> &points,
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
BallGridFunction<dim>::
evaluate_1(const ValueVector<GridPoint> &points,
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
BallGridFunction<dim>::
evaluate_2(const ValueVector<GridPoint> &points,
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

template<int dim>
void
BallGridFunction<dim>::
print_info(LogStream &out) const
{
  out.begin_item("BallGridFunction<" + std::to_string(dim) +">");
  out.end_item();

  out << "Name: " << this->name_ << std::endl;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

template<int dim>
SphereGridFunction<dim>::
SphereGridFunction(const SharedPtrConstnessHandler<GridType> &grid)
  :
  parent_t::FormulaGridFunction(grid)
{}



template<int dim>
auto
SphereGridFunction<dim>::
create(std::shared_ptr<GridType> grid) -> std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid)));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;

}

template<int dim>
auto
SphereGridFunction<dim>::
const_create(std::shared_ptr<const GridType> grid) -> std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<GridType>(grid)));
}

#ifdef MESH_REFINEMENT
template<int dim>
void
SphereGridFunction<dim>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
  const Grid<dim> &old_grid)
{
  this->grid_function_previous_refinement_ =
    self_t::const_create(
      this->get_grid()->get_grid_pre_refinement());
}
#endif // MESH_REFINEMENT



template<int dim>
template<int order>
auto
SphereGridFunction<dim>::get_aux_vals(const ValueVector<GridPoint> &points) const
{
  SafeSTLArray<SafeSTLArray<SafeSTLVector<SafeSTLArray<double, range> >, order>, 2> val_table;
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
    sin_val[0][qp][0] = R;
    for (int i = 1; i < range; ++i)
    {
      sin_val[0][qp][i]   = sin(points[qp][i-1]);
      cos_val[0][qp][i-1] = cos(points[qp][i-1]);
    }
    cos_val[0][qp][range - 1] = 1;

    for (int der = 1; der < order; ++der)
    {
      auto res = std::div(der,2);
      sin_val[der][qp][0] = R;
      for (int i = 1; i < range; ++i)
      {
        sin_val[der][qp][i] =
          std::pow(-1, res.quot) *
          (res.rem == 0? sin_val[0][qp][i]: cos_val[0][qp][i-1]);
        cos_val[der][qp][i-1] = -sin_val[der-1][qp][i];
      }
      cos_val[der][qp][range - 1] = 1.;
    }
  }
  return val_table;
}



template<int dim>
auto
SphereGridFunction<dim>::
evaluate_0(const ValueVector<GridPoint> &points,
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
    for (int i = 0; i < range; ++i)
    {
      y *= s[qp][i];
      x[i] = y * c[qp][i];
    }
  }
}



template<int dim>
auto
SphereGridFunction<dim>::
evaluate_1(const ValueVector<GridPoint> &points,
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

    for (int i = 0; i < range-1; ++i)
    {
      for (int j = 1; j < i+2; ++j)
      {
        double djy = 1.;
        for (int k = 1; k < i+1; ++k)
          djy *= k!=j ? s[qp][k] : s_p[qp][k];
        grad[j-1][i] = djy * (i+1!=j ? c[qp][i] : c_p[qp][i]);
      }
    }

    const int i = range-1;
    for (int j = 1; j < dim+1 ; ++j)
    {
      double djy = 1.;
      for (int k = 1; k < i+1; ++k)
        djy *= k!=j ? s[qp][k] : s_p[qp][k];
      grad[j-1][i] = djy;
    }
  }
}



template<int dim>
auto
SphereGridFunction<dim>::
evaluate_2(const ValueVector<GridPoint> &points,
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
    for (int i = 0; i < range-1; ++i)
    {
      for (int j = 1; j < i+2; ++j)
      {
        for (int k = 1; k < j+1; ++k)
        {
          double d2jy = 1.;
          for (int l = 1; l < i+1; ++l)
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

    const int i = range-1;
    for (int j = 1; j < dim+1; ++j)
      for (int k = 1; k < j+1; ++k)
      {
        double d2jy = 1.;
        for (int l = 1; l < i+1; ++l)
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


    for (int i = 0; i < range; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < j; ++k)
        {
          hessian[k][j][i] = hessian[j][k][i];
        }
  }
}

template<int dim>
void
SphereGridFunction<dim>::
print_info(LogStream &out) const
{
  out.begin_item("SphereGridFunction<" + std::to_string(dim) +">");
  out.end_item();

  out << "Name: " << this->name_ << std::endl;
}

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------


CylindricalAnnulusGridFunction::
CylindricalAnnulusGridFunction(
  const SharedPtrConstnessHandler<GridType> &grid,
  const Real r0,
  const Real r1,
  const Real h0,
  const Real h1,
  const Real theta0,
  const Real theta1)
  :
  parent_t::FormulaGridFunction(grid),
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



auto
CylindricalAnnulusGridFunction::
create(std::shared_ptr<GridType> grid,
       const Real r0,
       const Real r1,
       const Real h0,
       const Real h1,
       const Real theta0,
       const Real theta1) ->  std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid),r0,r1,h0,h1,theta0,theta1));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;
}


auto
CylindricalAnnulusGridFunction::
const_create(std::shared_ptr<const GridType> grid,
             const Real r0,
             const Real r1,
             const Real h0,
             const Real h1,
             const Real theta0,
             const Real theta1) -> std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(new self_t(
    SharedPtrConstnessHandler<GridType>(grid),r0,r1,h0,h1,theta0,theta1));
}

#ifdef MESH_REFINEMENT
void
CylindricalAnnulusGridFunction::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,3> &knots_to_insert,
  const Grid<3> &old_grid)
{
  this->grid_function_previous_refinement_ =
    self_t::const_create(
      this->get_grid()->get_grid_pre_refinement(),r0_,r1_,h0_,h1_,theta0_,theta1_);
}
#endif // MESH_REFINEMENT


auto
CylindricalAnnulusGridFunction::
evaluate_0(const ValueVector<GridPoint> &points,
           ValueVector<Value> &values) const -> void
{

  const int n_points = points.size();

  for (int qp = 0; qp < n_points; ++qp)
  {
    auto &F = values[qp];
    const auto &pt = points[qp];

    const Real theta = pt[0];
    const Real r     = pt[1];
    const Real z     = pt[2];

    F[0] = (dR_ * r + r0_) * cos(dT_ * theta);
    F[1] = (dR_ * r + r0_) * sin(dT_ * theta);
    F[2] = h0_ + z * dH_;
  }
}




auto
CylindricalAnnulusGridFunction::
evaluate_1(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<1>> &values) const -> void
{

  const int n_points = points.size();

  for (int qp = 0; qp < n_points; ++qp)
  {
    auto &dF = values[qp];
    const auto &pt = points[qp];

    const Real theta = pt[0];
    const Real r     = pt[1];

    const auto s_dt_theta = sin(dT_ * theta);
    const auto c_dt_theta = cos(dT_ * theta);

    dF[0][0] = - dT_ * (dR_ * r + r0_) * s_dt_theta;
    dF[0][1] =   dT_ * (dR_ * r + r0_) * c_dt_theta;
    dF[0][2] = 0.0;

    dF[1][0] = dR_ * c_dt_theta;
    dF[1][1] = dR_ * s_dt_theta;
    dF[1][2] = 0.0;

    dF[2][0] = 0.0;
    dF[2][1] = 0.0;
    dF[2][2] = dH_;
  }
}




auto
CylindricalAnnulusGridFunction::
evaluate_2(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<2>> &values) const -> void
{
  const int n_points = points.size();

  for (int qp = 0; qp < n_points; ++qp)
  {
    auto &d2F = values[qp];
    const auto &pt = points[qp];

    const Real theta = pt[0];
    const Real r     = pt[1];

    const auto s_dt_theta = sin(dT_ * theta);
    const auto c_dt_theta = cos(dT_ * theta);

    d2F[0][0][0] = - dT_ * dT_ * (dR_ * r + r0_) * c_dt_theta;
    d2F[0][0][1] = - dT_ * dT_ * (dR_ * r + r0_) * s_dt_theta;
    d2F[0][0][2] = 0.0;

    d2F[1][0][0] = -dT_ * dR_ * s_dt_theta;
    d2F[1][0][1] =  dT_ * dR_ * c_dt_theta;
    d2F[1][0][2] = 0.0;

    d2F[2][0][0] = 0.0;
    d2F[2][0][1] = 0.0;
    d2F[2][0][2] = 0.0;


    d2F[0][1][0] = - dT_ * dR_ * s_dt_theta;
    d2F[0][1][1] =   dT_ * dR_ * c_dt_theta;
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


void
CylindricalAnnulusGridFunction::
print_info(LogStream &out) const
{
  using std::endl;
  out.begin_item("CylindricalAnnulusGridFunction");
  out << "r0 = " << r0_ << endl;
  out << "r1 = " << r1_ << endl;
  out << "h0 = " << h0_ << endl;
  out << "h1 = " << h1_ << endl;
  out << "theta0 = " << theta0_ << endl;
  out << "theta1 = " << theta1_ << endl;
  out.end_item();

  out << "Name: " << this->name_ << std::endl;
}

//------------------------------------------------------------------------------
} // of namespace functions.

IGA_NAMESPACE_CLOSE

#include <igatools/functions/grid_function_lib.inst>

