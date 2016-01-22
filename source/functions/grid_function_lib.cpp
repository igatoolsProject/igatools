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
ConstantGridFunction(const SharedPtrConstnessHandler<GridType> &grid,
                     const Value &b)
  :
  parent_t(grid),
  b_(b)
{}



template<int dim, int range>
auto
ConstantGridFunction<dim,range>::
create(const std::shared_ptr<GridType> &grid,
       const Value &b) ->  std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid),b));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;

}

template<int dim, int range>
auto
ConstantGridFunction<dim,range>::
const_create(const std::shared_ptr<const GridType> &grid,
             const Value &b) ->  std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid), b));
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
LinearGridFunction(const SharedPtrConstnessHandler<GridType> &grid,
                   const Derivative<1> &A,
                   const Value &b)
  :
  parent_t(grid),
  A_(A),
  b_(b)
{}



template<int dim, int range>
auto
LinearGridFunction<dim,range>::
create(const std::shared_ptr<GridType> &grid,
       const Derivative<1> &A,
       const Value &b) ->  std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid),A,b));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;

}

template<int dim, int range>
auto
LinearGridFunction<dim,range>::
const_create(const std::shared_ptr<const GridType> &grid,
             const Derivative<1> &A,
             const Value &b) ->  std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid), A, b));
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
get_A() const -> const Derivative<1> &
{
  return A_;
}



template<int dim, int range>
auto
LinearGridFunction<dim,range>::
get_b() const -> const Value &
{
  return b_;
}


//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
template<int dim>
IdentityGridFunction<dim>::
IdentityGridFunction(const SharedPtrConstnessHandler<GridType> &grid)
  :
  parent_t(grid)
{}



template<int dim>
auto
IdentityGridFunction<dim>::
create(const std::shared_ptr<GridType> &grid)
->  std::shared_ptr<self_t>
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
IdentityGridFunction<dim>::
const_create(const std::shared_ptr<const GridType> &grid)
->  std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid)));
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
  const SharedPtrConstnessHandler<GridType> &grid)
  :
  parent_t::FormulaGridFunction(grid)
{}



auto
CylindricalAnnulusGridFunction::
create(const std::shared_ptr<GridType> &grid) ->  std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid)));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;
}


auto
CylindricalAnnulusGridFunction::
const_create(const std::shared_ptr<const GridType> &grid) -> std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(new self_t(
    SharedPtrConstnessHandler<GridType>(grid)));
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
      this->get_grid()->get_grid_pre_refinement());
}
#endif // MESH_REFINEMENT


void
CylindricalAnnulusGridFunction::
evaluate_0(const ValueVector<GridPoint> &points,
           ValueVector<Value> &values) const
{
  const int n_points = points.size();

  // x[0] = r
  // x[1] = theta
  // x[2] = zeta

  for (int qp = 0; qp < n_points; ++qp)
  {
    auto &F = values[qp];
    const auto &x = points[qp];

    F[0] = x[0] * cos(x[1]);
    F[1] = x[0] * sin(x[1]);
    F[2] = x[2];
  }
}




void
CylindricalAnnulusGridFunction::
evaluate_1(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<1>> &values) const
{
  const int n_points = points.size();

  // x[0] = r
  // x[1] = theta
  // x[2] = zeta

  for (int qp = 0; qp < n_points; ++qp)
  {
    auto &dF = values[qp];
    const auto &x = points[qp];

    const Real s = sin(x[1]);
    const Real c = cos(x[1]);

    dF[0][0] = c;
    dF[0][1] = s;
    dF[0][2] = 0.0;

    dF[1][0] = - x[0] * s;
    dF[1][1] =   x[0] * c;
    dF[1][2] = 0.0;

    dF[2][0] = 0.0;
    dF[2][1] = 0.0;
    dF[2][2] = 1.0;
  }
}




void
CylindricalAnnulusGridFunction::
evaluate_2(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<2>> &values) const
{
  const int n_points = points.size();


  // x[0] = r
  // x[1] = theta
  // x[2] = zeta
  SafeSTLArray<Real,3> x;
  SafeSTLArray<Real,3> dx;

  for (int qp = 0; qp < n_points; ++qp)
  {
    auto &d2F = values[qp];
    const auto &x = points[qp];

    const Real s = sin(x[1]);
    const Real c = cos(x[1]);

    d2F[0][0][0] = 0.0;
    d2F[0][0][1] = 0.0;
    d2F[0][0][2] = 0.0;

    d2F[0][1][0] = - s;
    d2F[0][1][1] =   c;
    d2F[0][1][2] = 0.0;

    d2F[0][2][0] = 0.0;
    d2F[0][2][1] = 0.0;
    d2F[0][2][2] = 0.0;

    d2F[1][0][0] = - s;
    d2F[1][0][1] =   c;
    d2F[1][0][2] = 0.0;

    d2F[1][1][0] = - x[0] * c;
    d2F[1][1][1] = - x[0] * s;
    d2F[1][1][2] = 0.0;

    d2F[1][2][0] = 0.0;
    d2F[1][2][1] = 0.0;
    d2F[1][2][2] = 0.0;

    d2F[2][0][0] = 0.0;
    d2F[2][0][1] = 0.0;
    d2F[2][0][2] = 0.0;

    d2F[2][1][0] = 0.0;
    d2F[2][1][1] = 0.0;
    d2F[2][1][2] = 0.0;

    d2F[2][2][0] = 0.0;
    d2F[2][2][1] = 0.0;
    d2F[2][2][2] = 0.0;
  }
}

auto
CylindricalAnnulusGridFunction::
evaluate_preimage(const ValueVector<Value> &phys_points) const
-> ValueVector<GridPoint>
{
  const int n_pts = phys_points.get_num_points();
  Assert(n_pts > 0,ExcEmptyObject());

  ValueVector<GridPoint> param_points(n_pts);

  const auto bounding_box = this->get_grid()->get_bounding_box();

  for (int pt = 0 ; pt < n_pts ; ++pt)
  {
    const auto &phys_pt = phys_points[pt];
    auto &param_pt = param_points[pt];

    const auto &x = phys_pt[0];
    const auto &y = phys_pt[1];
    const auto &z = phys_pt[2];

    param_pt[0] = sqrt(x*x + y*y);
    param_pt[1] = std::atan2(y,x);
    param_pt[2] = z;


    AssertThrow(bounding_box.is_point_inside(param_pt),
    ExcMessage("The pre-image of the point " + std::to_string(pt) + " is not in the parametric domain."));
  }
  return param_points;
}

void
CylindricalAnnulusGridFunction::
print_info(LogStream &out) const
{
  const auto box = this->get_grid()->get_bounding_box();

  using std::endl;
  out.begin_item("CylindricalAnnulusGridFunction");
  out << "r0 = " << box[0][0] << endl;
  out << "r1 = " << box[0][1] << endl;
  out << "h0 = " << box[2][0] << endl;
  out << "h1 = " << box[2][1] << endl;
  out << "theta0 = " << box[2][0] << endl;
  out << "theta1 = " << box[2][1] << endl;
  out.end_item();

  out << "Name: " << this->name_ << std::endl;
}

//------------------------------------------------------------------------------





//------------------------------------------------------------------------------
template <int range>
TriangleGridFunction<range>::
TriangleGridFunction(
  const SharedPtrConstnessHandler<GridType> &grid,
  const Points<range> &vertex_0,
  const Points<range> &vertex_1,
  const Points<range> &vertex_2)
  :
  parent_t(grid),
  vertices_{vertex_0,vertex_1,vertex_2}
{}


template <int range>
auto
TriangleGridFunction<range>::
create(const std::shared_ptr<GridType> &grid,
       const Points<range> &vertex_0,
       const Points<range> &vertex_1,
       const Points<range> &vertex_2) ->  std::shared_ptr<self_t>
{
  auto func = std::shared_ptr<self_t>(
    new self_t(SharedPtrConstnessHandler<GridType>(grid),vertex_0,vertex_1,vertex_2));

#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif

  return func;
}


template <int range>
auto
TriangleGridFunction<range>::
const_create(const std::shared_ptr<const GridType> &grid,
             const Points<range> &vertex_0,
             const Points<range> &vertex_1,
             const Points<range> &vertex_2) ->  std::shared_ptr<const self_t>
{
  return std::shared_ptr<self_t>(new self_t(
    SharedPtrConstnessHandler<GridType>(grid),vertex_0,vertex_1,vertex_2));
}


#ifdef MESH_REFINEMENT
template <int range>
void
TriangleGridFunction<range>::
rebuild_after_insert_knots(
  const SafeSTLArray<SafeSTLVector<Real>,2> &knots_to_insert,
  const Grid<2> &old_grid)
{
  this->grid_function_previous_refinement_ =
    self_t::const_create(
      this->get_grid()->get_grid_pre_refinement(),vertices_[0],vertices_[1],vertices_[2]);
}
#endif // MESH_REFINEMENT



template <int range>
void
TriangleGridFunction<range>::
evaluate_0(const ValueVector<GridPoint> &points,
           ValueVector<Value> &values) const
{
  const auto &knots = this->get_grid()->get_knots();

  const auto h0 = knots[0]->back() - knots[0]->front();
  const auto h1 = knots[1]->back() - knots[1]->front();

  const int n_points = points.size();

  Points<range> A = vertices_[0];
  Points<range> B = (vertices_[1] - vertices_[0]) / h0;
  Points<range> C = (vertices_[2] - vertices_[0]) / h1;
  Points<range> D = (vertices_[2] - A - h0 * B - h1 * C) / (h0 * h1);

  for (int qp = 0; qp < n_points; ++qp)
  {
    auto &f = values[qp];
    const auto &pt = points[qp];

    const Real &u = pt[0];
    const Real &v = pt[1];

    f = A + u * B + v * C + u * v * D;
  }
}


template <int range>
void
TriangleGridFunction<range>::
evaluate_1(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<1>> &values) const
{
  const auto &knots = this->get_grid()->get_knots();

  const auto h0 = knots[0]->back() - knots[0]->front();
  const auto h1 = knots[1]->back() - knots[1]->front();

  const int n_points = points.size();

  Points<range> A = vertices_[0];
  Points<range> B = (vertices_[1] - vertices_[0]) / h0;
  Points<range> C = (vertices_[2] - vertices_[0]) / h1;
  Points<range> D = (vertices_[2] - A - h0 * B - h1 * C) / (h0 * h1);

  for (int qp = 0; qp < n_points; ++qp)
  {
    auto &df = values[qp];
    const auto &pt = points[qp];

    const Real &u = pt[0];
    const Real &v = pt[1];

    df[0] = B + v * D;
    df[1] = C + u * D;
  }
}

template <int range>
void
TriangleGridFunction<range>::
evaluate_2(const ValueVector<GridPoint> &points,
           ValueVector<Derivative<2>> &values) const
{
  const auto &knots = this->get_grid()->get_knots();

  const auto h0 = knots[0]->back() - knots[0]->front();
  const auto h1 = knots[1]->back() - knots[1]->front();


  Points<range> A = vertices_[0];
  Points<range> B = (vertices_[1] - vertices_[0]) / h0;
  Points<range> C = (vertices_[2] - vertices_[0]) / h1;
  Points<range> D = (vertices_[2] - A - h0 * B - h1 * C) / (h0 * h1);

  for (auto &d2f : values)
  {
    d2f[0][0] = 0.0;
    d2f[1][0] = D;
    d2f[0][1] = D;
    d2f[1][1] = 0.0;
  }
}

template <int range>
void
TriangleGridFunction<range>::
print_info(LogStream &out) const
{
  using std::endl;
  out.begin_item("TriangleGridFunction");
  out << "Vertices: " << vertices_ << endl;
  out.end_item();

  out << "Name: " << this->name_ << std::endl;
}



//------------------------------------------------------------------------------


} // of namespace functions.

IGA_NAMESPACE_CLOSE

#include <igatools/functions/grid_function_lib.inst>

