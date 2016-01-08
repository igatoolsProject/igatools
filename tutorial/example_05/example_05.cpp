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

#include <igatools/basis_functions/bspline.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/nurbs.h>
#include <igatools/geometry/grid_function_lib.h>
#include <igatools/functions/ig_function.h>
#include <igatools/io/writer.h>

#include <igatools/geometry/grid_element.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/base/quadrature_lib.h>

#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/linear_algebra/epetra_matrix.h>
#include <igatools/linear_algebra/epetra_vector.h>
#include <igatools/linear_algebra/epetra_solver.h>

// [include]
#include <igatools/functions/formula_function.h>
#include <igatools/functions/function_lib.h>
#include <igatools/basis_functions/space_tools.h>
// [include]

using namespace iga;
using namespace std;
using namespace EpetraTools;

LogStream out;

// [custom_declaration]
template<int dim, int codim=0, int range=1, int rank=1>
class CustomFunction : public FormulaFunction<dim,codim,range,rank>
{
private:
  using self_t   = CustomFunction<dim,codim,range,rank>;
  using typename Function<dim,codim,range,rank>::DomainType;
  using typename FormulaFunction<dim,codim,range,rank>::Value;
  using typename FormulaFunction<dim,codim,range,rank>::Point;
  template <int order>
  using Derivative = typename FormulaFunction<dim,codim,range,rank>::template Derivative<order>;
  Value(*funct_D0)(const Point);
// [custom_declaration]

// [custom_constructor]
  CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain)
    : FormulaFunction<dim,codim,range,rank>(domain,"") {};

  CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain,
                 Value(*f_D0)(const Point))
    : FormulaFunction<dim,codim,range,rank>(domain,""), funct_D0(f_D0) {};
// [custom_constructor]

// [custom_evaluator]
  void evaluate_0(const ValueVector<Point> &points, ValueVector<Value> &values) const
  {
    auto point = points.begin();
    for (auto &val : values)
    {
      val = funct_D0(*point);
      ++point;
    }
  };
// [custom_evaluator]

// [custom_not_implemented]
  void evaluate_1(const ValueVector<Point> &points, ValueVector<Derivative<1>> &values) const
  {
    std::cout << "evaluation of first derivatives not implemented!" << std::endl;
  };
  void evaluate_2(const ValueVector<Point> &points, ValueVector<Derivative<2>> &values) const
  {
    std::cout << "evaluation of second derivatives not implemented!" << std::endl;
  };
// [custom_not_implemented]

// [custom_public]
public:
  static std::shared_ptr<const self_t> const_create(const std::shared_ptr<const DomainType> &domain,
                                                    Value(*f_D0)(const Point))
  {
    return std::shared_ptr<const self_t>(new self_t(SharedPtrConstnessHandler<DomainType>(domain),f_D0));
  };
// [custom_public]
  void print_info(LogStream &out) const
  {
    std::cout << "C makes it easy to shoot yourself in the foot." << std::endl;
    std::cout << "C++ makes it harder, but when you do, it blows away your whole leg." << std::endl;
  };
};

shared_ptr<const Domain<2>> quarter_annulus(const Size nel)
{
  using numbers::PI;
  BBox<2> box;
  box[0] = {{1.0,2.0}};
  box[1] = {{0.0,PI/2}};
  auto grid = Grid<2>::const_create(box,nel+1);
  auto geom_funct = grid_functions::BallGridFunction<2>::const_create(grid);
  return Domain<2>::const_create(geom_funct);
}

// [poisson_declaration]
template<int dim>
class PoissonProblem
{
private:
  shared_ptr<const SplineSpace<dim>>        ref_space;
  shared_ptr<const BSpline<dim>>            ref_basis;
  shared_ptr<const PhysicalSpaceBasis<dim>> phy_basis;
  shared_ptr<const QGauss<dim>>             quad;
  shared_ptr<const QGauss<dim-1>>           face_quad;
  shared_ptr<const Function<dim>>           source_term;
  std::map<Index,shared_ptr<const Function<dim-1,1,1>>> dirichlet_cond;
  shared_ptr<Matrix> mat;
  shared_ptr<Vector> rhs;
  shared_ptr<Vector> sol;
// [poisson_declaration]

// [poisson_constructor]
public:
  PoissonProblem(const shared_ptr<const Domain<dim>> domain, const Index deg,
                 const shared_ptr<const Function<dim>> source,
                 const std::map<Index,shared_ptr<const Function<dim-1,1,1>>> dirichlet)
  {
    source_term    = source;
    dirichlet_cond = dirichlet;
    auto grid = domain->get_grid_function()->get_grid();
    ref_space = SplineSpace<dim>::const_create(deg,grid);
    ref_basis = BSpline<dim>::const_create(ref_space);
    phy_basis = PhysicalSpaceBasis<dim>::const_create(ref_basis,domain);
    quad      = QGauss<dim>::const_create(deg+1);
    face_quad = QGauss<dim-1>::const_create(deg+1);
    mat = create_matrix(*phy_basis,DofProperties::active,Epetra_SerialComm());
    rhs = create_vector(mat->RangeMap());
    sol = create_vector(mat->DomainMap());
  };
// [poisson_constructor]

// [poisson_methods]
  void assemble();
  void solve();
  void save();
  Real error(shared_ptr<const Function<dim>> exact_solution);
};
// [poisson_methods]

template<int dim>
void PoissonProblem<dim>::assemble()
{

  auto basis_el      = phy_basis->begin();
  auto basis_el_end  = phy_basis->end();
  auto basis_handler = phy_basis->create_cache_handler();
  auto flag = space_element::Flags::value |
              space_element::Flags::gradient |
              space_element::Flags::w_measure;
  basis_handler->set_element_flags(flag);
  basis_handler->init_element_cache(basis_el,quad);

  auto funct_el      = source_term->begin();
  auto funct_handler = source_term->create_cache_handler();
  funct_handler->set_element_flags(function_element::Flags::D0);
  funct_handler->init_cache(funct_el,quad);

// [poisson_assemble]
  for (; basis_el!=basis_el_end; ++basis_el, ++funct_el)
  {
    basis_handler->fill_element_cache(basis_el);
    funct_handler->fill_element_cache(funct_el);

    DenseMatrix loc_mat = basis_el->integrate_gradu_gradv();
    auto funct  = funct_el->get_element_values_D0();
    DenseVector loc_rhs = basis_el->integrate_u_func(funct);

    const auto loc_dofs = basis_el->get_local_to_global();
    mat->add_block(loc_dofs, loc_dofs,loc_mat);
    rhs->add_block(loc_dofs, loc_rhs);
  }
  mat->FillComplete();
// [poisson_assemble]

// [poisson_dirichlet]
  std::map<Index,Real> dirichlet_vals;
  space_tools::project_boundary_values(dirichlet_cond,*phy_basis,face_quad,dirichlet_vals);
  dof_tools::apply_boundary_values(dirichlet_vals,*mat,*rhs,*sol);
}
// [poisson_dirichlet]

template<int dim>
void PoissonProblem<dim>::solve()
{
  auto solver = create_solver(*mat,*sol,*rhs);
  solver->solve();
}

template<int dim>
void PoissonProblem<dim>::save()
{
  auto domain = phy_basis->get_physical_domain();
  Writer<dim> writer(domain,5);
  auto solution = IgFunction<dim,0,1,1>::const_create(phy_basis, *sol);
  writer.add_field(*solution, "solution");
  string filename = "problem_" + to_string(dim) + "d" ;
  writer.save(filename);
}

template<int dim>
Real PoissonProblem<dim>::error(shared_ptr<const Function<dim>> exact_solution)
{

  auto domain         = phy_basis->get_physical_domain();
  auto domain_el      = domain->begin();
  auto domain_el_end  = domain->end();
  auto domain_handler = domain->create_cache_handler();
  domain_handler->set_element_flags(domain_element::Flags::w_measure);
  domain_handler->init_element_cache(domain_el,quad);

  auto discr_solution    = IgFunction<dim,0,1,1>::const_create(phy_basis, *sol);
  auto discr_sol_el      = discr_solution->begin();
  auto discr_sol_handler = discr_solution->create_cache_handler();
  discr_sol_handler->set_element_flags(function_element::Flags::D0);
  discr_sol_handler->init_cache(discr_sol_el,quad);

  auto exact_sol_el      = exact_solution->begin();
  auto exact_sol_handler = exact_solution->create_cache_handler();
  exact_sol_handler->set_element_flags(function_element::Flags::D0);
  exact_sol_handler->init_cache(exact_sol_el,quad);

  Real error = 0.0;
  for (; domain_el!=domain_el_end; ++domain_el, ++discr_sol_el, ++exact_sol_el)
  {
    domain_handler->fill_element_cache(domain_el);
    discr_sol_handler->fill_element_cache(discr_sol_el);
    exact_sol_handler->fill_element_cache(exact_sol_el);

    auto w_meas    = domain_el->get_element_w_measures();
    auto discr_val = discr_sol_el->get_element_values_D0();
    auto exact_val = exact_sol_el->get_element_values_D0();
    auto num_quad  = w_meas.get_num_points();

    for (int q=0; q<num_quad; q++)
    {
      error += std::pow(discr_val[q][0]-exact_val[q][0],2) * w_meas[q];
    }
  }
  return std::sqrt(error);
}

// [functions_u]]
#define SQR(x,y) (((x)*(x))+((y)*(y)))
using numbers::PI;
// exact solution
Values<2,1,1> u(Points<2> x)
{
  Values<2,1,1> y;
  y  = 0.1*sin(2.0*PI*SQR(x[0],x[1])) + std::sqrt(SQR(x[0],x[1]));
  return y;
}
// [functions_u]

// [functions_f]
Values<2,1,1> f(Points<2> x)
{
  Values<2,1,1> y;
  auto y1 = 1.6*PI*PI * SQR(x[0],x[1]) * sin(2.0*PI*SQR(x[0],x[1]));
  auto y2 = 0.8*PI * cos(2.0*PI*SQR(x[0],x[1]));
  auto y3 = std::pow(SQR(x[0],x[1]),-0.5);
  y = y1-y2-y3;
  return y;
}
// [functions_f]

// [functions_g2]
Values<1,1,1> g2(Points<2> x)
{
  Values<1,1,1> y;
  y  = 0.1*sin(2.0*PI*SQR(x[0],x[1]))+x[0];
  return y;
}
// [functions_g2]

// [functions_g3]
Values<1,1,1> g3(Points<2> x)
{
  Values<1,1,1> y;
  y  = 0.1*sin(2.0*PI*SQR(x[0],x[1]))+x[1];
  return y;
}
// [functions_g3]

// [main]
int main()
{
  const int nel = 16;
  const int deg = 2;

  auto annulus = quarter_annulus(nel);
  auto source  = CustomFunction<2>::const_create(annulus,f);
  auto exact_solution = CustomFunction<2>::const_create(annulus,u);
// [main]

// [main_pre_dirichlet]
  using SubGridElemMap = typename Grid<2>::template SubGridMap<1>;
  SubGridElemMap sub_grid_elem_map;
  auto grid = annulus->get_grid_function()->get_grid();
  std::map<Index,shared_ptr<const Function<1,1,1>>> dirichlet;
// [main_pre_dirichlet]

// [main_dirichlet_cond]
  const auto sub_grid_0    = grid->template get_sub_grid<1>(0,sub_grid_elem_map);
  const auto sub_annulus_0 = annulus->template get_sub_domain<1>(0,sub_grid_elem_map,sub_grid_0);
  auto g_0 = functions::ConstantFunction<1,1,1>::const_create(sub_annulus_0, {1.0});
  dirichlet[0] = dynamic_pointer_cast<const Function<1,1,1>>(g_0);
// [main_dirichlet_cond]

  const auto sub_grid_1    = grid->template get_sub_grid<1>(1,sub_grid_elem_map);
  const auto sub_annulus_1 = annulus->template get_sub_domain<1>(1,sub_grid_elem_map,sub_grid_1);
  auto g_1 = functions::ConstantFunction<1,1,1>::const_create(sub_annulus_1, {2.0});
  dirichlet[1] = dynamic_pointer_cast<const Function<1,1,1>>(g_1);

  const auto sub_grid_2    = grid->template get_sub_grid<1>(2,sub_grid_elem_map);
  const auto sub_annulus_2 = annulus->template get_sub_domain<1>(2,sub_grid_elem_map,sub_grid_2);
  auto g_2 = CustomFunction<1,1,1>::const_create(sub_annulus_2,g2);
// [main_dirichlet_cond_custom]
  dirichlet[2] = dynamic_pointer_cast<const Function<1,1,1>>(g_2);
// [main_dirichlet_cond_custom]

  const auto sub_grid_3    = grid->template get_sub_grid<1>(3,sub_grid_elem_map);
  const auto sub_annulus_3 = annulus->template get_sub_domain<1>(3,sub_grid_elem_map,sub_grid_3);
  auto g_3 = CustomFunction<1,1,1>::const_create(sub_annulus_3,g3);
  dirichlet[3] = dynamic_pointer_cast<const Function<1,1,1>>(g_3);

  PoissonProblem<2> problem_2d(annulus,deg,source,dirichlet);
  problem_2d.assemble();
  problem_2d.solve();
  out << "The L2-error of the Poisson problem is: ";
  out << problem_2d.error(exact_solution) << endl;
  problem_2d.save();

  return 0;
}