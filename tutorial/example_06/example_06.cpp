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
#include <igatools/functions/grid_function_lib.h>
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

#include <igatools/functions/formula_function.h>
#include <igatools/functions/function_lib.h>
#include <igatools/basis_functions/space_tools.h>

using namespace iga;
using namespace std;
using namespace EpetraTools;

LogStream out;

template<int dim, int codim=0, int range=1, int rank=1>
class CustomFunction : public FormulaFunction<dim,codim,range,rank> {
private:
  using self_t   = CustomFunction<dim,codim,range,rank>;
  using typename Function<dim,codim,range,rank>::DomainType;
  using typename FormulaFunction<dim,codim,range,rank>::Value;
  using typename FormulaFunction<dim,codim,range,rank>::Point;
  template <int order>
  using Derivative = typename FormulaFunction<dim,codim,range,rank>::template Derivative<order>;
  Value (*funct_D0)(const Point);
  CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain)
    : FormulaFunction<dim,codim,range,rank>(domain) {};

  CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain,
                 Value (*f_D0)(const Point))
    : FormulaFunction<dim,codim,range,rank>(domain), funct_D0(f_D0) {};
  // evaluators
  void evaluate_0(const ValueVector<Point> &points, ValueVector<Value> &values) const {
    auto point = points.begin();
    for (auto &val : values ) {
      val = funct_D0(*point);
      ++point;
    }
  };
  void evaluate_1(const ValueVector<Point> &points, ValueVector<Derivative<1>> &values) const {
    std::cout << "evaluation of first derivatives not implemented!" << std::endl;
  };
  void evaluate_2(const ValueVector<Point> &points, ValueVector<Derivative<2>> &values) const {
    std::cout << "evaluation of second derivatives not implemented!" << std::endl;
  };  
public:           
  static std::shared_ptr<const self_t> const_create(const std::shared_ptr<const DomainType> &domain,
                                                    Value (*f_D0)(const Point)) {
    return std::shared_ptr<const self_t>(new self_t(SharedPtrConstnessHandler<DomainType>(domain),f_D0));
  };
  void print_info(LogStream &out) const {
    std::cout << "Sometimes it pays to stay in bed on Monday, rather than" <<;
    std::cout << "spending the rest of the week debugging Monday's code." << std::endl;
  };
};

template<int dim>
class PoissonProblem {
  private:
    shared_ptr<const SplineSpace<dim>>        ref_space;
    shared_ptr<const BSpline<dim>>            ref_basis;
    shared_ptr<const PhysicalSpaceBasis<dim>> phy_basis;
    shared_ptr<const QGauss<dim>>   quad;
    shared_ptr<const QGauss<dim-1>> face_quad;
    shared_ptr<const Function<dim>> source_term;
    std::map<Index,shared_ptr<const Function<dim-1,1,1>>> dirichlet_cond;
// [neumann_data]
    std::map<Index,shared_ptr<const Function<dim-1,1,1>>> neumann_cond;
// [neumann_data]
    shared_ptr<Matrix> mat;
    shared_ptr<Vector> rhs;
    shared_ptr<Vector> sol;
// [constructor]
  public:
    PoissonProblem(const shared_ptr<const Domain<dim>> domain, const Index deg,
                   const shared_ptr<const Function<dim>> source,
                   const std::map<Index,shared_ptr<const Function<dim-1,1,1>>> dirichlet,
                   const std::map<Index,shared_ptr<const Function<dim-1,1,1>>> neumann) {
      source_term    = source;
      dirichlet_cond = dirichlet;
      neumann_cond   = neumann;
// [constructor]
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
    void assemble();
    void solve();
    void save();
    Real error(const shared_ptr<const Function<dim>> exact_solution);
};

template<int dim>
void PoissonProblem<dim>::assemble() {

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

  for (; basis_el!=basis_el_end; ++basis_el, ++funct_el) {
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

// [neu_loop_init]
  for (auto it=neumann_cond.begin(); it!=neumann_cond.end(); ++it) {

    auto face = it->first;
    auto funct = it->second;
// [neu_loop_init]

// [sub_grid_init]
    using SubGridElemMap = typename Grid<dim>::template SubGridMap<dim-1>;
    SubGridElemMap sub_grid_elem_map;
    auto grid     = phy_basis->get_domain()->get_grid_function()->get_grid();
    auto sub_grid = grid->template get_sub_grid<dim-1>(face,sub_grid_elem_map);
// [sub_grid_init]

// [standard_caches]
    auto sub_domain    = funct->get_domain();
    auto s_dom_el      = sub_domain->begin();
    auto s_dom_el_end  = sub_domain->end();
    auto s_dom_handler = sub_domain->create_cache_handler();
    s_dom_handler->set_element_flags(domain_element::Flags::w_measure);
    s_dom_handler->init_element_cache(s_dom_el,face_quad);

    auto funct_el      = funct->begin();
    auto funct_handler = funct->create_cache_handler();
    funct_handler->set_element_flags(function_element::Flags::D0);
    funct_handler->init_cache(funct_el,face_quad);
// [standard_caches]

// [basis_handler_init]
    auto basis_el = phy_basis->begin();
    basis_handler->template set_flags<dim-1>(space_element::Flags::value);
    basis_handler->init_face_cache(basis_el,face_quad);
// [basis_handler_init]

// [neu_el_loop_init]
    for (; s_dom_el!=s_dom_el_end; ++s_dom_el, ++funct_el) {
      s_dom_handler->fill_element_cache(s_dom_el);
      funct_handler->fill_element_cache(funct_el);
// [neu_el_loop_init]


// [get_func_elem_id]
      auto funct_elem_id = funct_el->get_index();
// [get_func_elem_id]
// [get_basis_el_id]
      auto basis_el_id = sub_grid_elem_map.at(funct_elem_id);
// [get_basis_el_id]
// [move_basis_el]
      basis_el->move_to(basis_el_id);
// [move_basis_el]
// [feel_da_cash]
      basis_handler->fill_face_cache(basis_el,face);
// [feel_da_cash]
      
// [get_data]
      auto vals   = basis_el->template get_basis_data<space_element::_Value,dim-1>(face);
      auto w_meas = s_dom_el->get_element_w_measures();
      auto funct  = funct_el->get_element_values_D0();
// [get_data]

// [computational_kernel]
      auto num_quad  = vals.get_num_points();
      auto num_basis = vals.get_num_functions();
      DenseVector loc_rhs(num_basis); loc_rhs = 0.0;
      for (int i=0; i<num_basis; i++) {
        const auto &vals_i = vals.get_function_view(i);
        for (int q=0; q<num_quad; q++) {
          loc_rhs(i) += funct[q][0] * vals_i[q][0] * w_meas[q];
        }
      }

      const auto loc_dofs = basis_el->get_local_to_global();
      rhs->add_block(loc_dofs, loc_rhs);
    }
  }
// [computational_kernel]
  
  auto dof_distribution = phy_basis->get_dof_distribution();
  std::map<Index,Real> dirichlet_vals;
  space_tools::project_boundary_values(dirichlet_cond,*phy_basis,face_quad,dirichlet_vals);
  dof_tools::apply_boundary_values(dirichlet_vals,*mat,*rhs,*sol);
}

template<int dim>
void PoissonProblem<dim>::solve() {
  auto solver = create_solver(*mat,*sol,*rhs);
  solver->solve();
}

template<int dim>
Real PoissonProblem<dim>::error(shared_ptr<const Function<dim>> exact_solution) {

  auto domain         = phy_basis->get_domain();
  auto domain_el      = domain->begin();
  auto domain_el_end  = domain->end();
  auto domain_handler = domain->create_cache_handler();
  domain_handler->set_element_flags(domain_element::Flags::w_measure);
  domain_handler->init_element_cache(domain_el,quad);

  auto discr_solution    = IgFunction<dim,0,1,1>::const_create(phy_basis, *sol);
  auto discr_sol_el      = discr_solution->begin();
  auto discr_sol_handler = discr_solution->create_cache_handler();
  discr_sol_handler->set_element_flags(function_element::Flags::D0 |
                                       function_element::Flags::w_measure);
  discr_sol_handler->init_cache(discr_sol_el,quad);

  auto exact_sol_el      = exact_solution->begin();
  auto exact_sol_handler = exact_solution->create_cache_handler();
  exact_sol_handler->set_element_flags(function_element::Flags::D0);
  exact_sol_handler->init_cache(exact_sol_el,quad);

  Real error = 0.0;
  for (; domain_el!=domain_el_end; ++domain_el, ++discr_sol_el, ++exact_sol_el) {
    domain_handler->fill_element_cache(domain_el);
    discr_sol_handler->fill_element_cache(discr_sol_el);
    exact_sol_handler->fill_element_cache(exact_sol_el);

    auto w_meas    = domain_el->get_element_w_measures();
    auto discr_val = discr_sol_el->get_element_values_D0();
    auto exact_val = exact_sol_el->get_element_values_D0();
    auto num_quad  = w_meas.get_num_points();

    for (int q=0; q<num_quad; q++) {
      const auto diff = discr_val[q]-exact_val[q];
      error += diff.norm_square() * w_meas[q];
    }
  }
  return std::sqrt(error);
}
 
template<int dim>
void PoissonProblem<dim>::save() {
  auto domain = phy_basis->get_domain();
  Writer<dim> writer(domain,5);
  auto discrete_solution = IgFunction<dim,0,1,1>::const_create(phy_basis,*sol);
  writer.add_field(*discrete_solution, "solution");
  string filename = "problem_" + to_string(dim) + "d" ;
  writer.save(filename);
}


shared_ptr<const Domain<2>> quarter_annulus(const Size nel) {
  using numbers::PI;
  BBox<2> box;
  box[0] = {{1.0,2.0}};
  box[1] = {{0.0,PI/2}};
  auto grid = Grid<2>::const_create(box,nel+1);
  auto geom_funct = grid_functions::BallGridFunction<2>::const_create(grid);
  return Domain<2>::const_create(geom_funct);
}

// [exact_solution]
Values<2,1,1> u(Points<2> x) {
  Values<2,1,1> y;
  y = exp(x[0])*sin(x[1]);
  return y;
}
// [exact_solution]
// [dirichlet_conds]
Values<1,1,1> gi(Points<2> x) {
  Values<1,1,1> y;
  y = exp(x[0])*sin(x[1]);
  return y;
}
// [dirichlet_conds]
// [neumann_conds]
Values<1,1,1> h3(Points<2> x) {
  Values<1,1,1> y;
  y = -sin(x[1]);
  return y;
}
// [neumann_conds]

int main()
{
  const int nel = 8;
  const int deg = 2;

  auto domain = quarter_annulus(nel);
  auto source  = functions::ConstantFunction<2,0,1,1>::const_create(domain,{0.0});

// [impose_boundary_conds]
  using SubGridElemMap = typename Grid<2>::template SubGridMap<1>;
  SubGridElemMap sub_grid_elem_map;
  auto grid = domain->get_grid_function()->get_grid();
  std::map<Index,shared_ptr<const Function<1,1,1,1>>> dirichlet;
  std::map<Index,shared_ptr<const Function<1,1,1,1>>> neumann;

  for (int face=0; face<3; face++) {
    const auto sub_grid    = grid->template get_sub_grid<1>(face,sub_grid_elem_map);
    const auto sub_annulus = domain->template get_sub_domain<1>(face,sub_grid_elem_map,sub_grid);
    auto g_i = CustomFunction<1,1,1>::const_create(sub_annulus,gi);
    dirichlet[face] = dynamic_pointer_cast<const Function<1,1,1>>(g_i);
  }
  const auto sub_grid_3    = grid->template get_sub_grid<1>(3,sub_grid_elem_map);
  const auto sub_annulus_3 = domain->template get_sub_domain<1>(3,sub_grid_elem_map,sub_grid_3);
  auto g_3 = CustomFunction<1,1,1>::const_create(sub_annulus_3,h3);
  neumann[3] = dynamic_pointer_cast<const Function<1,1,1>>(g_3);
// [impose_boundary_conds]

  PoissonProblem<2> problem(domain,deg,source,dirichlet,neumann);
  problem.assemble();
  problem.solve();
  auto exact_solution = CustomFunction<2>::const_create(domain,u);
  printf("The L^2-error is :  %e\n",problem.error(exact_solution));
  problem.save();

  return 0;
}