
#include <igatools/basis_functions/bspline.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/nurbs.h>
#include <igatools/functions/grid_function_lib.h>
#include <igatools/functions/ig_function.h>
//#include <igatools/io/writer.h>

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


#include <igatools/io/objects_container_xml_writer.h>
#include <igatools/base/objects_container.h>

using namespace iga;
using namespace std;
using namespace EpetraTools;
LogStream out;

// -------------------------------------------------------
//   custom function
// -------------------------------------------------------
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
    std::cout << "C makes it easy to shoot yourself in the foot." << std::endl;
    std::cout << "C++ makes it harder, but when you do, it blows away your whole leg." << std::endl;
  };
};

// -------------------------------------------------------
//   class definition
// -------------------------------------------------------
template<int dim>
class LinearElasticity {
  private:
    // spaces
    shared_ptr<const SplineSpace<dim,dim,1>>          ref_space;
    shared_ptr<const BSpline<dim,dim,1>>              ref_basis;
    shared_ptr<const PhysicalBasis<dim,dim,1,0>> phy_basis;
    // quadrature rules
    shared_ptr<const QGauss<dim>>   quad;
    shared_ptr<const QGauss<dim-1>> face_quad;
    // coefficients
    Real lambda;
    Real mu;
    // initial data
    shared_ptr<const Function<dim,0,dim,1>> source_term;
    std::map<Index,shared_ptr<const Function<dim-1,1,dim,1>>> dirichlet_cond;
    std::map<Index,shared_ptr<const Function<dim-1,1,dim,1>>> neumann_cond;
    // linear system
    shared_ptr<Matrix> mat;
    shared_ptr<Vector> rhs;
    shared_ptr<Vector> sol;

    using IgFunc_t = IgFunction<dim,0,dim,1>;
    
  public:
    // constructor
    LinearElasticity(const shared_ptr<const Domain<dim>> domain, const Index deg,
                     const Real l, const Real m,
                     const shared_ptr<const Function<dim,0,dim,1>> source,
                     const std::map<Index,shared_ptr<const Function<dim-1,1,dim,1>>> dirichlet,
                     const std::map<Index,shared_ptr<const Function<dim-1,1,dim,1>>> neumann) {
      lambda = l;
      mu     = m;
      source_term    = source;
      dirichlet_cond = dirichlet;
      neumann_cond   = neumann;
      auto grid = domain->get_grid_function()->get_grid();
      ref_space = SplineSpace<dim,dim,1>::const_create(deg,grid);
      ref_basis = BSpline<dim,dim,1>::const_create(ref_space);
      phy_basis = PhysicalBasis<dim,dim,1,0>::const_create(ref_basis,domain);
      quad      = QGauss<dim>::const_create(deg+1);
      face_quad = QGauss<dim-1>::const_create(deg+1);
      mat = create_matrix(*phy_basis,DofProperties::active,Epetra_SerialComm());
      rhs = create_vector(mat->RangeMap());
      sol = create_vector(mat->DomainMap());
    };
    // methods
    void assemble();
    void solve();
    shared_ptr<const IgFunc_t> get_solution();
    void save();
    void save(shared_ptr<const Function<dim,0,dim,1>> exact_solution);
    void save_plugin();
    void save_plugin(shared_ptr<const Function<dim,0,dim,1>> exact_solution);
    Real error(shared_ptr<const Function<dim,0,dim,1>> exact_solution);
};

// -------------------------------------------------------
//   assemble
// -------------------------------------------------------
template<int dim>
void LinearElasticity<dim>::assemble() {

  auto basis_el      = phy_basis->begin();
  auto basis_el_end  = phy_basis->end();
  auto basis_handler = phy_basis->create_cache_handler();
  auto flag = space_element::Flags::value |
              space_element::Flags::gradient |
              space_element::Flags::divergence |
              space_element::Flags::w_measure;
  basis_handler->set_element_flags(flag);
  basis_handler->init_element_cache(basis_el,quad);

  auto funct_el      = source_term->begin();
  auto funct_handler = source_term->create_cache_handler();
  funct_handler->set_element_flags(function_element::Flags::D0);
  funct_handler->init_cache(funct_el,quad);

  auto num_basis  = basis_el->get_num_basis();
  auto num_points = quad->get_num_points();
  auto l = lambda;
  auto m = 2.0*mu;
  for (; basis_el!=basis_el_end; ++basis_el, ++funct_el) {
    // feel da cache
    basis_handler->fill_element_cache(basis_el);
    funct_handler->fill_element_cache(funct_el);

    // fetching data
    auto values = basis_el->get_element_values();
    auto grads  = basis_el->get_element_gradients();
    auto divs   = basis_el->get_element_divergences();
    auto w_meas = basis_el->get_element_w_measures();

    DenseMatrix loc_mat(num_basis,num_basis); loc_mat = 0.0;
    // element loop
    for (int i=0; i<num_basis; i++) {
      const auto &grd_i = grads.get_function_view(i);
      const auto &div_i = divs.get_function_view(i);
      for (int j=0; j<num_basis; j++) {
        const auto &grd_j = grads.get_function_view(j);
        const auto &div_j = divs.get_function_view(j);
        Real part_1 = 0.0; Real part_2 = 0.0;
        for (int q=0; q<num_points; q++) {
          part_1 += div_i[q][0] * div_j[q][0] * w_meas[q];
          part_2 += scalar_product(symmetric_tensor(grd_i[q]),symmetric_tensor(grd_j[q])) * w_meas[q];
          //part_2 += scalar_product(grd_i[q],grd_j[q]) * w_meas[q];
        }
        loc_mat(i,j) = l*part_1 + m*part_2;
      }
    }

    auto funct  = funct_el->get_element_values_D0();
    DenseVector loc_rhs = basis_el->integrate_u_func(funct);

    const auto loc_dofs = basis_el->get_local_to_global();
    mat->add_block(loc_dofs, loc_dofs,loc_mat);
    rhs->add_block(loc_dofs, loc_rhs);
  }
  mat->FillComplete();

  std::map<Index,Real> dirichlet_vals;
  space_tools::project_boundary_values(dirichlet_cond,*phy_basis,face_quad,dirichlet_vals);
  dof_tools::apply_boundary_values(dirichlet_vals,*mat,*rhs,*sol);
}

// -------------------------------------------------------
//   solve
// -------------------------------------------------------
template<int dim>
void LinearElasticity<dim>::solve() {
  auto solver = create_solver(*mat,*sol,*rhs);
  solver->solve();
}

// -------------------------------------------------------
//   L^2-error
// -------------------------------------------------------
template<int dim>
Real LinearElasticity<dim>::error(shared_ptr<const Function<dim,0,dim,1>> exact_solution) {

  auto domain         = phy_basis->get_domain();
  auto domain_el      = domain->begin();
  auto domain_el_end  = domain->end();
  auto domain_handler = domain->create_cache_handler();
  domain_handler->set_element_flags(domain_element::Flags::w_measure);
  domain_handler->init_element_cache(domain_el,quad);

  auto discr_solution    = IgFunction<dim,0,dim,1>::const_create(phy_basis, *sol);
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

// -------------------------------------------------------
//   vtk save
// -------------------------------------------------------
template<int dim>
void LinearElasticity<dim>::save() {
  auto domain = phy_basis->get_domain();
  Writer<dim> writer(domain,5);
  auto solution = IgFunction<dim,0,dim,1>::const_create(phy_basis,*sol);
  writer.add_field(*solution, "solution");
  string filename = "problem_" + to_string(dim) + "d" ;
  writer.save(filename);
}
// -------------------------------------------------------
//   vtk save with exact solution
// -------------------------------------------------------
template<int dim>
void LinearElasticity<dim>::save(shared_ptr<const Function<dim,0,dim,1>> exact_solution) {
  auto domain = phy_basis->get_domain();
  Writer<dim> writer(domain,5);
  auto discrete_solution = IgFunction<dim,0,dim,1>::const_create(phy_basis,*sol);
  writer.add_field(*discrete_solution, "discrete solution");
  writer.add_field(*exact_solution, "exact solution");
  string filename = "problem_" + to_string(dim) + "d" ;
  writer.save(filename);
}

// -------------------------------------------------------
//   Pablo's get_solution
// -------------------------------------------------------
template<int dim>
auto LinearElasticity<dim>::get_solution() -> std::shared_ptr<const IgFunc_t> {
  return IgFunc_t::const_create(phy_basis, *sol);
}
// -------------------------------------------------------
//   ParaView plugin save
// -------------------------------------------------------
template<int dim>
void LinearElasticity<dim>::save_plugin() {
  const auto solution = this->get_solution();
  string filename = "problem_plugin_" + to_string(dim) + "d" ;
  const auto solution_non_const = std::const_pointer_cast<IgFunc_t>(solution);
  solution_non_const->set_name("solution");
#ifdef SERIALIZATION
  const auto objs_container = ObjectsContainer::create();
  objs_container->insert_const_object<Function<dim,0,dim,1>>(solution);
  {
    std::ofstream xml_ostream(filename + ".iga");
    OArchive xml_out(xml_ostream);
    xml_out << *objs_container;
  }
#else
#ifdef XML_IO
  const auto objs_container = ObjectsContainer::create();
  objs_container->insert_const_object<Function<dim,0,dim,1>>(solution);
  ObjectsContainerXMLWriter::write(filename + ".iga", objs_container);
#else
  save();
#endif
#endif
}
// -------------------------------------------------------
//   ParaView plugin save with exact solution
// -------------------------------------------------------
template<int dim>
void LinearElasticity<dim>::save_plugin(shared_ptr<const Function<dim,0,dim,1>> exact_solution) {
  const auto solution = this->get_solution();
  string filename = "problem_plugin_" + to_string(dim) + "d" ;
  const auto solution_non_const = std::const_pointer_cast<IgFunc_t>(solution);
  solution_non_const->set_name("solution");
#ifdef SERIALIZATION
  printf(" -- serialization on ---\n");
  const auto objs_container = ObjectsContainer::create();
  objs_container->insert_const_object<Function<dim,0,dim,1>>(solution);
  objs_container->insert_const_object<Function<dim,0,dim,1>>(exact_solution);
  {
    std::ofstream xml_ostream(filename + ".iga");
    OArchive xml_out(xml_ostream);
    xml_out << *objs_container;
  }
#else
#ifdef XML_IO
  printf(" -- XML i/o on ---\n");
  const auto objs_container = ObjectsContainer::create();
  objs_container->insert_const_object<Function<dim,0,dim,1>>(solution);
  //auto exsol = dynamic_pointer_cast<const Function<dim,0,dim,1>>(exact_solution);
  objs_container->insert_const_object<Function<dim,0,dim,1>>(exact_solution);
  ObjectsContainerXMLWriter::write(filename + ".iga", objs_container);
#else
  save(exact_solution);
#endif
#endif
}// */

// -------------------------------------------------------
//   useful domains
// -------------------------------------------------------
shared_ptr<const Domain<2>> quarter_annulus(const Size nel) {
  using numbers::PI;
  BBox<2> box;
  box[0] = {{1.0,2.0}};
  box[1] = {{0.0,PI/2}};
  auto grid = Grid<2>::const_create(box,nel+1);
  auto geom_funct = grid_functions::BallGridFunction<2>::const_create(grid);
  return Domain<2>::const_create(geom_funct);
}

template<int dim>
shared_ptr<const Domain<dim>> hypercube(const Size nel) {
  auto grid = Grid<dim>::const_create(nel+1);
  auto idty = grid_functions::IdentityGridFunction<dim>::const_create(grid);
  return Domain<dim>::const_create(idty);
}


// -------------------------------------------------------
//   problem functions
// -------------------------------------------------------
#define LAMBDA 0.57692
#define MU     0.38462
using numbers::PI;
template<int dim>
Values<dim,dim,1> u(Points<dim> x) {
  Values<dim,dim,1> y;
  y[0] = 0.1;
  for (int d=0; d<dim; d++) y[0] *= sin(PI*x[d]);
  for (int d=1; d<dim; d++) y[d]  = 0.0;
  return y;
}
template<int dim>
Values<dim,dim,1> f(Points<dim> x) {
  Values<dim,dim,1> y;
  y[0] = 0.1*PI*PI*(LAMBDA + (dim+1)*MU);
  for (int d=0; d<dim; d++) y[0] *= sin(PI*x[d]);
  Real coef = -0.1*PI*PI*(LAMBDA + MU);
  for (int d=1; d<dim; d++) {
    y[d] = coef * cos(PI*x[0]) * cos(PI*x[d]);
    for (int e=1;   e<d;   e++) y[d] *= sin(PI*x[e]);
    for (int e=d+1; e<dim; e++) y[d] *= sin(PI*x[e]);
  }
  return y;
}

// -------------------------------------------------------
//   main body
// -------------------------------------------------------
int main()
{

  // problem dimension
  const int dim = 2;
  // problem domain
  const Size nel = 4;
  // discretization degrees
  const Index deg = 2;

//   Real e0=0.0, h0=0.0, h1=0.0;
//   for (int deg=1; deg<=4; deg++) {
//   printf("\n");
//   for (int nel=2; nel<=8; nel+=2) {
  
  auto domain = hypercube<dim>(nel);

  // source term
  auto source_term  = CustomFunction<dim,0,dim,1>::const_create(domain,f);
  // coefficients
  Real lambda = LAMBDA;
  Real mu     = MU;
  // homogeneous Dirichlet boundary conditions
  using SubGridElemMap = typename Grid<dim>::template SubGridMap<dim-1>;
  SubGridElemMap sub_grid_elem_map;
  auto grid = domain->get_grid_function()->get_grid();

  std::map<Index,shared_ptr<const Function<dim-1,1,dim,1>>> dirichlet;
  Values<dim,dim,1> zeros; for (int d=0; d<dim; d++) zeros[d] = 0.0;

  for (int face=0; face<2*dim; face++) {
    const auto sub_grid    = grid->template get_sub_grid<dim-1>(face,sub_grid_elem_map);
    const auto sub_annulus = domain->template get_sub_domain<dim-1>(face,sub_grid_elem_map,sub_grid);
    auto g = functions::ConstantFunction<dim-1,1,dim,1>::const_create(sub_annulus,zeros);
    dirichlet[face] = dynamic_pointer_cast<const Function<dim-1,1,dim,1>>(g);
  }

  std::map<Index,shared_ptr<const Function<dim-1,1,dim,1>>> neumann;
  
  // creating problems
  auto problem = LinearElasticity<dim>(domain,deg,lambda,mu,source_term,dirichlet,neumann);
  problem.assemble();
  problem.solve();
  auto exact_solution = CustomFunction<dim,0,dim,1>::const_create(domain,u);

  auto e1  = problem.error(exact_solution);
  printf("The L^2-error is : %1.5e\n",e1);
  problem.save(exact_solution);
  problem.save_plugin(exact_solution);

//   h1  = 1.0/nel;
//   auto eoc = log(e1/e0)/log(h1/h0);
//   printf("deg %d nel %2d:    error = %1.5e\teoc = %1.5f\n",deg,nel,e1,eoc);
//   e0=e1;
//   h0=h1;
//   }
//   }

  return 0;
}