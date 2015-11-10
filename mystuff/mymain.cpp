// headers for spline spaces
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
// headers for quadrature
#include <igatools/base/quadrature_lib.h>
// headers for linear algebra objects
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
// headers for system solving
#include <igatools/linear_algebra/epetra_solver.h>
#include <igatools/linear_algebra/dof_tools.h>
// headers for function representations
#include <igatools/functions/function_lib.h>
#include <igatools/geometry/grid_function_lib.h>
// headers for output
#include <igatools/io/writer.h>
#include <igatools/base/logstream.h>
// finally! my stuff!
#include "my_formula_grid_function.h"
#include "custom_grid_function.h"

using namespace iga;
using namespace std;
using namespace EpetraTools;
LogStream out;

// ----------------------------------------------------------------------------
//    MY SUPERCOOL CLASS TURBO SPECIAL CONTAINER
// ----------------------------------------------------------------------------
template<int dim_>
class MySpace {

  private:
    using self_t = MySpace<dim_>; // creating the alias
    static int count;

    // constructors
  public:
    MySpace() = delete;
    MySpace(const Size nel, const Size deg);
    MySpace(const TensorSize<dim_> nel, const TensorIndex<dim_> deg);
    static const int dim = dim_;

  public:
      // space data
    shared_ptr<const Grid<dim_>>        grid;
    shared_ptr<const SplineSpace<dim_>> space;
    shared_ptr<const BSpline<dim_>>     basis;
      // linear system data
    shared_ptr<Matrix> mat;
    shared_ptr<Vector> rhs;
    shared_ptr<Vector> sol;

      // methods: creators
    static std::shared_ptr<self_t>
    create(const Size nel,const Size deg);

    static std::shared_ptr<const self_t> 
    const_create(const Size nel, const Size deg);

      // methods:
    void how_are_you_doin() const;
    void print_system() const;
    void assemble(bool print) const;
    void solve() const;
};
template<int dim_>
int MySpace<dim_>::count=0;

// ----------------------------------------------------------------------------
// constructors: the simple one
template<int dim_>
MySpace<dim_>::MySpace(const Size nel, const Size deg)
  : grid  {Grid<dim_>::const_create(nel+1)}
  , space {SplineSpace<dim_>::const_create(deg,grid)}
  , basis {BSpline<dim_>::const_create(space)}
  , mat   {create_matrix(*basis,DofProperties::active,Epetra_SerialComm())}
  , rhs   {create_vector(mat->RangeMap())}
  , sol   {create_vector(mat->DomainMap())}
{count++;}
// constructors: the cool one
template<int dim_>
MySpace<dim_>::MySpace(const TensorSize<dim_> nel, const TensorIndex<dim_> deg) {
  TensorSize<dim_> nknt;
  for (int idim=0; idim<dim_; idim++) nknt[idim]=nel[idim]+1;
  grid  = Grid<dim_>::const_create(nknt);
  space = SplineSpace<dim_>::const_create(deg,grid);
  basis = BSpline<dim_>::const_create(space);
  mat   = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
  rhs   = create_vector(mat->RangeMap());
  sol   = create_vector(mat->DomainMap());
  count++;
}

// ----------------------------------------------------------------------------
// creators: the simple ones in const/non const versions
template<int dim_>
auto MySpace<dim_>::create(const Size nel, const Size deg) -> shared_ptr<self_t> {
  return shared_ptr<self_t>(new self_t(nel,deg));
} 

template<int dim_>
auto MySpace<dim_>::const_create(const Size nel, const Size deg) -> shared_ptr<const self_t> {
  return create(nel,deg);
}

// ----------------------------------------------------------------------------
// methods: check that everybody is fine
template<int dim_>
void MySpace<dim_>::how_are_you_doin() const {
  out << "regardless copyright infringements, this is the MySpace class instance number " << count << ":" << endl;
  out << "        elements: " << grid->get_num_elements() << " = " << grid->get_num_intervals() << endl;
  out << "         degrees: " << space->get_degree_table()[0] << endl;
  TensorSize<dim_> nbf;
  for (int idim=0; idim<dim_; idim++) nbf[idim]=space->get_num_basis(0,idim);
  out << " basis functions: " << space->get_num_basis() << " = " << nbf << endl;
  out << "   system matrix: " << mat->NumGlobalRows() << " x " << mat->NumGlobalCols() << endl;
}
template<int dim_>
void MySpace<dim_>::print_system() const {
  if (mat->Filled()) {
    auto size = mat->NumMyRows();
    for (int irow=0; irow<size; irow++) {  
      // extracting the global row
      double *val; int *ind; int nnz;
      int ierr = mat->ExtractMyRowView(irow,nnz,val,ind); if (ierr!=0) cout << "I hate Trilinos!" << endl;
      double *bval;
      ierr = rhs->ExtractView(&bval); if (ierr!=0) cout << "I hate Trilinos!" << endl;
      int inz=0;
      cout << " |";
      for (int icol=0; icol<size; icol++) {
        if (icol==ind[inz]) {
          double en = val[inz];
          if (fabs(en)<1e-10) printf("   .o ");
          else if (en>0)      printf("  %1.2f",en);
          else                printf(" %1.2f",en);
          inz++;
        }
        else printf("   .  ");
      }
      if (irow==0)
        printf(" | = |");
      else
        printf(" |   |");
      double en=bval[irow];
      if (fabs(en)<1e-10) printf("   .o ");
      else if (en>0)      printf("  %1.2f",en);
      else                printf(" %1.2f",en);
      printf(" |\n");
    }
    cout << endl;
  }
  else printf("system is not filled yet!\n");
}


// methods: system matrix and right hand side vector assemble
template<int dim_>
void MySpace<dim_>::assemble(bool print) const {

  // starting the cache handler for the basis functions:
  auto basis_handler = basis->create_cache_handler();
  auto basis_el = basis->begin();
  const auto basis_eld = basis->end();
  // setting the flags
  using Flags = space_element::Flags;
  auto flag = Flags::value | Flags::gradient | Flags::w_measure;
  basis_handler->set_element_flags(flag);
  // setting the quarature rule
  TensorSize<dim> nqn;
  typename SplineSpace<dim_>::DegreeTable deg = space->get_degree_table();
  for (int idim=0; idim<dim_; idim++) nqn[idim] = deg[0][idim]+1;
  auto quad = QGauss<dim_>::create(nqn);
  auto Nqn = quad->get_num_points();
  basis_handler->init_element_cache(basis_el,quad);

  // starting the cache handler for the (constant) function f:
  typename Function<dim_,0,1,1>::Value f_val {5.0};
  const auto f = grid_functions::CustomGridFunction<dim_,1>::const_create(grid);
  auto funct_handler = f->create_cache_handler();
  auto funct_el = f->begin();
  funct_handler->template set_flags<dim_>(grid_function_element::Flags::D0);
  funct_handler->init_cache(funct_el,quad);

  // retrieving the last datum and then starting the loop
  for (int iel=0; basis_el!=basis_eld; ++basis_el) {
    basis_handler->fill_element_cache(basis_el);
    funct_handler->fill_element_cache(funct_el);
    // preparing some stuff: creating the local matrices
    auto Nbf = basis_el->get_num_basis(DofProperties::active);
    DenseMatrix loc_mat(Nbf,Nbf); loc_mat=0.0;
    DenseVector loc_rhs(Nbf);     loc_rhs=0.0;
    // gathering the requested data
    auto values = basis_el->get_element_values();
    auto grads  = basis_el->get_element_gradients();
    auto w_meas = basis_el->get_element_w_measures();
    auto f_vals = funct_el->get_element_values_D0();
    // finally, the loop
    for (int ibf=0; ibf<Nbf; ibf++) {
      // loop for the stiffness local matrix
      const auto &Dbfi = grads.get_function_view(ibf); // view for the i-th basis function gradient
      for (int jbf=0; jbf<Nbf; jbf++) {
        const auto &Dbfj = grads.get_function_view(jbf); // view for the j-th basis function gradient
        for (int iqn=0; iqn<Nqn; iqn++) {
	        loc_mat(ibf,jbf) += scalar_product(Dbfi[iqn],Dbfj[iqn]) * w_meas[iqn];
	      }
      }
      // loop for the right hand side local vector
      const auto &bfi = values.get_function_view(ibf);
      for (int iqn=0; iqn<Nqn; iqn++) {
        loc_rhs(ibf) += scalar_product(bfi[iqn],f_vals[iqn]) * w_meas[iqn]; 
      }
    }
    iel++;
    // spatashing element matrix into the global matrix
    const auto loc_dofs = basis_el->get_local_to_global(DofProperties::active);
    mat->add_block(loc_dofs,loc_dofs,loc_mat);
    rhs->add_block(loc_dofs,loc_rhs);
  }
  mat->FillComplete();
  if (print) print_system();

  using space_tools::project_boundary_values;
  using dof_tools::apply_boundary_values;
  //using RefSpace = ReferenceSpaceBasis<dim>;
  //using Function = Function<dim_,0,1,1>;

  // applying boundary condition
  /*typename Function<dim_,0,1,1>::Value g_val {0.0};
  const auto g = grid_functions::ConstantGridFunction<dim_,1>::const_create(grid,{0.0});
  auto face_quad = QGauss<dim_-1>::create(nqn[0]);
  const set<boundary_id> dir_id {0};
  std::map<Index,Real> values;
  project_boundary_values<dim_,1>(*g,*basis,face_quad,dir_id,values);*/
  const set<boundary_id> dir_id {0};
  auto bdr_dofs = space_tools::get_boundary_dofs<ReferenceSpaceBasis<dim_>>(basis,dir_id);
  std::map<Index,Real> bdr_vals;

  //int ind; double val=0.0;
  for (set<Index>::iterator it=bdr_dofs.begin(); it!=bdr_dofs.end(); it++) {
    bdr_vals[*it]=0.0;
    //rhs->ReplaceMyValues(1,&val,&ind);
  }

  //cout << bdr_vals << endl;
  apply_boundary_values(bdr_vals,*mat,*rhs,*sol);
  if (print) print_system();
}

template<int dim_>
void MySpace<dim_>::solve() const {
  auto solver = create_solver(*mat,*sol,*rhs);
  solver->solve();
}

// ----------------------------------------------------------------------------
//   MY CUSTOM FUNCTION
// ----------------------------------------------------------------------------
template<int dim>
double test_function(Points<dim> pts) {
  return 711;
}

int main() {

  const int dim  = 2;
  const int nel  = 4; 
  const int deg  = 2;

  // cool constructor for everything: grid, space, basis
  //TensorSize<dim>  mnel; for (int idim=0; idim<dim; idim++) mnel[idim]=nel+idim;
  //TensorIndex<dim> mdeg; for (int idim=0; idim<dim; idim++) mdeg[idim]=deg+idim;
  // testing the simple constructor
  /*MySpace<dim> space1(nel,deg);
  space1.how_are_you_doin();
  // testing the cool constructor
  MySpace<dim> space2(mnel,mdeg);
  space2.how_are_you_doin();
  // testing the simple creator
  auto space3 = MySpace<dim>::create(nel,deg+2);
  space3->how_are_you_doin();*/
  // testing the simple constant creator
  //auto space4 = MySpace<dim>::const_create(nel,deg);
  //space4->how_are_you_doin();
  //space4->assemble(false);
  //space4->solve();

  /*auto test_class = grid_functions::CustomGridFunction<dim,1>::const_create(space4->grid);
  auto quad = QGauss<dim>::create(2);
  auto funct_handler = test_class->create_cache_handler();
  auto funct_el = test_class->GridFunction<dim,1>::begin();
  const auto funct_eld = test_class->GridFunction<dim,1>::end();
  funct_handler->template set_flags<dim>(grid_function_element::Flags::D0);
  funct_handler->template set_flags<dim>(grid_function_element::Flags::D1);
  funct_handler->init_element_cache(funct_el,quad);
  for (int iel=0; funct_el!=funct_eld; ++funct_el) {
    funct_handler->fill_element_cache(funct_el);
    auto f_vals = funct_el->get_element_values_D0();
    auto f_grad = funct_el->get_element_values_D1();
    std::cout << "evaluation of element " << iel << std::endl;
    f_vals.print_info(out);
    out << endl;
    f_grad.print_info(out);
    out << endl;
    iel++;
  }*/

  auto grid = Grid<dim>::create(nel);
  auto test_class = grid_functions::CustomGridFunction<dim,1>::create(grid);
  //test_class->grid_functions::CustomGridFunction<dim,1>::funct[0] = &test_function;
  test_class->funct[0]=&test_function;

  return 0;
}
