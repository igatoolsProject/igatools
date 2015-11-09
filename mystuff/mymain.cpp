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
    void assemble() const;
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

// methods: system matrix and right hand side vector assemble
template<int dim_>
void MySpace<dim_>::assemble() const {

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
  //using functions::ConstantFunction;
  //using ConstFunction = ConstantFunction<dim_,0,dim_,1>;
  auto id = grid_functions::IdentityGridFunction<dim_>::const_create(grid);
  auto f = functions::ConstantFunction<dim_,0,1,1>::create(grid,f_val); 


  // retrieving the last datum and then starting the loop
  for (int iel=0; basis_el!=basis_eld; ++basis_el) {
    basis_handler->fill_element_cache(basis_el);
    // preparing some stuff: creating the local matrices
    auto Nbf = basis_el->get_num_basis(DofProperties::active);
    DenseMatrix loc_mat(Nbf,Nbf); loc_mat=0.0;
    DenseVector loc_rhs(Nbf);     loc_rhs=0.0;
    // gathering the requested data
    auto values = basis_el->get_element_values();
    auto grads  = basis_el->get_element_gradients();
    auto w_meas = basis_el->get_element_w_measures();
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
        loc_rhs(ibf) += bfi[ibf][0] * w_meas[ibf]; 
      }
    }
    // plotting the computed stuff
    /*out << "element matrix " << basis_el->get_index() << " = " << iel << endl;
    //out << lA << endl;
    for (int ibf=0; ibf<Nbf; ibf++) {
      for (int jbf=0; jbf<Nbf; jbf++) {
        auto en = loc_mat(ibf,jbf);
        if (fabs(en)<1e-10) printf("   o.o  ");
        else if (en>0)      printf("  %1.4f",en);
        else                printf(" %1.4f",en);
      }
      cout << endl;
    }
    cout << endl;*/
    iel++;
    // spatashing element matrix into the global matrix
    const auto loc_dofs = basis_el->get_local_to_global(DofProperties::active);
    mat->add_block(loc_dofs,loc_dofs,loc_mat);
    rhs->add_block(loc_dofs,loc_rhs);
  }
  mat->FillComplete();

  auto size = mat->NumMyRows();
  for (int irow=0; irow<size; irow++) {  
    // extracting the global row
    double *val; int *ind; int nnz;
    const int ierr = mat->ExtractMyRowView(irow,nnz,val,ind); if (ierr!=0) cout << "I hate Trilinos!" << endl;
    int inz=0;
    for (int icol=0; icol<size; icol++) {
      if (icol==ind[inz]) {
        double en = val[inz];
        if (fabs(en)<1e-10) printf("  o.o ");
        else if (en>0)      printf("  %1.2f",en);
        else                printf(" %1.2f",en);
        inz++;
      }
      else printf("   .  ");
    }
    cout << endl;
  }
  cout << endl;

  // apply boundary condition


}

int main() {

  const int dim  = 2;
  const int nel  = 2; 
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
  auto space4 = MySpace<dim>::const_create(nel,deg);
  space4->how_are_you_doin();
  space4->assemble();

  return 0;
}
