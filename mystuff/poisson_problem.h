//#include<igatools/linear_algebra/epetra.h>

#include "AztecOO_config.h"
#include "AztecOO.h"
#include "AztecOO_ConditionNumber.h"
#include "Ifpack_PointRelaxation.h"

template<int dim>
class Geometry {
  public:
  TensorSize<dim>  nel;
  TensorIndex<dim> deg;
  IgCoefficients   coefs;
  IgCoefficients   weights;
};

template<int dim_>
class PoissonProblem {

  private:
    using self_t = PoissonProblem<dim_>; // creating the alias
    static int count;

    // constructors
  public:
    PoissonProblem() = delete;
    PoissonProblem(const Size nel, const Size deg);
    PoissonProblem(const TensorSize<dim_> nel, const TensorIndex<dim_> deg);
    PoissonProblem(const Size nel, const Size deg, const Geometry<dim_> &geom);

    static const int dim = dim_;

  public:
      // space data
    shared_ptr<const Grid<dim_>>           grid;
    shared_ptr<const SplineSpace<dim_>>    space;
    shared_ptr<const Basis<dim_,0,1,1>>    basis;
    shared_ptr<const QGauss<dim_>>         quad;
    shared_ptr<const Domain<dim_>>       domain;
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
    //using Funct = typename grid_functions::FormulaGridFunction<dim_,1>::CustomGridFunction<dim_,1>;
    void assemble(std::shared_ptr<const GridFunction<dim_,1>> f) const;
    void solve() const;
    void custom_solve(int &it1, double &cond1, int &it2, double &cond2) const;
    Real l2_error(std::shared_ptr<const GridFunction<dim_,1>> u) const;
};
template<int dim_>
int PoissonProblem<dim_>::count=0;

// ----------------------------------------------------------------------------
//   CONSTRUCTORS
// ----------------------------------------------------------------------------
template<int dim_> // constructor for reference domain problem
PoissonProblem<dim_>::PoissonProblem(const Size nel, const Size deg)
  : grid  {Grid<dim_>::const_create(nel+1)}
  , space {SplineSpace<dim_>::const_create(deg,grid)}
  , basis {BSpline<dim_>::const_create(space)}
  , quad  {QGauss<dim_>::const_create(deg+1)}
  , mat   {create_matrix(*basis,DofProperties::active,Epetra_SerialComm())}
  , rhs   {create_vector(mat->RangeMap())}
  , sol   {create_vector(mat->DomainMap())}
{count++;}

template<int dim_> // constructor for reference domain problem
PoissonProblem<dim_>::PoissonProblem(const TensorSize<dim_> nel, const TensorIndex<dim_> deg) {
  TensorSize<dim_> nknt, nqn;
  for (int idim=0; idim<dim_; idim++) {
    nknt[idim] = nel[idim]+1;
    nqn[idim]  = deg[idim]+1;
  }
  grid  = Grid<dim_>::const_create(nknt);
  space = SplineSpace<dim_>::const_create(deg,grid);
  basis = BSpline<dim_>::const_create(space);
  quad  = QGauss<dim_>::const_create(nqn);
  mat   = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
  rhs   = create_vector(mat->RangeMap());
  sol   = create_vector(mat->DomainMap());
  count++;
}

template<int dim_> // constructor for pysical domain problem
PoissonProblem<dim_>::PoissonProblem(const Size nel, const Size deg, const Geometry<dim_> &geom) {
  // CREATING THE PHYSICAL DOMAIN
  // geometry grid
  TensorIndex<dim_> nknt;
  TensorSize<dim_>  nqn;
  TensorSize<dim_>  neel;
  for (int idim=0; idim<dim_; idim++) {
    nknt[idim] = geom.nel[idim]+1;
    nqn[idim]  = deg+1;
    neel[idim] = geom.nel[idim]+1;
  }
  grid         = Grid<dim_>::const_create(neel);
  // B-spline vector field for the geometry function
  auto vect_space   = SplineSpace<dim,dim>::const_create(geom.deg,grid);
  auto vect_bspline = BSpline<dim,dim>::const_create(vect_space);
  // B-spline scalar field for the weight function
  auto scal_space   = SplineSpace<dim,1>::const_create(geom.deg,grid);
  auto scal_bspline = BSpline<dim,1>::const_create(scal_space);
  // the weight function
  auto weight_funct = IgGridFunction<dim,1>::const_create(scal_bspline,geom.weights);
  // NURBS vector field for the geometry function
  auto vect_nurbs   = NURBS<dim,dim>::const_create(vect_bspline,weight_funct);
  // geometry function 
  auto geom_funct   = IgGridFunction<dim,dim>::const_create(vect_bspline,geom.coefs);
  // finally, the domain
  domain   = Domain<dim>::const_create(geom_funct);
  // the basis functions in the physical domain
  basis    = PhysicalSpaceBasis<dim_>::const_create(scal_bspline,domain);
  // refine everything
  //grid->refine();
  // qudrature rule, linear system 
  quad  = QGauss<dim_>::const_create(nqn);
  mat   = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
  rhs   = create_vector(mat->RangeMap());
  sol   = create_vector(mat->DomainMap());
}

// ----------------------------------------------------------------------------
//   CREATORS
// ----------------------------------------------------------------------------
template<int dim_> // creator: the non const simple one
auto PoissonProblem<dim_>::create(const Size nel, const Size deg) -> shared_ptr<self_t> {
  return shared_ptr<self_t>(new self_t(nel,deg));
}

template<int dim_> // creator: the const simple one
auto PoissonProblem<dim_>::const_create(const Size nel, const Size deg) -> shared_ptr<const self_t> {
  return create(nel,deg);
}

/*template<int dim_> // creator
auto PoissonProblem<dim_>::const_create(const Size nel, const Size deg, const Geometry<dim_> &geom) -> shared_ptr<self_t> {
  return shared_ptr<self_t>(new self_t(nel,deg));
}

template<int dim_> // creator
auto PoissonProblem<dim_>::const_create(const Size nel, const Size deg, const Geometry<dim_> &geom) -> shared_ptr<const self_t> {
  return create(nel,deg,geom);
}*/

// ----------------------------------------------------------------------------
//   METHODS
// ----------------------------------------------------------------------------
template<int dim_> // assemble the system
void PoissonProblem<dim_>::assemble(std::shared_ptr<const GridFunction<dim_,1>> source_term) const {

  // starting the cache handler for the basis functions:
  auto basis_handler = basis->create_cache_handler();
  auto basis_el = basis->begin();
  const auto basis_eld = basis->end();
  // setting the flags
  using Flags = space_element::Flags;
  auto flag = Flags::value | Flags::gradient | Flags::w_measure;
  basis_handler->set_element_flags(flag);
  // setting the quarature rule
  auto Nqn = quad->get_num_points();
  basis_handler->init_element_cache(basis_el,quad);

  // starting the cache handler for the (constant) function f:
  auto funct_handler = source_term->create_cache_handler();
  auto funct_el = source_term->begin();
  //funct_handler->template set_flags<dim_>(function_element::Flags::D0);
  funct_handler->set_element_flags(grid_function_element::Flags::D0);
  funct_handler->init_cache(funct_el,quad);

  // retrieving the last datum and then starting the loop
  for (int iel=0; basis_el!=basis_eld; ++basis_el, ++funct_el) {
    basis_handler->fill_element_cache(basis_el);
    //funct_handler->fill_element_cache(funct_el);
    //funct_handler->template fill_cache<dim_>(funct_el,0);
    // preparing some stuff: creating the local matrices
    auto Nbf = basis_el->get_num_basis(DofProperties::active);
    DenseMatrix loc_mat(Nbf,Nbf); loc_mat=0.0;
    DenseVector loc_rhs(Nbf);     loc_rhs=0.0;
    // gathering the requested data
    auto values = basis_el->get_element_values();
    auto grads  = basis_el->get_element_gradients();
    auto w_meas = basis_el->get_element_w_measures();
    //using _D0 = typename function_element::template _D<0>;
    //auto f_vals = funct_el->template get_values_from_cache<_D0,dim_>(0);
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
      //const auto &bfi = values.get_function_view(ibf); // view for the i-th basis function
      //for (int iqn=0; iqn<Nqn; iqn++) {
      //  loc_rhs(ibf) += scalar_product(bfi[iqn],f_vals[iqn]) * w_meas[iqn]; 
      //}
    }
    iel++;
    // spatashing element matrix into the global matrix
    const auto loc_dofs = basis_el->get_local_to_global(DofProperties::active);
    mat->add_block(loc_dofs,loc_dofs,loc_mat);
    //rhs->add_block(loc_dofs,loc_rhs);
  }
  mat->FillComplete();

  // applying the boundary conditions
  using space_tools::project_boundary_values;
  using dof_tools::apply_boundary_values;
  const set<boundary_id> dir_id {0};
  auto bdr_dofs = space_tools::get_boundary_dofs<Basis<dim_,0,1,1>>(basis,dir_id);
  std::map<Index,Real> bdr_vals;

  for (set<Index>::iterator it=bdr_dofs.begin(); it!=bdr_dofs.end(); it++) {
    bdr_vals[*it]=0.0;
  }
  rhs->Random();
  apply_boundary_values(bdr_vals,*mat,*rhs,*sol);

}

template<int dim_> // solver for the linear system
void PoissonProblem<dim_>::solve() const {
  auto solver = create_solver(*mat,*sol,*rhs);
  solver->solve();
}


using OP = Epetra_Operator;
using MV = Epetra_MultiVector;
using SolverPtr = Teuchos::RCP<Belos::SolverManager<double, MV, OP> >;
template<int dim_> // custom siuppacool solver
void PoissonProblem<dim_>::custom_solve(int &it1, double &cond, int &it2, double &cond2) const {

  double cond1;

  // setting up the problem
  Epetra_LinearProblem problem1(&*mat,&*sol,&*rhs);
  AztecOO solver1(problem1);
  // setting up the solver
  solver1.SetAztecOption(AZ_solver,  AZ_cg);
  solver1.SetAztecOption(AZ_kspace,  100);
  solver1.SetAztecOption(AZ_precond, AZ_none);
  solver1.SetAztecOption(AZ_output,  AZ_none);
  // solve, for god's sake! SOLVE!
  solver1.Iterate(100, 1.0E-7);
  // extracting info
  it1 = solver1.NumIters();

  // condition number estimate
  AztecOOConditionNumber condest1;
  condest1.initialize(*mat);
  //AztecOOConditionNumber::SolverType solver_type = CG_;
  condest1.computeConditionNumber(100,1.0E-7);
  cond1 = condest1.getConditionNumber();

  // try this new trick!
  Epetra_Vector diag(mat->DomainMap());
  Epetra_Vector diag_reciprocal(mat->DomainMap());
  mat->ExtractDiagonalCopy(diag);
  //double *val;
  //diag.ExtractView(&val);
  for (int idof=0; idof<diag.MyLength(); idof++) {
    diag_reciprocal[idof]=sqrt(1.0/diag[idof]);
  }
  mat->LeftScale(diag_reciprocal);
  mat->RightScale(diag_reciprocal);

  // setting up the problem
  Epetra_Vector sol2(mat->DomainMap());
  Epetra_LinearProblem problem2(&*mat,&sol2,&*rhs);
  AztecOO solver2(problem2);
  // setting up the solver
  solver2.SetAztecOption(AZ_solver,  AZ_cg);
  solver2.SetAztecOption(AZ_kspace,  100);
  solver2.SetAztecOption(AZ_precond, AZ_none);
  solver2.SetAztecOption(AZ_output,  AZ_none);
  // solve, for god's sake! SOLVE!
  solver2.Iterate(100, 1.0E-7);
  // extracting info
  it2 = solver2.NumIters();

  // condition number estimate
  AztecOOConditionNumber cond_estimator;
  cond_estimator.initialize(*mat);
  cond_estimator.computeConditionNumber(100,1.0E-7);
  cond = cond_estimator.getConditionNumber();

}

template<int dim_> // compute the L2 error given the exact solution
Real PoissonProblem<dim_>::l2_error(std::shared_ptr<const GridFunction<dim_,1>> u_ex) const {
//                                 Real &l2_err, bool compute_h1, Real &h1_err = 0.0) const {

  // creating the discrete solution
  auto u_ds = IgGridFunction<dim_,1>::const_create(basis,*sol);
  
  // cache for the discrete solution
  auto u_ds_handler = u_ds->create_cache_handler();
  auto u_ds_el = u_ds->begin();
  u_ds_handler->set_element_flags(grid_function_element::Flags::D0);
  u_ds_handler->init_element_cache(u_ds_el,quad);

  // cache for the exact solution
  auto u_ex_handler = u_ex->create_cache_handler();
  auto u_ex_el = u_ex->begin();
  u_ex_handler->set_element_flags(grid_function_element::Flags::D0);
  u_ex_handler->init_element_cache(u_ex_el,quad);

  // cache for the grid
  auto grid_handler = grid->create_cache_handler();
  auto grid_el = grid->begin();
  const auto grid_eld = grid->end();
  grid_handler->set_element_flags(grid_element::Flags::weight);
  grid_handler->init_element_cache(grid_el,quad);  

  const auto nqn = quad->get_num_points();
  Real err = 0.0;
  for (; grid_el != grid_eld; ++grid_el, ++u_ex_el, ++u_ds_el) {
    u_ds_handler->fill_element_cache(u_ds_el);
    u_ex_handler->fill_element_cache(u_ex_el);
    grid_handler->fill_element_cache(grid_el);

    // retrieving data
    auto uh = u_ds_el->get_element_values_D0();
    auto ue = u_ex_el->get_element_values_D0();
    auto qw = grid_el->get_element_weights();

    // loopin' on the quadrature points
    for (int iqn=0; iqn<nqn; iqn++) {
      Real diff = uh[iqn][0][0] - ue[iqn][0][0];
      err += diff * diff * qw[iqn];
    }
  }

  return sqrt(err);
}

// ----------------------------------------------------------------------------
//   INFO PRINTER
// ----------------------------------------------------------------------------
template<int dim_> // problem overview
void PoissonProblem<dim_>::how_are_you_doin() const {
  out << endl;
  out << "        elements: " << grid->get_num_elements() << " = " << grid->get_num_intervals() << endl;
  out << "         degrees: " << space->get_degree_table()[0] << endl;
  TensorSize<dim_> nbf;
  for (int idim=0; idim<dim_; idim++) nbf[idim]=space->get_num_basis(0,idim);
  out << " basis functions: " << space->get_num_basis() << " = " << nbf << endl;
  out << "   system matrix: " << mat->NumGlobalRows() << " x " << mat->NumGlobalCols() << endl;
}

template<int dim_> // system matrix printer
void PoissonProblem<dim_>::print_system() const {
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
