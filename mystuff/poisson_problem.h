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
    static const int dim = dim_;

  public:
      // space data
    shared_ptr<const Grid<dim_>>        grid;
    shared_ptr<const SplineSpace<dim_>> space;
    shared_ptr<const BSpline<dim_>>     basis;
    shared_ptr<const QGauss<dim_>>      quad;
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
    Real l2_error(std::shared_ptr<const GridFunction<dim_,1>> u) const;
};
template<int dim_>
int PoissonProblem<dim_>::count=0;

// ----------------------------------------------------------------------------
// constructors: the simple one
template<int dim_>
PoissonProblem<dim_>::PoissonProblem(const Size nel, const Size deg)
  : grid  {Grid<dim_>::const_create(nel+1)}
  , space {SplineSpace<dim_>::const_create(deg,grid)}
  , basis {BSpline<dim_>::const_create(space)}
  , quad  {QGauss<dim_>::const_create(deg+1)}
  , mat   {create_matrix(*basis,DofProperties::active,Epetra_SerialComm())}
  , rhs   {create_vector(mat->RangeMap())}
  , sol   {create_vector(mat->DomainMap())}
{count++;}
// constructors: the cool one
template<int dim_>
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

// ----------------------------------------------------------------------------
// creators: the simple ones in const/non const versions
template<int dim_>
auto PoissonProblem<dim_>::create(const Size nel, const Size deg) -> shared_ptr<self_t> {
  return shared_ptr<self_t>(new self_t(nel,deg));
} 

template<int dim_>
auto PoissonProblem<dim_>::const_create(const Size nel, const Size deg) -> shared_ptr<const self_t> {
  return create(nel,deg);
}

// ----------------------------------------------------------------------------
// methods: check that everybody is fine
template<int dim_>
void PoissonProblem<dim_>::how_are_you_doin() const {
  out << "regardless copyright infringements, this is the PoissonProblem class instance number " << count << ":" << endl;
  out << "        elements: " << grid->get_num_elements() << " = " << grid->get_num_intervals() << endl;
  out << "         degrees: " << space->get_degree_table()[0] << endl;
  TensorSize<dim_> nbf;
  for (int idim=0; idim<dim_; idim++) nbf[idim]=space->get_num_basis(0,idim);
  out << " basis functions: " << space->get_num_basis() << " = " << nbf << endl;
  out << "   system matrix: " << mat->NumGlobalRows() << " x " << mat->NumGlobalCols() << endl;
}
template<int dim_>
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


// methods: system matrix and right hand side vector assemble
template<int dim_>
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
  //TensorSize<dim> nqn;
  //typename SplineSpace<dim_>::DegreeTable deg = space->get_degree_table();
  //for (int idim=0; idim<dim_; idim++) nqn[idim] = deg[0][idim]+1;
  //auto quad = QGauss<dim_>::create(nqn);
  auto Nqn = quad->get_num_points();
  basis_handler->init_element_cache(basis_el,quad);

  // starting the cache handler for the (constant) function f:
  //typename Function<dim_,0,1,1>::Value f_val {5.0};
  //const auto ff = grid_functions::CustomGridFunction<dim_,1>::const_create(grid);
  auto funct_handler = source_term->create_cache_handler();
  auto funct_el = source_term->begin();
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
  //if (print) print_system();

  using space_tools::project_boundary_values;
  using dof_tools::apply_boundary_values;
  //using RefSpace = ReferenceSpaceBasis<dim>;
  //using Function = Function<dim_,0,1,1>;

  // applying boundary condition
  //typename Function<dim_,0,1,1>::Value g_val {0.0};
  //const auto g = grid_functions::ConstantGridFunction<dim_,1>::const_create(grid,{0.0});
  //auto face_quad = QGauss<dim_-1>::create(nqn[0]);
  //const set<boundary_id> dir_id {0};
  //std::map<Index,Real> values;
  //project_boundary_values<dim_,1>(*g,*basis,face_quad,dir_id,values);
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
  //if (print) print_system();
}

template<int dim_>
void PoissonProblem<dim_>::solve() const {
  auto solver = create_solver(*mat,*sol,*rhs);
  solver->solve();
}

template<int dim_>
Real PoissonProblem<dim_>::l2_error(std::shared_ptr<const GridFunction<dim_,1>> u_ex) const {

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
  for (; grid_el != grid_eld; ++grid_el) {
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
