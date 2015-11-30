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
  void load(const char* fname) {
    FILE *fp;
    int check;
    fp=fopen(fname,"r");
    if (fp==0) cout << "cannot open the .nurbs file!" <<endl;
    else {
      // checking dimension
      fscanf(fp,"%d",&check);
      if (check!=dim) cout << "wrong geometry dimension!" << endl;
      else {
        // reading degrees
        for (int idim=0; idim<dim; idim++)
          fscanf(fp,"%d",&deg[idim]);
        // reading elements
        for (int idim=0; idim<dim; idim++)
          fscanf(fp,"%d",&nel[idim]);
        // reading number of control points
        int ncp;
        fscanf(fp,"%d",&ncp);
        // reading control points
        double data;
        for (int icp=0; icp<ncp*dim; icp++) {
          fscanf(fp,"%lf",&data);
          coefs[icp]=data;
        }
        // reading weights
        for (int icp=0; icp<ncp; icp++) {
          fscanf(fp,"%lf",&data);
          weights[icp]=data;
        }
      }
    }
  }
};

template<int dim>
class ElasticityProblem {

    // constructors
  public:
    ElasticityProblem() = delete;
    ElasticityProblem(const TensorSize<dim> nel, 
                      const TensorIndex<dim> deg,
                      const Geometry<dim> &geom);
    ElasticityProblem(const TensorSize<dim> nel,
                      const TensorIndex<dim> deg);
  private:
      // space data
    shared_ptr<Grid<dim>>                   grid;
    shared_ptr<Domain<dim>>                 domain;
    shared_ptr<Domain<dim>>                 deformed_domain;
    shared_ptr<BSpline<dim,dim,1>>          basis;
    shared_ptr<const QGauss<dim>>           elem_quad;
    shared_ptr<const QGauss<dim-1>>         face_quad;
      // linear system data
    shared_ptr<Matrix> mat;
    shared_ptr<Vector> rhs;
    shared_ptr<Vector> sol;
    bool is_grid;

  public:
    shared_ptr<Grid<dim>> get_grid() const {
      return grid;
    }
    void assemble(const Real lambda, const Real mu, shared_ptr<const FormulaGridFunction<dim,dim>> source_term) const;
    void solve() const;
    void check() const;
    void deform() const;
    void output() const;
    void output(shared_ptr<const FormulaGridFunction<dim,dim>> exact_solution) const;
};

// ----------------------------------------------------------------------------
//   CONSTRUCTOR
// ----------------------------------------------------------------------------
template<int dim>
ElasticityProblem<dim>::ElasticityProblem(const TensorSize<dim> nel,
                                          const TensorIndex<dim> deg,
                                          const Geometry<dim> &geom) {
  // CREATING THE PHYSICAL DOMAIN
  // computing the number of knots required
  TensorSize<dim> num_knots;
  for (int idim=0; idim<dim; idim++) {
    num_knots[idim] = geom.nel[idim]+1;
  }
  // base underlying grid for the geometry function
  grid = Grid<dim>::create(num_knots);
  // B-spline vector field for the geometry function
  auto vect_space   = SplineSpace<dim,dim>::create(geom.deg, grid);
  auto vect_bspline = BSpline<dim,dim>::create(vect_space);
  // B-spline scalar field for the weight function
  auto scal_space   = SplineSpace<dim,1>::create(geom.deg, grid);
  auto scal_bspline = BSpline<dim,1>::create(scal_space);
  // the weight function
  auto weight_funct = IgGridFunction<dim,1>::create(scal_bspline, geom.weights);
  // NURBS vector field for the geometry function
  auto vect_nurbs   = NURBS<dim,dim>::create(vect_bspline, weight_funct);
  // geometry function
  auto geom_funct   = IgGridFunction<dim,dim>::create(vect_nurbs, geom.coefs);
  // domain Omega
  domain = Domain<dim>::create(geom_funct);
  // h-refining everything
  SafeSTLArray<bool,dim> directions;
  SafeSTLArray<Size,dim> subdivisions;
  for (int idim=0; idim<dim; idim++) {
    directions[idim] = true;
    subdivisions[idim] = nel[idim];
  }
  grid->refine_directions(directions, subdivisions);

  // CREATING THE DISCRETE SPACE IN PHYSICAL DOMAIN
  // creating the reference space with custom degrees
  auto reference_space   = SplineSpace<dim,dim>::create(deg, grid);
  auto reference_bspline = BSpline<dim,dim>::create(reference_space);
  // at last, the B-spline basis functions space in physical domain
  //using grid_functions::IdentityGridFunction<dim>;
  //auto trivial_domain = Domain<dim>::create(grid_functions::IdentityGridFunction<dim>::create(grid));
  basis = PhysicalSpaceBasis<dim,dim>::create(reference_bspline, domain);

  // CREATING THE QUADRATURE RULES
  TensorSize<dim> elem_num_quad;
  for (int idim=0; idim<dim; idim++) {
    elem_num_quad[idim] = deg[idim]+1;
  }
  elem_quad = QGauss<dim>::const_create(elem_num_quad);
  TensorSize<dim-1> face_num_quad;
  for (int idim=0; idim<dim-1; idim++) {
    face_num_quad[idim] = deg[idim]+1;
  }
  face_quad = QGauss<dim-1>::const_create(face_num_quad);

  // the basis functions in the physical domain
  //basis    = PhysicalSpaceBasis<dim>::create(scal_bspline,domain);
  // refine everything
  //grid->print_info(out);
  //out << endl;
  // qudrature rule, linear system 
  //quad  = QGauss<dim>::const_create(nqn);
  // CREATING THE LINEAR SYSTEM
  mat   = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
  rhs   = create_vector(mat->RangeMap());
  sol   = create_vector(mat->DomainMap());
  is_grid=false;
}

template<int dim>
ElasticityProblem<dim>::ElasticityProblem(const TensorSize<dim> nel,
                                          const TensorIndex<dim> deg) {
  TensorSize<dim> num_knots;
  for (int idim=0; idim<dim; idim++) {
    num_knots[idim] = nel[idim]+1;
  }
  grid = Grid<dim>::create(num_knots);
  auto space = SplineSpace<dim,dim>::create(deg,grid);
  basis = BSpline<dim,dim>::create(space);
  TensorSize<dim> elem_num_quad;
  for (int idim=0; idim<dim; idim++) {
    elem_num_quad[idim] = deg[idim]+1;
  }
  elem_quad = QGauss<dim>::const_create(elem_num_quad);
  is_grid=true;
  mat   = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
  rhs   = create_vector(mat->RangeMap());
  sol   = create_vector(mat->DomainMap());
  auto id_funct = grid_functions::IdentityGridFunction<dim>::create(grid);
  domain = Domain<dim>::create(id_funct);
}

// ----------------------------------------------------------------------------
//   METHODS
// ----------------------------------------------------------------------------
template<int dim> // assemble the system
void ElasticityProblem<dim>::assemble(const Real lambdaa,
                                      const Real muu,
                                      shared_ptr<const FormulaGridFunction<dim,dim>> source_term) const {

  //out << "Number of elements: " << grid->get_num_intervals() << endl;
  //auto d = basis->get_spline_space()->get_degree_table();
  //out << "Degrees: " << d[0] << " " << d[1] << endl;
  //out << "Number of basis functions: " << basis->get_num_basis() << endl;

  // starting the cache handler for the basis functions:
  auto basis_handler = basis->create_cache_handler();
  auto basis_el      = basis->begin();
  const auto basis_eld = basis->end();
  // setting the flags
  using Flags = space_element::Flags;
  auto flag = Flags::value | Flags::gradient | Flags::w_measure | Flags::divergence;
  basis_handler->set_element_flags(flag);
  // setting the quarature rule
  auto Nqn = elem_quad->get_num_points();
  basis_handler->init_element_cache(basis_el,elem_quad);

  // ATTEMPT 1: GridFunction
  //auto source_term   = grid_functions::ConstantGridFunction<dim,dim>::const_create(grid,{0.0,0.0,-0.0001});
  auto funct_handler = source_term->create_cache_handler();
  auto funct_el      = source_term->begin();
  funct_handler->set_element_flags(grid_function_element::Flags::D0);
  funct_handler->init_cache(funct_el,elem_quad);

  // ATTEMPT 2: Function
  /*auto source_term = functions::ConstantFunction<dim,0,dim>::const_create(domain,{0.0,0.0,-1.0});
  auto funct_handler = source_term->create_cache_handler();
  auto funct_el = source_term->begin();
  funct_handler->set_element_flags(function_element::Flags::D0);
  funct_handler->init_cache(funct_el,elem_quad);//*/

  // retrieving the last datum and then starting the loop
  const auto m = 2.0 * mu;
  for (int iel=0; basis_el!=basis_eld; iel++, ++basis_el, ++funct_el) {
    basis_handler->fill_element_cache(basis_el);
    funct_handler->fill_element_cache(funct_el);
    //funct_handler->template fill_cache<dim>(funct_el,0);

    // preparing some stuff: creating the local matrices
    auto Nbf = basis_el->get_num_basis(DofProperties::active);
    DenseMatrix loc_mat(Nbf,Nbf); loc_mat=0.0;
    DenseVector loc_rhs(Nbf);     loc_rhs=0.0;

    // gathering the requested data
    auto values = basis_el->get_element_values();
    auto grads  = basis_el->get_element_gradients();
    auto divs   = basis_el->get_element_divergences();
    auto w_meas = basis_el->get_element_w_measures();
    auto &fval  = funct_el->get_element_values_D0();

    //out << "Element " << iel << endl;
    //for (int iqn=0; iqn<Nqn; iqn++) {
    //  out << fval[iqn] << endl;
    //}

    // precomputing epsilon(v) = 0.5 * (\grad(v) + \grad(v)^T)
    using Der = typename SpaceElement<dim,0,dim,1>::template Derivative<1>;
    ValueTable<Der> epsils(Nbf,Nqn);
    for (int ibf=0; ibf<Nbf; ibf++) {
      auto epsil = epsils.get_function_view(ibf);
      auto grad  = grads.get_function_view(ibf);
      for (int iqn=0; iqn<Nqn; iqn++) {
        //epsil[iqn] = grad[iqn] + transpose(co_tensor(grad[iqn]));
        epsil[iqn] = symmetric_tensor(grad[iqn]);
      }
    }
    
    // the basis function loop
    for (int ibf=0; ibf<Nbf; ibf++) { 
      // loop for the stiffness local matrix
      const auto &div_i = divs.get_function_view(ibf); // view for the i-th basis function divergence
      const auto &eps_i = epsils.get_function_view(ibf); 
      for (int jbf=0; jbf<Nbf; jbf++) {
        const auto &div_j = divs.get_function_view(jbf); // view for the j-th basis function divergence
        const auto &eps_j = epsils.get_function_view(jbf);
        double part_1 = 0.0, part_2 = 0.0;
        for (int iqn=0; iqn<Nqn; iqn++) {
          // PART 1: assembling  lambda \int div(v_i)*div(v_j)
          //part_1 += scalar_product(div_i[iqn], div_j[iqn]) * w_meas[iqn];
          // PART 2: assembling  2mu \int eps(v_i):eps(v_j)
          part_2 += scalar_product(eps_i[iqn], eps_j[iqn]) * w_meas[iqn];
        }
        //loc_mat(ibf,jbf) = lambda*part_1 + mu*part_2;
        loc_mat(ibf,jbf) = m * part_2;
      }
      // loop for the right hand side local vector
      //const auto &val_i = values.get_function_view(ibf); // view for the i-th basis function
      //for (int iqn=0; iqn<Nqn; iqn++) {
        // PART 3: assembling \int v_i*1
        //loc_rhs(ibf) += scalar_product(val_i[iqn][0],fval[iqn]) * w_meas[iqn];
      //}
    }
    auto loc_mat_bis = m * basis_el->integrate_gradu_gradv();
    
    /*for (int ibf=0; ibf<Nbf; ibf++) {
      for (int jbf=0; jbf<Nbf; jbf++) {
        if (fabs(loc_mat(ibf,jbf)-loc_mat_bis(ibf,jbf))>1e-10) {
          cout << "houston, we have a problem with basis function ";
          printf("%2d,%2d:  %1.3e != %1.3e\n",ibf,jbf,loc_mat(ibf,jbf),loc_mat_bis(ibf,jbf));
	      }
      }
    }// */

    loc_rhs = basis_el->integrate_u_func(fval);

    // spatashing element matrix into the global matrix
    const auto loc_dofs = basis_el->get_local_to_global(DofProperties::active);
    mat->add_block(loc_dofs,loc_dofs,loc_mat_bis);
    rhs->add_block(loc_dofs,loc_rhs);
  }
  mat->FillComplete();//*/

  // applying the boundary conditions
  using space_tools::project_boundary_values;
  using dof_tools::apply_boundary_values;


    // ATTEMPT 1
  std::map<Index,Real> bdr_vals;
  std::set<Index> faces = {0};
  auto dof_distribution = basis->get_dof_distribution(); 
  Topology<dim-1> sub_elem_topology;
  for (set<Index>::iterator iface=faces.begin(); iface!=faces.end(); iface++) {
    auto bdr_dofs = dof_distribution->get_boundary_dofs(*iface,sub_elem_topology);
    for (set<Index>::iterator ival=bdr_dofs.begin(); ival!=bdr_dofs.end(); ival++) {
      bdr_vals[*ival]=0.0;
     // cout << "changed boundary value " << *ival << endl;
    }
  //bdr_dofs = dof_distribution->get_boundary_dofs(1,sub_elem_topology);
  //for (set<Index>::iterator it=bdr_dofs.begin(); it!=bdr_dofs.end(); it++) {
    //bdr_vals[*it]=0.0;
    //cout << "changed boundary value " << *it << endl;
  //}
  }
  apply_boundary_values(bdr_vals,*mat,*rhs,*sol);//*/

    // ATTEMPT 2
  /*using SubGridElemMap = typename Grid<dim>::template SubGridMap<dim-1>;
  using FuncAtBndry = Function<dim-1,1,dim,1>;
  SafeSTLMap<int,std::shared_ptr<const FuncAtBndry>> bndry_funcs;
  for (int face_id=0; face_id<1; ++face_id) {
    SubGridElemMap sub_grid_elem_map;
    const auto sub_grid = grid->template get_sub_grid<dim-1>(face_id,sub_grid_elem_map);
    //const auto sub_grid = domain->get_grid_function()->get_grid()->template get_sub_grid<dim-1>(face_id,sub_grid_elem_map);
    const auto sub_domain = domain->template get_sub_domain<dim-1>(face_id,sub_grid_elem_map,sub_grid);
    bndry_funcs[face_id] = functions::ConstantFunction<dim-1,1,dim,1>::const_create(sub_domain, {0.,0.});
  }
  SafeSTLMap<Index, Real> bndry_values;
  project_boundary_values(bndry_funcs,*basis,face_quad,bndry_values);
  apply_boundary_values(bndry_values, *mat, *rhs, *sol);//*/

  //rhs->print_info(out);
  //out << rhs << endl;
}

template<int dim>
void ElasticityProblem<dim>::check() const {

  // dirichlet face
  for (int iface=0; iface<6; iface++) {
    cout << endl << "  FACE " << iface << endl;
    auto dof_distribution = basis->get_dof_distribution();
    Topology<dim-1> sub_elem_topology;
    auto bdr_dofs = dof_distribution->get_boundary_dofs(iface,sub_elem_topology);   
    set<Index>::iterator it=bdr_dofs.begin();
    for (int idim=0; idim<dim; idim++) {
      for (int irow=0; irow<6; irow++) {
        for (int icol=0; icol<6; icol++) {
          if ((*sol)[*it]<0)
            printf(" %1.6f",(*sol)[*it]);
          else
            printf("  %1.6f",(*sol)[*it]);
          ++it;
        }
        printf("\n");
      }
      printf("\n");
    }
    cout << endl;
  }

  // other face (hihihi)
  //cout << "other face" << endl;
  //dof_distribution = basis->get_dof_distribution();
  //bdr_dofs = dof_distribution->get_boundary_dofs(1,sub_elem_topology);   
  //for (set<Index>::iterator it=bdr_dofs.begin(); it!=bdr_dofs.end(); it++) {
  //  out << (*sol)[*it];
  //}


}


template<int dim> // solver for the linear system
void ElasticityProblem<dim>::solve() const {
  //auto solver = create_solver(*mat,*sol,*rhs);
  //solver->solve();

  // setting up the problem
  Epetra_LinearProblem problem(&*mat,&*sol,&*rhs);
  AztecOO solver(problem);
  // setting up the solver
  solver.SetAztecOption(AZ_solver,  AZ_gmres);
  solver.SetAztecOption(AZ_precond, AZ_none);
  solver.SetAztecOption(AZ_output,  AZ_none);
  // solve, for god's sake! SOLVE!
  solver.Iterate(1000, 1.0E-7);
}

/*void ElasticityProblem<dim>::deform() const {

  auto geom_funct = domain->get_grid_function();
  auto coefs      = geom_funct->get_coefficients();
  auto vect_nurbs = geom_funct->get_basis();
  IgCoefficiencts deformed_coefs;
  for (int ien=0; ien<sol->MyLength(); ien++) {
    deformed_coefs[ien] = coefs[ien] + sol[ien];
  }  
}//*/


/*using OP = Epetra_Operator;
using MV = Epetra_MultiVector;
using SolverPtr = Teuchos::RCP<Belos::SolverManager<double, MV, OP> >;
template<int dim_> // custom siuppacool solver
void ElasticityProblem<dim_>::custom_solve(int &it1, double &cond1, int &it2, double &cond2) const {

  //double cond1;

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
  cond2 = cond_estimator.getConditionNumber();

}

template<int dim_> // compute the L2 error given the exact solution
Real ElasticityProblem<dim_>::l2_error(std::shared_ptr<const GridFunction<dim_,1>> u_ex) const {
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
}*/


template<int dim>
void ElasticityProblem<dim>::output() const
{
  const int num_plot_pts = 10;
  //auto domain = basis->get_physical_domain();
  Writer<dim> writer(domain, num_plot_pts);
  //using IgFunc = IgFunction<dim,0,1,1>;
  if (is_grid) {
    auto solution = IgGridFunction<dim,dim>::create(basis, *sol);
    writer.template add_field<dim>(*solution, "solution");
  }
  else {
    //auto solution_function = IgFunction<dim,0,dim,1>::create(basis, *sol);
    //writer.template add_field<dim,1>(*solution_function, "solution");
  }
  string filename = "plot_" + to_string(dim) + "d" ;
  writer.save(filename);
}

template<int dim>
void ElasticityProblem<dim>::output(shared_ptr<const FormulaGridFunction<dim,dim>> exact_solution) const
{
  const int num_plot_pts = 10;
  //auto domain = basis->get_physical_domain();
  Writer<dim> writer(domain, num_plot_pts);
  //using IgFunc = IgFunction<dim,0,1,1>;
  if (is_grid) {
    auto solution = IgGridFunction<dim,dim>::create(basis, *sol);
    writer.template add_field<dim>(*solution, "solution");
    writer.template add_field<dim>(*exact_solution,"exact solution");
  }
  else {
    //auto solution_function = IgFunction<dim,0,dim,1>::create(basis, *sol);
    //writer.template add_field<dim,1>(*solution_function, "solution");
  }
  string filename = "plot_" + to_string(dim) + "d" ;
  writer.save(filename);
}
