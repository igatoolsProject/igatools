#include "AztecOO_config.h"
#include "AztecOO.h"

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
    shared_ptr<Basis<dim,0,dim,1>>          basis;
    shared_ptr<const QGauss<dim>>           quad;
      // linear system data
    shared_ptr<Matrix> mat;
    shared_ptr<Vector> rhs;
    shared_ptr<Vector> sol;
    void complete_construction();

  public:
    shared_ptr<Grid<dim>> get_grid() const {
      return grid;
    }
    void assemble(const Real lambda, const Real mu, shared_ptr<const FormulaGridFunction<dim,dim>> source_term) const;
    void solve() const;
    Real l2_error(shared_ptr<const GridFunction<dim,dim>> uex) const;
    void check() const;
    // outputs
    void output() const;
    void output(shared_ptr<const FormulaGridFunction<dim,dim>> exact_solution) const;
    void save_system();
};

template<int dim, int range, int rank=1>
std::set<Index> get_boundary_dofs_comp(shared_ptr<const DofDistribution<dim,range,rank>> dof_distribution,
                                       int face, int comp) {

  // extracting dof distribution
  int first_dof = 0;
  for (int icomp=0; icomp<comp; icomp++) {
    first_dof += dof_distribution->get_num_dofs_comp(icomp);
  }
  // hand-made creation of dof distriubtion of component comp
  std::set<Index> comp_dofs;
  for (int idof=0; idof<dof_distribution->get_num_dofs_comp(comp); idof++) {
    comp_dofs.insert(first_dof + idof);
  }
  // all boundary dofs on the face
  Topology<dim-1> sub_elem_topology;
  auto bdr_dofs = dof_distribution->get_boundary_dofs(face,sub_elem_topology);

  // creating bounadry dofs
  std::set<Index> bdr_dofs_comp;
  std::set_intersection(comp_dofs.begin(), comp_dofs.end(),
                        bdr_dofs.begin(), bdr_dofs.end(),
                        std::insert_iterator<std::set<Index>>(bdr_dofs_comp,bdr_dofs_comp.begin()));

  /*cout << "FACE " << face << " COMP " << comp << endl;
  //auto all_dofs = dof_distribution->get_global_dofs();
  //cout << "All dofs:" << endl;
  //for (auto it=all_dofs.begin(); it!=all_dofs.end(); ++it) {
  //  cout << "  " << *it;
  //}
  //cout << endl;
  cout << "Boundary dofs of face " << face << ":" << endl;
  for (auto it=bdr_dofs.begin(); it!=bdr_dofs.end(); ++it) {
    cout << "  " << *it;
  }
  cout << endl;
  cout << "Component " << comp << " dofs:" << endl;
  for (auto it=comp_dofs.begin(); it!=comp_dofs.end(); ++it) {
    cout << "  " << *it;
  }
  cout << endl;
  cout << "Intersection dofs:" << endl;
  for (auto it=bdr_dofs_comp.begin(); it!=bdr_dofs_comp.end(); ++it) {
    cout << "  " << *it;
  }
  cout << endl;// */

  printf("face %d, comp %d:",face,comp);
  for (auto it=bdr_dofs_comp.begin(); it!=bdr_dofs_comp.end(); ++it) {
    printf("  %2d",*it);
  }
  printf("\n");// */

  return bdr_dofs_comp;
}

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
  // underlying grid for the geometry function
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
  basis = PhysicalSpaceBasis<dim,dim>::create(reference_bspline, domain);

  // COMPLETING THE CONSTRUCTION
  complete_construction();

  // CREATING THE QUADRATURE RULES
  /*TensorSize<dim> num_quad;
  for (int idim=0; idim<dim; idim++) {
    num_quad[idim] = deg[idim]+1;
  }
  quad = QGauss<dim>::const_create(num_quad);

  // CREATING THE LINEAR SYSTEM
  mat   = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
  rhs   = create_vector(mat->RangeMap());
  sol   = create_vector(mat->DomainMap());*/
}

template<int dim>
ElasticityProblem<dim>::ElasticityProblem(const TensorSize<dim> nel,
                                          const TensorIndex<dim> deg) {
  // CREATING THE TRIVIAL PHYSICAL DOMAIN
  // computing the number of knots required
  TensorSize<dim> num_knots;
  for (int idim=0; idim<dim; idim++) {
    num_knots[idim] = nel[idim]+1;
  }
  // underlying grid as physical space
  grid = Grid<dim>::create(num_knots);

  // CREATING THE DISCRETE SPACE IN PARAMETRIC DOMAIN
  // spline space for the basis
  auto space = SplineSpace<dim,dim>::create(deg,grid);
  // B-spline space as basis space
  basis = BSpline<dim,dim>::create(space);
  // trivial domain defined as imaage of the identity vector field
  auto id_funct = grid_functions::IdentityGridFunction<dim>::create(grid);
  domain = Domain<dim>::create(id_funct);

  // COMPLETING THE CONSTRUCTION
  complete_construction();
  
  /*TensorSize<dim> num_quad;
  for (int idim=0; idim<dim; idim++) {
    num_quad[idim] = deg[idim]+1;
  }
  quad = QGauss<dim>::const_create(num_quad);
  mat   = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
  rhs   = create_vector(mat->RangeMap());
  sol   = create_vector(mat->DomainMap());
  auto id_funct = grid_functions::IdentityGridFunction<dim>::create(grid);
  domain = Domain<dim>::create(id_funct);*/
}

template<int dim>
void ElasticityProblem<dim>::complete_construction() {
  // getting the degrees
  auto deg = basis->get_spline_space()->get_degree_table();
  TensorSize<dim> num_quad;
  for (int idim=0; idim<dim; idim++) {
    num_quad[idim] = deg[idim][0]+1;
  }
  quad = QGauss<dim>::const_create(num_quad);
  mat   = create_matrix(*basis,DofProperties::active,Epetra_SerialComm());
  rhs   = create_vector(mat->RangeMap());
  sol   = create_vector(mat->DomainMap());
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
  auto Nqn = quad->get_num_points();
  basis_handler->init_element_cache(basis_el,quad);

  // ATTEMPT 1: GridFunction
  //auto source_term   = grid_functions::ConstantGridFunction<dim,dim>::const_create(grid,{0.0,0.0,-0.0001});
  auto funct_handler = source_term->create_cache_handler();
  auto funct_el      = source_term->begin();
  funct_handler->set_element_flags(grid_function_element::Flags::D0);
  funct_handler->init_cache(funct_el,quad);

  // ATTEMPT 2: Function
  /*auto source_term = functions::ConstantFunction<dim,0,dim>::const_create(domain,{0.0,0.0,-1.0});
  auto funct_handler = source_term->create_cache_handler();
  auto funct_el = source_term->begin();
  funct_handler->set_element_flags(function_element::Flags::D0);
  funct_handler->init_cache(funct_el,quad);//*/

  // retrieving the last datum and then starting the loop
  const auto l = lambda;
  const auto m = 2.0 * mu;
  for (int iel=0; basis_el!=basis_eld; iel++, ++basis_el, ++funct_el) {
    basis_handler->fill_element_cache(basis_el);
    funct_handler->fill_element_cache(funct_el);
    //funct_handler->template fill_cache<dim>(funct_el,0);

    // preparing some stuff: creating the local matrices
    auto Nbf = basis_el->get_num_basis(DofProperties::active);
    DenseMatrix loc_mat(Nbf,Nbf); loc_mat=0.0;
    //DenseMatrix loc_mat_bis(Nbf,Nbf); loc_mat_bis=0.0;
    DenseVector loc_rhs(Nbf);     loc_rhs=0.0;

    // gathering the requested data
    auto values = basis_el->get_element_values();
    auto grads  = basis_el->get_element_gradients();
    auto divs   = basis_el->get_element_divergences();
    auto w_meas = basis_el->get_element_w_measures();
    auto fval   = funct_el->get_element_values_D0();

    //out << "Element " << iel << endl;
    //for (int iqn=0; iqn<Nqn; iqn++) {
    //  out << fval[iqn] << endl;
    //}

    // precomputing epsilon(v) = 0.5 * (\grad(v) + \grad(v)^T)
    /*using Der = typename SpaceElement<dim,0,dim,1>::template Derivative<1>;
    ValueTable<Der> epsils(Nbf,Nqn);
    for (int ibf=0; ibf<Nbf; ibf++) {
      auto epsil = epsils.get_function_view(ibf);
      auto grad  = grads.get_function_view(ibf);
      for (int iqn=0; iqn<Nqn; iqn++) {
        //epsil[iqn] = grad[iqn] + transpose(co_tensor(grad[iqn]));
        epsil[iqn] = symmetric_tensor(grad[iqn]);
      }
    }// */

    // the basis function loop
    for (int ibf=0; ibf<Nbf; ibf++) { 
      // loop for the stiffness local matrix
      const auto &div_i = divs.get_function_view(ibf); // view for the i-th basis function divergence
      const auto &grd_i = grads.get_function_view(ibf); // view for the i-th basis function divergence
      //const auto &eps_i = epsils.get_function_view(ibf); 
      for (int jbf=0; jbf<Nbf; jbf++) {
        const auto &div_j = divs.get_function_view(jbf); // view for the j-th basis function divergence
        const auto &grd_j = grads.get_function_view(jbf); // view for the i-th basis function divergence
        //const auto &eps_j = epsils.get_function_view(jbf);
        double part_1 = 0.0, part_2 = 0.0;
        for (int iqn=0; iqn<Nqn; iqn++) {
          // PART 1: assembling  lambda \int div(v_i)*div(v_j)
          part_1 += scalar_product(div_i[iqn], div_j[iqn]) * w_meas[iqn];
          // PART 2: assembling  2mu \int eps(v_i):eps(v_j)
          //test += scalar_product(eps_i[iqn], eps_j[iqn]) * w_meas[iqn];
          part_2 += scalar_product(grd_i[iqn], grd_j[iqn]) * w_meas[iqn];
          //part_2 += scalar_product(grd_i[iqn],symmetric_tensor(grd_j[iqn])) * w_meas[iqn];
        }
        loc_mat(ibf,jbf) = l * part_1 + m * part_2;
        //printf("eps = %1.5f\tgrad = %1.5f\tepsbis = %1.5f\n",part_2,test,test2);
      }
      // loop for the right hand side local vector
      //const auto &val_i = values.get_function_view(ibf); // view for the i-th basis function
      //for (int iqn=0; iqn<Nqn; iqn++) {
        // PART 3: assembling \int v_i*1
        //loc_rhs(ibf) += scalar_product(val_i[iqn][0],fval[iqn]) * w_meas[iqn];
      //}
    }// */

    /*for (int ibf=0; ibf<Nbf; ibf++) {
      for (int jbf=0; jbf<Nbf; jbf++) {
        if (fabs(loc_mat(ibf,jbf)-loc_mat_bis(ibf,jbf))>1e-10) {
          cout << "houston, we have a problem with basis function ";
          printf("%2d,%2d:  %1.3e != %1.3e\n",ibf,jbf,loc_mat(ibf,jbf),loc_mat_bis(ibf,jbf));
	      }
      }
    }// */

    //loc_mat += m * basis_el->integrate_gradu_gradv();
    loc_rhs = basis_el->integrate_u_func(fval);

    // spatashing element matrix into the global matrix
    const auto loc_dofs = basis_el->get_local_to_global(DofProperties::active);
    mat->add_block(loc_dofs,loc_dofs,loc_mat);
    rhs->add_block(loc_dofs,loc_rhs);
  }
  mat->FillComplete();//*/

  // applying the boundary conditions
  using space_tools::project_boundary_values;
  using dof_tools::apply_boundary_values;
  auto dof_distribution = basis->get_dof_distribution();

  /*std::map<Index,Real> bdr_vals_1;
  auto bdr_dofs_1 = get_boundary_dofs_comp<dim,dim>(dof_distribution,0,0); 
  for (set<Index>::iterator ival=bdr_dofs_1.begin(); ival!=bdr_dofs_1.end(); ival++) {
      bdr_vals_1[*ival] = 0.0;
  }
  apply_boundary_values(bdr_vals_1,*mat,*rhs,*sol);

  std::map<Index,Real> bdr_vals_2;
  auto bdr_dofs_2 = get_boundary_dofs_comp<dim,dim>(dof_distribution,2,1); 
  for (set<Index>::iterator ival=bdr_dofs_2.begin(); ival!=bdr_dofs_2.end(); ival++) {
      bdr_vals_2[*ival] = 0.0;
  }
  apply_boundary_values(bdr_vals_2,*mat,*rhs,*sol);

  std::map<Index,Real> bdr_vals_3;
  auto bdr_dofs_3 = get_boundary_dofs_comp<dim,dim>(dof_distribution,4,2); 
  for (set<Index>::iterator ival=bdr_dofs_3.begin(); ival!=bdr_dofs_3.end(); ival++) {
      bdr_vals_3[*ival] = 0.0;
  }
  apply_boundary_values(bdr_vals_3,*mat,*rhs,*sol);//*/


  std::map<Index,Real> bdr_vals;
  Topology<dim-1> sub_elem_topology;
  auto bdr_dofs = dof_distribution->get_boundary_dofs(0,sub_elem_topology);
  for (set<Index>::iterator it=bdr_dofs.begin(); it!=bdr_dofs.end(); it++) {
    bdr_vals[*it]=0.0;
    //cout << "changed boundary value " << *it << endl;
  }
  apply_boundary_values(bdr_vals,*mat,*rhs,*sol);//*/

/*  auto dof_distribution = basis->get_dof_distribution(); 
  Topology<dim-1> sub_elem_topology;
  for (set<Index>::iterator iface=faces.begin(); iface!=faces.end(); iface++) {
    auto bdr_dofs = dof_distribution->get_boundary_dofs(*iface,sub_elem_topology);
    for (set<Index>::iterator ival=bdr_dofs.begin(); ival!=bdr_dofs.end(); ival++) {
      bdr_vals[*ival] = 0.0;
     // cout << "changed boundary value " << *ival << endl;
    }
  //bdr_dofs = dof_distribution->get_boundary_dofs(1,sub_elem_topology);
  //for (set<Index>::iterator it=bdr_dofs.begin(); it!=bdr_dofs.end(); it++) {
    //bdr_vals[*it]=0.0;
    //cout << "changed boundary value " << *it << endl;
  //}
  }
  apply_boundary_values(bdr_vals,*mat,*rhs,*sol);//*/

}

template<int dim> // solver for the linear system
void ElasticityProblem<dim>::solve() const {
  //auto solver = create_solver(*mat,*sol,*rhs);
  //solver->solve();

  // setting up the problem
  Epetra_LinearProblem problem(&*mat,&*sol,&*rhs);
  AztecOO solver(problem);
  // setting up the solver
  solver.SetAztecOption(AZ_solver,  AZ_cg);
  solver.SetAztecOption(AZ_precond, AZ_none);
  solver.SetAztecOption(AZ_output,  AZ_none);
  // solve, for god's sake! SOLVE!
  solver.Iterate(1000, 1.0E-7);
  //sol->print_info(out);
  //out << endl;
}

template<int dim> // compute the L2 error given the exact solution
Real ElasticityProblem<dim>::l2_error(std::shared_ptr<const GridFunction<dim,dim>> u_ex) const {
//                                 Real &l2_err, bool compute_h1, Real &h1_err = 0.0) const {

  // creating the discrete solution
  auto u_ds = IgFunction<dim,0,dim,1>::const_create(basis,*sol);
  
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
  for (int iel=0; grid_el != grid_eld; ++grid_el, ++u_ex_el, ++u_ds_el, ++iel) {
    u_ds_handler->fill_element_cache(u_ds_el);
    u_ex_handler->fill_element_cache(u_ex_el);
    grid_handler->fill_element_cache(grid_el);

    // retrieving data
    auto uh = u_ds_el->get_element_values_D0();
    auto ue = u_ex_el->get_element_values_D0();
    auto qw = grid_el->get_element_weights();

    // loopin' on the quadrature points
    //cout << " ELEMENT " << iel << endl;
    for (int iqn=0; iqn<nqn; iqn++) {
      //auto diff = uh[iqn] - ue[iqn];
      uh[iqn] -= ue[iqn];
      //cout << "difference: " << diff.norm_square() << endl;
      err += uh[iqn].norm_square() * qw[iqn];
    }
  }

  return sqrt(err);
}// */


/*template<int dim>
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
}// */

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
}// */

template<int dim>
void ElasticityProblem<dim>::save_system() {

  FILE *matfp, *solfp, *rhsfp;
  auto nel = grid->get_num_intervals()[0];
  auto deg = basis->get_spline_space()->get_degree_table()[0][0];

  string path    = "../../../workspace/results/";
  string discr = to_string(dim) + "d_" + to_string(nel) + "_" + to_string(deg) + ".m";
  string matname = path + "mat_" + discr;
  string solname = path + "sol_" + discr;
  string rhsname = path + "rhs_" + discr;
  matfp = fopen(matname.c_str(),"w");
  solfp = fopen(solname.c_str(),"w");
  rhsfp = fopen(rhsname.c_str(),"w");

  fprintf(solfp,"sol = [\n");
  fprintf(rhsfp,"rhs = [\n");
  for (int ien=0; ien<sol->MyLength(); ien++) {
    fprintf(solfp,"%1.16f\n",(*sol)[ien]);
    fprintf(rhsfp,"%1.16f\n",(*rhs)[ien]);
  }
  fprintf(solfp,"];\n");
  fprintf(rhsfp,"];\n");

  //fprintf(matfp,"zzz = zeros(%d,3);\n",mat->NumGlobalNonzeros());
  fprintf(matfp,"zzz = [\n");
  for (int irow=0; irow<mat->NumMyRows(); irow++) {
    int *ind, nnz;
    double *val;
    mat->ExtractMyRowView(irow,nnz,val,ind);
    for (int ien=0; ien<nnz; ien++) {
      //if (val[ien]!=0.0)
        fprintf(matfp,"%d %d %1.64f\n",irow+1,ind[ien]+1,val[ien]);
    }
  }
  fprintf(matfp,"];\n");
  fprintf(matfp,"mat = spconvert(zzz);\n");
  fprintf(matfp,"clear zzz;\n");

  fclose(matfp);
  fclose(solfp);
  fclose(rhsfp);
}
