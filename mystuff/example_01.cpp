// headers for spline spaces
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/physical_space_basis.h>
// headers for quadrature
#include <igatools/base/quadrature_lib.h>
// headers for linear algebra objects
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
// headers for system solving
//#include <igatools/linear_algebra/epetra_solver.h>
#include <igatools/linear_algebra/dof_tools.h>
// headers for function representations
#include <igatools/functions/function_lib.h>
#include <igatools/functions/function_element.h>
#include <igatools/geometry/grid_function_lib.h>
// headers for output
#include <igatools/io/writer.h>
#include <igatools/base/logstream.h>
// finally! my stuff!
//#include "my_formula_grid_function.h"
#include "custom_grid_function.h"
#include "custom_function.h"
#include "geometry.h"

// headers for Trilinos stuff
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>

using namespace iga;
using namespace std;
using namespace EpetraTools;
LogStream out;

const Real lambda = 0.57692;
const Real mu     = 0.38462;
#include "elasticity_problem.h"

// ----------------------------------------------------------------------------
//   MY CUSTOM FUNCTION
// ----------------------------------------------------------------------------
#define PI 3.141592653589793
template<int dim>
Values<dim,dim,1> u(Points<dim> pts)
{
  Values<dim,dim,1> x;
  for (int idim=0; idim<dim-1; idim++) x[idim] = 0.0;
  x[dim-1] = 0.1 * (cos(2.0*PI*pts[0]) - 1.0);
  //for (int idim=0; idim<dim; idim++) {
  //  x *= sin( 4.0 * pts[idim] * PI);
  //}
  return x;
}
template<int dim>
Values<dim,dim,1> f(Points<dim> pts)
{
  Values<dim,dim,1> x;
  for (int idim=0; idim<dim-1; idim++) x[idim] = 0.0;
  x[dim-1] = mu*PI*PI*0.8 * cos(2.0*PI*pts[0]);// * cos(2.0*PI*pts[1]);
  //for (int idim=0; idim<dim; idim++) {
  //  x *= sin(4.0 * pts[idim] * PI);
  //}
  return x;
}

Real eoc(Real e1, Real e2, Real h1, Real h2)
{
  return (log(e2)-log(e1))/(log(h2)-log(h1));
}

// ----------------------------------------------------------------------------
//   MAIN
// ----------------------------------------------------------------------------
int main()
{

  // problem dimension
  const int dim = 2;

  // geometry definition
  Geometry<dim> geometry;
  if (dim==2) geometry.load("../../../workspace/geometries/ring.nurbs");
  if (dim==3) geometry.load("../../../workspace/geometries/hose.nurbs");

  /*geometry.nel   = {1,1};
  geometry.deg   = {1,2};
  // weights
  geometry.weights[0] = 1.0;
  geometry.weights[1] = 1.0;
  geometry.weights[2] = sqrt(2.0)/2.0;
  geometry.weights[3] = sqrt(2.0)/2.0;
  geometry.weights[4] = 1.0;
  geometry.weights[5] = 1.0;
  // control points
  geometry.coefs[ 0] = 1.0;  geometry.coefs[ 6] = 0.0;
  geometry.coefs[ 1] = 2.0;  geometry.coefs[ 7] = 0.0;
  geometry.coefs[ 2] = 1.0;  geometry.coefs[ 8] = 1.0;
  geometry.coefs[ 3] = 2.0;  geometry.coefs[ 9] = 2.0;
  geometry.coefs[ 4] = 0.0;  geometry.coefs[10] = 1.0;
  geometry.coefs[ 5] = 0.0;  geometry.coefs[11] = 2.0;// */

  // linear elasticity problem creation
  /*std::set<int> nels = {2,4,6,8,10,12,14,16};
  std::set<int> degs = {1,2,3};
  //std::set<int> nels = {4};
  //std::set<int> degs = {2};
  std::map<int,Real> e1, e2;
  double h1, h2;
  printf("\n\n                     deg 1");
  printf("                           deg 2");
  printf("                           deg 3\n");
  //printf("                           deg 4");
  //printf("                           deg 5\n");
  for (set<int>::iterator iel=nels.begin(); iel!=nels.end(); iel++) {
    h2=1.0/(*iel);
    printf("nel %2d:   ",*iel);
    for (set<int>::iterator ideg=degs.begin(); ideg!=degs.end(); ideg++) {
      // assigning discretization setup
      TensorSize<dim>  nel;
      TensorIndex<dim> deg;
      for (int idim=0; idim<dim; idim++) {
        nel[idim]=(*iel);
        deg[idim]=(*ideg);
      }
      // creating the problem
      auto problem =  ElasticityProblem<dim>(nel,deg);
      auto grid = problem.get_grid();
      using namespace grid_functions;
      auto source_term = CustomGridFunction<dim,dim>::const_create(grid,f);
      problem.assemble(lambda,mu,source_term);
      problem.solve();
      auto exact_solution = CustomGridFunction<dim,dim>::const_create(grid,u);
      e2[*ideg] = problem.l2_error(exact_solution);
      printf("\terr = %1.2e\teoc = %1.2f",e2[*ideg],eoc(e1[*ideg],e2[*ideg],h1,h2));
      //problem.output(exact_solution);
      e1[*ideg] = e2[*ideg];
    }
    h1 = h2;
    printf("\n");
  }// */

  TensorSize<dim>  nel;
  TensorIndex<dim> deg;
  for (int idim=0; idim<dim; idim++)
  {
    nel[idim]=2;
    deg[idim]=1;
  }
  auto problem = ElasticityProblem<dim>(nel,deg);
  using namespace grid_functions;
  //auto source_term = ConstantGridFunction<dim,dim>::const_create(problem.get_grid(),{1.0,0.0});
  //auto source_term = ConstantGridFunction<dim,dim>::const_create(problem.get_grid(),{1.0,0.0,0.0});
  auto source_term = CustomGridFunction<dim,dim>::const_create(problem.get_grid(),f);
  problem.assemble(lambda,mu,source_term);
  problem.solve();
  //problem.output();// */



  return 0;
}


