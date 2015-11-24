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

// headers for Trilinos stuff
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>

using namespace iga;
using namespace std;
using namespace EpetraTools;
LogStream out;

#include "poisson_problem.h"

// ----------------------------------------------------------------------------
//   MY CUSTOM FUNCTION
// ----------------------------------------------------------------------------
#define PI 3.141592653589793
template<int dim> // exact solution
Values<dim,1,1> exact_solution(Points<dim> pts) {
  Values<dim,1,1> x = {1.0};
  for (int idim=0; idim<dim; idim++) {
    x *= sin( 4.0 * pts[idim] * PI);
  }
  return x;
}

template<int dim> // source term
Values<dim,1,1> source_term(Points<dim> pts) {
  Values<dim,1,1> x = {16.0 * dim * PI * PI};
  for (int idim=0; idim<dim; idim++) {
    x *= sin(4.0 * pts[idim] * PI);
  }
  return x;
}

// ----------------------------------------------------------------------------
//   MAIN
// ----------------------------------------------------------------------------
int main() {

  // problem dimension
  const int dim = 2;

  // geometry definition
  Geometry<dim> geometry;
  geometry.nel   = {1,1};
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
  geometry.coefs[ 5] = 0.0;  geometry.coefs[11] = 2.0;

  // linear elasticity problem creation
  TensorSize<dim>  nel = {4,4};
  TensorIndex<dim> deg = {3,3};
  auto problem =  ElasticityProblem<dim>(nel,deg,geometry);

  problem.output();

  return 0;
}


