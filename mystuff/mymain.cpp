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
//#include "my_formula_grid_function.h"
#include "custom_grid_function.h"

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
    x *= sin( 2.0 * pts[idim] * PI);
  }
  return x;
}
//template<int dim> // exact solution gradient

template<int dim> // source term
Values<dim,1,1> source_term(Points<dim> pts) {
  Values<dim,1,1> x = {4.0 * dim * PI * PI};
  for (int idim=0; idim<dim; idim++) {
    x *= sin(2.0 * pts[idim] * PI);
  }
  //printf(" evaluation complete:  f(%f,%f) = %f\n",pts[0],pts[1],x);
  return x;
}


// ----------------------------------------------------------------------------
//   MAIN
// ----------------------------------------------------------------------------
int main() {

  const int dim  = 2;
  const int nel  = 32;
  const int deg  = 2;

  Real err1, err2, eoc, h1, h2;
  
  using CustomFunct = grid_functions::CustomGridFunction<dim,1>;
  /*for (int iel=4; iel>nel; iel+=4) {
    auto problem = PoissonProblem<dim>::create(iel,deg);
    auto f = CustomFunct::const_create(problem->grid,&source_term);
    problem->assemble(f);
    problem->solve();
    auto u = CustomFunct::const_create(problem->grid,&exact_solution);
    err2 = problem->l2_error(u);
    h2   = (double)1.0/iel;
    printf("deg = %d, nel = %2d:\t l2_err = %le",deg,iel,err2);
    if (iel==4)
      printf("  eoc = -\n");
    else {
      eoc = (log(err2)-log(err1))/(log(h2)-log(h1));
      printf("  eoc = %1.5f\n",eoc);
    }
    err1 = err2;
    h1 = h2; 
  }*/

  auto problem = PoissonProblem<dim>::create(8,deg);
  auto f = CustomFunct::const_create(problem->grid,&source_term);
  problem->assemble(f);
  problem->solve();
  auto u = CustomFunct::const_create(problem->grid,&exact_solution);
  err2 = problem->l2_error(u);

  auto solution = IgGridFunction<dim,1>::const_create(problem->basis,*problem->sol);
  const int npt = 10;
  Writer<dim> output(problem->grid,npt);
  output.template add_field(*solution,"solution");
  output.template add_field(*u,"exact solution");
  output.template add_field(*f,"source term");
  output.save("plot");

  return 0;
}
