#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

// class definition
template<int dim>
class PoissonAssembling {
  public:
    // constructors
    PoissonAssembling(const int nknt, const int deg);
    // the method that does everything
    void assemble();
  private:
    shared_ptr<const Grid<dim>>    grid;
    shared_ptr<const BSpline<dim>> space;
};

// constructors
template <int dim>
PoissonAssembling<dim>::PoissonAssembling(const int nknt, const int deg)
  :
  grid  {Grid<dim>::const_create(nknt)},
  space {BSpline<dim>::const_create(SplineSpace<dim>::const_create(deg,grid))}
{}

template <int dim>
void PoissonAssembling<dim>::assemble() {
  cout << "the assembling will be done, someday... I slightly promise!" << endl;
  cout << "cheking that private members are correctely created:" << endl;
  cout << "\t grid has " << grid->get_num_elements() << " elements" << endl;
  cout << "\tspace has " << space->get_num_basis() << " basis functions" << endl;

  // starting up the machinery for the cache
  auto       el  = space->begin();
  const auto eld = space->end();
  auto cache_handler = space->create_cache_handler();

  auto flag = space_element::Flags::value |
              space_element::Flags::gradient;
  cache_handler->set_element_flags(flag);

  auto quad = QGauss<dim>::create(space->get_degree)

}

int main() {

  const int dim  = 2;
  const int nknt = 5;
  const int deg  = 2;
  cout << "creating a " << pow((nknt-1),dim) << " elements grid with ";
  cout << pow((nknt+deg-1),2) << " basis functions!" << endl;

  PoissonAssembling<dim> prob(nknt,deg);
  prob.assemble();

  return 0;
}



