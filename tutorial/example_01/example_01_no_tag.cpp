
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/io/writer.h>

using namespace iga;
using namespace std;

int main()
{
    const int dim = 2;

    const int n_knots = 3;
    cout << "Creating a " << dim << " dimensional cartesian grid" << endl;
    auto grid = CartesianGrid<dim>::create(n_knots);
    cout << "Number of elements: ";
    cout << grid->get_num_elements() << endl;

    Writer<dim> output(grid);
    output.save("grid");

    const int degree = 2;
    cout << "Creating a spline space of degree " << degree << endl;
    auto space = BSplineSpace<dim>::create(grid, degree);
    cout << "Number of basis functions: ";
    cout << space->get_num_basis() << endl;

    return 0;
}



