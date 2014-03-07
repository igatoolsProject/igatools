
/*
 * Example for looping on the elements of a grid-like container
 * and accessing information that do not require the use of cache
 */
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

template <int dim>
void loop_on_grid()
{
    out << "Traversing the elements of a " << dim;
    out << "-dimensional grid." << endl;
    const int n_knots = 3;
    auto grid = CartesianGrid<dim>::create(n_knots);
    for (auto elem : *grid)
    {
        out << "The center of element: " << elem.get_flat_index();
        out << " is: "<< elem.center() << endl;
    }
    out << endl;
}


template <int dim>
void loop_on_space()
{
    out << "Traversing the elements of a " << dim;
    out << "-dimensional B-spline space." << endl;
    const int n_knots = 3;
    auto grid = CartesianGrid<dim>::create(n_knots);
    const int degree = 2;
    auto space = BSplineSpace<dim>::create(grid, degree);

    for (auto elem : *space)
    {
        out << "Element: " << elem.get_flat_index();
        out << " has global basis: " << elem.get_local_to_global() << endl;
    }
    out << endl;
}


int main()
{
    loop_on_grid<1>();
    loop_on_grid<2>();
    loop_on_grid<3>();

    out << endl;

    loop_on_space<1>();
    loop_on_space<2>();
    loop_on_space<3>();

    return 0;
}



