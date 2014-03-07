
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

template <int dim>
void loop_on_grid_with_cache()
{
    out << "Traversing the elements of a " << dim << "-dimensional grid." << endl;
    const int n_knots = 3;
    auto grid = CartesianGrid<dim>::create(n_knots);

    auto elem = grid->begin();
    const auto elem_end = grid->end();
    QGauss<dim> quad(2);
    ValueFlags fill_flag = ValueFlags::w_measure;
    elem->init_values(fill_flag, quad);

    for (; elem != elem_end; ++elem)
    {
        elem->fill_values();
        out << "The center of element: " << elem->get_flat_index();
        out << " is: "<< elem->center() << endl;

        auto w_meas = elem->get_w_measures();
        out << "The weighted measure is: ";
        w_meas.print_info(out);
        out << endl;
    }
    out << endl;
}


template <int dim>
void loop_on_space_with_cache()
{
    out << "Traversing the elements of a " << dim;
    out << "-dimensional B-spline space." << endl;
    const int n_knots = 3;
    auto grid = CartesianGrid<dim>::create(n_knots);
    const int degree = 2;
    auto space = BSplineSpace<dim>::create(grid, degree);

    auto elem = space->begin();
    const auto elem_end = space->end();
    elem->init_values(ValueFlags::value, QGauss<dim>(1));
    for (; elem != elem_end; ++elem)
    {
        elem->fill_values();
        out << "Element: " << elem->get_flat_index();
        out << " has global basis: " << elem->get_local_to_global() << endl;
        elem->get_basis_values().print_info(out);
        out<< endl;
    }
    out << endl;
}



int main()
{

    loop_on_grid_with_cache<1>();
    loop_on_grid_with_cache<2>();
    loop_on_grid_with_cache<3>();

    loop_on_space_with_cache<1>();
    loop_on_space_with_cache<2>();
    loop_on_space_with_cache<3>();

    return 0;
}



