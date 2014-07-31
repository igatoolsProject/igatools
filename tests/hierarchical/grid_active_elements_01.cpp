/*
 *  Test for active grid elements
 *
 *  author: pauletti
 *  date: Aug 29, 2014
 *
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>

template<int dim>
void
test()
{
    const int n_knots = 4;
    auto grid = CartesianGrid<dim>::create(n_knots);
    grid->print_info(out);

    for (auto elem : *grid)
        out << elem.get_flat_index() << endl;

    for (auto elem : *grid)
    {
        if (elem.get_flat_index() % 2 == 0)
            elem.set_active(false);
    }

    grid->print_info(out);
    for (auto elem : *grid)
        out << elem.get_flat_index() << endl;
}




int main ()
{

    out.depth_console(10);

    test<0>();
    test<1>();
    test<2>();
    test<3>();

    return 0;
}
