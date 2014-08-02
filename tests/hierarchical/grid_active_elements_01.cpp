/*
 *  Test for active grid elements
 *
 *  author: pauletti
 *  date: Aug 29, 2014
 *
 */

#include "../tests.h"

#include <igatools/geometry/cartesian_grid.h>


template<dim>
void
test()
{
    const int n_knots = 4;
    auto grid = CartesianGrid<dim>::create(n_knots);
    grid.print_info(out);

}




int main ()
{

    out.depth_console(10);

    test<2>();

    return 0;
}
