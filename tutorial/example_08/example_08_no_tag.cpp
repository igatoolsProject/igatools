
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/basis_functions/physical_space.h>

#include <igatools/io/writer.h>

using namespace iga;
using namespace std;
using numbers::PI;

template<int dim>
void physical_space(const int deg)
{
    using RefSpace = BSplineSpace<dim>;
    using PushFw   = PushForward<Transformation::h_grad, dim>;
    using Space    = PhysicalSpace<RefSpace, PushFw>;

    BBox<dim> box;
    box[0] = {{0.5,1}};
    for (int i=1; i<dim; ++i)
        box[i] = {{PI/4,PI/2}};

    const int n_knots = 3;
    auto grid = CartesianGrid<dim>::create(box, n_knots);
    auto ref_space = RefSpace::create(grid, deg);

    auto map  = BallMapping<dim>::create(grid);
    auto push_fordward = PushFw::create(map);
    auto space = Space::create(ref_space, push_fordward);

    const int n_plot_points = 2;
    Writer<dim> writer(map, n_plot_points);
    string filename = "ball_geometry-" + to_string(dim) + "d" ;
    writer.save(filename);
}

int main()
{
    physical_space<2>(1);


    return  0;
}
