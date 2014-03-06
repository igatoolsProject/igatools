
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/io/writer.h>
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

template <int dim>
void plot_basis(const int deg)
{
    const int n_knots = deg + 2;
    auto grid  = CartesianGrid<dim>::create(n_knots);
    auto space = BSplineSpace<dim>::create(grid, deg);

    const int n_basis = space->get_num_basis();
    Vector coeffs(n_basis);

    TensorIndex<dim> basis_t_index(deg);
    auto basis_index = space->tensor_to_flat(basis_t_index);
    coeffs(basis_index) = 1.;

    out << "Coefficient vector of: " << basis_index << "-th basis" << endl;
    coeffs.print(out);
    out << endl;

    out << "Saving basis plot" << endl;
    const int n_plot_points = 5;
    Writer<dim> output(grid, n_plot_points);

    string field_name = "basis " + to_string(basis_index);
    output.add_field(space, coeffs, field_name);

    string file_name = "bspline_basis-" + to_string(dim) + "d";
    output.save(file_name);

}


int main()
{
    const int deg = 2;

    plot_basis<1>(deg);
    plot_basis<2>(deg);
    plot_basis<3>(deg);

    return 0;
}



