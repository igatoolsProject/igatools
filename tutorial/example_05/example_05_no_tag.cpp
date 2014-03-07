


#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

template <int dim>
class PoissonPreparation
{
public:
    PoissonPreparation(const int n_knots,  const int deg);
    void local_assemble();

private:
    shared_ptr<CartesianGrid<dim>>  grid;
    shared_ptr<BSplineSpace<dim>>   space;
};


template <int dim>
PoissonPreparation<dim>::PoissonPreparation(const int n_knots,  const int deg)
    :
    grid {CartesianGrid<dim>::create(n_knots)},
     space {BSplineSpace<dim>::create(grid, deg)}
{}


template <int dim>
void  PoissonPreparation<dim>::local_assemble()
{
    out << "Assembling local contributions for the " << dim;
    out << "-dimensional laplace problem." << endl;

    const int n_basis = space->get_num_basis_per_element();
    DenseMatrix loc_mat(n_basis, n_basis);
    DenseVector loc_rhs(n_basis);

    const QGauss<dim> quad(2);
    ValueFlags fill_flags = ValueFlags::value|
                            ValueFlags::gradient|
                            ValueFlags::w_measure;
    auto elem = space->begin();
    const auto elem_end = space->end();
    elem->init_values(fill_flags, quad);

    const int n_qp = quad.get_num_points();
    for (; elem != elem_end; ++elem)
    {
        loc_mat = 0.;
        loc_rhs = 0.;

        elem->fill_values();
        auto values = elem->get_basis_values();
        auto grads  = elem->get_basis_gradients();
        auto w_meas = elem->get_w_measures();

        for (int i=0; i<n_basis; ++i)
        {
            auto grd_phi_i = grads.get_function_view(i);
            for (int j=0; j<n_basis; ++j)
            {
                auto grd_phi_j = grads.get_function_view(j);
                for (int qp=0; qp<n_qp; ++qp)
                    loc_mat(i,j) +=
                        scalar_product(grd_phi_i[qp], grd_phi_j[qp])
                        * w_meas[qp];
            }
            auto phi_i = values.get_function_view(i);
            for (int qp=0; qp<n_qp; ++qp)
                loc_rhs(i) += phi_i[qp][0] * w_meas[qp];
        }

        out << "Element matrix:" << endl;
        out << loc_mat << endl;
        out << "Element vector:" << endl;
        out << loc_rhs << endl;
    }
}



int main()
{
    const int deg = 1;
    const int n_knots = 2;
    PoissonPreparation<2> problem_2d(n_knots, deg);
    problem_2d.local_assemble();

    return 0;
}



