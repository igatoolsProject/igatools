//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
//
// This file is part of the igatools library.
//
// The igatools library is free software: you can use it, redistribute
// it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-+--------------------------------------------------------------------



#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
// [new includes]
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
// [new includes]
#include <igatools/base/logstream.h>

using namespace iga;
using namespace std;

LogStream out;

// [class declaration]
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
// [class declaration]


// [constructor]
template <int dim>
PoissonPreparation<dim>::PoissonPreparation(const int n_knots,  const int deg)
    :
    grid {CartesianGrid<dim>::create(n_knots)},
     space {BSplineSpace<dim>::create(deg, grid)}
{}
// [constructor]


// [assemble function]
template <int dim>
void  PoissonPreparation<dim>::local_assemble()
{
    out << "Assembling local contributions for the " << dim;
    out << "-dimensional laplace problem." << endl;
    // [assemble function]


    // [iterate as before]
    const QGauss<dim> quad(2);
    ValueFlags fill_flags = ValueFlags::value|
                            ValueFlags::gradient|
                            ValueFlags::w_measure;
    auto elem = space->begin();
    const auto elem_end = space->end();
    elem->init_cache(fill_flags, quad);

    const int n_qp = quad.get_num_points();
    for (; elem != elem_end; ++elem)
    {
        // [iterate as before]

        // [local matrix]
        const int n_basis = elem->get_num_basis();
        DenseMatrix loc_mat(n_basis, n_basis);
        DenseVector loc_rhs(n_basis);
        // [local matrix]

        // [get the values]
        elem->fill_cache();
        auto values = elem->get_basis_values();
        auto grads  = elem->get_basis_gradients();
        auto w_meas = elem->get_w_measures();
        // [get the values]

        // [assemble]
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
        // [assemble]

        // [print info]
        out << "Element matrix:" << endl;
        out << loc_mat << endl;
        out << "Element vector:" << endl;
        out << loc_rhs << endl;
        // [print info]
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



