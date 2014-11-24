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

#include <igatools/geometry/mapping_element.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <boost/numeric/ublas/operation.hpp>
IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
auto
MappingElement<dim_, codim_>::
get_external_normals() const -> ValueVector<Points<space_dim> >
{
	Assert(codim==1, ExcNotImplemented());
	ValueVector<Points<space_dim> > res;
	const auto &DF = this->template get_values<1, dim>(0);
	const auto n_points = DF.get_num_points();

	res.resize(n_points);
	for (int i = 0; i< n_points; ++i)
	{
		res[i] = cross_product<dim, codim>(DF[i]);
		res[i] /= res[i].norm();
	}

	return res;
}



template<int dim_, int codim_>
auto
MappingElement<dim_, codim_>::
get_principal_curvatures() const -> ValueVector<vector<Real>>
{
	Assert(codim==1, ExcNotImplemented());

	const auto &D2_F  = this->template get_values<2, dim>(0);
	const auto normal = this->get_external_normals();

	const auto n_points = D2_F.get_num_points();
	ValueVector<vector<Real>> res(n_points);
	const auto G_inv = compute_inv_first_fundamental_form();
	DenseMatrix A(dim, dim);
	for (int pt = 0; pt < n_points; ++pt)
	{
		const auto B = unroll_to_matrix(G_inv[pt]);
		for (int i = 0; i<dim; ++i)
			for (int j = 0; j<dim; ++j)
				A(i,j) = -scalar_product(D2_F[pt][i][j], normal[pt]);
		DenseMatrix C(dim, dim);
		boost::numeric::ublas::axpy_prod(B, A, C, true);

		res[pt] = C.eigen_values();
	}

	return res;
}



//    template<int sub_dim>
//    ValueVector<Points<space_dim> >
//    get_boundary_normals(const int s_id) const
//    {
//        Assert(dim==sub_dim+1, ExcNotImplemented());
//        ValueVector<Points<space_dim>> res;
//        const auto &DF_inv = get_inverse_values<1, sub_dim>(s_id);
//        const auto n_hat  = this->get_grid()->template get_boundary_normals<sub_dim>(s_id)[0];
//
//        const auto n_points = DF_inv.get_num_points();
//        res.resize(n_points);
//        for (int i = 0; i< n_points; ++i)
//        {
//            const auto DF_inv_t = co_tensor(transpose(DF_inv[i]));
//            res[i] = action(DF_inv_t, n_hat);
//            res[i] /= res[i].norm();
//        }
//        return res;
//    }

IGA_NAMESPACE_CLOSE
