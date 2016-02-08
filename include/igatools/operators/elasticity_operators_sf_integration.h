//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#ifndef ELASTICITY_OPERATORS_SF_INTEGRATION_H_
#define ELASTICITY_OPERATORS_SF_INTEGRATION_H_

#include <igatools/base/config.h>
#include <igatools/operators/elliptic_operators_sf_integration.h>
#include <igatools/operators/integrator_sum_factorization.h>
//#include <igatools/utils/multi_array_utils.h>

//#include <vector>

IGA_NAMESPACE_OPEN

template <class PhysSpace>
class ElasticityOperatorsSFIntegration : public EllipticOperatorsSFIntegration<PhysSpace,PhysSpace>
{
public:

    /** Type for the element accessor of the physical space. */
    using Elem = typename PhysSpace::ElementAccessor;

    static const int       dim = PhysSpace::dim;
    static const int space_dim = PhysSpace::space_dim;


    void eval_operator_Bt_D_B(
        const Elem &elem,
        const std::array< std::array<vector<iga::DenseMatrix>,dim>,dim> &D_hat,
        DenseMatrix &local_stiffness) const;

private:
};



template <class PhysSpace>
inline
void
ElasticityOperatorsSFIntegration<PhysSpace>::
eval_operator_Bt_D_B(
    const Elem &elem,
    const std::array< std::array<vector<iga::DenseMatrix>,dim>,dim> &D_hat,
    DenseMatrix &local_stiffness) const
{
    const auto &phys_space = elem.get_physical_space();
    const auto   &ref_space = phys_space->get_reference_space();

    Assert(ref_space->has_tensor_product_structure(),ExcInvalidState());
//  Assert(ref_space->is_range_homogeneous(),ExcInvalidState());


    //--------------------------------------------------------------------------
    // getting the number of basis along each coordinate direction

    const Index comp = 0; // only scalar spaces for the moment

    TensorIndex<dim> degree = ref_space->get_degree()[comp];
    TensorSize<dim> n_basis_elem(degree + 1);

    Assert(n_basis_elem.flat_size() == elem.get_num_basis()/3,
           ExcDimensionMismatch(n_basis_elem.flat_size(),elem.get_num_basis()/3));

    const auto weight_basis = MultiArrayUtils<dim>::compute_weight(n_basis_elem);

    const TensorSize<dim> n_basis = n_basis_elem;
    const Index n_basis_comp_flat = n_basis.flat_size();
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // getting the 1D values -- begin
    std::array< ValueTable<Real>,dim> phi_1D;
    std::array< ValueTable<Real>,dim> grad_phi_1D;

    const auto &ref_elem_accessor = elem.get_ref_space_accessor().get_bspline_accessor();

    const auto &quad_scheme =ref_elem_accessor.get_quad_points();

    const auto phi_1D_table =
        ref_elem_accessor.evaluate_univariate_derivatives_at_points(0,quad_scheme);

    const auto grad_phi_1D_table =
        ref_elem_accessor.evaluate_univariate_derivatives_at_points(1,quad_scheme);

    phi_1D = phi_1D_table[0]; // only valid for scalar or isoparametric spaces
    grad_phi_1D = grad_phi_1D_table[0]; // only valid for scalar or isoparametric spaces

    // getting the 1D values -- end
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Checking dimensions of D_hat
#ifndef NDEBUG
    const auto n_points_tmp = D_hat[0][0].size();
    for (const auto &it0 : D_hat)
    {
        for (const auto &it1 : it0)
        {
            Assert(it1.size() == n_points_tmp,
                   ExcDimensionMismatch(it1.size(), n_points_tmp));
            for (const auto &it2 : it1)
            {
                Assert(it2.size1() == dim,
                       ExcDimensionMismatch(it2.size1(), dim));
                Assert(it2.size2() == dim,
                       ExcDimensionMismatch(it2.size2(), dim));
            }
        }
    }
#endif
    //--------------------------------------------------------------------------





    //--------------------------------------------------------------------------
    std::array< ValueTable<Real>,dim> phi_hat_alpha_i_1D;
    std::array< ValueTable<Real>,dim> phi_hat_beta_j_1D;

    const auto &quad_weights = quad_scheme.get_weights();

    const std::array<Real,dim> length_element_edge =
        elem.get_ref_space_accessor().get_coordinate_lengths();

    TensorSize<dim> n_points_1D = quad_scheme.get_num_points_direction();

    TensorSize<3> tensor_size_C0;
    tensor_size_C0[0] = n_points_1D.flat_size(); // theta size
    tensor_size_C0[1] = 1; // alpha size
    tensor_size_C0[2] = 1; // beta size

    DynamicMultiArray<Real,3> C0(tensor_size_C0);
    const Size n_entries = tensor_size_C0.flat_size();


    IntegratorSumFactorization<dim> integrate_sf;

//    local_stiffness.resize(dim*n_basis_comp_flat,dim*n_basis_comp_flat);
//    local_stiffness.clear();

    DenseMatrix loc_stiffness_tmp(n_basis_comp_flat,n_basis_comp_flat);
    for (int i = 0 ; i < dim ; ++i)
    {
        for (int k = 0 ; k < dim ; ++k)
        {
            if (k != i)
                phi_hat_alpha_i_1D[k] = phi_1D[k];
            else
                phi_hat_alpha_i_1D[k] = grad_phi_1D[k];
        }

        for (int j = 0 ; j < dim ; ++j)
        {
            for (int k = 0 ; k < dim ; ++k)
            {
                if (k != j)
                    phi_hat_beta_j_1D[k] = phi_1D[k];
                else
                    phi_hat_beta_j_1D[k] = grad_phi_1D[k];
            }


            //----------------------------------------------------
            // precalculation of the J[i](theta_i,alpha_i,beta_i) terms
            // (i.e. the weigths[theta_i] * phi_trial[alpha_i] * phi_test[beta_i] )
            const auto J = this->evaluate_w_phi1Dtrial_phi1Dtest(
                               phi_hat_alpha_i_1D,
                               phi_hat_beta_j_1D,
                               quad_weights,
                               length_element_edge);
            //----------------------------------------------------




            //----------------------------------------------------
            const auto &D_hat_pts = D_hat[i][j];

//            const auto n_points = D_hat_pts.size();
            Assert(D_hat_pts.size() == n_points_1D.flat_size(),
                   ExcDimensionMismatch(D_hat_pts.size(),n_points_1D.flat_size()));


            Assert(n_entries == D_hat_pts.size(),
                   ExcDimensionMismatch(n_entries,D_hat_pts.size()));
            //----------------------------------------------------

            for (int l = 0 ; l < dim ; ++l)
            {
                for (int m = 0 ; m < dim ; ++m)
                {

                    for (Index entry_id = 0 ; entry_id < n_entries ; ++entry_id)
//                         C0[entry_id] = D_hat_pts[entry_id][l][m];
                        C0[entry_id] = D_hat_pts[entry_id](l, m);


                    //----------------------------------------------------
                    loc_stiffness_tmp.clear();

                    integrate_sf(false, //non symmetric
                                 n_points_1D,
                                 n_basis,
                                 n_basis,
                                 J,
                                 C0,
                                 loc_stiffness_tmp);


                    for (int alpha = 0 ; alpha < n_basis_comp_flat ; ++alpha)
                        for (int beta = 0 ; beta < n_basis_comp_flat ; ++beta)
                            local_stiffness(alpha + l*n_basis_comp_flat, beta + m*n_basis_comp_flat) +=
                                loc_stiffness_tmp(alpha,beta);
                    //----------------------------------------------------

                } // end loop m
            } // end loop l
        }
    }
    //--------------------------------------------------------------------------

}




IGA_NAMESPACE_CLOSE



#endif // #ifndef ELASTICITY_OPERATORS_SF_INTEGRATION_H_
