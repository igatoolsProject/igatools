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

#ifndef ELLIPTIC_OPERATORS_SF_INTEGRATION_H_
#define ELLIPTIC_OPERATORS_SF_INTEGRATION_H_

#include <igatools/base/config.h>
#include <igatools/operators/elliptic_operators.h>
#include <igatools/operators/integrator_sum_factorization.h>
#include <igatools/utils/multi_array_utils.h>
#include <igatools/basis_functions/bspline_element.h>

#include <vector>

IGA_NAMESPACE_OPEN


#define TIME_PROFILING


/**
 * @brief Class containing the methods for the evaluation of some elliptic operators
 * (mass matrix, stiffness matrix, etc.) on a Bezier element
 * using the <em>sum-factorization quadrature approach</em>.
 *
 * All public methods take as input two PhysicalSpaceElementAccessor object
 * (one for the test space and the other for the trial space) and has an output argument that is
 * the local matrix relative to the elliptic operator that is evaluated.
 *
 *
 * @author M. Martinelli
 * @date 16 Apr 2014
 */

template <int dim_,int range_,int rank_>
class EllipticOperatorsSFIntegrationBSpline
  : public EllipticOperators<dim_,0,range_,rank_>
{
public:
  using base_t = EllipticOperators<dim_,0,range_,rank_>;
  using self_t = EllipticOperatorsSFIntegrationBSpline<dim_,range_,rank_>;


  using base_t::dim;
  using base_t::space_dim;

  using Clock = typename base_t::Clock;
  using TimePoint = typename base_t::TimePoint;
  using Duration = typename base_t::Duration;


  /** Type for the element accessor of the <em>test</em> physical space. */
  using ElemTest = BSplineElement<dim_,range_,rank_>;

  /** Type for the element accessor of the <em>trial</em> physical space. */
  using ElemTrial = ElemTest;

  /** The constructors and destructor are inherithed from the base class. */
  using base_t::base_t;



  static const int n_components = BSpline<dim_,range_,rank_>::n_components;



  /** @name Assignment operators */
  ///@{
  /** Copy assignment operator. */
  self_t &
  operator=(const self_t &in) = delete;


  /** Move assignment operator. */
  self_t &
  operator=(self_t &&in) = delete;
  ///@}


  /**
   * This function evaluates the local (i.e. element-based) matrix \f$ A_e \f$
   * for witch its entries are
   * \f[
        (A_e)_{ij} = \int_{\Omega_e} \phi^{e,\text{test}}_i(x) C(x) \phi^{e,\text{trial}}_j(x) \; d \Omega.
     \f]
   * The matrix \f$ A_e \f$ is commonly referred as <em>local mass-matrix</em>.
   */
  template <int sdim>
  void eval_operator_u_v(
    const ElemTest &elem_test,
    const ElemTrial &elem_trial,
    const ValueVector<Real> &c,
    const int s_id,
    DenseMatrix &operator_u_v) const;

  /**
   * This function evaluates the local (i.e. element-based) matrix \f$ A_e \f$
   * for witch its entries are
   * \f[
        (A_e)_{ij} = \int_{\Omega_e} \sum_{r=1}^{sp\_dim} \sum_{s=1}^{sp\_dim}
        \bigl( \nabla \phi^{e,\text{test}}_i \bigr)_r
        \, C_{sr}(x) \,
        \bigl( \nabla \phi^{e,\text{trial}}_j \bigr)_s \; d \Omega.
     \f]
   * The matrix \f$ A_e \f$ is commonly referred as <em>local stiffness-matrix</em>.
   */
  template <int sdim>
  void eval_operator_gradu_gradv(
    const ElemTest &elem_test,
    const ElemTrial &elem_trial,
    const SafeSTLArray<ValueVector<Real>,dim_> &coeffs,
    const int s_id,
    DenseMatrix &operator_gradu_gradv) const;


  /**
   * This function evaluates the local (i.e. element-based) vector \f$ f_e \f$
   * for witch its entries are
   * \f[
        (f_e)_{i} = \int_{\Omega_e} \phi^{e,\text{test}}_i
        f(x)  \; d \Omega.
     \f]
   */
  void eval_operator_rhs_v(
    const ElemTest &elem_test,
    const ValueVector<typename BSpline<dim_,range_,rank_>::Value> &f,
    DenseVector &operator_rhs_v) const;


  /**
   * This function evaluates the local (i.e. element-based) matrix \f$ A_e \f$
   * for witch its entries are
   * \f[
        (A_e)_{ij} = \int_{\Omega_e} \sum_{s=1}^{sp\_dim}
        \phi^{e,\text{test}}_i
        \, \beta_{s}(x) \,
        \bigl( \nabla \phi^{e,\text{trial}}_j \bigr)_s \; d \Omega
        = \int_{\Omega_e}
        \phi^{e,\text{test}}_i
        \, \vec{\beta}(x) \, \cdot \,
        \nabla \phi^{e,\text{trial}}_j \; d \Omega .
     \f]
   */
  void eval_operator_gradu_v(
    const ElemTest &elem_test,
    const ElemTrial &elem_trial,
    const ValueVector<typename BSpline<dim_,range_,rank_>::Gradient> &beta,
    DenseMatrix &operator_gradu_v) const;

protected:

  /**
   * Returns the quadrature weights multiplied by the one-dimensional basis
   * for test and trial space, as needed by the integration using the
   * sum factorization technique.
   * \f[ J[i]_{\theta_i,\alpha_i,\beta_i} =
     w_{i,\theta_i}
     \phi^{\text{test}}_{\beta_i}(x_{\theta_i})
     \phi^{\text{trial}}_{\alpha_i}(x_{\theta_i})
     \f]
   * where \f$ w_{i,\theta_i} \f$ is the quadrature weight
   * relative to the the point \f$x_{\theta_i}\f$
   * along the \f$i\f$-th direction.
   */
  SafeSTLArray<DynamicMultiArray<Real,3>,dim_>
  evaluate_w_phi1Dtrial_phi1Dtest(
    const SafeSTLArray<ValueTable<Real>,dim_> &phi_1D_test,
    const SafeSTLArray<ValueTable<Real>,dim_> &phi_1D_trial,
    const TensorProductArray<dim_> &quad_weights,
    const SafeSTLArray<Real,dim_> &length_element_edge) const;


  template <int sdim>
  void
  integrate_add_operator_general_order(
    const ElemTest &elem_test,
    const int comp_test,
    const int row_id_begin,
    const int row_id_last,
    const TensorIndex<dim_> &deriv_order_test,
    const ElemTrial &elem_trial,
    const int comp_trial,
    const int col_id_begin,
    const int col_id_last,
    const TensorIndex<dim_> &deriv_order_trial,
    const ValueVector<Real> &c,
    const int s_id,
    DenseMatrix &op) const;
};











template <int dim_,int range_,int rank_>
inline
auto
EllipticOperatorsSFIntegrationBSpline<dim_,range_,rank_>::
evaluate_w_phi1Dtrial_phi1Dtest(
  const SafeSTLArray<ValueTable<Real>,dim_> &phi_1D_test,
  const SafeSTLArray<ValueTable<Real>,dim_> &phi_1D_trial,
  const TensorProductArray<dim_> &quad_weights,
  const SafeSTLArray<Real,dim_> &length_element_edge) const
-> SafeSTLArray<DynamicMultiArray<Real,3>,dim_>
{
  SafeSTLArray<DynamicMultiArray<Real,3>,dim_> moments;

  for (int dir = 0 ; dir < dim ; ++dir)
  {
    const auto &phi_test  = phi_1D_test [dir];
    const auto &phi_trial = phi_1D_trial[dir];
    const auto &w = quad_weights.get_data_direction(dir);

    const Size n_basis_test  = phi_test.get_num_functions();
    const Size n_basis_trial = phi_trial.get_num_functions();
    const Size n_pts = w.size();

    TensorSize<3> moments1D_tensor_size;
    moments1D_tensor_size[0] = n_pts;
    moments1D_tensor_size[1] = n_basis_trial;
    moments1D_tensor_size[2] = n_basis_test;

    auto &moments1D = moments[dir];
    moments1D.resize(moments1D_tensor_size);

    Assert(phi_test.get_num_points() == n_pts,
    ExcDimensionMismatch(phi_test.get_num_points(),n_pts));
    Assert(phi_trial.get_num_points() == n_pts,
    ExcDimensionMismatch(phi_trial.get_num_points(),n_pts));


    SafeSTLVector<Real> w_times_edge_length(n_pts);

    const Real edge_length = length_element_edge[dir];
    for (int jpt = 0 ; jpt < n_pts ; ++jpt)
      w_times_edge_length[jpt] = w[jpt] * edge_length;


    Index flat_id_I = 0 ;
    for (Index f_id_test = 0 ; f_id_test < n_basis_test ; ++f_id_test)
    {
      const auto phi_1D_test = phi_test.get_function_view(f_id_test);

      for (Index f_id_trial = 0 ; f_id_trial < n_basis_trial ; ++f_id_trial)
      {
        const auto phi_1D_trial = phi_trial.get_function_view(f_id_trial);

        for (int pt = 0 ; pt < n_pts ; ++pt)
        {
          moments1D[flat_id_I++] =
          w_times_edge_length[pt] * phi_1D_test[pt] * phi_1D_trial[pt];
        } // end loop pt
      } // end loop mu1
    } // end loop mu2
  } // end loop dir

  return moments;
}



template <int dim_,int range_,int rank_>
template <int sdim>
inline
void
EllipticOperatorsSFIntegrationBSpline<dim_,range_,rank_>::
integrate_add_operator_general_order(
  const ElemTest &elem_test,
  const int comp_test,
  const int row_id_begin,
  const int row_id_last,
  const TensorIndex<dim_> &deriv_order_test,
  const ElemTrial &elem_trial,
  const int comp_trial,
  const int col_id_begin,
  const int col_id_last,
  const TensorIndex<dim_> &deriv_order_trial,
  const ValueVector<Real> &coeffs_test_trial,
  const int s_id,
  DenseMatrix &op) const
{


#ifdef TIME_PROFILING
  const auto start_initialization = Clock::now();
#endif //#ifdef TIME_PROFILING

  //--------------------------------------------------------------------------
  const auto test_basis = elem_test.get_bspline_basis();
  const auto trial_basis = elem_trial.get_bspline_basis();

  //TODO: only the symmetric case is tested. In the non symmetric case, we need to check that
  // the physical space iterators have the same grid, map, reference space, index, etc.
  const bool same_basis = (test_basis == trial_basis);
  Assert(same_basis,ExcNotImplemented());
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  const auto &grid_elem_test = elem_test.get_grid_element();
  const auto quad_elem_test = grid_elem_test.template get_quad<sdim>();
  Assert(quad_elem_test != nullptr,ExcNullPtr());


  const auto &grid_elem_trial = elem_trial.get_grid_element();
  const auto quad_elem_trial = grid_elem_trial.template get_quad<sdim>();
  Assert(quad_elem_trial != nullptr,ExcNullPtr());


  // checks that the elements on the grid are the same
  Assert(test_basis->get_grid() == trial_basis->get_grid(),
         ExcMessage("Different grids between test space and trial space."));


  Assert(grid_elem_test.get_index() == grid_elem_trial.get_index(),
         ExcMessage("Test and trial elements have different indices."));
  const auto &grid_elem = grid_elem_test;

  Assert(quad_elem_test == quad_elem_trial,
         ExcMessage("Test and trial elements have different quadrature schemes."));

  const auto quad_scheme = extend_sub_elem_quad<sdim,dim_>(*quad_elem_test,s_id);

//  const auto quad_scheme = quad_elem_test;
  Assert(quad_scheme.is_tensor_product(),
         ExcMessage("The quadrature scheme has not the tensor-product structure."))
  const auto points_t_size = quad_scheme.get_num_coords_direction();
  //--------------------------------------------------------------------------




  //--------------------------------------------------------------------------
  // getting the number of basis along each coordinate direction of the current scalar component of the test space
//  const int n_rows = elem_test.get_num_basis_comp(comp_test);
  const auto basis_t_size_elem_test = elem_test.get_num_splines_1D(comp_test);
  Assert(basis_t_size_elem_test.flat_size() == elem_test.get_num_basis_comp(comp_test),
         ExcDimensionMismatch(basis_t_size_elem_test.flat_size(),elem_test.get_num_basis_comp(comp_test)));
  //--------------------------------------------------------------------------



  //--------------------------------------------------------------------------
  // getting the number of basis along each coordinate direction of the current scalar component of the trial space
//  const int n_cols = elem_trial.get_num_basis_comp(comp_trial);
  const auto basis_t_size_elem_trial = elem_trial.get_num_splines_1D(comp_trial);
  Assert(basis_t_size_elem_trial.flat_size() == elem_trial.get_num_basis_comp(comp_trial),
         ExcDimensionMismatch(basis_t_size_elem_trial.flat_size(),elem_trial.get_num_basis_comp(comp_trial)));
  //--------------------------------------------------------------------------






  //--------------------------------------------------------------------------
  //TODO (martinelli, Jan 29, 2016): up to now the cost of copying the matrix entries for the
  // symmetric part is equal to compute the entries of the non-symmetric case.

  const bool is_symmetric = same_basis &&
                            (comp_test == comp_trial) &&
                            (deriv_order_test == deriv_order_trial) &&
                            false;
  //--------------------------------------------------------------------------





  //--------------------------------------------------------------------------
  // getting the 1D values for the test and trial space -- begin
  SafeSTLArray<ValueTable<Real>,dim> phi_1D_test;
  SafeSTLArray<ValueTable<Real>,dim> phi_1D_trial;

  const auto &phi_1D_table_test  = elem_test.template get_splines1D_table(dim_,0); // the 0 index is because sdim == dim
  const auto &phi_1D_table_trial = elem_test.template get_splines1D_table(dim_,0); // the 0 index is because sdim == dim

  const auto &phi_1D_comp_test  = phi_1D_table_test [comp_test];
  const auto &phi_1D_comp_trial = phi_1D_table_trial[comp_trial];

  for (int i = 0 ; i < dim ; ++i)
  {
    const int n_pts_1D = points_t_size[i];

    const auto &v_test = phi_1D_comp_test[i].get_derivative(deriv_order_test[i]);
    Assert(v_test.get_num_rows() == basis_t_size_elem_test [i],
           ExcDimensionMismatch(v_test.get_num_rows(),basis_t_size_elem_test[i]));
    Assert(v_test.get_num_cols() == n_pts_1D,
           ExcDimensionMismatch(v_test.get_num_cols(),n_pts_1D));

    const auto &v_trial = phi_1D_comp_trial[i].get_derivative(deriv_order_trial[i]);
    Assert(v_trial.get_num_rows() == basis_t_size_elem_trial[i],
           ExcDimensionMismatch(v_trial.get_num_rows(),basis_t_size_elem_trial[i]));
    Assert(v_trial.get_num_cols() == n_pts_1D,
           ExcDimensionMismatch(v_trial.get_num_cols(),n_pts_1D));

    phi_1D_test [i].resize(basis_t_size_elem_test [i],n_pts_1D);
    phi_1D_trial[i].resize(basis_t_size_elem_trial[i],n_pts_1D);

    for (int fn = 0 ; fn < basis_t_size_elem_test[i] ; ++fn)
    {
      auto phi_1D_test_fn = phi_1D_test[i].get_function_view(fn);
      for (int pt = 0 ; pt < n_pts_1D ; ++pt)
        phi_1D_test_fn[pt] = v_test(fn,pt); // only valid for scalar spaces
    }

    for (int fn = 0 ; fn < basis_t_size_elem_trial[i] ; ++fn)
    {
      auto phi_1D_trial_fn = phi_1D_trial[i].get_function_view(fn);
      for (int pt = 0 ; pt < n_pts_1D ; ++pt)
        phi_1D_trial_fn[pt] = v_trial(fn,pt); // only valid for scalar spaces
    }
  }
  // getting the 1D values for the test and trial space -- end
  //--------------------------------------------------------------------------




#ifdef TIME_PROFILING
  this->elapsed_time_initialization_ += Clock::now() - start_initialization;
#endif //#ifdef TIME_PROFILING
  //--------------------------------------------------------------------------





  //----------------------------------------------------
  // precalculation of the J[i](theta_i,alpha_i,beta_i) terms
  // (i.e. the weigths[theta_i] * phi_trial[alpha_i] * phi_test[beta_i] )
#ifdef TIME_PROFILING
  const auto start_compute_phi1Dtest_phi1Dtrial = Clock::now();
#endif //#ifdef TIME_PROFILING



  const auto &l_tmp = grid_elem.template get_side_lengths<sdim>(s_id);
  const auto &sub_element = UnitElement<dim>::template get_elem<sdim>(s_id);
  const auto &active_directions = sub_element.active_directions;

  SafeSTLArray<Real,dim> length_element_edges(1.0);
  for (int i = 0 ; i < sdim ; ++i)
  {
    const auto active_dir = active_directions[i];
    length_element_edges[active_dir] = l_tmp[i];
  }


  const auto w_phi1Dtrial_phi1Dtest = evaluate_w_phi1Dtrial_phi1Dtest(
                                        phi_1D_test,
                                        phi_1D_trial,
                                        quad_scheme.get_weights_1d(),
                                        length_element_edges);

#ifdef TIME_PROFILING
  this->elapsed_time_compute_phi1Dtest_phi1Dtrial_ +=
    Clock::now() - start_compute_phi1Dtest_phi1Dtrial;
#endif //#ifdef TIME_PROFILING
  //----------------------------------------------------




  //----------------------------------------------------
  // Assembly of the local mass matrix using sum-factorization -- begin
#ifdef TIME_PROFILING
  const auto start_sum_factorization = Clock::now();
#endif //#ifdef TIME_PROFILING

  TensorSize<3> tensor_size_C0;
  tensor_size_C0[0] = points_t_size.flat_size(); // theta size
  tensor_size_C0[1] = 1; // alpha size
  tensor_size_C0[2] = 1; // beta size

  DynamicMultiArray<Real,3> C0(tensor_size_C0);
  const Size n_entries = tensor_size_C0.flat_size();


  Assert(coeffs_test_trial.size() == quad_scheme.get_num_points(),
         ExcDimensionMismatch(coeffs_test_trial.size(),quad_scheme.get_num_points()));
  Assert(n_entries == coeffs_test_trial.flat_size(),
         ExcDimensionMismatch(n_entries,coeffs_test_trial.flat_size()));
  for (Index entry_id = 0 ; entry_id < n_entries ; ++entry_id)
    C0[entry_id] = (coeffs_test_trial)[entry_id];


  IntegratorSumFactorization<dim> integrate_sf;
  integrate_sf(is_symmetric,
               points_t_size,
               basis_t_size_elem_test,
               basis_t_size_elem_trial,
               w_phi1Dtrial_phi1Dtest,
               C0,
               row_id_begin, row_id_last,
               col_id_begin, col_id_last,
               op);


#ifdef TIME_PROFILING
  this->elapsed_time_sum_factorization_ += Clock::now() - start_sum_factorization;
#endif //#ifdef TIME_PROFILING
  // Assembly of the local mass matrix using sum-factorization -- end
  //----------------------------------------------------

//  return op;
}


template <int dim_,int range_,int rank_>
template <int sdim>
inline
void
EllipticOperatorsSFIntegrationBSpline<dim_,range_,rank_>::
eval_operator_u_v(
  const ElemTest &elem_test,
  const ElemTrial &elem_trial,
  const ValueVector<Real> &coeffs,
  const int s_id,
  DenseMatrix &operator_u_v) const
{
  // TODO (martinelli, Jan 25, 2016): for the moment, only the vector case is treated
  Assert(rank_ == 1,ExcDimensionMismatch(rank_,1));


#ifdef TIME_PROFILING
  this->elapsed_time_initialization_.zero();
  this->elapsed_time_compute_phi1Dtest_phi1Dtrial_.zero();
  this->elapsed_time_sum_factorization_.zero();

#endif //#ifdef TIME_PROFILING


  const TensorIndex<dim_> deriv_order_test(0);
  const TensorIndex<dim_> deriv_order_trial(0);

  int row_id_begin = 0;
  int row_id_last = 0;
  for (int comp_test = 0 ; comp_test < n_components ; ++comp_test)
  {
    // getting the number of basis along each coordinate direction of the current scalar component of the test space
    const int n_rows = elem_test.get_num_basis_comp(comp_test);
    row_id_last = row_id_begin + n_rows - 1;


    int col_id_begin = 0;
    int col_id_last = 0;
    for (int comp_trial = 0 ; comp_trial < n_components ; ++comp_trial)
    {
      // getting the number of basis along each coordinate direction of the current scalar component of the trial space
      const int n_cols = elem_trial.get_num_basis_comp(comp_trial);
      col_id_last = col_id_begin + n_cols - 1;

      if (comp_test == comp_trial)
      {
        this->integrate_add_operator_general_order<sdim>(
          elem_test,
          comp_test,
          row_id_begin,row_id_last,
          deriv_order_test,
          elem_trial,
          comp_trial,
          col_id_begin,col_id_last,
          deriv_order_trial,
          coeffs,
          s_id,
          operator_u_v);
      }
      col_id_begin = col_id_last + 1;
    } // end loop comp_trial

    row_id_begin = row_id_last + 1;
  } // end loop comp_test


#ifdef TIME_PROFILING
  std::cout << "Elapsed_seconds initialization = " << this->elapsed_time_initialization_.count() << std::endl;
  std::cout << "Elapsed seconds w * phi1d_trial * phi1d_test = "
            << this->elapsed_time_compute_phi1Dtest_phi1Dtrial_.count() << std::endl;

  std::cout << "Elapsed seconds sum-factorization = " << this->elapsed_time_sum_factorization_.count() << std::endl;

  const Duration elapsed_time_assemble = this->elapsed_time_sum_factorization_ +
                                         this->elapsed_time_compute_phi1Dtest_phi1Dtrial_ +
                                         this->elapsed_time_initialization_ ;
  std::cout << "Elapsed seconds assemblying = " << elapsed_time_assemble.count() << std::endl;



  std::cout << std::endl;
#endif //#ifdef TIME_PROFILING
  // Assembly of the local mass matrix using sum-factorization -- end
  //----------------------------------------------------

}



template <int dim_,int range_,int rank_>
template<int sdim>
inline
void
EllipticOperatorsSFIntegrationBSpline<dim_,range_,rank_>::
eval_operator_gradu_gradv(
  const ElemTest &elem_test,
  const ElemTrial &elem_trial,
  const SafeSTLArray<ValueVector<Real>,dim_> &coeffs,
  const int s_id,
  DenseMatrix &operator_gradu_gradv) const
{
  // TODO (martinelli, Jan 25, 2016): for the moment, only the scalar case is treated
//  Assert(range_ == 1,ExcDimensionMismatch(range_,1));
  Assert(rank_ == 1,ExcDimensionMismatch(rank_,1));


#ifdef TIME_PROFILING
  this->elapsed_time_initialization_.zero();
  this->elapsed_time_compute_phi1Dtest_phi1Dtrial_.zero();
  this->elapsed_time_sum_factorization_.zero();

#endif //#ifdef TIME_PROFILING



  int row_id_begin = 0;
  int row_id_last = 0;
  for (int comp_test = 0 ; comp_test < n_components ; ++comp_test)
  {
    // getting the number of basis along each coordinate direction of the current scalar component of the test space
    const int n_rows = elem_test.get_num_basis_comp(comp_test);
    row_id_last = row_id_begin + n_rows - 1;


    int col_id_begin = 0;
    int col_id_last = 0;
    for (int comp_trial = 0 ; comp_trial < n_components ; ++comp_trial)
    {
      // getting the number of basis along each coordinate direction of the current scalar component of the trial space
      const int n_cols = elem_trial.get_num_basis_comp(comp_trial);
      col_id_last = col_id_begin + n_cols - 1;

      if (comp_test == comp_trial)
      {
        for (int k = 0 ; k < dim_ ; ++k)
        {
          TensorIndex<dim_> deriv_order_test(0);
          TensorIndex<dim_> deriv_order_trial(0);

          deriv_order_test[k] = 1;
          deriv_order_trial[k] = 1;


          this->integrate_add_operator_general_order<sdim>(
            elem_test,
            comp_test,
            row_id_begin,row_id_last,
            deriv_order_test,
            elem_trial,
            comp_trial,
            col_id_begin,col_id_last,
            deriv_order_trial,
            coeffs[k],
            s_id,
            operator_gradu_gradv);
        } // end loop k
      } // end if (comp_test == comp_trial)

      col_id_begin = col_id_last + 1;
    } // end loop comp_trial

    row_id_begin = row_id_last + 1;
  } // end loop comp_test


#ifdef TIME_PROFILING
  std::cout << "Elapsed_seconds initialization = " << this->elapsed_time_initialization_.count() << std::endl;
  std::cout << "Elapsed seconds w * phi1d_trial * phi1d_test = "
            << this->elapsed_time_compute_phi1Dtest_phi1Dtrial_.count() << std::endl;

  std::cout << "Elapsed seconds sum-factorization = " << this->elapsed_time_sum_factorization_.count() << std::endl;

  const Duration elapsed_time_assemble = this->elapsed_time_sum_factorization_ +
                                         this->elapsed_time_compute_phi1Dtest_phi1Dtrial_ +
                                         this->elapsed_time_initialization_ ;
  std::cout << "Elapsed seconds assemblying = " << elapsed_time_assemble.count() << std::endl;



  std::cout << std::endl;
#endif //#ifdef TIME_PROFILING
  // Assembly of the local mass matrix using sum-factorization -- end
  //----------------------------------------------------


#if 0

  //----------------------------------------------------
  // Assembly of the local stiffness matrix using sum-factorization -- begin
#ifdef TIME_PROFILING
  const TimePoint start_assembly_stiffness_matrix = Clock::now();

  //--------------------------------------------------------------------------
  const auto start_initialization = Clock::now();
#endif //#ifdef TIME_PROFILING




  //--------------------------------------------------------------------------
  // getting the number of basis along each coordinate direction for the test and trial space

  const Index comp = 0; // only scalar spaces for the moment

  // test space -- begin
  TensorIndex<dim> degree_test = elem_test.get_physical_space()->get_reference_space()->get_degree()[comp];
  TensorSize<dim> n_basis_elem_test(degree_test + 1);

  Assert(n_basis_elem_test.flat_size()==elem_test.get_num_basis(),
         ExcDimensionMismatch(n_basis_elem_test.flat_size(),elem_test.get_num_basis()));

  const auto weight_basis_test = MultiArrayUtils<dim>::compute_weight(n_basis_elem_test);

  const TensorSize<dim> n_basis_test = n_basis_elem_test;
  // test space -- end


  // trial space -- begin
  TensorIndex<dim> degree_trial = elem_trial.get_physical_space()->get_reference_space()->get_degree()[comp];
  TensorSize<dim> n_basis_elem_trial(degree_trial + 1);

  Assert(n_basis_elem_trial.flat_size()==elem_trial.get_num_basis(),
         ExcDimensionMismatch(n_basis_elem_trial.flat_size(),elem_trial.get_num_basis()));

  const auto weight_basis_trial = MultiArrayUtils<dim>::compute_weight(n_basis_elem_trial);

  const TensorSize<dim> n_basis_trial = n_basis_elem_trial;
  // trial space -- end

//    const Size n_basis = n_basis_elem.flat_size();
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  // getting the 1D values for the test space -- begin
  std::array< ValueTable<Real>,dim> phi_1D_test;
  std::array< ValueTable<Real>,dim> grad_phi_1D_test;
  {
    const auto &ref_elem_accessor = elem_test.get_ref_space_accessor().get_bspline_accessor();

    const auto &quad_points = ref_elem_accessor.get_quad_points();

    const auto phi_1D_test_table =
      ref_elem_accessor.evaluate_univariate_derivatives_at_points(0,quad_points);

    const auto grad_phi_1D_test_table =
      ref_elem_accessor.evaluate_univariate_derivatives_at_points(1,quad_points);

    phi_1D_test = phi_1D_test_table[0]; // only valid for scalar spaces
    grad_phi_1D_test = grad_phi_1D_test_table[0];  // only valid for scalar spaces
  }
  // getting the 1D values for the test space -- end
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  // getting the 1D values for the trial space -- begin
  std::array< ValueTable<Real>,dim> phi_1D_trial;
  std::array< ValueTable<Real>,dim> grad_phi_1D_trial;
  {
    const auto &ref_elem_accessor = elem_trial.get_ref_space_accessor().get_bspline_accessor();

    const auto &quad_points = ref_elem_accessor.get_quad_points();

    const auto phi_1D_trial_table =
      ref_elem_accessor.evaluate_univariate_derivatives_at_points(0,quad_points);

    const auto grad_phi_1D_trial_table =
      ref_elem_accessor.evaluate_univariate_derivatives_at_points(1,quad_points);

    phi_1D_trial = phi_1D_trial_table[0]; // only valid for scalar spaces
    grad_phi_1D_trial = grad_phi_1D_trial_table[0];  // only valid for scalar spaces
  }
  // getting the 1D values for the trial space -- end
  //--------------------------------------------------------------------------

#ifdef TIME_PROFILING
  const auto end_initialization = Clock::now();
  const Duration elapsed_time_initialization = end_initialization - start_initialization;
  std::cout << "Elapsed_seconds initialization = " << elapsed_time_initialization.count() << std::endl;
#endif //#ifdef TIME_PROFILING
  //--------------------------------------------------------------------------






  //----------------------------------------------------
  // Coefficient evaluation phase -- begin
#ifdef TIME_PROFILING
  const TimePoint start_coefficient_evaluation = Clock::now();
#endif //#ifdef TIME_PROFILING


  // checks that the mapping used in the test space and in the trial space is the same
  Assert(elem_test.get_physical_space()->get_push_forward()->get_mapping() ==
         elem_trial.get_physical_space()->get_push_forward()->get_mapping(),
         ExcMessage("Test and trial spaces must have the same mapping (and the same grid)!"));


  // checks that the elements on the grid are the same
  Assert(static_cast<const GridElement<dim> &>(elem_test.get_ref_space_accessor()) ==
         static_cast<const GridElement<dim> &>(elem_trial.get_ref_space_accessor()),
         ExcMessage("Different elements for test space and trial space."));


  // performs the evaluation of the function DF^{-1} * C * DF^{-T} * det(DF) at the quadrature points
  const auto &det_DF = elem_test.get_measures() ;

  const auto &invDF = elem_test.get_push_forward_accessor().get_inv_gradients();

  Size n_points = coeffs.size();
  Assert(det_DF.size() == n_points,
         ExcDimensionMismatch(det_DF.size(), n_points));
  Assert(invDF.size() == n_points,
         ExcDimensionMismatch(invDF.size(), n_points));


  TensorSize<dim> n_points_1D = elem_test.get_ref_space_accessor().get_quad_points().get_num_points_direction();
  Assert(n_points_1D.flat_size() == n_points,
         ExcDimensionMismatch(n_points_1D.flat_size(),n_points));

//    LogStream out;
  std::vector<TMatrix<dim,dim>> C_hat(n_points);
  for (Index ipt = 0 ; ipt < n_points ; ++ipt)
  {
    TMatrix<dim,dim> &C_hat_ipt = C_hat[ipt];

    const TMatrix<space_dim,space_dim> C_ipt = coeffs[ipt] * det_DF[ipt];
    const auto &invDF_ipt = invDF[ipt];
//        out << "C at point " << ipt << "= " << C_ipt << endl;

    for (Index r = 0 ; r < dim ; ++r)
      for (Index s = 0 ; s < dim ; ++s)
        for (Index k = 0 ; k < space_dim ; ++k)
          for (Index l = 0 ; l < space_dim ; ++l)
            C_hat_ipt[r][s] += invDF_ipt[k][r](0) * C_ipt[k][l] * invDF_ipt[l][s](0);


//        out << "C_hat at point " << ipt << "= " << C_hat_ipt << endl;
  }


#ifdef TIME_PROFILING
  const TimePoint end_coefficient_evaluation = Clock::now();
  const Duration elapsed_time_coefficient_evaluation =
    end_coefficient_evaluation - start_coefficient_evaluation;
  std::cout << "Elapsed seconds coefficient evaluation stiffness = "
            << elapsed_time_coefficient_evaluation.count() << std::endl;
#endif //#ifdef TIME_PROFILING
  // Coefficient evaluation phase -- end
  //----------------------------------------------------






  //----------------------------------------------------
  // Assembly of the local stiffness matrix using sum-factorization -- begin
#ifdef TIME_PROFILING
  Duration elapsed_time_compute_phi1Dtest_phi1Dtrial;

  const auto start_sum_factorization = Clock::now();
#endif //#ifdef TIME_PROFILING

  TensorSize<3> tensor_size_C0;
  tensor_size_C0[0] = n_points_1D.flat_size(); // theta size
  tensor_size_C0[1] = 1; // alpha size
  tensor_size_C0[2] = 1; // beta size

  DynamicMultiArray<Real,3> C0(tensor_size_C0);
  const Size n_entries = tensor_size_C0.flat_size();

  DenseMatrix operator_gradu_gradv_tmp(n_basis_test.flat_size(),n_basis_trial.flat_size());
  operator_gradu_gradv.clear();

  std::array<DynamicMultiArray<Real,3>,dim> J;

  std::array< ValueTable<Real>,dim> trial_1D;
  std::array< ValueTable<Real>,dim> test_1D;
  for (Index k = 0 ; k < dim ; ++k)
  {
    for (Index l = 0 ; l < dim ; ++l)
    {
      for (Index i = 0 ; i < dim ; ++i)
      {
        if (i == k && i == l)
        {
          trial_1D[i] = grad_phi_1D_trial[i];
          test_1D[i] = grad_phi_1D_test [i];
        }
        else if (i == k && i != l)
        {
          trial_1D[i] = grad_phi_1D_trial[i];
          test_1D[i] =      phi_1D_test [i];
        }
        else if (i != k && i == l)
        {
          trial_1D[i] =      phi_1D_trial[i];
          test_1D[i] = grad_phi_1D_test [i];
        }
        else if (i != k && i != l)
        {
          trial_1D[i] = phi_1D_trial[i];
          test_1D[i] = phi_1D_test [i];
        }
      } // end loop i


      //----------------------------------------------------
      // precalculation of the J[i](theta_i,alpha_i,beta_i) terms
      // (i.e. the weigths[theta_i] * phi_trial[alpha_i] * phi_test[beta_i] )
#ifdef TIME_PROFILING
      const auto start_compute_phi1Dtest_phi1Dtrial = Clock::now();
#endif //#ifdef TIME_PROFILING

      const std::array<Real,dim> length_element_edge =
        elem_test.get_ref_space_accessor().get_coordinate_lengths();

      const auto J = evaluate_w_phi1Dtrial_phi1Dtest(
                       test_1D,
                       trial_1D,
                       elem_test.get_ref_space_accessor().get_quad_points().get_weights(),
                       length_element_edge);

#ifdef TIME_PROFILING
      const auto end_compute_phi1Dtest_phi1Dtrial = Clock::now();
      elapsed_time_compute_phi1Dtest_phi1Dtrial +=
        end_compute_phi1Dtest_phi1Dtrial- start_compute_phi1Dtest_phi1Dtrial;
#endif //#ifdef TIME_PROFILING
      //----------------------------------------------------



      Assert(n_entries == C_hat.size(),
             ExcDimensionMismatch(n_entries,C_hat.size()));
      for (Index entry_id = 0 ; entry_id < n_entries ; ++entry_id)
        C0[entry_id] = C_hat[entry_id][k][l];

      operator_gradu_gradv_tmp.clear();


      IntegratorSumFactorization<dim> integrate_sf;
      integrate_sf(false, //non symmetric
                   n_points_1D,
                   n_basis_trial,
                   n_basis_test,
                   J,
                   C0,
                   operator_gradu_gradv_tmp);

      operator_gradu_gradv += operator_gradu_gradv_tmp;

    }

  }



#ifdef TIME_PROFILING
  const auto end_sum_factorization = Clock::now();

  std::cout << "Elapsed seconds w * trial * test = "
            << elapsed_time_compute_phi1Dtest_phi1Dtrial.count() << std::endl;

  Duration elapsed_time_sum_factorization = end_sum_factorization - start_sum_factorization;
  std::cout << "Elapsed seconds sum-factorization = " << elapsed_time_sum_factorization.count() << std::endl;
  // Assembly of the local stiffness matrix using sum-factorization -- end
  //----------------------------------------------------


  const Duration elapsed_time_assemble = elapsed_time_sum_factorization +
                                         elapsed_time_coefficient_evaluation +
                                         elapsed_time_initialization ;
  std::cout << "Elapsed seconds assemblying = " << elapsed_time_assemble.count() << std::endl;


  const TimePoint end_assembly_stiffness_matrix = Clock::now();

  const_cast<Duration &>(this->elapsed_time_operator_gradu_gradv_)
    = end_assembly_stiffness_matrix - start_assembly_stiffness_matrix;
  std::cout << "Elapsed seconds operator gradu_gradv sum-factorization= "
            << this->elapsed_time_operator_gradu_gradv_.count() << std::endl;
#endif //#ifdef TIME_PROFILING
  // Assembly of the local stiffness matrix using sum-factorization -- end
  //----------------------------------------------------

#endif
//  Assert(false,ExcNotImplemented());
//  AssertThrow(false,ExcNotImplemented());
}


template <int dim_,int range_,int rank_>
inline
void
EllipticOperatorsSFIntegrationBSpline<dim_,range_,rank_>::
eval_operator_rhs_v(
  const ElemTest &elem_test,
  const ValueVector<typename BSpline<dim_,range_,rank_>::Value> &f,
  DenseVector &operator_rhs_v) const
{
  //TODO: (martinelli 22 Sep 2014): this function is not implemented using sum_factorization. Fix it!
  const Size n_basis_test  = elem_test .get_num_basis();

  const auto &phi_test  = elem_test.get_basis_values();
  const auto &w_meas  = elem_test.get_w_measures();


  Assert(operator_rhs_v.size() == n_basis_test,
         ExcDimensionMismatch(operator_rhs_v.size(),n_basis_test));


  const Size n_qp = f.get_num_points();
  Assert(n_qp == phi_test.get_num_points(),ExcDimensionMismatch(n_qp,phi_test.get_num_points()));

  std::vector<Real> f_times_w_meas(n_qp);
  for (int qp = 0; qp < n_qp; ++qp)
    f_times_w_meas[qp] = f[qp](0) * w_meas[qp];


  operator_rhs_v.clear();
  for (int i = 0; i < n_basis_test; ++i)
  {
    const auto phi_i = phi_test.get_function_view(i);
    for (int qp = 0; qp < n_qp; ++qp)
      operator_rhs_v(i) += phi_i[qp](0) * f_times_w_meas[qp];
  }
}


template <int dim_,int range_,int rank_>
inline
void
EllipticOperatorsSFIntegrationBSpline<dim_,range_,rank_>::
eval_operator_gradu_v(
  const ElemTest &elem_test,
  const ElemTrial &elem_trial,
  const ValueVector<typename BSpline<dim_,range_,rank_>::Gradient> &beta,
  DenseMatrix &operator_gradu_v) const
{

  //----------------------------------------------------
  // Assembly of the local convection matrix using sum-factorization -- begin
#ifdef TIME_PROFILING
  const TimePoint start_assembly_convection_matrix = Clock::now();

  //--------------------------------------------------------------------------
  const auto start_initialization = Clock::now();
#endif //#ifdef TIME_PROFILING




  //--------------------------------------------------------------------------
  // getting the number of basis along each coordinate direction for the test and trial space

  const Index comp = 0; // only scalar spaces for the moment

  // test space -- begin
  TensorIndex<dim> degree_test = elem_test.get_physical_space()->get_reference_space()->get_degree()[comp];
  TensorSize<dim> n_basis_elem_test(degree_test + 1);

  Assert(n_basis_elem_test.flat_size()==elem_test.get_num_basis(),
         ExcDimensionMismatch(n_basis_elem_test.flat_size(),elem_test.get_num_basis()));

  const auto weight_basis_test = MultiArrayUtils<dim>::compute_weight(n_basis_elem_test);

  const TensorSize<dim> n_basis_test = n_basis_elem_test;
  // test space -- end


  // trial space -- begin
  TensorIndex<dim> degree_trial = elem_trial.get_physical_space()->get_reference_space()->get_degree()[comp];
  TensorSize<dim> n_basis_elem_trial(degree_trial + 1);

  Assert(n_basis_elem_trial.flat_size()==elem_trial.get_num_basis(),
         ExcDimensionMismatch(n_basis_elem_trial.flat_size(),elem_trial.get_num_basis()));

  const auto weight_basis_trial = MultiArrayUtils<dim>::compute_weight(n_basis_elem_trial);

  const TensorSize<dim> n_basis_trial = n_basis_elem_trial;
  // trial space -- end

//    const Size n_basis = n_basis_elem.flat_size();
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  // getting the 1D values for the test space -- begin
  std::array< ValueTable<Real>,dim> phi_1D_test;
  {
    const auto &ref_elem_accessor = elem_test.get_ref_space_accessor().get_bspline_accessor();
    const auto &quad_points = ref_elem_accessor.get_quad_points();

    const auto phi_1D_test_table =
      ref_elem_accessor.evaluate_univariate_derivatives_at_points(0,quad_points);

    phi_1D_test = phi_1D_test_table[0]; // only valid for scalar spaces
  }
  // getting the 1D values for the test space -- end
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  // getting the 1D values for the trial space -- begin
  std::array< ValueTable<Real>,dim> phi_1D_trial;
  std::array< ValueTable<Real>,dim> grad_phi_1D_trial;
  {
    const auto &ref_elem_accessor = elem_trial.get_ref_space_accessor().get_bspline_accessor();

    const auto &quad_points = ref_elem_accessor.get_quad_points();

    const auto phi_1D_trial_table =
      ref_elem_accessor.evaluate_univariate_derivatives_at_points(0,quad_points);

    const auto grad_phi_1D_trial_table =
      ref_elem_accessor.evaluate_univariate_derivatives_at_points(1,quad_points);

    phi_1D_trial = phi_1D_trial_table[0]; // only valid for scalar spaces
    grad_phi_1D_trial = grad_phi_1D_trial_table[0];  // only valid for scalar spaces
  }
  // getting the 1D values for the trial space -- end
  //--------------------------------------------------------------------------

#ifdef TIME_PROFILING
  const auto end_initialization = Clock::now();
  const Duration elapsed_time_initialization = end_initialization - start_initialization;
  std::cout << "Elapsed_seconds initialization = " << elapsed_time_initialization.count() << std::endl;
#endif //#ifdef TIME_PROFILING
  //--------------------------------------------------------------------------



  //----------------------------------------------------
  // Coefficient evaluation phase -- begin
#ifdef TIME_PROFILING
  const TimePoint start_coefficient_evaluation = Clock::now();
#endif //#ifdef TIME_PROFILING


  // checks that the mapping used in the test space and in the trial space is the same
  Assert(elem_test.get_physical_space()->get_push_forward()->get_mapping() ==
         elem_trial.get_physical_space()->get_push_forward()->get_mapping(),
         ExcMessage("Test and trial spaces must have the same mapping (and the same grid)!"));


  // checks that the elements on the grid are the same
  Assert(static_cast<const GridElement<dim> &>(elem_test.get_ref_space_accessor()) ==
         static_cast<const GridElement<dim> &>(elem_trial.get_ref_space_accessor()),
         ExcMessage("Different elements for test space and trial space."));


  // performs the evaluation of the function DF^{-1} * beta * det(DF) at the quadrature points
  const auto &det_DF = elem_test.get_measures() ;

  const auto &invDF = elem_test.get_push_forward_accessor().get_inv_gradients();

  Size n_points = beta.size();
  Assert(det_DF.size() == n_points,
         ExcDimensionMismatch(det_DF.size(), n_points));
  Assert(invDF.size() == n_points,
         ExcDimensionMismatch(invDF.size(), n_points));


  TensorSize<dim> n_points_1D = elem_test.get_ref_space_accessor().get_quad_points().get_num_points_direction();
  Assert(n_points_1D.flat_size() == n_points,
         ExcDimensionMismatch(n_points_1D.flat_size(),n_points));

//    LogStream out;
  ValueVector<typename BSpline<dim_,range_,rank_>::Gradient> vel_hat(n_points);
  for (Index ipt = 0 ; ipt < n_points ; ++ipt)
  {
    auto &vel_hat_ipt = vel_hat[ipt];

    const auto vel_ipt = beta[ipt] * det_DF[ipt];
    const auto &invDF_ipt = invDF[ipt];

    for (Index r = 0 ; r < dim ; ++r)
      for (Index k = 0 ; k < space_dim ; ++k)
        vel_hat_ipt[r](0) += invDF_ipt[k][r](0) * vel_ipt[k](0);
  }


#ifdef TIME_PROFILING
  const TimePoint end_coefficient_evaluation = Clock::now();
  const Duration elapsed_time_coefficient_evaluation =
    end_coefficient_evaluation - start_coefficient_evaluation;
  std::cout << "Elapsed seconds coefficient evaluation stiffness = "
            << elapsed_time_coefficient_evaluation.count() << std::endl;
#endif //#ifdef TIME_PROFILING
  // Coefficient evaluation phase -- end
  //----------------------------------------------------






  //----------------------------------------------------
  // Assembly of the local stiffness matrix using sum-factorization -- begin
#ifdef TIME_PROFILING
  Duration elapsed_time_compute_phi1Dtest_phi1Dtrial;

  const auto start_sum_factorization = Clock::now();
#endif //#ifdef TIME_PROFILING

  TensorSize<3> tensor_size_C0;
  tensor_size_C0[0] = n_points_1D.flat_size(); // theta size
  tensor_size_C0[1] = 1; // alpha size
  tensor_size_C0[2] = 1; // beta size

  DynamicMultiArray<Real,3> C0(tensor_size_C0);
  const Size n_entries = tensor_size_C0.flat_size();

  DenseMatrix operator_gradu_v_tmp(n_basis_test.flat_size(),n_basis_trial.flat_size());
  operator_gradu_v.clear();

  std::array<DynamicMultiArray<Real,3>,dim> J;

  std::array< ValueTable<Real>,dim> trial_1D;
  std::array< ValueTable<Real>,dim> test_1D;
  for (Index l = 0 ; l < dim ; ++l)
  {
    for (Index i = 0 ; i < dim ; ++i)
    {
      test_1D[i] = phi_1D_test[i];

      if (i == l)
        trial_1D[i] = grad_phi_1D_trial[i];
      else // (i != l)
        trial_1D[i] = phi_1D_trial[i];
    } // end loop i


    //----------------------------------------------------
    // precalculation of the J[i](theta_i,alpha_i,beta_i) terms
    // (i.e. the weigths[theta_i] * phi_trial[alpha_i] * phi_test[beta_i] )
#ifdef TIME_PROFILING
    const auto start_compute_phi1Dtest_phi1Dtrial = Clock::now();
#endif //#ifdef TIME_PROFILING

    const std::array<Real,dim> length_element_edge =
      elem_test.get_ref_space_accessor().get_coordinate_lengths();

    const auto J = evaluate_w_phi1Dtrial_phi1Dtest(
                     test_1D,
                     trial_1D,
                     elem_test.get_ref_space_accessor().get_quad_points().get_weights(),
                     length_element_edge);

#ifdef TIME_PROFILING
    const auto end_compute_phi1Dtest_phi1Dtrial = Clock::now();
    elapsed_time_compute_phi1Dtest_phi1Dtrial +=
      end_compute_phi1Dtest_phi1Dtrial- start_compute_phi1Dtest_phi1Dtrial;
#endif //#ifdef TIME_PROFILING
    //----------------------------------------------------



    Assert(n_entries == vel_hat.size(),
           ExcDimensionMismatch(n_entries,vel_hat.size()));
    for (Index entry_id = 0 ; entry_id < n_entries ; ++entry_id)
      C0[entry_id] = vel_hat[entry_id][l](0);

    operator_gradu_v_tmp.clear();


    IntegratorSumFactorization<dim> integrate_sf;
    integrate_sf(false, //non symmetric
                 n_points_1D,
                 n_basis_trial,
                 n_basis_test,
                 J,
                 C0,
                 operator_gradu_v_tmp);

    operator_gradu_v += operator_gradu_v_tmp;

  }



#ifdef TIME_PROFILING
  const auto end_sum_factorization = Clock::now();

  std::cout << "Elapsed seconds w * trial * test = "
            << elapsed_time_compute_phi1Dtest_phi1Dtrial.count() << std::endl;

  Duration elapsed_time_sum_factorization = end_sum_factorization - start_sum_factorization;
  std::cout << "Elapsed seconds sum-factorization = " << elapsed_time_sum_factorization.count() << std::endl;
  // Assembly of the local convection matrix using sum-factorization -- end
  //----------------------------------------------------


  const Duration elapsed_time_assemble = elapsed_time_sum_factorization +
                                         elapsed_time_coefficient_evaluation +
                                         elapsed_time_initialization ;
  std::cout << "Elapsed seconds assemblying = " << elapsed_time_assemble.count() << std::endl;


  const TimePoint end_assembly_convection_matrix = Clock::now();

  const_cast<Duration &>(this->elapsed_time_operator_gradu_v_)
    = end_assembly_convection_matrix - start_assembly_convection_matrix;
  std::cout << "Elapsed seconds operator gradu_v sum-factorization= "
            << this->elapsed_time_operator_gradu_v_.count() << std::endl;
#endif //#ifdef TIME_PROFILING
  // Assembly of the local convection matrix using sum-factorization -- end
  //----------------------------------------------------


//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
}





IGA_NAMESPACE_CLOSE


#endif // #ifndef ELLIPTIC_OPERATORS_SF_INTEGRATION_H_
