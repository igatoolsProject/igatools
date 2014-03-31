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
// [old includes]
#include <igatools/base/function_lib.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/linear_solver.h>
#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/io/writer.h>
#include <igatools/utils/value_table.h>
#include <igatools/utils/multi_array_utils.h>
#include <igatools/contrib/table_handler.h>

#include <boost/math/special_functions/binomial.hpp>

#include <cblas.h>

#include <chrono>
// [old includes]

// [unqualified names]
using namespace iga;
using namespace std;
using namespace std::chrono;

using functions::ConstantFunction;
using space_tools::project_boundary_values;
using dof_tools::apply_boundary_values;
using dof_tools::get_sparsity_pattern;
using numbers::PI;
// [unqualified names]


// [Problem class]
template<int dim,class DerivedClass>
class PoissonProblem
{
public:
    PoissonProblem(const TensorSize<dim> &n_knots, const int deg);
    void run();

private:
//    void assemble();
    void solve();
    void output();
    // [Problem class]

    // [type aliases]
protected:
    using RefSpace  = BSplineSpace<dim>;
    using PushFw    = PushForward<Transformation::h_grad, dim>;
    using Space     = PhysicalSpace<RefSpace, PushFw>;
//    using Value     = typename Function<dim>::Value;
    // [type aliases]

    shared_ptr<Mapping<dim>> map;
    shared_ptr<Space>        space;

    const Quadrature<dim>   elem_quad;
    const Quadrature<dim-1> face_quad;

    const boundary_id dir_id = 0;

    std::shared_ptr<Matrix> matrix;
    std::shared_ptr<Vector> rhs;
    std::shared_ptr<Vector> solution;
};



template<int dim,class DerivedClass>
PoissonProblem<dim,DerivedClass>::
PoissonProblem(const TensorSize<dim> &n_knots, const int deg)
    :
    elem_quad(QGauss<dim>(deg+1)),
    face_quad(QGauss<dim-1>(deg+1))
{
    BBox<dim> box;
    for (int i=0 ; i < dim ; ++i)
    {
        box[i][0] = 0.0;
        box[i][1] = 1.0;
    }
    box[0][1] = 0.5;

    auto grid = CartesianGrid<dim>::create(box,n_knots);
    auto ref_space = RefSpace::create(grid, deg);
//    map       = BallMapping<dim>::create(grid);
    map       = IdentityMapping<dim,0>::create(grid);
    space     = Space::create(ref_space, PushFw::create(map));

    const auto n_basis = space->get_num_basis();
    matrix   = Matrix::create(get_sparsity_pattern(*space));
    rhs      = Vector::create(n_basis);
    solution = Vector::create(n_basis);
}






template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
solve()
{
    LinearSolver solver(LinearSolver::Type::CG);
    solver.solve(*matrix, *rhs, *solution);
}



template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
output()
{
    const int n_plot_points = 2;
    Writer<dim> writer(map, n_plot_points);

    writer.add_field(space, *solution, "solution");
    string filename = "poisson_problem-" + to_string(dim) + "d" ;
    writer.save(filename);
}



template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
run()
{
    static_cast<DerivedClass &>(*this).assemble();
    solve();
    output();
}



template<int dim>
class PoissonProblemStandardIntegration :
    public PoissonProblem< dim, PoissonProblemStandardIntegration<dim> >
{
public:
    // inheriting the constructors from PoissonProblem
    using PoissonProblem< dim, PoissonProblemStandardIntegration<dim> >::PoissonProblem;

    void assemble();

private:
    using base_t = PoissonProblem< dim, PoissonProblemStandardIntegration<dim> >;
    using typename base_t::Space;
};


template<int dim>
void
PoissonProblemStandardIntegration<dim>::
assemble()
{
    const int n_basis = this->space->get_num_basis_per_element();
    DenseMatrix loc_mat(n_basis, n_basis);
    DenseVector loc_rhs(n_basis);
    vector<Index> loc_dofs(n_basis);

    const int n_qp = this->elem_quad.get_num_points();
    ConstantFunction<dim> f({0.5});
    vector< typename Function<dim>::ValueType > f_values(n_qp);

    auto elem = this->space->begin();
    const auto elem_end = this->space->end();
    ValueFlags fill_flags = ValueFlags::value | ValueFlags::gradient |
                            ValueFlags::w_measure | ValueFlags::point;
    elem->init_values(fill_flags, this->elem_quad);

    for (; elem != elem_end; ++elem)
    {
        loc_mat.clear();
        loc_rhs.clear();
        elem->fill_values();

        auto points  = elem->get_points();
        auto phi     = elem->get_basis_values();
        auto grd_phi = elem->get_basis_gradients();
        auto w_meas  = elem->get_w_measures();

        f.evaluate(points, f_values);

        for (int i = 0; i < n_basis; ++i)
        {
            auto grd_phi_i = grd_phi.get_function_view(i);
            for (int j = 0; j < n_basis; ++j)
            {
                auto grd_phi_j = grd_phi.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                    loc_mat(i,j) +=
                        scalar_product(grd_phi_i[qp], grd_phi_j[qp])
                        * w_meas[qp];
            }

            auto phi_i = phi.get_function_view(i);
            for (int qp = 0; qp < n_qp; ++qp)
                loc_rhs(i) += scalar_product(phi_i[qp], f_values[qp])
                              * w_meas[qp];
        }

        loc_dofs = elem->get_local_to_global();
        this->matrix->add_block(loc_dofs, loc_dofs, loc_mat);
        this->rhs->add_block(loc_dofs, loc_rhs);
    }

    this->matrix->fill_complete();

    // [dirichlet constraint]
    ConstantFunction<dim> g({0.0});
    std::map<Index, Real> values;
    project_boundary_values<Space>(g, this->space, this->face_quad, this->dir_id, values);
    apply_boundary_values(values, *this->matrix, *this->rhs, *this->solution);
    // [dirichlet constraint]
}




template<int dim>
class PoissonProblemSumFactorization :
    public PoissonProblem< dim, PoissonProblemSumFactorization<dim> >
{
public:

    /** @name Constructors & destructor. */
    ///@{
    PoissonProblemSumFactorization() = delete ;

    /**
     * Constructor.
     * @param[in] n_knots Number of knots along each direction.
     * @param[in] space_deg Polynomial degree for the solution space.
     * @param[in] proj_deg Polynomial degree for the projection space.
     */
    PoissonProblemSumFactorization(const TensorSize<dim> &n_knots,
                                   const int space_deg,const int proj_deg);

    /** Copy constructor. Not allowed to be used. */
    PoissonProblemSumFactorization(const PoissonProblemSumFactorization<dim> &in) = delete;

    /** Move constructor. Not allowed to be used. */
    PoissonProblemSumFactorization(PoissonProblemSumFactorization<dim> &&in) = delete;

    /** Destructor. */
    ~PoissonProblemSumFactorization() = default;
    ///@}


    Real get_elapsed_time_projection() const
    {
        return elapsed_time_projection_.count();
    }

    Real get_elapsed_time_assembly_mass_matrix() const
    {
        return elapsed_time_assembly_mass_matrix_.count();
    }

    Real get_elapsed_time_assembly_mass_matrix_old() const
    {
        return elapsed_time_assembly_mass_matrix_old_.count();
    }


    void assemble();

private:
    using base_t = PoissonProblem< dim, PoissonProblemSumFactorization<dim> >;

    int space_deg_ = 0;
    int proj_deg_ = 0;

    /** Number of one-dimensional Bernstein's polynomials used in the projection */
    Size n_bernst_1D_ = 0 ;

    /** Number of d-dimensional Bernstein's polynomials used in the projection */
    Size n_basis_proj_ = 0;

    /** Number of points for the one-dimensional quadrature scheme used for projection. */
    Size n_quad_proj_1D_ = 0;


    /** Index weight used to convert tensor-to-flat id (and viceversa) for the Bernstein polynomials.*/
    TensorIndex<dim> bernst_tensor_weight_;

    /** Quadrature scheme used in the projection phase. */
    QGauss<dim> quad_proj_;



    chrono::duration<Real> elapsed_time_projection_;
    chrono::duration<Real> elapsed_time_assembly_mass_matrix_;

    chrono::duration<Real> elapsed_time_assembly_mass_matrix_old_;
    chrono::duration<Real> elapsed_time_compute_I_;

    DenseMatrix eval_matrix_proj() const;

    /** Mass matrix of the Bernstein polynomials used in the projection, in [0,1]^dim.*/
    DenseMatrix B_proj_;

    /** Inverse of the mass matrix of the Bernstein polynomials used in the projection, in [0,1]^dim.*/
    DenseMatrix inv_B_proj_;

};

template<int dim>
PoissonProblemSumFactorization<dim>::
PoissonProblemSumFactorization(const TensorSize<dim> &n_knots,
                               const int space_deg,const int proj_deg)
    :
    base_t(n_knots,space_deg),
    space_deg_(space_deg),
    proj_deg_(proj_deg),
    n_bernst_1D_(proj_deg_+1),
    n_basis_proj_(pow(n_bernst_1D_,dim)),
    n_quad_proj_1D_(proj_deg_+1),
    bernst_tensor_weight_(MultiArrayUtils<dim>::compute_weight(TensorSize<dim>(n_bernst_1D_))),
    quad_proj_(QGauss<dim>(TensorSize<dim>(n_bernst_1D_)))
{
    using boost::math::binomial_coefficient;

    //------------------------------------
    // evaluation of the Bernstein's mass matrix
    DenseMatrix B_proj_1D(n_bernst_1D_,n_bernst_1D_);

    const int q = proj_deg_;
    for (int i = 0 ; i < n_bernst_1D_ ; ++i)
    {
        const Real tmp = binomial_coefficient<double>(q,i) / (2.0 * q + 1.0);

        for (int j = 0 ; j < n_bernst_1D_ ; ++j)
        {
            B_proj_1D(i,j) =
                tmp / binomial_coefficient<double>(2*q,i+j) * binomial_coefficient<double>(q,j) ;
        }
    }


    using MAUtils = MultiArrayUtils<dim>;
    B_proj_.resize(n_basis_proj_,n_basis_proj_);


    for (Size row = 0 ; row < n_basis_proj_ ; ++row)
    {
        const auto row_tensor_id = MAUtils::flat_to_tensor_index(row,bernst_tensor_weight_);

        for (Size col = 0 ; col < n_basis_proj_ ; ++col)
        {
            const auto col_tensor_id = MAUtils::flat_to_tensor_index(col,bernst_tensor_weight_);

            Real tmp = 1.0;
            for (int k = 0 ; k < dim ; ++k)
                tmp *= B_proj_1D(row_tensor_id[k],col_tensor_id[k]);

            B_proj_(row,col) = tmp;
        }
    }


    inv_B_proj_ = B_proj_.inverse();
    //------------------------------------

    /*
    LogStream out;
    out << "Bernstein "<< dim <<"D mass matrix = " << B_proj_ << endl;
    out << "inverse of Bernstein "<< dim <<"D mass matrix = " << inv_B_proj_ << endl;
    //*/
}



template <int dim,int k=0>
class IntegratorSumFacRHS
{
public:
    DynamicMultiArray<Real,k> operator()(
        const DynamicMultiArray<Real,dim> &Gamma,
        const array<ValueTable<typename Function<1>::ValueType>,dim> &omega_B_1d,
        const TensorIndex<dim> &bernstein_tensor_id)
    {

        IntegratorSumFacRHS<dim,k+1> integrate ;
        DynamicMultiArray<Real,k+1> C1 = integrate(Gamma,omega_B_1d,bernstein_tensor_id);
        const auto tensor_size_C1 = C1.tensor_size();

        TensorSize<k> tensor_size_C;
        for (int i = 0 ; i < k ; ++i)
        {
            Assert(tensor_size_C1[i] == tensor_size_C1[i+1],
                   ExcDimensionMismatch(tensor_size_C1[i],tensor_size_C1[i+1]));
            tensor_size_C[i] = tensor_size_C1[i];
        }

        DynamicMultiArray<Real,k> C(tensor_size_C) ;

        TensorIndex<k+1> tensor_id_C1;

        const Size n_entries_C = C.flat_size();

        const Index eta_k1 = bernstein_tensor_id[k] ;

        const auto &w_B_k_view = omega_B_1d[k].get_function_view(eta_k1);
        vector<Real> w_B_k;
        for (const auto & w_B_k_value : w_B_k_view)
            w_B_k.emplace_back(w_B_k_value(0));

        const Size n_pts = w_B_k.size();
        for (Index flat_id_C = 0 ; flat_id_C < n_entries_C ; ++flat_id_C)
        {
            TensorIndex<k> tensor_id_C = C.flat_to_tensor(flat_id_C);

            for (int i = 0 ; i < k ; ++i)
                tensor_id_C1[i] = tensor_id_C[i];

            Real sum = 0.0;
            for (int ipt = 0 ; ipt < n_pts ; ++ipt) // quadrature along the last tensor index
            {
                tensor_id_C1[k] = ipt ;
                sum += C1(tensor_id_C1) * w_B_k[ipt];
            }
            C(tensor_id_C) = sum;

//          cout << "sum=" << sum << endl;
        }

        /*
        LogStream out ;
        out << "C="<< C << endl;
        cout << "foo<"<<k<<">" << endl;
        //*/
        return C;
    }
};

template <int dim>
class IntegratorSumFacRHS<dim,dim>
{
public:
    DynamicMultiArray<Real,dim> operator()(
        const DynamicMultiArray<Real,dim> &Gamma,
        const array<ValueTable<typename Function<1>::ValueType>,dim> &omega_B_1d,
        const TensorIndex<dim> &bernstein_tensor_id)
    {
        return Gamma;
    }
};




template <int dim, int r=dim>
class MassMatrixIntegrator
{
public:
	DynamicMultiArray<Real,3>
	operator()(
			const TensorSize<dim> &tensor_size_alphabeta,
			const TensorSize<dim> &tensor_size_theta,
			const array<DynamicMultiArray<Real,3>,dim> &J,
			const DynamicMultiArray<Real,3> &C_k_1) const
	{
	const int k = dim-r+1;

    // (alpha_1,...alpha_{k-1})
    TensorSize<k-1> tensor_size_alphabeta_k_1;
    for (int i = 0 ; i < k-1 ; ++i)
        tensor_size_alphabeta_k_1[i] = tensor_size_alphabeta[i];


    // (alpha_1,...alpha_k)
    TensorSize<k> tensor_size_alphabeta_k;
    for (int i = 0 ; i < k ; ++i)
        tensor_size_alphabeta_k  [i] = tensor_size_alphabeta[i];


    // (theta_{k+1},...theta_dim)
    TensorSize<dim-k> tensor_size_theta_k_1;
    for (int i = 0 ; i < dim-k ; ++i)
        tensor_size_theta_k_1[i] = tensor_size_theta[i+k];


    // (theta_k,...theta_dim)
    TensorSize<dim-k+1> tensor_size_theta_k;
    for (int i = 0 ; i <= dim-k ; ++i)
        tensor_size_theta_k[i] = tensor_size_theta[i+k-1];





//    TensorSize<3> tensor_size_C_k_1 = C_k_1.tensor_size();
    Assert(C_k_1.tensor_size()[0] == tensor_size_theta_k.flat_size(),
           ExcDimensionMismatch(C_k_1.tensor_size()[0],tensor_size_theta_k.flat_size()));
    Assert(C_k_1.tensor_size()[1] == tensor_size_alphabeta_k_1.flat_size(),
           ExcDimensionMismatch(C_k_1.tensor_size()[1],tensor_size_alphabeta_k_1.flat_size()));
    Assert(C_k_1.tensor_size()[2] == tensor_size_alphabeta_k_1.flat_size(),
           ExcDimensionMismatch(C_k_1.tensor_size()[2],tensor_size_alphabeta_k_1.flat_size()));


    const auto &J_k = J[k-1];
//    TensorSize<3> tensor_size_J_k = J_k.tensor_size();
    Assert(J_k.tensor_size()[0] == tensor_size_theta[k-1],
           ExcDimensionMismatch(J_k.tensor_size()[0],tensor_size_theta[k-1]));
    Assert(J_k.tensor_size()[1] == tensor_size_alphabeta[k-1],
           ExcDimensionMismatch(J_k.tensor_size()[1],tensor_size_alphabeta[k-1]));
    Assert(J_k.tensor_size()[2] == tensor_size_alphabeta[k-1],
           ExcDimensionMismatch(J_k.tensor_size()[2],tensor_size_alphabeta[k-1]));




    TensorIndex<3> tensor_id_J_k;
    TensorIndex<3> tensor_id_C_k_1;
    TensorIndex<3> tensor_id_C_k;

    const Size size_flat_theta_k_1 = (dim-k>0)?tensor_size_theta_k_1.flat_size():1;
    const Size size_flat_alpha_k_1 = tensor_size_alphabeta_k_1.flat_size();
    const Size size_flat_beta_k_1  = size_flat_alpha_k_1;


    TensorSize<3> tensor_size_C_k;
    tensor_size_C_k[0] = size_flat_theta_k_1;
    tensor_size_C_k[1] = tensor_size_alphabeta_k.flat_size();
    tensor_size_C_k[2] = tensor_size_alphabeta_k.flat_size();
    DynamicMultiArray<Real,3> C_k(tensor_size_C_k);


#define OPTIMIZED

#ifndef OPTIMIZED

    tensor_id_C_k[2] = 0;
    for (Index flat_beta_k_1 = 0 ; flat_beta_k_1 < size_flat_beta_k_1 ; ++flat_beta_k_1)
    {
        tensor_id_C_k_1[2] = flat_beta_k_1;


        for (int beta_k = 0 ; beta_k < tensor_size_alphabeta[k-1] ; ++beta_k)
        {
            tensor_id_J_k[2] = beta_k;
            tensor_id_C_k[1] = 0 ;

            for (Index flat_alpha_k_1 = 0 ; flat_alpha_k_1 < size_flat_alpha_k_1 ; ++flat_alpha_k_1)
            {
                tensor_id_C_k_1[1] = flat_alpha_k_1;

                for (int alpha_k = 0 ; alpha_k < tensor_size_alphabeta[k-1] ; ++alpha_k)
                {
                    tensor_id_J_k[1] = alpha_k;

                    tensor_id_C_k_1[0] = 0 ;
                    for (Index flat_theta_k_1 = 0 ; flat_theta_k_1 < size_flat_theta_k_1 ; ++flat_theta_k_1)
                    {
                        tensor_id_C_k[0] = flat_theta_k_1;

                        Real sum = 0.0;
                        for (int theta_k = 0 ; theta_k < tensor_size_theta[k-1] ; ++theta_k)
                        {

                            tensor_id_J_k[0] = theta_k;
                            sum += C_k_1(tensor_id_C_k_1) * J_k(tensor_id_J_k);

                            tensor_id_C_k_1[0]++;

                        } // end loop theta_k

                        C_k(tensor_id_C_k) = sum;

                    } //end loop flat_theta_k_1

                    tensor_id_C_k[1]++;

                } // end loop alpha_k
            } // end loop flat_alpha_k_1

            tensor_id_C_k[2]++;
        } // end loop beta_k

    } // end loop flat_beta_k_1

#else
    // OPTIMIZED branch

    tensor_id_C_k[2] = 0;
    for (Index flat_beta_k_1 = 0 ; flat_beta_k_1 < size_flat_beta_k_1 ; ++flat_beta_k_1)
    {
        tensor_id_C_k_1[2] = flat_beta_k_1;


        for (int beta_k = 0 ; beta_k < tensor_size_alphabeta[k-1] ; ++beta_k)
        {
            tensor_id_J_k[2] = beta_k;
            tensor_id_C_k[1] = 0 ;

            for (Index flat_alpha_k_1 = 0 ; flat_alpha_k_1 < size_flat_alpha_k_1 ; ++flat_alpha_k_1)
            {
                tensor_id_C_k_1[1] = flat_alpha_k_1;

                for (int alpha_k = 0 ; alpha_k < tensor_size_alphabeta[k-1] ; ++alpha_k)
                {
                    tensor_id_J_k[1] = alpha_k;

                    tensor_id_C_k_1[0] = 0 ;
                    const Real *C_k_1_ptr = &C_k_1(tensor_id_C_k_1);
                    for (Index flat_theta_k_1 = 0 ; flat_theta_k_1 < size_flat_theta_k_1 ; ++flat_theta_k_1)
                    {
                        tensor_id_C_k[0] = flat_theta_k_1;


                        tensor_id_J_k[0] = 0;
                        const Real *J_k_ptr = &J_k(tensor_id_J_k);
                        /*
                        for ( ; J_k_ptr != J_k_ptr_end ; ++J_k_ptr, ++C_k_1_ptr)
//                        for (int theta_k = 0 ; theta_k < tensor_size_theta[k-1] ; ++theta_k, ++J_k_ptr)
                        {

//                            tensor_id_J_k[0] = theta_k;
                            sum += (*C_k_1_ptr) * (*J_k_ptr);

//                            tensor_id_C_k_1[0]++;

                        } // end loop theta_k
//*/
                        C_k(tensor_id_C_k) = cblas_ddot(tensor_size_theta[k-1], C_k_1_ptr, 1, J_k_ptr, 1);

                    } //end loop flat_theta_k_1

                    tensor_id_C_k[1]++;

                } // end loop alpha_k
            } // end loop flat_alpha_k_1

            tensor_id_C_k[2]++;
        } // end loop beta_k

    } // end loop flat_beta_k_1


#endif


    MassMatrixIntegrator<dim,r-1> mass_matrix_integrator;
    return mass_matrix_integrator(tensor_size_alphabeta,tensor_size_theta,J,C_k);
}
};


template <int dim>
class MassMatrixIntegrator<dim,0>
{
public:
	DynamicMultiArray<Real,3> operator()(
			const TensorSize<dim> &tensor_size_alphabeta,
			const TensorSize<dim> &tensor_size_theta,
			const array<DynamicMultiArray<Real,3>,dim> &J,
			const DynamicMultiArray<Real,3> &C) const
	{
		return C;
	}
};

template<int dim>
void
PoissonProblemSumFactorization<dim>::
assemble()
{
//    using RefSpaceAccessor = typename base_t::RefSpace::ElementAccessor;

    /*
        string filename;
    #ifndef OPTIMIZED
        filename = "matrix_diff_orig.txt";
    #else
        filename = "matrix_diff_optim.txt";
    #endif

        ofstream out(filename);
        //*/

    LogStream out;
    //*/
    using MAUtils = MultiArrayUtils<dim>;


    const int n_basis = this->space->get_num_basis_per_element();
    DenseMatrix loc_mat(n_basis, n_basis);
    DenseVector loc_rhs(n_basis);
    vector<Index> loc_dofs(n_basis);

    ConstantFunction<dim> f({1.0});



    //-----------------------------------------------------------------
    /**
     * Initialization of the container for the values of the function
     * that must be projected.
     */
    using ValueType = typename Function<dim>::ValueType;
    vector<ValueType> f_values_proj(quad_proj_.get_num_points());
    //-----------------------------------------------------------------




    auto elem = this->space->begin();
    const auto elem_end = this->space->end();
    ValueFlags fill_flags = ValueFlags::value |
                            ValueFlags::gradient |
                            ValueFlags::measure |
                            ValueFlags::w_measure |
                            ValueFlags::point;

    elem->init_values(fill_flags, this->elem_quad);



    //-----------------------------------------------------------------
    std::array<DenseMatrix,dim> B_proj_1D;
    for (int i = 0 ; i < dim ; ++i)
    {
        // values of the one-dimensional Bernstein basis at the projection pts
        B_proj_1D[i] = BernsteinBasis::evaluate(
                           proj_deg_,
                           quad_proj_.get_points().get_data_direction(i));
        Assert(n_bernst_1D_ == B_proj_1D[i].size1(),
               ExcDimensionMismatch(n_bernst_1D_, B_proj_1D[i].size1()));
        Assert(quad_proj_.get_num_points_direction()[i] == B_proj_1D[i].size2(),
               ExcDimensionMismatch(quad_proj_.get_num_points_direction()[i], B_proj_1D[i].size2()));
    }


    using ValueType1D = Function<1>::ValueType;
    array<ValueTable<ValueType1D>,dim> w_B_proj_1D;

    const auto &ref_elem_accessor = elem->get_ref_space_accessor();
    array<Real,dim> element_edge_size = ref_elem_accessor.get_coordinate_lengths();

    for (int i = 0 ; i < dim ; ++i)
    {
        const Size n_pts = quad_proj_.get_num_points_direction()[i];
        w_B_proj_1D[i].resize(n_bernst_1D_,n_pts);

        const auto &B_i = B_proj_1D[i];

        const vector<Real> w_i = quad_proj_.get_weights().get_data_direction(i);

        for (Index ifn = 0 ; ifn < n_bernst_1D_ ; ++ifn)
        {
            auto w_B_ifn = w_B_proj_1D[i].get_function_view(ifn);

            for (Index jpt = 0 ; jpt < n_pts ; ++jpt)
                w_B_ifn[jpt] = w_i[jpt] * B_i(ifn,jpt) * element_edge_size[i];
        }
    }
    //-----------------------------------------------------------------




    IntegratorSumFacRHS<dim> integrate_rhs;

    TensorIndex<dim> bernst_tensor_id;

    DenseVector integral_rhs(n_basis_proj_);
    DenseVector k_rhs(n_basis_proj_);

    // number of points along each directin for the quadrature scheme.
    TensorSize<dim> n_quad_points = this->elem_quad.get_num_points_direction();

    TensorSize<dim> n_basis_elem(space_deg_+1);



    using Clock = chrono::high_resolution_clock;
    using TimePoint = chrono::time_point<Clock>;

    TimePoint start_projection;
    TimePoint   end_projection;

    TimePoint start_assembly_mass_matrix;
    TimePoint   end_assembly_mass_matrix;

    TimePoint start_assembly_mass_matrix_old;
    TimePoint   end_assembly_mass_matrix_old;

    DenseMatrix loc_mass_matrix_sf(n_basis,n_basis);

    const auto weight_basis = MultiArrayUtils<dim>::compute_weight(n_basis_elem);

    for (; elem != elem_end; ++elem)
    {


        loc_rhs.clear();
        //*/
        elem->fill_values();


        //--------------------------------------------------------------------------
        const auto &ref_elem_accessor = elem->get_ref_space_accessor();

        const Index comp = 0;
        const auto &scalar_evaluators = ref_elem_accessor.get_scalar_evaluators()(comp);

        array< ValueTable<ValueType1D>,dim>  phi_1D;
        array< ValueTable<ValueType1D>,dim> Dphi_1D;
        for (int i = 0 ; i < dim ; ++i)
        {
            phi_1D[i].resize(n_basis_elem[i],n_quad_points[i]);
            Dphi_1D[i].resize(n_basis_elem[i],n_quad_points[i]);
        }

        const Size n_basis = n_basis_elem.flat_size();
        for (Index flat_fn_id = 0 ; flat_fn_id < n_basis ; ++flat_fn_id)
        {
        	const TensorIndex<dim> tensor_fn_id = MultiArrayUtils<dim>::flat_to_tensor_index(flat_fn_id,weight_basis);
        	const auto bspline_evaluator = scalar_evaluators(tensor_fn_id);

            const auto &bspline1D_values = bspline_evaluator->get_derivative_components_view(0);
            const auto &bspline1D_derivatives = bspline_evaluator->get_derivative_components_view(1);

            for (int i = 0 ; i < dim ; ++i)
            {
                auto  phi_1D_ifn =  phi_1D[i].get_function_view(tensor_fn_id[i]);
                auto Dphi_1D_ifn = Dphi_1D[i].get_function_view(tensor_fn_id[i]);

                const auto & bsp_val = bspline1D_values[i];
                const auto & bsp_der = bspline1D_derivatives[i];

                for (int jpt = 0 ; jpt < n_quad_points[i] ; ++jpt)
                {
                     phi_1D_ifn[jpt] = bsp_val(jpt);
                    Dphi_1D_ifn[jpt] = bsp_der(jpt);
                }
            }
        }
        //--------------------------------------------------------------------------






//        auto points  = elem->get_points();
        auto phi     = elem->get_basis_values();
//        auto grd_phi = elem->get_basis_gradients();
        auto w_meas  = elem->get_w_measures();



        auto elem_measure = ref_elem_accessor.get_measure();

        f.evaluate(quad_proj_.get_points().get_flat_cartesian_product(), f_values_proj);





        //----------------------------------------------------
        // Projection phase -- begin
        start_projection = Clock::now();
        DynamicMultiArray<Real,dim> c;
        c = DynamicMultiArray<Real,dim>(TensorSize<dim>(n_quad_proj_1D_));
        for (Index i = 0 ; i < c.flat_size() ; ++i)
        {
            c(i) = f_values_proj[i](0) ;
        }

        for (Index bernst_flat_id = 0 ; bernst_flat_id < n_basis_proj_ ; ++bernst_flat_id)
        {
            bernst_tensor_id = MAUtils::flat_to_tensor_index(bernst_flat_id,bernst_tensor_weight_);

            integral_rhs(bernst_flat_id) =
                (integrate_rhs(c,w_B_proj_1D,bernst_tensor_id)).get_data()[0];
        }

        // coefficient of the rhs projection with the Bersntein basis
        k_rhs = boost::numeric::ublas::prod(inv_B_proj_, integral_rhs) ;

        end_projection = Clock::now();
        elapsed_time_projection_ += end_projection - start_projection;

        // Projection phase -- end
        //----------------------------------------------------




        //----------------------------------------------------
        // Assembly of the local mass matrix using sum-factorization -- begin
        start_assembly_mass_matrix = Clock::now();




        //----------------------------------------------------
        // precalculation of the I[i](lambda,mu1,mu2) terms
        const auto start_compute_I = Clock::now();

        array<DynamicMultiArray<Real,3>,dim> I_container;
        for (int dir = 0 ; dir < dim ; ++dir)
        {
            TensorSize<3> tensor_size_I;
            tensor_size_I[0] = n_bernst_1D_;
            tensor_size_I[1] = n_basis_elem[dir];
            tensor_size_I[2] = n_basis_elem[dir];

            I_container[dir].resize(tensor_size_I);

            auto &I = I_container[dir];

            const auto &phi_dir = phi_1D[dir];

            const auto &w_B_dir = w_B_proj_1D[dir];

            const Size n_pts = n_quad_points[dir];
//            LogStream out ;

            vector<Real> phi_mu1_mu2(n_pts);

            Index flat_id_I = 0 ;
            for (int mu2 = 0 ; mu2 < n_basis_elem[dir] ; ++mu2)
            {
                const auto phi_1D_mu2 = phi_dir.get_function_view(mu2);

                for (int mu1 = 0 ; mu1 < n_basis_elem[dir] ; ++mu1)
                {
                    const auto phi_1D_mu1 = phi_dir.get_function_view(mu1);

                    for (int jpt = 0 ; jpt < n_pts ; ++jpt)
                        phi_mu1_mu2[jpt] = phi_1D_mu2[jpt](0) * phi_1D_mu1[jpt](0);

                    for (int lambda = 0 ; lambda < n_bernst_1D_ ; ++lambda)
                    {
                        const auto w_B_lambda = w_B_dir.get_function_view(lambda);

                        Real sum = 0.0;
                        for (int jpt = 0 ; jpt < n_pts ; ++jpt)
                            sum += w_B_lambda[jpt](0) * phi_mu1_mu2[jpt];

                        I(flat_id_I++) = sum;
                    }
                }
            }
            /*
                        out << "w_B[" << dir << "]=" << endl;
                        w_B_proj_1D[dir].print_info(out);

                        out << "phi[" << dir << "]=" << endl;
                        phi_1D[dir].print_info(out);

                        //*/
        }
        const auto end_compute_I = Clock::now();
        elapsed_time_compute_I_ = end_compute_I - start_compute_I;
        //----------------------------------------------------





        TensorSize<dim> k_tensor_size;
        for (int i = 0 ; i <dim ; ++i)
            k_tensor_size[i] = n_bernst_1D_;

        DynamicMultiArray<Real,dim> K(k_tensor_size);
        AssertThrow(K.flat_size() == k_rhs.size(),
                    ExcDimensionMismatch(K.flat_size(),k_rhs.size()));

        const Size flat_size = K.flat_size();
        for (int flat_id = 0 ; flat_id < flat_size ; ++ flat_id)
            K(flat_id) = k_rhs(flat_id) / elem_measure ;

        loc_mass_matrix_sf.clear();

            const int k = 1;
            TensorSize<dim> tensor_size_alphabeta;
            TensorSize<dim> tensor_size_theta;
            for (int i = 0 ; i < dim ; ++i)
            {
                tensor_size_alphabeta[i] = n_basis_elem[i];
                tensor_size_theta[i] = n_bernst_1D_;
            }


            Size size_flat_theta_k_1 = 1;
            for (int i = k ; i < dim ; ++i)
                size_flat_theta_k_1 *= n_bernst_1D_;

            TensorSize<3> tensor_size_C_0;
            tensor_size_C_0[0] = tensor_size_theta.flat_size(); // theta size
            tensor_size_C_0[1] = 1; // alpha size
            tensor_size_C_0[2] = 1; // beta size
            DynamicMultiArray<Real,3> C_0(tensor_size_C_0); // C0


            for (Index flat_id = 0 ; flat_id < tensor_size_C_0[0] ; ++flat_id)
                C_0(flat_id) = K(flat_id);

//            const DynamicMultiArray<Real,3> &J_k = ;  // J1

            MassMatrixIntegrator<dim> integrate_mass_matrix;
            DynamicMultiArray<Real,3> C_ab = integrate_mass_matrix(tensor_size_alphabeta,tensor_size_theta,I_container,C_0);

            Assert(C_ab.tensor_size()(1) == n_basis,ExcDimensionMismatch(C_ab.tensor_size()(1),n_basis));
            Assert(C_ab.tensor_size()(2) == n_basis,ExcDimensionMismatch(C_ab.tensor_size()(2),n_basis));


            Index flat_id = 0 ;
            for (int alpha_flat_id = 0 ; alpha_flat_id < n_basis ; ++alpha_flat_id)
                for (int beta_flat_id = 0 ; beta_flat_id < n_basis ; ++beta_flat_id)
                    loc_mass_matrix_sf(alpha_flat_id,beta_flat_id) = C_ab(flat_id++);



/*
        for (int alpha_flat_id = 0 ; alpha_flat_id < n_basis ; ++alpha_flat_id)
            for (int beta_flat_id = 0 ; beta_flat_id < alpha_flat_id ; ++beta_flat_id)
                loc_mass_matrix_sf(alpha_flat_id,beta_flat_id) = loc_mass_matrix_sf(beta_flat_id,alpha_flat_id);
//*/
        //        out << "Local Mass Matrix = " << local_mass_matrix << endl;

        end_assembly_mass_matrix = Clock::now();

        elapsed_time_assembly_mass_matrix_ += end_assembly_mass_matrix - start_assembly_mass_matrix;

        // Assembly of the local mass matrix using sum-factorization -- end
        //----------------------------------------------------




        const int n_qp = this->elem_quad.get_num_points();

        start_assembly_mass_matrix_old = Clock::now();
        loc_mat.clear();

        for (int i = 0; i < n_basis; ++i)
        {
            const auto phi_i = phi.get_function_view(i);
            for (int j = i; j < n_basis; ++j)
            {
                const auto phi_j = phi.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                    loc_mat(i,j) += phi_i[qp](0) * phi_j[qp](0) * w_meas[qp];
            }
        }
        for (int i = 0; i < n_basis; ++i)
            for (int j = 0; j < i; ++j)
                loc_mat(i,j) = loc_mat(j,i);

        end_assembly_mass_matrix_old = Clock::now();
        elapsed_time_assembly_mass_matrix_old_ += end_assembly_mass_matrix_old - start_assembly_mass_matrix_old;
        //*/
        /*
        loc_dofs = elem->get_local_to_global();
        this->matrix->add_block(loc_dofs, loc_dofs, loc_mat);
        this->rhs->add_block(loc_dofs, loc_rhs);
        //*/

//        out<< "Local mass matrix sum-factorization=" << loc_mass_matrix_sf << endl << endl;
//        out<< "Local mass matrix original=" << loc_mat << endl << endl;
//        out<< "mass matrix difference=" << loc_mat - loc_mass_matrix_sf << endl << endl;

    }

    this->matrix->fill_complete();
//*/

    out << "Dim=" << dim << "         space_deg=" << space_deg_ << "         proj_deg=" << proj_deg_ << endl;
    out << "Elapsed seconds projection           = " << elapsed_time_projection_.count() << endl;
    out << "Elapsed seconds I computation        = " << elapsed_time_compute_I_.count() << endl;
    out << "Elapsed seconds assembly mass matrix = " << elapsed_time_assembly_mass_matrix_.count() << endl;
    out << "Elapsed seconds assembly mass matrix old = " << elapsed_time_assembly_mass_matrix_old_.count() << endl;
    out << endl;



    // AssertThrow(false,ExcNotImplemented());
}


template <int dim>
void
do_test()
{
    const int n_knots = 2;

    TableHandler elapsed_time_table;

    string time_proj = "Time projection";
    string time_mass_sum_fac = "Time mass-matrix sum_fac";
    string time_mass_orig = "Time mass-matrix orig";

    int degree_min = 3;
    int degree_max = 3;
    for (int degree = degree_min ; degree <= degree_max ; ++degree)
    {
        const int space_deg = degree;
        const int  proj_deg = degree;

        PoissonProblemSumFactorization<dim> poisson_sf(TensorSize<dim>(n_knots), space_deg, proj_deg);
        poisson_sf.run();
        //*/
        elapsed_time_table.add_value("Space degree",space_deg);
        elapsed_time_table.add_value("Projection degree",proj_deg);
        elapsed_time_table.add_value(time_proj,poisson_sf.get_elapsed_time_projection());
        elapsed_time_table.add_value(time_mass_sum_fac,poisson_sf.get_elapsed_time_assembly_mass_matrix());
        elapsed_time_table.add_value(time_mass_orig,poisson_sf.get_elapsed_time_assembly_mass_matrix_old());

    }

    elapsed_time_table.set_precision(time_proj,10);
    elapsed_time_table.set_scientific(time_proj,true);

    elapsed_time_table.set_precision(time_mass_sum_fac,10);
    elapsed_time_table.set_scientific(time_mass_sum_fac,true);

    elapsed_time_table.set_precision(time_mass_orig,10);
    elapsed_time_table.set_scientific(time_mass_orig,true);

    ofstream elapsed_time_file("sum_factorization_time_"+to_string(dim)+"D.txt");
    elapsed_time_table.write_text(elapsed_time_file);

}


int main()
{

    do_test<3>();
//*/
    return  0;
}
