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
// [old includes]

// [unqualified names]
using namespace iga;
using namespace std;
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
    box[0] = {{0.5,1}};
    for (int i=1; i<dim; ++i)
        box[i] = {{PI/4,PI/2}};

    auto grid = CartesianGrid<dim>::create(box, n_knots);
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




    void assemble();

private:
    using base_t = PoissonProblem< dim, PoissonProblemSumFactorization<dim> >;

    int space_deg_ = 0;
    int proj_deg_ = 0;


    DenseMatrix eval_matrix_proj() const;
};

template<int dim>
PoissonProblemSumFactorization<dim>::
PoissonProblemSumFactorization(const TensorSize<dim> &n_knots,
                               const int space_deg,const int proj_deg)
    :
    base_t(n_knots,space_deg),
    space_deg_(space_deg),
    proj_deg_(proj_deg)
{}

template<int dim>
DenseMatrix
PoissonProblemSumFactorization<dim>::
eval_matrix_proj() const
{
    const int q = proj_deg_;

    DenseMatrix matrix_proj(q+1,q+1);

    for (int row = 0 ; row <= q ; ++row)
    {
        const Real tmp = constexpr_binomial_coefficient(q,row) / (2.0 *q + 1.0) ;
        for (int col = 0 ; col <= row ; ++col)
        {
            matrix_proj(row,col) = tmp / constexpr_binomial_coefficient(2*q,row+col) *
                                   constexpr_binomial_coefficient(q,col) ;
        }
    }

    for (int row = 0 ; row <= q ; ++row)
        for (int col = row+1 ; col <= q ; ++col)
            matrix_proj(row,col) = matrix_proj(col,row) ;

    return matrix_proj;
}


template<int k>
MultiArray<Real,k-1> sum_factorization_rhs_step(
    const MultiArray<Real,k> &Gamma_k,
    const vector<Real> &omega_B_etai)
{

    const auto tensor_size = Gamma_k.get_sizes();

#ifndef NDEBUG
    for (int i = 1 ; i < k ; ++i)
    {
        Assert(tensor_size[i] == tensor_size[0],
               ExcMessage("Gamma_k computed with different number of points along the rank directions."));
    }
#endif

    const Size n_pts = omega_B_etai.size();
//  sum_factorization_rhs_quad<dim,dim-1>(Gamma_k,omega_B_etai)
//  Gamma_k[][][][][i]
}


template <int dim,int k=0>
class IntegratorSumFacRHS
{
public:
    DynamicMultiArray<Real,k> operator()(
        const DynamicMultiArray<Real,dim> &Gamma,
        const ValueTable<typename Function<1>::ValueType> &omega_B_1d,
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

        const Size n_pts = omega_B_1d.get_num_points();
        const auto &B = omega_B_1d.get_function_view(eta_k1);

        for (Index flat_id_C = 0 ; flat_id_C < n_entries_C ; ++flat_id_C)
        {
            TensorIndex<k> tensor_id_C = C.flat_to_tensor(flat_id_C);

            for (int i = 0 ; i < k ; ++i)
                tensor_id_C1[i] = tensor_id_C[i];

            Real sum = 0.0;
            for (int ipt = 0 ; ipt < n_pts ; ++ipt) // quadrature along the last tensor index
            {
                tensor_id_C1[k] = ipt ;
                sum += C1(tensor_id_C1) * B[ipt](0);
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
        const ValueTable<typename Function<1>::ValueType> &omega_B_1d,
        const TensorIndex<dim> &bernstein_tensor_id)
    {
        return Gamma;
    }
};


template <int dim,int k>
class IntegratorSumFacMass
{
public:
    using ValueTable1D = ValueTable<typename Function<1>::ValueType> ;

    Real
    operator()(const DynamicMultiArray<Real,k> &lambda_k,
               const ValueTable1D &omega_B,
               const array<ValueTable1D,dim> &phi,
               const TensorIndex<dim> &alpha_tensor_id,
               const TensorIndex<dim> &beta_tensor_id)
    {
        TensorSize<k> tensor_size_lambda_k = lambda_k.tensor_size();

        TensorSize<k-1> tensor_size_lambda_k_1;
        for (int i = 0 ; i < k-1 ; ++i)
            tensor_size_lambda_k_1(i) = tensor_size_lambda_k(i);

        DynamicMultiArray<Real,k-1> lambda_k_1(tensor_size_lambda_k_1);


        const Size n_bernst = omega_B.get_num_functions();
        const Size n_quad = omega_B.get_num_points();

        Assert(n_bernst == tensor_size_lambda_k(0),
               ExcDimensionMismatch(n_bernst,tensor_size_lambda_k(0)));


        TensorIndex<k> tensor_id_lambda_k;
        TensorIndex<k-1> tensor_id_lambda_k_1;
        const int dir = dim-k;

        Assert(n_quad == phi[dir].get_num_points(),
               ExcDimensionMismatch(n_quad,phi[dir].get_num_points()));
        const auto phi_alpha = phi[dir].get_function_view(alpha_tensor_id[dir]);
        const auto phi_beta  = phi[dir].get_function_view(beta_tensor_id[dir]);

        Index flat_id_entry_lambda_k_1 = 0;
        for (auto & value_entry_lambda_k_1 : lambda_k_1)
        {
            const auto tensor_id_lambda_k_1 = lambda_k_1.flat_to_tensor(flat_id_entry_lambda_k_1);

            for (int i = dim-1 ; i >= 1 ; --i)
                tensor_id_lambda_k[i] = tensor_id_lambda_k_1[i-1];

            Real sum = 0.0;
            for (int ifn = 0 ; ifn < n_bernst ; ++ifn)
            {
                const auto omega_B_theta = omega_B.get_function_view(ifn);

                Real sum_quad_pt = 0.0;
                for (int jpt = 0 ; jpt < n_quad ; ++jpt)
                    sum_quad_pt += omega_B_theta[jpt](0) * phi_alpha[jpt](0) * phi_beta[jpt](0);

                tensor_id_lambda_k[0] = ifn;
                sum += lambda_k(tensor_id_lambda_k) * sum_quad_pt;
            }

            value_entry_lambda_k_1 = sum;

            flat_id_entry_lambda_k_1++;
        }


        IntegratorSumFacMass<dim,k-1> integrate;

        return integrate(lambda_k_1,omega_B,phi,alpha_tensor_id,beta_tensor_id);
    };
};

template <int dim>
class IntegratorSumFacMass<dim,1>
{
public:
    using ValueTable1D = ValueTable<typename Function<1>::ValueType> ;

    Real
    operator()(const DynamicMultiArray<Real,1> &lambda_1,
               const ValueTable1D &omega_B,
               const array<ValueTable1D,dim> phi,
               const TensorIndex<dim> &alpha_tensor_id,
               const TensorIndex<dim> &beta_tensor_id)
    {
        const Size size_lambda_1 = lambda_1.flat_size();

        const Size n_bernst = omega_B.get_num_functions();
        const Size n_quad = omega_B.get_num_points();

        Assert(n_bernst == size_lambda_1,
               ExcDimensionMismatch(n_bernst,size_lambda_1));


        const int dir = dim-1;

        Assert(n_quad == phi[dir].get_num_points(),
               ExcDimensionMismatch(n_quad,phi[dir].get_num_points()));
        const auto phi_alpha = phi[dir].get_function_view(alpha_tensor_id[dir]);
        const auto phi_beta  = phi[dir].get_function_view(beta_tensor_id[dir]);

        Real lambda_0 = 0.0;
        for (int ifn = 0 ; ifn < n_bernst ; ++ifn)
        {
            const auto omega_B_theta = omega_B.get_function_view(ifn);

            Real sum_quad_pt = 0.0;
            for (int jpt = 0 ; jpt < n_quad ; ++jpt)
                sum_quad_pt += omega_B_theta[jpt](0) * phi_alpha[jpt](0) * phi_beta[jpt](0);

            lambda_0 += lambda_1(ifn) * sum_quad_pt;
        }

        return lambda_0;
    };
};



template<int dim>
void
PoissonProblemSumFactorization<dim>::
assemble()
{
//    using std::cout;
    LogStream out ;

    DenseMatrix matrix_proj = eval_matrix_proj();

    DenseMatrix inv_matrix_proj = matrix_proj.inverse();

//    cout << "    proj matrix: " <<     matrix_proj <<endl;
//    cout << "inv proj matrix: " << inv_matrix_proj <<endl;


    const int n_basis = this->space->get_num_basis_per_element();
    DenseMatrix loc_mat(n_basis, n_basis);
    DenseVector loc_rhs(n_basis);
    vector<Index> loc_dofs(n_basis);

//    const int n_qp = this->elem_quad.get_num_points();
    ConstantFunction<dim> f({0.5});

    using Val = typename Function<dim>::ValueType;

    auto elem = this->space->begin();
    const auto elem_end = this->space->end();
    ValueFlags fill_flags = ValueFlags::value | ValueFlags::gradient | ValueFlags::w_measure | ValueFlags::point;
    elem->init_values(fill_flags, this->elem_quad);

    // quadrature scheme for the projection on Bernstein basis
    const Size n_quad_proj = proj_deg_+1;
    QGauss<dim> quad_proj(n_quad_proj);
    QGauss<1> quad_proj_1d(n_quad_proj);
    vector<Val> f_values_proj(quad_proj.get_num_points());

    const vector<Real> x_1d = quad_proj_1d.get_points().get_data_direction(0);
    const vector<Real> w_1d = quad_proj_1d.get_weights().get_data_direction(0);

    // values of the one-dimensional Bernstein basis
    const DenseMatrix B_1d = BernsteinBasis::evaluate(proj_deg_,x_1d);
    const Size n_bernstein_1d = B_1d.size1();
    Assert(n_bernstein_1d == proj_deg_+1, ExcDimensionMismatch(n_bernstein_1d, proj_deg_+1));



    //------------------------------------
    // evaluation of the Bernstein's mass matrix
    DenseMatrix Bernstein_1d_mass_matrix(n_bernstein_1d,n_bernstein_1d);
    const int q = proj_deg_;
    for (int i = 0 ; i < n_bernstein_1d ; ++i)
    {
        const Real tmp_i = constexpr_binomial_coefficient(q,i) / (2.0 * q + 1.0);
        for (int j = 0 ; j < n_bernstein_1d ; ++j)
        {
            Bernstein_1d_mass_matrix(i,j) =
                tmp_i * constexpr_binomial_coefficient(q,j) /
                constexpr_binomial_coefficient(2*q,i+j);
        }
    }
    out << "Bernstein 1D mass matrix = " << Bernstein_1d_mass_matrix << endl;


    Size n_rows_B_mass_matrix = pow(n_bernstein_1d,dim);
    DenseMatrix Bernstein_mass_matrix(n_rows_B_mass_matrix,n_rows_B_mass_matrix);

    TensorSize<dim> bernstein_tensor_size(n_bernstein_1d);
    TensorIndex<dim> bernstein_tensor_weight = MultiArrayUtils<dim>::compute_weight(bernstein_tensor_size);
    for (Size row = 0 ; row < n_rows_B_mass_matrix ; ++row)
    {
        TensorIndex<dim> row_tensor_id = MultiArrayUtils<dim>::flat_to_tensor_index(row,bernstein_tensor_weight);

        for (Size col = 0 ; col < n_rows_B_mass_matrix ; ++col)
        {
            TensorIndex<dim> col_tensor_id = MultiArrayUtils<dim>::flat_to_tensor_index(col,bernstein_tensor_weight);

            Real tmp = 1.0;
            for (int k = 0 ; k < dim ; ++k)
                tmp *= Bernstein_1d_mass_matrix(row_tensor_id[k],col_tensor_id[k]);

            Bernstein_mass_matrix(row,col) = tmp;
        }
    }
    out << "Bernstein "<< dim <<"D mass matrix = " << Bernstein_mass_matrix << endl;


    DenseMatrix inverse_Bernstein_mass_matrix = Bernstein_mass_matrix.inverse();
    out << "inverse of Bernstein "<< dim <<"D mass matrix = " << inverse_Bernstein_mass_matrix << endl;
    //------------------------------------




    using ValueType1D = Function<1>::ValueType;
    ValueTable<ValueType1D> omega_B_1d(n_bernstein_1d,n_quad_proj);
    for (Index ifunc = 0 ; ifunc < n_bernstein_1d ; ++ifunc)
    {
        auto omega_B_ifunc =  omega_B_1d.get_function_view(ifunc);
        for (Index ipt = 0 ; ipt < n_quad_proj ; ++ipt)
            omega_B_ifunc[ipt] = w_1d[ipt] * B_1d(ifunc,ipt);
    }

    IntegratorSumFacRHS<dim> integrate_rhs;

    TensorIndex<dim> bernst_tensor_id;

    const Size bernstein_flat_size = bernstein_tensor_size.flat_size();
    DenseVector integral_rhs(bernstein_flat_size);
    DenseVector k_rhs(bernstein_flat_size);

    // number of points along each directin for the quadrature scheme.
    TensorSize<dim> n_quad_points = this->elem_quad.get_num_points_direction();

    TensorSize<dim> n_basis_elem(space_deg_+1);

    for (; elem != elem_end; ++elem)
    {
        loc_mat.clear();
        loc_rhs.clear();
        elem->fill_values();


        //--------------------------------------------------------------------------
        // getting the 1D basis functions values
        const auto &splines1d_direction = elem->elem_univariate_values_(0);

        array< ValueTable<ValueType1D>,dim>  phi_1D;
        array< ValueTable<ValueType1D>,dim> Dphi_1D;
        for (int i = 0 ; i < dim ; ++i)
        {
            const auto &splines1d = *(splines1d_direction[i]);

            phi_1D[i].resize(n_basis_elem[i],n_quad_points[i]);
            Dphi_1D[i].resize(n_basis_elem[i],n_quad_points[i]);

            for (int ifn = 0 ; ifn < n_basis_elem[i] ; ++ifn)
            {
                auto  phi_1D_ifn =  phi_1D[i].get_function_view(ifn);
                auto Dphi_1D_ifn = Dphi_1D[i].get_function_view(ifn);

                for (int jpt = 0 ; jpt < n_quad_points[i] ; ++jpt)
                {
                    phi_1D_ifn[jpt] = splines1d[0](ifn,jpt);
                    Dphi_1D_ifn[jpt] = splines1d[1](ifn,jpt);
                }
            }
        }
        //--------------------------------------------------------------------------






        auto points  = elem->get_points();
//        auto phi     = elem->get_basis_values();
//        auto grd_phi = elem->get_basis_gradients();
        auto w_meas  = elem->get_w_measures();

        f.evaluate(quad_proj.get_points().get_flat_cartesian_product(), f_values_proj);

        DynamicMultiArray<Real,dim> c(n_quad_proj);
        for (Index i = 0 ; i < c.flat_size() ; ++i)
        {
            c(i) = f_values_proj[i](0) ;
        }

        for (Index bernst_flat_id = 0 ; bernst_flat_id < bernstein_flat_size ; ++bernst_flat_id)
        {
            bernst_tensor_id = MultiArrayUtils<dim>::flat_to_tensor_index(bernst_flat_id,bernstein_tensor_weight);

            integral_rhs(bernst_flat_id) = (integrate_rhs(c,omega_B_1d,bernst_tensor_id)).get_data()[0];
        }
        out << "integral rhs = " << integral_rhs << endl;

        // coefficient of the rhs projection with the Bersntein basis
        k_rhs = boost::numeric::ublas::prod(inverse_Bernstein_mass_matrix, integral_rhs) ;
        out << "k rhs = " << k_rhs << endl;

        if (dim==3)
        {
            TensorSize<dim> k_tensor_size;
            for (int i = 0 ; i <dim ; ++i)
                k_tensor_size[i] = n_bernstein_1d;

            DynamicMultiArray<Real,dim> K(k_tensor_size);
            AssertThrow(K.flat_size() == k_rhs.size(),
                        ExcDimensionMismatch(K.flat_size(),k_rhs.size()));

            const Size flat_size = K.flat_size();
            for (int flat_id = 0 ; flat_id < flat_size ; ++ flat_id)
                K(flat_id) = k_rhs(flat_id);


            Size n_flat_alpha = elem->get_num_basis();
            Size n_flat_beta = elem->get_num_basis();

            const auto weight_alpha = MultiArrayUtils<dim>::compute_weight(n_basis_elem);
            const auto weight_beta  = MultiArrayUtils<dim>::compute_weight(n_basis_elem);

            DenseMatrix local_mass_matrix(n_flat_alpha,n_flat_beta);

            IntegratorSumFacMass<dim,dim> integrate;

            for (int alpha_flat_id = 0 ; alpha_flat_id < n_flat_alpha ; ++alpha_flat_id)
            {
                TensorIndex<dim> alpha_tensor_id =
                    MultiArrayUtils<dim>::flat_to_tensor_index(alpha_flat_id,weight_alpha);

                for (int beta_flat_id = 0 ; beta_flat_id < n_flat_beta ; ++beta_flat_id)
                {
                    TensorIndex<dim> beta_tensor_id =
                        MultiArrayUtils<dim>::flat_to_tensor_index(beta_flat_id,weight_beta);

                    /*
                    int dir;


                    DynamicMultiArray<Real,3> Lambda_3 = K;
                    TensorSize<3> tensor_size_Lambda_3 = Lambda_3.tensor_size();

                    TensorSize<2> tensor_size_Lambda_2;
                    tensor_size_Lambda_2(0) = tensor_size_Lambda_3(1);
                    tensor_size_Lambda_2(1) = tensor_size_Lambda_3(2);
                    DynamicMultiArray<Real,2> Lambda_2(tensor_size_Lambda_2);

                    TensorIndex<3> tensor_id_lambda_3;
                    TensorIndex<2> tensor_id_lambda_2;
                    dir = 0;
                    for (int i2 = 0 ; i2 < n_bernstein_1d ; ++i2)
                    {
                        tensor_id_lambda_3[2] = i2;
                        tensor_id_lambda_2[1] = i2;
                        for (int i1 = 0 ; i1 < n_bernstein_1d ; ++i1)
                        {
                            tensor_id_lambda_3[1] = i1;
                            tensor_id_lambda_2[0] = i1;

                            Real sum = 0.0;
                            for (int i0 = 0 ; i0 < n_bernstein_1d ; ++i0)
                            {
                                tensor_id_lambda_3[0] = i0;

                                const auto omega_B_theta =  omega_B_1d.get_function_view(i0);

                                const auto phi_alpha =  phi_1D[dir].get_function_view(alpha_tensor_id[dir]);
                                const auto phi_beta =  phi_1D[dir].get_function_view(beta_tensor_id[dir]);

                                Real sum_quad_pt = 0.0;
                                for (int gamma = 0 ; gamma < q ; ++gamma)
                                {
                                    sum_quad_pt += omega_B_theta[gamma](0) * phi_alpha[gamma](0) * phi_beta[gamma](0);
                                }

                                sum += Lambda_3(tensor_id_lambda_3) * sum_quad_pt;
                            }
                            Lambda_2(tensor_id_lambda_2) = sum;
                        }
                    }


                    TensorSize<1> tensor_size_Lambda_1;
                    tensor_size_Lambda_1(0) = tensor_size_Lambda_2(1);
                    DynamicMultiArray<Real,1> Lambda_1(tensor_size_Lambda_1);

                    TensorIndex<1> tensor_id_lambda_1;
                    dir = 1;
                    for (int i1 = 0 ; i1 < n_bernstein_1d ; ++i1)
                    {
                        tensor_id_lambda_2[1] = i1;
                        tensor_id_lambda_1[0] = i1;

                        Real sum = 0.0;
                        for (int i0 = 0 ; i0 < n_bernstein_1d ; ++i0)
                        {
                            tensor_id_lambda_2[0] = i0;

                            const auto omega_B_theta =  omega_B_1d.get_function_view(i0);

                            const auto phi_alpha =  phi_1D[dir].get_function_view(alpha_tensor_id[dir]);
                            const auto phi_beta =  phi_1D[dir].get_function_view(beta_tensor_id[dir]);

                            Real sum_quad_pt = 0.0;
                            for (int gamma = 0 ; gamma < q ; ++gamma)
                            {
                                sum_quad_pt += omega_B_theta[gamma](0) * phi_alpha[gamma](0) * phi_beta[gamma](0);
                            }

                            sum += Lambda_2(tensor_id_lambda_2) * sum_quad_pt;
                        }
                        Lambda_1(tensor_id_lambda_1) = sum;
                    }


                    Real Lambda_0;
                    Real sum = 0.0;
                    dir = 2;
                    for (int i0 = 0 ; i0 < n_bernstein_1d ; ++i0)
                    {
                        tensor_id_lambda_1[0] = i0;

                        const auto omega_B_theta =  omega_B_1d.get_function_view(i0);

                        const auto phi_alpha =  phi_1D[dir].get_function_view(alpha_tensor_id[dir]);
                        const auto phi_beta =  phi_1D[dir].get_function_view(beta_tensor_id[dir]);

                        Real sum_quad_pt = 0.0;
                        for (int gamma = 0 ; gamma < q ; ++gamma)
                        {
                            sum_quad_pt += omega_B_theta[gamma](0) * phi_alpha[gamma](0) * phi_beta[gamma](0);
                        }

                        sum += Lambda_1(tensor_id_lambda_1) * sum_quad_pt;
                    }
                    Lambda_0 = sum;

                    local_mass_matrix(alpha_flat_id,beta_flat_id) = Lambda_0;
                    //*/

                    local_mass_matrix(alpha_flat_id,beta_flat_id) = integrate(K,omega_B_1d,phi_1D,alpha_tensor_id,beta_tensor_id);

                }
            }
            out << "Local Mass Matrix = " << local_mass_matrix << endl;
        }
        else
        {
            AssertThrow(false,ExcNotImplemented());
        }

        /*
                for (int i = 0; i < n_basis; ++i)
                {
                    auto grd_phi_i = grd_phi.get_function(i);
                    for (int j = 0; j < n_basis; ++j)
                    {
                        auto grd_phi_j = grd_phi.get_function(j);
                        for (int qp = 0; qp < n_qp; ++qp)
                            loc_mat(i,j) +=
                                scalar_product(grd_phi_i[qp], grd_phi_j[qp])
                                * w_meas[qp];
                    }

                    auto phi_i = phi.get_function(i);
                    for (int qp = 0; qp < n_qp; ++qp)
                        loc_rhs(i) += scalar_product(phi_i[qp], f_values[qp])
                                      * w_meas[qp];
                }
        //*/
        loc_dofs = elem->get_local_to_global();
        this->matrix->add_block(loc_dofs, loc_dofs, loc_mat);
        this->rhs->add_block(loc_dofs, loc_rhs);
    }

    this->matrix->fill_complete();






    // AssertThrow(false,ExcNotImplemented());
}


int main()
{
    const int n_knots = 2;
    const int space_deg = 2;
    const int  proj_deg = 2;
    /*
        PoissonProblemStandardIntegration<1> poisson_1d({n_knots}, space_deg);
        poisson_1d.run();
    //*/

//    PoissonProblemStandardIntegration<2> poisson_2d({n_knots, n_knots}, space_deg);
//    poisson_2d.run();

//    PoissonProblemStandardIntegration<3> poisson_3d({n_knots, n_knots, n_knots}, space_deg);
//    poisson_3d.run();

    /*
        PoissonProblemSumFactorization<1> poisson_1d_sf({n_knots}, space_deg, proj_deg);
        poisson_1d_sf.run();

        PoissonProblemSumFactorization<2> poisson_2d_sf({n_knots,n_knots}, space_deg, proj_deg);
        poisson_2d_sf.run();
    //*/
    PoissonProblemSumFactorization<3> poisson_3d_sf({n_knots,n_knots,n_knots}, space_deg, proj_deg);
    poisson_3d_sf.run();
//*/
    return  0;
}
