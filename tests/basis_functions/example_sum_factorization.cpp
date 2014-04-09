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

#include <numeric>

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



template <class PhysSpaceTest,class PhysSpaceTrial>
class EllipticOperators
{
public:
	/** Type for the element accessor of the <em>test</em> physical space. */
	using ElemTest = typename PhysSpaceTest::ElementAccessor;

	/** Type for the element accessor of the <em>trial</em> physical space. */
	using ElemTrial = typename PhysSpaceTrial::ElementAccessor;

	static const int dim = PhysSpaceTest::dim;

	/**
	 * Default constructor.
	 * In Debug mode, it checks if the template arguments are consistent.
	 */
	EllipticOperators();

	virtual void eval_operator_u_v(
			const ElemTest &elem_test,
			const ElemTrial &elem_trial,
			const ValueVector<Real> &c,
			DenseMatrix &operator_u_v) const = 0;


	/** Returns true if the space for the test functions and the trial functions is the same. */
	bool test_same_space(const ElemTest &elem_test,const ElemTrial &elem_trial) const ;
};

template <class PhysSpaceTest,class PhysSpaceTrial>
inline
EllipticOperators<PhysSpaceTest,PhysSpaceTrial>::
EllipticOperators()
{
    //-----------------------------------------------------------------
    Assert(PhysSpaceTest::dim == PhysSpaceTrial::dim,
    		ExcDimensionMismatch(PhysSpaceTest::dim,PhysSpaceTrial::dim));
    Assert(PhysSpaceTest::space_dim == PhysSpaceTrial::space_dim,
    		ExcDimensionMismatch(PhysSpaceTest::space_dim,PhysSpaceTrial::space_dim));
    Assert(PhysSpaceTest::range == PhysSpaceTrial::range,
    		ExcDimensionMismatch(PhysSpaceTest::range,PhysSpaceTrial::range));
    Assert(PhysSpaceTest::rank == PhysSpaceTrial::rank,
    		ExcDimensionMismatch(PhysSpaceTest::rank,PhysSpaceTrial::rank));

//    const int dim = PhysSpaceTest::dim;
//    const int space_dim = PhysSpace::space_dim;
//    const int range = PhysSpace::range;
//    const int rank = PhysSpace::rank;

    Assert(PhysSpaceTest::range == 1,ExcDimensionMismatch(PhysSpaceTest::range,1));
    Assert(PhysSpaceTest::rank == 1,ExcDimensionMismatch(PhysSpaceTest::rank,1));
    //-----------------------------------------------------------------
}


template <class PhysSpaceTest,class PhysSpaceTrial>
bool
EllipticOperators<PhysSpaceTest,PhysSpaceTrial>::
test_same_space(const ElemTest &elem_test,const ElemTrial &elem_trial) const
{
    //--------------------------------------------------------------------------
    // checks that the mapping used in the test space and in the trial space is the same
    Assert(elem_test.get_physical_space()->get_push_forward()->get_mapping() ==
    		elem_trial.get_physical_space()->get_push_forward()->get_mapping(),
    		ExcMessage("Test and trial spaces must have the same mapping (and the same grid)!"));
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // checks that the elements on the grid are the same
    using GridElem = CartesianGridElementAccessor<dim>;
    Assert(static_cast<const GridElem &>(elem_test.get_ref_space_accessor()) ==
    	   static_cast<const GridElem &>(elem_trial.get_ref_space_accessor()),
    	   ExcMessage("Different elements for test space and trial space."));
    //--------------------------------------------------------------------------

    // the test is true only if the element accessors reference the same memory location
    //TODO(MM 08 apr 2014): this is a really crude/raw test. Maybe abetter test would be
    // a comparison of iterator index, grid, knots values and multiplicities of the underlying
    // reference space and also push-forwardand mapping.
    return (&elem_test==&elem_trial)?true:false;
}


template <class PhysSpaceTest,class PhysSpaceTrial>
class EllipticOperatorsSumFactorizationIntegration :
		public EllipticOperators<PhysSpaceTest,PhysSpaceTrial>
{
public:
	/** Type for the element accessor of the <em>test</em> physical space. */
	using ElemTest = typename PhysSpaceTest::ElementAccessor;

	/** Type for the element accessor of the <em>trial</em> physical space. */
	using ElemTrial = typename PhysSpaceTrial::ElementAccessor;

	static const int dim = PhysSpaceTest::dim;

    using Clock = chrono::high_resolution_clock;
    using TimePoint = chrono::time_point<Clock>;
    using Duration = chrono::duration<Real>;

    /**
     * Constructor. Builds an object that implements the method for computing several
     * elliptic operators (in the sense of local element-contribution) using the
     * sum factorization technique.
     * The input arguments of this constructor specifies the polynomial degree and the
     * quadrature scheme used to project the non-tensor product functions in the different
     * operator, using the Bernstein's basis.
     */
    EllipticOperatorsSumFactorizationIntegration(
    		const TensorSize<dim> &projection_degree,
    		const Quadrature<dim> &projection_quadrature);



	virtual void eval_operator_u_v(
			const ElemTest &elem_test,
			const ElemTrial &elem_trial,
			const ValueVector<Real> &c,
			DenseMatrix &operator_u_v) const override;

private:

	/**
	 * @brief Performs the L2 projection on a function evaluated on points inside a single
	 * Bezier element, using as basis of the projection space the Benrstein's polynomials.
	 * @param func_to_proj_at_pts Values of the function to project at the quadrature points.
	 */
	DynamicMultiArray<Real,dim>
	projection_l2_bernstein_basis(const ValueVector<Real> &func_to_proj_at_pts) const;

	/**
	 * @brief Returns the moments for the mass matrix, as needed by the integration using the
	 * sum factorization technique.
	 * \f$ m_{i,\theta_i,\alpha_i,\beta_i} =
	   \int_{I_i} B_{\theta_i}(x_i) u_{\beta_i}(x_i) v_{\alpha_i}(x_i) \; dx_i \, .\f$
	 *
	 */
	array<DynamicMultiArray<Real,3>,dim>
	evaluate_moments_op_u_v(
		const array<ValueTable<Function<1>::ValueType>,dim> &phi_1D_test,
	    const array<ValueTable<Function<1>::ValueType>,dim> &phi_1D_trial,
	    const array<Real,dim> &length_element_edge) const;


    /**
     * Degree (along each coordinate direction) of the Bernstein's polynomial used
     * for the construction of the projection matrix B_proj_.
     */
    TensorSize<dim> proj_degree_;

    /** Quadrature scheme used in the projection phase to built the projection matrix B_proj_. */
    Quadrature<dim> proj_quad_;


    /**
     * Number of one-dimensional Bernstein's polynomials (along each coordinate direction)
     * used to build the projection matrix B_proj_.
     */
    TensorSize<dim> n_bernst_1D_;

    /**
     * Number of multi-dimensional Bernstein's polynomial used to define the basis of the projection space.
     */
    Size n_basis_proj_;


	/**
	 * One-dimensional quadrature weights multiplied by the one-dimensional
	 * Bernstein's basis values used to build the rhs term in the projection.
	 * @note All the quantities refers to the unitary hypercube [0,1]^{dim}.
	 */
    array<ValueTable<Real>,dim> w_times_B_proj_1D_;

    /** Mass matrix of the Bernstein's polynomials used in the projection, in [0,1]^dim.*/
    DenseMatrix B_proj_;

    /** Inverse of the mass matrix of the Bernstein polynomials used in the projection, in [0,1]^dim.*/
    DenseMatrix inv_B_proj_;


	Duration elapsed_time_operator_u_v_;
};


template <class PhysSpaceTest,class PhysSpaceTrial>
class EllipticOperatorsStandardIntegration :
		public EllipticOperators<PhysSpaceTest,PhysSpaceTrial>
{
public:
	/** Type for the element accessor of the <em>test</em> physical space. */
	using ElemTest = typename PhysSpaceTest::ElementAccessor;

	/** Type for the element accessor of the <em>trial</em> physical space. */
	using ElemTrial = typename PhysSpaceTrial::ElementAccessor;


	EllipticOperatorsStandardIntegration();

	virtual void eval_operator_u_v(
			const ElemTest &elem_test,
			const ElemTrial &elem_trial,
			const ValueVector<Real> &c,
			DenseMatrix &operator_u_v) const override;
};

template <class PhysSpaceTest,class PhysSpaceTrial>
inline
EllipticOperatorsStandardIntegration<PhysSpaceTest,PhysSpaceTrial>::
EllipticOperatorsStandardIntegration()
:
EllipticOperators<PhysSpaceTest,PhysSpaceTrial>()
{}


template <class PhysSpaceTest,class PhysSpaceTrial>
inline
void
EllipticOperatorsStandardIntegration<PhysSpaceTest,PhysSpaceTrial>::
eval_operator_u_v(
		const ElemTest &elem_test,
		const ElemTrial &elem_trial,
		const ValueVector<Real> &coeffs,
		DenseMatrix &operator_u_v) const
{
	//TODO: only the symmetric case is tested. In the non symmetric case, we need to check that
	// the physical space iterators have the same grid, map, reference space, index, etc.
	Assert(&elem_test == &elem_trial,ExcNotImplemented());

	const bool is_symmetric = this->test_same_space(elem_test,elem_trial);

	const Size n_basis_test  = elem_test .get_num_basis();
	const Size n_basis_trial = elem_trial.get_num_basis();

    const auto &phi_test  = elem_test.get_basis_values();
    const auto &phi_trial = elem_trial.get_basis_values();
    const auto &w_meas  = elem_test.get_w_measures();


    Assert(operator_u_v.get_num_rows() == n_basis_test,
    		ExcDimensionMismatch(operator_u_v.get_num_rows(),n_basis_test));
    Assert(operator_u_v.get_num_cols() == n_basis_trial,
    		ExcDimensionMismatch(operator_u_v.get_num_cols(),n_basis_trial));


    const Size n_qp = coeffs.size();
    Assert(n_qp == phi_test.get_num_points(),ExcDimensionMismatch(n_qp,phi_test.get_num_points()));
    Assert(n_qp == phi_trial.get_num_points(),ExcDimensionMismatch(n_qp,phi_trial.get_num_points()));

    vector<Real> coeffs_times_w_meas(n_qp);
	for (int qp = 0; qp < n_qp; ++qp)
		coeffs_times_w_meas[qp] = coeffs[qp] * w_meas[qp];


    operator_u_v.clear();

    if (!is_symmetric)
    {
    	for (int i = 0; i < n_basis_test; ++i)
    	{
    		const auto phi_i = phi_test.get_function_view(i);
    		for (int j = 0; j < n_basis_trial; ++j)
    		{
    			const auto phi_j = phi_trial.get_function_view(j);
    			for (int qp = 0; qp < n_qp; ++qp)
    				operator_u_v(i,j) += phi_j[qp](0) * ( phi_i[qp](0) * coeffs_times_w_meas[qp] );
    		}
    	}

    } // end if (!is_symmetric)
    else
    {
    	for (int i = 0; i < n_basis_test; ++i)
    	{
    		const auto phi_i = phi_test.get_function_view(i);
    		for (int j = i; j < n_basis_trial; ++j)
    		{
    			const auto phi_j = phi_trial.get_function_view(j);
    			for (int qp = 0; qp < n_qp; ++qp)
    				operator_u_v(i,j) += phi_j[qp](0) * ( phi_i[qp](0) * coeffs_times_w_meas[qp] );
    		}
    	}
    	for (int i = 0; i < n_basis_test; ++i)
    		for (int j = 0; j < i; ++j)
    			operator_u_v(i,j) = operator_u_v(j,i);
    } // end if (is_symmetric)
}


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
    using RefSpace  = BSplineSpace<dim>;
    using PushFw    = PushForward<Transformation::h_grad, dim>;
    using Space     = PhysicalSpace<RefSpace, PushFw>;

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

    Quadrature<dim> quad_proj_;

    chrono::duration<Real> elapsed_time_projection_;
    chrono::duration<Real> elapsed_time_assembly_mass_matrix_;

    chrono::duration<Real> elapsed_time_assembly_mass_matrix_old_;
    chrono::duration<Real> elapsed_time_compute_I_;

    DenseMatrix eval_matrix_proj() const;

};

template<int dim>
PoissonProblemSumFactorization<dim>::
PoissonProblemSumFactorization(const TensorSize<dim> &n_knots,
                               const int space_deg,const int proj_deg)
    :
    base_t(n_knots,space_deg),
    space_deg_(space_deg),
    proj_deg_(proj_deg),
    quad_proj_(QGauss<dim>(TensorSize<dim>(proj_deg_+1)))
{}


template <int dim,int k=dim>
class IntegratorSumFacRHS_ver2
{
public:
    Real operator()(
        const DynamicMultiArray<Real,k> &Cpre,
        const array<ValueTable<Real>,dim> &omega_B_1d,
        const TensorIndex<dim> &t_id_bernst)
    {
    	const int dir = dim-k;
    	const auto & w_times_B = omega_B_1d[dir];

    	const TensorSize<k> t_size_Cpre = Cpre.tensor_size();
    	for ( int i = 0 ; i < k ; ++i)
    	{
    		Assert(t_size_Cpre(i)==omega_B_1d[dim-k+i].get_num_points(),
    				ExcDimensionMismatch(t_size_Cpre(i),omega_B_1d[dim-k+i].get_num_points()));
    	}

    	TensorSize<k-1> t_size_Cpost;
    	for ( int i = 0 ; i < k-1 ; ++i)
    	{
    		t_size_Cpost(i) = t_size_Cpre(i+1);
    	}
        DynamicMultiArray<Real,k-1> Cpost(t_size_Cpost);

        const Size n_quad_pts = t_size_Cpre(0);
        const Size f_size_Cpost = t_size_Cpost.flat_size();

        Index fn_id = t_id_bernst(dir);
        const auto w_times_B_fn = w_times_B.get_function_view(fn_id);
        Index f_id_Cpre = 0;
        for (Index f_id_Cpost = 0 ; f_id_Cpost < f_size_Cpost ; f_id_Cpost++)
        {
        	Real sum = 0.0;
        	for (Index pt = 0 ; pt < n_quad_pts ; ++pt)
        	{
        		sum += w_times_B_fn[pt] * Cpre(f_id_Cpre);
        		++f_id_Cpre;
        	}
        	Cpost(f_id_Cpost) = sum;
        }

        IntegratorSumFacRHS_ver2<dim,k-1> integrate_rhs;
        return integrate_rhs(Cpost,omega_B_1d,t_id_bernst);
    }
};


template <int dim>
class IntegratorSumFacRHS_ver2<dim,1>
{
public:
    Real operator()(
        const DynamicMultiArray<Real,1> &Cpre,
        const array<ValueTable<Real>,dim> &omega_B_1d,
        const TensorIndex<dim> &t_id_bernst)
    {
    	const auto & w_times_B = omega_B_1d[dim-1];

    	const Size f_size_Cpre = Cpre.flat_size();
    	Assert(f_size_Cpre == omega_B_1d[dim-1].get_num_points(),
    			ExcDimensionMismatch(f_size_Cpre,omega_B_1d[dim-1].get_num_points()));

        const Size n_quad_pts = f_size_Cpre;

        const Index fn_id = t_id_bernst(dim-1);
        const auto w_times_B_fn = w_times_B.get_function_view(fn_id);

        Real Cpost = 0.0;
        for (Index pt = 0 ; pt < n_quad_pts ; ++pt)
        	Cpost += w_times_B_fn[pt] * Cpre(pt);

    	return Cpost;
    }
};





template <int dim,int k=dim>
class IntegratorSumFacRHS_ver3
{
public:
    void operator()(
        const DynamicMultiArray<Real,dim> &Cpre,
        const array<ValueTable<Real>,dim> &omega_B_1d,
        DenseVector &integral_rhs)
    {
#ifndef NDEBUG
    	int n_basis = 1;
    	for (int i = 0 ; i < dim ; ++i)
    		n_basis *= omega_B_1d[i].get_num_functions();
    	Assert(integral_rhs.size() == n_basis,
    			ExcDimensionMismatch(integral_rhs.size(),n_basis));
#endif

    	using std::cout;
    	using std::endl;

    	const int dir = dim-k;

    	const TensorSize<dim> t_size_Cpre = Cpre.tensor_size();
    	TensorSize<dim> t_size_Cpost;

    	TensorSize<dim-k+1> t_size_alpha_post;
    	TensorSize<k-1> t_size_theta_post;
    	for (Index i = 0 ; i < dim-k ; ++i)
    	{
    		t_size_alpha_post(i) = t_size_Cpre(i);
    		t_size_Cpost(i) = t_size_alpha_post(i);
    	}
    	t_size_alpha_post(dim-k) = omega_B_1d[dim-k].get_num_functions();
		t_size_Cpost(dim-k) = t_size_alpha_post(dim-k);

    	for (Index i = dim-k+1 ; i < dim ; ++i)
    	{
    		t_size_theta_post(i-dim+k-1) = t_size_Cpre(i);
    		t_size_Cpost(i) = t_size_Cpre(i);
    	}
    	DynamicMultiArray<Real,dim> Cpost(t_size_Cpost);

    	const TensorIndex<dim-k+1> wgt_alpha_post = MultiArrayUtils<dim-k+1>::compute_weight(t_size_alpha_post);
    	const TensorIndex<k-1> wgt_theta_post = MultiArrayUtils<k-1>::compute_weight(t_size_theta_post);

    	const Size f_size_theta_post = t_size_theta_post.flat_size();
    	const Size f_size_alpha_post = t_size_alpha_post.flat_size();
		TensorIndex<dim> t_id_Cpre;

		const Index n_quad_pts = t_size_Cpre(dir);
    	const auto & w_times_B = omega_B_1d[dir];

    	Index f_id_Cpost = 0;
    	for (Index f_id_theta_post = 0 ; f_id_theta_post < f_size_theta_post ; ++f_id_theta_post)
    	{
    		const TensorIndex<k-1> t_id_theta_post =
    				MultiArrayUtils<k-1>::flat_to_tensor_index(f_id_theta_post,wgt_theta_post);

        	for (Index f_id_alpha_post = 0 ; f_id_alpha_post < f_size_alpha_post ; ++f_id_alpha_post)
        	{
        		const TensorIndex<dim-k+1> t_id_alpha_post =
        				MultiArrayUtils<dim-k+1>::flat_to_tensor_index(f_id_alpha_post,wgt_alpha_post);

        		const Index alpha = t_id_alpha_post(dim-k);
        		const auto w_times_B_alpha = w_times_B.get_function_view(alpha);


        		for (int i = 0 ; i < dir ;++ i)
        			t_id_Cpre(i) = t_id_alpha_post(i);
        		for (int i = dir+1 ; i < dim ;++ i)
        			t_id_Cpre(i) = t_id_theta_post(i-dir-1);

        		Real sum = 0.0;
        		for (int pt = 0 ; pt < n_quad_pts ; ++pt)
        		{
        			t_id_Cpre(dir) = pt;
        			sum += w_times_B_alpha[pt] * Cpre(t_id_Cpre);
        		}
        		Cpost(f_id_Cpost) = sum;
        		f_id_Cpost++;
        	} // end loop f_id_alpha_post

    	} // end loop f_id_theta_post

        IntegratorSumFacRHS_ver3<dim,k-1> integrate_rhs;
        integrate_rhs(Cpost,omega_B_1d,integral_rhs);
    }
};



template <int dim>
class IntegratorSumFacRHS_ver3<dim,1>
{
public:
    void operator()(
        const DynamicMultiArray<Real,dim> &Cpre,
        const array<ValueTable<Real>,dim> &omega_B_1d,
        DenseVector &integral_rhs)
    {
#ifndef NDEBUG
    	int n_basis = 1;
    	for (int i = 0 ; i < dim ; ++i)
    		n_basis *= omega_B_1d[i].get_num_functions();
    	Assert(integral_rhs.size() == n_basis,
    			ExcDimensionMismatch(integral_rhs.size(),n_basis));
#endif
    	const int dir = dim-1;

    	const TensorSize<dim> t_size_Cpre = Cpre.tensor_size();
    	TensorSize<dim> t_size_Cpost;

    	TensorSize<dim> t_size_alpha_post;
    	for (Index i = 0 ; i < dim-1 ; ++i)
    	{
    		t_size_alpha_post(i) = t_size_Cpre(i);
    		t_size_Cpost(i) = t_size_alpha_post(i);
    	}
    	t_size_alpha_post(dim-1) = omega_B_1d[dim-1].get_num_functions();
		t_size_Cpost(dim-1) = t_size_alpha_post(dim-1);


    	const TensorIndex<dim> wgt_alpha_post = MultiArrayUtils<dim>::compute_weight(t_size_alpha_post);

    	const Size f_size_alpha_post = t_size_alpha_post.flat_size();
		TensorIndex<dim> t_id_Cpre;

		const Index n_quad_pts = t_size_Cpre(dir);
    	const auto & w_times_B = omega_B_1d[dir];

    	Index f_id_Cpost = 0;
        for (Index f_id_alpha_post = 0 ; f_id_alpha_post < f_size_alpha_post ; ++f_id_alpha_post)
        {
        	const TensorIndex<dim> t_id_alpha_post =
        			MultiArrayUtils<dim>::flat_to_tensor_index(f_id_alpha_post,wgt_alpha_post);

        	const Index alpha = t_id_alpha_post(dim-1);
        	const auto w_times_B_alpha = w_times_B.get_function_view(alpha);


        	for (int i = 0 ; i < dir ;++ i)
        		t_id_Cpre(i) = t_id_alpha_post(i);

        	Real sum = 0.0;
        	for (int pt = 0 ; pt < n_quad_pts ; ++pt)
        	{
        		t_id_Cpre(dir) = pt;
        		sum += w_times_B_alpha[pt] * Cpre(t_id_Cpre);
        	}
        	integral_rhs(f_id_Cpost) = sum;
        	f_id_Cpost++;
        } // end loop f_id_alpha_post

    }
};



template <int dim, int r=dim>
class MassMatrixIntegrator
{
public:
    void
    operator()(
        const bool is_symmetric,
        const TensorSize<dim> &t_size_theta,
        const TensorSize<dim> &t_size_alpha,
        const TensorSize<dim> &t_size_beta,
        const array<DynamicMultiArray<Real,3>,dim> &J,
        const DynamicMultiArray<Real,3> &Cpre,
    	DenseMatrix &local_mass_matrix) const
    {
        const int k = dim-r+1;

        // (alpha_1,...alpha_{k-1})
        TensorSize<k-1> t_size_alpha_1_km1;
        // (beta_1,...beta_{k-1})
        TensorSize<k-1> t_size_beta_1_km1;
        for (int i = 0 ; i < k-1 ; ++i)
        {
            t_size_alpha_1_km1[i] = t_size_alpha[i];
            t_size_beta_1_km1 [i] = t_size_beta[i];
        }

        // (alpha_1,...alpha_k)
        TensorSize<k> t_size_alpha_1_k;
        // (beta_1,...beta_k)
        TensorSize<k> t_size_beta_1_k;
        for (int i = 0 ; i < k ; ++i)
        {
            t_size_alpha_1_k[i] = t_size_alpha[i];
            t_size_beta_1_k [i] = t_size_beta[i];
        }

        // (theta_{k+1},...theta_{dim})
        TensorSize<dim-k> t_size_theta_kp1_d;
        for (int i = 0 ; i < dim-k ; ++i)
            t_size_theta_kp1_d[i] = t_size_theta[i+k];


        // (theta_k,...theta_{dim})
        TensorSize<dim-k+1> t_size_theta_k_d;
        for (int i = 0 ; i <= dim-k ; ++i)
            t_size_theta_k_d[i] = t_size_theta[i+k-1];


        const Size f_size_theta_kp1_d = (dim-k>0)?t_size_theta_kp1_d.flat_size():1;
        const Size f_size_alpha_1_km1 = t_size_alpha_1_km1.flat_size();
        const Size f_size_beta_1_km1  = t_size_beta_1_km1.flat_size();


        const auto & Jk = J[k-1];
        TensorSize<3> t_size_Jk = Jk.tensor_size();
        Assert(t_size_Jk(0) == t_size_theta[k-1],ExcDimensionMismatch(t_size_Jk(0),t_size_theta[k-1]));
        Assert(t_size_Jk(1) == t_size_alpha[k-1],ExcDimensionMismatch(t_size_Jk(1),t_size_alpha[k-1]));
        Assert(t_size_Jk(2) == t_size_beta [k-1],ExcDimensionMismatch(t_size_Jk(2),t_size_beta [k-1]));
        TensorIndex<3> t_wgt_Jk = MultiArrayUtils<3>::compute_weight(t_size_Jk);


        const TensorSize<3> t_size_Cpre = Cpre.tensor_size();
        Assert(t_size_Cpre(0) == t_size_theta_k_d.flat_size(),
        		ExcDimensionMismatch(t_size_Cpre(0),t_size_theta_k_d.flat_size()));
        Assert(t_size_Cpre(1) == t_size_alpha_1_km1.flat_size(),
        		ExcDimensionMismatch(t_size_Cpre(1),t_size_alpha_1_km1.flat_size()));
        Assert(t_size_Cpre(2) == t_size_beta_1_km1.flat_size(),
        		ExcDimensionMismatch(t_size_Cpre(2),t_size_beta_1_km1.flat_size()));
        TensorIndex<3> t_wgt_Cpre = MultiArrayUtils<3>::compute_weight(t_size_Cpre);


        TensorSize<3> t_size_Cpost;
        t_size_Cpost[0] = f_size_theta_kp1_d;
        t_size_Cpost[1] = t_size_alpha_1_k.flat_size();
        t_size_Cpost[2] = t_size_beta_1_k.flat_size();
        TensorIndex<3> t_wgt_Cpost = MultiArrayUtils<3>::compute_weight(t_size_Cpost);
        DynamicMultiArray<Real,3> Cpost(t_size_Cpost);


        TensorIndex<3> tid_Jk;
        TensorIndex<3> tid_Cpre;
        TensorIndex<3> tid_Cpost;

        if (!is_symmetric)
        {

            tid_Cpost[2] = 0;
            tid_Jk[0] = 0;
            for (Index flat_beta_k_1 = 0 ; flat_beta_k_1 < f_size_beta_1_km1 ; ++flat_beta_k_1)
            {
                tid_Cpre[2] = flat_beta_k_1;

                for (int beta_k = 0 ; beta_k < t_size_beta[k-1] ; ++beta_k)
                {
                    tid_Jk[2] = beta_k;
                    tid_Cpost[1] = 0 ;

                    for (Index flat_alpha_k_1 = 0 ; flat_alpha_k_1 < f_size_alpha_1_km1 ; ++flat_alpha_k_1)
                    {
                        tid_Cpre[1] = flat_alpha_k_1;

                        for (int alpha_k = 0 ; alpha_k < t_size_alpha[k-1] ; ++alpha_k)
                        {
                            tid_Jk[1] = alpha_k;
                            tid_Cpre[0] = 0 ;
                            tid_Cpost[0] = 0 ;
                            for (Index fid_theta_kp1_d = 0 ; fid_theta_kp1_d < f_size_theta_kp1_d ; ++fid_theta_kp1_d)
                            {
                            	const Index f_id_Cpost = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpost,t_wgt_Cpost);

                            	const Index f_id_Cpre = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpre,t_wgt_Cpre);

                            	const Index f_id_Jk = MultiArrayUtils<3>::tensor_to_flat_index(tid_Jk,t_wgt_Jk);

                                Cpost(f_id_Cpost) = std::inner_product(
                                                  &Jk(f_id_Jk),
                                                  &Jk(f_id_Jk)+t_size_theta[k-1],
                                                  &Cpre(f_id_Cpre),
                                                  0.0);

                            	tid_Cpost[0]++;
                            	tid_Cpre [0] += t_size_theta[k-1];
                            } //end loop flat_theta_k_1

                            tid_Cpost[1]++;

                        } // end loop alpha_k
                    } // end loop flat_alpha_k_1

                    tid_Cpost[2]++;
                } // end loop beta_k

            } // end loop flat_beta_k_1


        } // end if(!is_symmetric)
        else
        {
            using MAUtils_k = MultiArrayUtils<k>;
            using MAUtils_km1 = MultiArrayUtils<k-1>;

            const Size f_size_alpha_1_k = t_size_alpha_1_k.flat_size();
            const TensorIndex<k> wgt_alpha_1_k = MAUtils_k::compute_weight(t_size_alpha_1_k);
            const TensorIndex<k-1> wgt_alpha_1_km1 = MAUtils_km1::compute_weight(t_size_alpha_1_km1);

            const Size f_size_beta_1_k = t_size_beta_1_k.flat_size();
            const TensorIndex<k> wgt_beta_1_k = MAUtils_k::compute_weight(t_size_beta_1_k);
            const TensorIndex<k-1> wgt_beta_1_km1 = MAUtils_km1::compute_weight(t_size_beta_1_km1);

            TensorIndex<k-1> tid_alpha_1_km1;
            TensorIndex<k-1> tid_beta_1_km1;
            for (Index fid_beta_1_k = 0 ; fid_beta_1_k < f_size_beta_1_k ; ++fid_beta_1_k)
            {
                const TensorIndex<k> tid_beta_1_k =
                    MAUtils_k::flat_to_tensor_index(fid_beta_1_k,wgt_beta_1_k);

                for (int i = 0 ; i < k-1 ; ++i)
                    tid_beta_1_km1(i) = tid_beta_1_k(i);

                const Index beta_k = tid_beta_1_k(k-1);

                const Index fid_beta_1_km1 =
                    (k>1)?MAUtils_km1::tensor_to_flat_index(tid_beta_1_km1,wgt_beta_1_km1):0;


                for (Index fid_alpha_1_k = fid_beta_1_k ; fid_alpha_1_k < f_size_alpha_1_k ; ++fid_alpha_1_k)
                {
                    const TensorIndex<k> tid_alpha_1_k =
                        MAUtils_k::flat_to_tensor_index(fid_alpha_1_k,wgt_alpha_1_k);
                    for (int i = 0 ; i < k-1 ; ++i)
                        tid_alpha_1_km1(i) = tid_alpha_1_k(i);

                    const Index alpha_k = tid_alpha_1_k(k-1);

                    const Index fid_alpha_1_km1 =
                        (k>1)?MAUtils_km1::tensor_to_flat_index(tid_alpha_1_km1,wgt_alpha_1_km1):0;

                    tid_Cpre[1] = max(fid_alpha_1_km1,fid_beta_1_km1);
                    tid_Cpre[2] = min(fid_alpha_1_km1,fid_beta_1_km1);

                    tid_Jk[1] = max(alpha_k,beta_k);
                    tid_Jk[2] = min(alpha_k,beta_k);

                    tid_Cpost[2] = fid_beta_1_k ;
                    tid_Cpost[1] = fid_alpha_1_k;

                    tid_Cpre [0] = 0;
                    tid_Cpost[0] = 0;
                    for (Index fid_theta_kp1_d = 0 ; fid_theta_kp1_d < f_size_theta_kp1_d ; ++fid_theta_kp1_d)
                    {
                    	const Index f_id_Cpost = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpost,t_wgt_Cpost);

                    	const Index f_id_Cpre = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpre,t_wgt_Cpre);

                    	const Index f_id_Jk = MultiArrayUtils<3>::tensor_to_flat_index(tid_Jk,t_wgt_Jk);

                        Cpost(f_id_Cpost) = std::inner_product(
                                          &Jk(f_id_Jk),
                                          &Jk(f_id_Jk)+t_size_theta[k-1],
                                          &Cpre(f_id_Cpre),
                                          0.0);


                    	tid_Cpost[0]++;
                    	tid_Cpre [0] += t_size_theta[k-1];
                    } //end loop flat_theta_k_1

                }// end loop fid_alpha_1_k

            } // end loop fid_beta_1_k
            //*/

        } // end if (symmetric)

        MassMatrixIntegrator<dim,r-1> mass_matrix_integrator;
        mass_matrix_integrator(
        		   is_symmetric,
                   t_size_theta,
                   t_size_alpha,
                   t_size_beta,
                   J,Cpost,
                   local_mass_matrix);
    }
};


template <int dim>
class MassMatrixIntegrator<dim,1>
{
public:
    void operator()(
        const bool is_symmetric,
        const TensorSize<dim> &t_size_theta,
        const TensorSize<dim> &t_size_alpha,
        const TensorSize<dim> &t_size_beta,
        const array<DynamicMultiArray<Real,3>,dim> &J,
        const DynamicMultiArray<Real,3> &Cpre,
    	DenseMatrix & local_mass_matrix) const
    {
        const int k = dim;

        // (alpha_1,...alpha_{k-1})
        TensorSize<k-1> t_size_alpha_1_km1;
        // (beta_1,...beta_{k-1})
        TensorSize<k-1> t_size_beta_1_km1;
        for (int i = 0 ; i < k-1 ; ++i)
        {
            t_size_alpha_1_km1[i] = t_size_alpha[i];
            t_size_beta_1_km1 [i] = t_size_beta[i];
        }

        // (alpha_1,...alpha_k)
        TensorSize<k> t_size_alpha_1_k;
        // (beta_1,...beta_k)
        TensorSize<k> t_size_beta_1_k;
        for (int i = 0 ; i < k ; ++i)
        {
            t_size_alpha_1_k[i] = t_size_alpha[i];
            t_size_beta_1_k [i] = t_size_beta[i];
        }

        // (theta_{k+1},...theta_{dim})
        TensorSize<dim-k> t_size_theta_kp1_d;
        for (int i = 0 ; i < dim-k ; ++i)
            t_size_theta_kp1_d[i] = t_size_theta[i+k];


        // (theta_k,...theta_{dim})
        TensorSize<dim-k+1> t_size_theta_k_d;
        for (int i = 0 ; i <= dim-k ; ++i)
            t_size_theta_k_d[i] = t_size_theta[i+k-1];


        const Size f_size_alpha_1_km1 = t_size_alpha_1_km1.flat_size();
        const Size f_size_beta_1_km1  = t_size_beta_1_km1.flat_size();


        const auto & Jk = J[k-1];
        TensorSize<3> t_size_Jk = Jk.tensor_size();
        Assert(t_size_Jk(0) == t_size_theta[k-1],ExcDimensionMismatch(t_size_Jk(0),t_size_theta[k-1]));
        Assert(t_size_Jk(1) == t_size_alpha[k-1],ExcDimensionMismatch(t_size_Jk(1),t_size_alpha[k-1]));
        Assert(t_size_Jk(2) == t_size_beta [k-1],ExcDimensionMismatch(t_size_Jk(2),t_size_beta [k-1]));
        TensorIndex<3> t_wgt_Jk = MultiArrayUtils<3>::compute_weight(t_size_Jk);


        const TensorSize<3> t_size_Cpre = Cpre.tensor_size();
        Assert(t_size_Cpre(0) == t_size_theta_k_d.flat_size(),
        		ExcDimensionMismatch(t_size_Cpre(0),t_size_theta_k_d.flat_size()));
        Assert(t_size_Cpre(1) == t_size_alpha_1_km1.flat_size(),
        		ExcDimensionMismatch(t_size_Cpre(1),t_size_alpha_1_km1.flat_size()));
        Assert(t_size_Cpre(2) == t_size_beta_1_km1.flat_size(),
        		ExcDimensionMismatch(t_size_Cpre(2),t_size_beta_1_km1.flat_size()));
        TensorIndex<3> t_wgt_Cpre = MultiArrayUtils<3>::compute_weight(t_size_Cpre);


        TensorIndex<3> tid_Jk;
        TensorIndex<3> tid_Cpre;
        TensorIndex<3> tid_Cpost;


        const Size f_size_alpha_1_k = t_size_alpha_1_k.flat_size();
        const Size f_size_beta_1_k = t_size_beta_1_k.flat_size();

        Assert(local_mass_matrix.get_num_rows() == f_size_beta_1_k,
        		ExcDimensionMismatch(local_mass_matrix.get_num_rows(),f_size_beta_1_k));
        Assert(local_mass_matrix.get_num_cols() == f_size_alpha_1_k,
        		ExcDimensionMismatch(local_mass_matrix.get_num_cols(),f_size_alpha_1_k));

        if (!is_symmetric)
        {

            tid_Jk[0] = 0;
            tid_Cpre[0] = 0 ;
            for (Index flat_beta_k_1 = 0 ; flat_beta_k_1 < f_size_beta_1_km1 ; ++flat_beta_k_1)
            {
                tid_Cpre[2] = flat_beta_k_1;

                for (int beta_k = 0 ; beta_k < t_size_beta[k-1] ; ++beta_k)
                {
                    tid_Jk[2] = beta_k;

                    const Index f_id_test = flat_beta_k_1 + beta_k;

                    for (Index flat_alpha_k_1 = 0 ; flat_alpha_k_1 < f_size_alpha_1_km1 ; ++flat_alpha_k_1)
                    {
                        tid_Cpre[1] = flat_alpha_k_1;

                        for (int alpha_k = 0 ; alpha_k < t_size_alpha[k-1] ; ++alpha_k)
                        {
                            tid_Jk[1] = alpha_k;

                            const Index f_id_trial = flat_alpha_k_1 + alpha_k;

                            const Index f_id_Cpre = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpre,t_wgt_Cpre);

                            const Index f_id_Jk = MultiArrayUtils<3>::tensor_to_flat_index(tid_Jk,t_wgt_Jk);

                            local_mass_matrix(f_id_test,f_id_trial) =
                            		std::inner_product(
                                              &Jk(f_id_Jk),
                                              &Jk(f_id_Jk)+t_size_theta[k-1],
                                              &Cpre(f_id_Cpre),
                                              0.0);
                        } // end loop alpha_k
                    } // end loop flat_alpha_k_1

                } // end loop beta_k

            } // end loop flat_beta_k_1


        } // end if(!is_symmetric)
        else
        {
            using MAUtils_k = MultiArrayUtils<k>;
            using MAUtils_km1 = MultiArrayUtils<k-1>;

            const TensorIndex<k> wgt_alpha_1_k = MAUtils_k::compute_weight(t_size_alpha_1_k);
            const TensorIndex<k-1> wgt_alpha_1_km1 = MAUtils_km1::compute_weight(t_size_alpha_1_km1);

            const TensorIndex<k> wgt_beta_1_k = MAUtils_k::compute_weight(t_size_beta_1_k);
            const TensorIndex<k-1> wgt_beta_1_km1 = MAUtils_km1::compute_weight(t_size_beta_1_km1);

            TensorIndex<k-1> tid_alpha_1_km1;
            TensorIndex<k-1> tid_beta_1_km1;
            tid_Cpre [0] = 0;
            for (Index fid_beta_1_k = 0 ; fid_beta_1_k < f_size_beta_1_k ; ++fid_beta_1_k)
            {
                const TensorIndex<k> tid_beta_1_k =
                    MAUtils_k::flat_to_tensor_index(fid_beta_1_k,wgt_beta_1_k);

                for (int i = 0 ; i < k-1 ; ++i)
                    tid_beta_1_km1(i) = tid_beta_1_k(i);

                const Index beta_k = tid_beta_1_k(k-1);

                const Index fid_beta_1_km1 =
                    (k>1)?MAUtils_km1::tensor_to_flat_index(tid_beta_1_km1,wgt_beta_1_km1):0;


                for (Index fid_alpha_1_k = fid_beta_1_k ; fid_alpha_1_k < f_size_alpha_1_k ; ++fid_alpha_1_k)
                {
                    const TensorIndex<k> tid_alpha_1_k =
                        MAUtils_k::flat_to_tensor_index(fid_alpha_1_k,wgt_alpha_1_k);
                    for (int i = 0 ; i < k-1 ; ++i)
                        tid_alpha_1_km1(i) = tid_alpha_1_k(i);

                    const Index alpha_k = tid_alpha_1_k(k-1);

                    const Index fid_alpha_1_km1 =
                        (k>1)?MAUtils_km1::tensor_to_flat_index(tid_alpha_1_km1,wgt_alpha_1_km1):0;

                    tid_Cpre[1] = max(fid_alpha_1_km1,fid_beta_1_km1);
                    tid_Cpre[2] = min(fid_alpha_1_km1,fid_beta_1_km1);

                    tid_Jk[1] = max(alpha_k,beta_k);
                    tid_Jk[2] = min(alpha_k,beta_k);

                    const Index f_id_Cpre = MultiArrayUtils<3>::tensor_to_flat_index(tid_Cpre,t_wgt_Cpre);

                    const Index f_id_Jk = MultiArrayUtils<3>::tensor_to_flat_index(tid_Jk,t_wgt_Jk);

                    local_mass_matrix(fid_beta_1_k,fid_alpha_1_k) =
                    		std::inner_product(
                                          &Jk(f_id_Jk),
                                          &Jk(f_id_Jk)+t_size_theta[k-1],
                                          &Cpre(f_id_Cpre),
                                          0.0);

                }// end loop fid_alpha_1_k

            } // end loop fid_beta_1_k
            //*/

            // here we copy the upper triangular part of the matrix on the lower triangular part
            for (int test_id = 0 ; test_id < f_size_beta_1_k ; ++test_id)
                for (int trial_id = 0; trial_id < test_id ; ++trial_id)
                    local_mass_matrix(test_id,trial_id) = local_mass_matrix(trial_id,test_id);

        } // end if (symmetric)
    }
};

#define SPECIALIZED
#ifdef SPECIALIZED
template <>
class MassMatrixIntegrator<1,1>
{
public:
    void operator()(
    	const bool is_symmetric,
        const TensorSize<1> &t_size_theta,
        const TensorSize<1> &t_size_alpha,
        const TensorSize<1> &t_size_beta,
        const array<DynamicMultiArray<Real,3>,1> &J,
        const DynamicMultiArray<Real,3> &C,
    	DenseMatrix & local_mass_matrix) const
    {
    	const Size n_basis_test = local_mass_matrix.get_num_rows();
    	const Size n_basis_trial = local_mass_matrix.get_num_cols();

        Assert(t_size_alpha.flat_size() == n_basis_trial,
               ExcDimensionMismatch(t_size_alpha.flat_size(),n_basis_trial));
        Assert(t_size_beta.flat_size() == n_basis_test,
               ExcDimensionMismatch(t_size_beta.flat_size(),n_basis_test));


        TensorIndex<3> t_id_J;
        TensorIndex<3> t_id_C;

        if (!is_symmetric)
        {
            //--------------------------------------------------------------
            t_id_C[1] = 0;
            t_id_C[2] = 0;

            for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
            {
                t_id_J[2] = beta_0;
                for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                {
                    t_id_J[1] = alpha_0;

                    Real sum = 0.0;
                    for (Index theta_0 = 0; theta_0 < t_size_theta[0] ; ++theta_0)
                    {
                        t_id_J[0] = theta_0;
                        t_id_C[0] = theta_0;
                        sum += C(t_id_C) * J[0](t_id_J);
                    }

                    local_mass_matrix(beta_0,alpha_0) = sum;
                }
            }
            //--------------------------------------------------------------
        } // end if (!is_symmetric)
        else
        {
            //--------------------------------------------------------------
            t_id_C[1] = 0;
            t_id_C[2] = 0;
            for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
            {
                t_id_J[2] = beta_0;
                for (Index alpha_0 = beta_0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                {
                    t_id_J[1] = alpha_0;

                    Real sum = 0.0;
                    for (Index theta_0 = 0; theta_0 < t_size_theta[0] ; ++theta_0)
                    {
                        t_id_J[0] = theta_0;
                        t_id_C[0] = theta_0;
                        sum += C(t_id_C) * J[0](t_id_J);
                    }

                    local_mass_matrix(beta_0,alpha_0) = sum;
                }
            }

            // here we copy the upper triangular part of the matrix on the lower triangular part
            for (int test_id = 0 ; test_id < n_basis_test ; ++test_id)
                for (int trial_id = 0; trial_id < test_id ; ++trial_id)
                    local_mass_matrix(test_id,trial_id) = local_mass_matrix(trial_id,test_id);

            //--------------------------------------------------------------

        } // end if (is_symmetric)
    }
};

template <>
class MassMatrixIntegrator<2,2>
{
public:
    void operator()(
        const bool is_symmetric,
        const TensorSize<2> &t_size_theta,
        const TensorSize<2> &t_size_alpha,
        const TensorSize<2> &t_size_beta,
        const array<DynamicMultiArray<Real,3>,2> &J,
        const DynamicMultiArray<Real,3> &C,
    	DenseMatrix & local_mass_matrix) const
    {
    	const Size n_basis_test = local_mass_matrix.get_num_rows();
    	const Size n_basis_trial = local_mass_matrix.get_num_cols();

        Assert(t_size_alpha.flat_size() == n_basis_trial,
               ExcDimensionMismatch(t_size_alpha.flat_size(),n_basis_trial));
        Assert(t_size_beta.flat_size() == n_basis_test,
               ExcDimensionMismatch(t_size_beta.flat_size(),n_basis_test));


        TensorIndex<3> t_id_J;
        TensorIndex<3> t_id_C;


        //--------------------------------------------------------------
        TensorIndex<3> t_id_C1;
        TensorSize<3> t_size_C1;
        t_size_C1[0] = t_size_theta[1];
        t_size_C1[1] = t_size_alpha[0];
        t_size_C1[2] = t_size_beta[0];
        DynamicMultiArray<Real,3> C1(t_size_C1);
        //--------------------------------------------------------------


        //--------------------------------------------------------------
        t_id_C[1] = 0;
        t_id_C[2] = 0;
        for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
        {
            t_id_J [2] = beta_0;
            t_id_C1[2] = beta_0;
            for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
            {
                t_id_J [1] = alpha_0;
                t_id_C1[1] = alpha_0;

                Index theta_0_1 = 0;
                for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1)
                {
                    Real sum = 0.0;
                    for (Index theta_0 = 0; theta_0 < t_size_theta[0] ; ++theta_0,++theta_0_1)
                    {
                        t_id_J[0] = theta_0;
                        t_id_C[0] = theta_0_1;
                        sum += C(t_id_C) * J[0](t_id_J);
                    }

                    t_id_C1[0] = theta_1;
                    C1(t_id_C1) = sum;
                } //end loop theta_1
            } //end loop alpha_0
        } // end loop beta_0
        //--------------------------------------------------------------


        if (!is_symmetric)
        {
            //--------------------------------------------------------------
            for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
            {
                t_id_J[2] = beta_1;
                for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
                {
                    Index beta_0_1 = beta_1*t_size_beta[0] + beta_0;

                    t_id_C1[2] = beta_0;

                    for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                    {
                        t_id_J[1] = alpha_1;
                        for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                        {
                            Index alpha_0_1 = alpha_1*t_size_alpha[0] + alpha_0;

                            t_id_C1[1] = alpha_0;

                            Real sum = 0.0;
                            for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1)
                            {
                                t_id_J [0] = theta_1;
                                t_id_C1[0] = theta_1;
                                sum += C1(t_id_C1) * J[1](t_id_J);
                            } // end loop theta_1

                            local_mass_matrix(beta_0_1,alpha_0_1) = sum;
                        } //end loop alpha_0
                    } //end loop alpha_1
                } // end loop beta_0
            } // end loop beta_1
            //--------------------------------------------------------------
        }//end if (!is_symmetric)
        else
        {
            //--------------------------------------------------------------
            for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
            {
                t_id_J[2] = beta_1;
                for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
                {
                    Index beta_0_1 = beta_1*t_size_beta[0] + beta_0;

                    t_id_C1[2] = beta_0;

                    for (Index alpha_1 = beta_1 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                    {
                        t_id_J[1] = alpha_1;
                        for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                        {
                            Index alpha_0_1 = alpha_1*t_size_alpha[0] + alpha_0;

                            t_id_C1[1] = alpha_0;

                            Real sum = 0.0;
                            for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1)
                            {
                                t_id_J [0] = theta_1;
                                t_id_C1[0] = theta_1;
                                sum += C1(t_id_C1) * J[1](t_id_J);
                            } // end loop theta_1

                            local_mass_matrix(beta_0_1,alpha_0_1) = sum;
                        } //end loop alpha_0
                    } //end loop alpha_1
                } // end loop beta_0
            } // end loop beta_1


            // here we copy the upper triangular part of the matrix on the lower triangular part
            for (int test_id = 0 ; test_id < n_basis_test ; ++test_id)
                for (int trial_id = 0; trial_id < test_id ; ++trial_id)
                    local_mass_matrix(test_id,trial_id) = local_mass_matrix(trial_id,test_id);

            //--------------------------------------------------------------
        }//end if (is_symmetric)
    }
};


template <>
class MassMatrixIntegrator<3,3>
{
public:
    void operator()(
        const bool is_symmetric,
        const TensorSize<3> &t_size_theta,
        const TensorSize<3> &t_size_alpha,
        const TensorSize<3> &t_size_beta,
        const array<DynamicMultiArray<Real,3>,3> &J,
        const DynamicMultiArray<Real,3> &C,
        DenseMatrix & local_mass_matrix) const
    {
    	const Size n_basis_test = local_mass_matrix.get_num_rows();
    	const Size n_basis_trial = local_mass_matrix.get_num_cols();

        Assert(t_size_alpha.flat_size() == n_basis_trial,
               ExcDimensionMismatch(t_size_alpha.flat_size(),n_basis_trial));
        Assert(t_size_beta.flat_size() == n_basis_test,
               ExcDimensionMismatch(t_size_beta.flat_size(),n_basis_test));



        TensorIndex<3> t_id_J;
        TensorIndex<3> t_id_C;

        //--------------------------------------------------------------
        TensorIndex<3> t_id_C1;
        TensorSize<3> t_size_C1;
        t_size_C1[0] = t_size_theta[1] * t_size_theta[2];
        t_size_C1[1] = t_size_alpha[0];
        t_size_C1[2] = t_size_beta[0];
        DynamicMultiArray<Real,3> C1(t_size_C1);
        //--------------------------------------------------------------


        //--------------------------------------------------------------
        TensorIndex<3> t_id_C2;
        TensorSize<3> t_size_C2;
        t_size_C2[0] = t_size_theta[2];
        t_size_C2[1] = t_size_alpha[0] * t_size_alpha[1];
        t_size_C2[2] = t_size_beta [0] * t_size_beta [1];
        DynamicMultiArray<Real,3> C2(t_size_C2);
        //--------------------------------------------------------------


        //--------------------------------------------------------------
        t_id_C[1] = 0;
        t_id_C[2] = 0;
        for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
        {
            t_id_J [2] = beta_0;
            t_id_C1[2] = beta_0;
            for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
            {
                t_id_J [1] = alpha_0;
                t_id_C1[1] = alpha_0;

                Index theta_1_2 = 0;
                Index theta_0_1_2 = 0;
                for (Index theta_2 = 0; theta_2 < t_size_theta[2] ; ++theta_2)
                {
                    for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1,++theta_1_2)
                    {
                        Real sum = 0.0;
                        for (Index theta_0 = 0; theta_0 < t_size_theta[0] ; ++theta_0,++theta_0_1_2)
                        {
                            t_id_J[0] = theta_0;
                            t_id_C[0] = theta_0_1_2;
                            sum += C(t_id_C) * J[0](t_id_J);
                        }

                        t_id_C1[0] = theta_1_2;
                        C1(t_id_C1) = sum;
                    } // end loop theta_1
                } // end loop theta_2
            } // end loop alpha_0
        } // end loop beta_0
        //--------------------------------------------------------------



        //--------------------------------------------------------------
        for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
        {
            t_id_J[2] = beta_1;
            for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
            {
                Index beta_0_1 = beta_1*t_size_beta[0] + beta_0;

                t_id_C1[2] = beta_0;
                t_id_C2[2] = beta_0_1;

                for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                {
                    t_id_J[1] = alpha_1;
                    for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                    {
                        Index alpha_0_1 = alpha_1*t_size_alpha[0] + alpha_0;

                        t_id_C1[1] = alpha_0;
                        t_id_C2[1] = alpha_0_1;

                        Index theta_1_2 = 0;
                        for (Index theta_2 = 0; theta_2 < t_size_theta[2] ; ++theta_2)
                        {
                            Real sum = 0.0;
                            for (Index theta_1 = 0; theta_1 < t_size_theta[1] ; ++theta_1,++theta_1_2)
                            {
                                t_id_J [0] = theta_1;
                                t_id_C1[0] = theta_1_2;
                                sum += C1(t_id_C1) * J[1](t_id_J);
                            } // end loop theta_1

                            t_id_C2[0] = theta_2;
                            C2(t_id_C2) = sum;
                        } // end loop theta_2
                    } //end loop alpha_0
                } //end loop alpha_1
            } // end loop beta_0
        } // end loop beta_1
        //--------------------------------------------------------------

        if (!is_symmetric)
        {

            //--------------------------------------------------------------
            for (Index beta_2 = 0 ; beta_2 < t_size_beta[2] ; ++beta_2)
            {
                t_id_J[2] = beta_2;
                for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
                {
                    for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0)
                    {
                        Index beta_0_1 = beta_1 * t_size_beta[0] + beta_0 ;
                        Index beta_0_1_2 = (beta_2*t_size_beta[1] + beta_1) * t_size_beta[0] + beta_0 ;

                        t_id_C2[2] = beta_0_1 ;

                        for (Index alpha_2 = 0 ; alpha_2 < t_size_alpha[2] ; ++alpha_2)
                        {
                            t_id_J[1] = alpha_2;
                            for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                            {
                                for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0)
                                {
                                    Index alpha_0_1 = alpha_1 * t_size_alpha[0] + alpha_0 ;
                                    Index alpha_0_1_2 = (alpha_2*t_size_alpha[1] + alpha_1) * t_size_alpha[0] + alpha_0 ;

                                    t_id_C2[1] = alpha_0_1 ;

                                    Real sum = 0.0;
                                    for (Index theta_2 = 0; theta_2 < t_size_theta[2] ; ++theta_2)
                                    {
                                        t_id_J [0] = theta_2;
                                        t_id_C2[0] = theta_2;
                                        sum += C2(t_id_C2) * J[2](t_id_J);
                                    } // end loop theta_1

                                    local_mass_matrix(beta_0_1_2,alpha_0_1_2) = sum;

                                } // end loop alpha_0
                            } // end loop alpha_1
                        } // end loop alpha_2
                    } // end loop beta_0
                } // end loop beta_1
            } // end loop beta_2
            //--------------------------------------------------------------
        } //end if (!is_symmetric)
        else
        {

            //--------------------------------------------------------------
            for (Index beta_2 = 0 ; beta_2 < t_size_beta[2] ; ++beta_2)
            {
                t_id_J[2] = beta_2;

                const Index b_tmp_2 = beta_2*t_size_beta[1];

                Index beta_0_1 = 0;
                for (Index beta_1 = 0 ; beta_1 < t_size_beta[1] ; ++beta_1)
                {
                    const Index b_tmp_1 = (b_tmp_2 + beta_1) * t_size_beta[0];

                    for (Index beta_0 = 0 ; beta_0 < t_size_beta[0] ; ++beta_0, ++beta_0_1)
                    {
                        Index beta_0_1_2 = b_tmp_1 + beta_0 ;

                        t_id_C2[2] = beta_0_1;

                        for (Index alpha_2 = beta_2 ; alpha_2 < t_size_alpha[2] ; ++alpha_2)
                        {
                            t_id_J[1] = alpha_2;

                            const Index a_tmp_2 = alpha_2*t_size_alpha[1];

                            Index alpha_0_1 = 0;
                            for (Index alpha_1 = 0 ; alpha_1 < t_size_alpha[1] ; ++alpha_1)
                            {
                                const Index a_tmp_1 = (a_tmp_2 + alpha_1) * t_size_alpha[0];

                                for (Index alpha_0 = 0 ; alpha_0 < t_size_alpha[0] ; ++alpha_0,++alpha_0_1)
                                {
                                    const Index alpha_0_1_2 = a_tmp_1 + alpha_0 ;

                                    t_id_C2[1] = alpha_0_1;

                                    Real sum = 0.0;
                                    for (Index theta_2 = 0; theta_2 < t_size_theta[2] ; ++theta_2)
                                    {
                                        t_id_J [0] = theta_2;
                                        t_id_C2[0] = theta_2;

                                        sum += C2(t_id_C2) * J[2](t_id_J);
                                    } // end loop theta_1

                                    local_mass_matrix(beta_0_1_2,alpha_0_1_2) = sum;
                                } // end loop alpha_0
                            } // end loop alpha_1
                        } // end loop alpha_2
                    } // end loop beta_0
                } // end loop beta_1
            } // end loop beta_2

            // here we copy the upper triangular part of the matrix on the lower triangular part
            for (int test_id = 0 ; test_id < n_basis_test ; ++test_id)
                for (int trial_id = 0; trial_id < test_id ; ++trial_id)
                    local_mass_matrix(test_id,trial_id) = local_mass_matrix(trial_id,test_id);

            //--------------------------------------------------------------
        } // end if (is_symmetric)
    }
};
#endif


template<class PhysSpaceTest, class PhysSpaceTrial>
inline
auto
EllipticOperatorsSumFactorizationIntegration<PhysSpaceTest,PhysSpaceTrial>::
projection_l2_bernstein_basis(
    const ValueVector<Real> &func_to_proj_at_pts) const -> DynamicMultiArray<Real,dim>
{
    Assert(inv_B_proj_.get_num_rows() == inv_B_proj_.get_num_cols(),
    		ExcDimensionMismatch(inv_B_proj_.get_num_rows(),inv_B_proj_.get_num_cols()));

    const Size n_basis = inv_B_proj_.get_num_rows();

    TensorSize<dim> n_basis_1D;
    TensorSize<dim> n_points_1D;
    for (int i = 0 ; i < dim ; ++i)
    {
        n_basis_1D(i) = w_times_B_proj_1D_[i].get_num_functions();
        n_points_1D(i) = w_times_B_proj_1D_[i].get_num_points();

        Assert(n_points_1D(i) == proj_quad_.get_num_points_direction()(i),
               ExcDimensionMismatch(n_points_1D(i),proj_quad_.get_num_points_direction()(i)));
    }

    Assert(n_basis_1D.flat_size() == n_basis,
           ExcDimensionMismatch(n_basis_1D.flat_size(),n_basis));

//    const TensorIndex<dim> basis_tensor_wgt =
//        MultiArrayUtils<dim>::compute_weight(n_basis_1D);

    const Size n_points = n_points_1D.flat_size();
    Assert(func_to_proj_at_pts.size() == n_points,
           ExcDimensionMismatch(func_to_proj_at_pts.size(),n_points));

    DynamicMultiArray<Real,dim> coeffs_to_project(n_points_1D);
    for (Index i = 0 ; i < n_points ; ++i)
    	coeffs_to_project(i) = func_to_proj_at_pts[i] ;


    DenseVector integral_rhs(n_basis);
/*
    IntegratorSumFacRHS_ver2<dim> integrate_rhs;
    for (Index basis_flat_id = 0 ; basis_flat_id < n_basis ; ++basis_flat_id)
    {
        const auto basis_tensor_id =
            MultiArrayUtils<dim>::flat_to_tensor_index(basis_flat_id,basis_tensor_wgt);

        integral_rhs(basis_flat_id) =
        		integrate_rhs(coeffs_to_project,w_times_B_proj_1D_,basis_tensor_id);
    }
//*/
    IntegratorSumFacRHS_ver3<dim> integrate_rhs;
    integrate_rhs(coeffs_to_project,w_times_B_proj_1D_,integral_rhs);


    // coefficient of the L2 projection using the tensor product basis
    DenseVector proj_coefs_tmp = boost::numeric::ublas::prod(inv_B_proj_, integral_rhs) ;


    DynamicMultiArray<Real,dim> proj_coefs(n_basis_1D);
    for (Index f_id_test = 0 ; f_id_test < n_basis ; ++f_id_test)
        proj_coefs(f_id_test) = proj_coefs_tmp(f_id_test);

    return proj_coefs;
}





template<class PhysSpaceTest, class PhysSpaceTrial>
inline
auto
EllipticOperatorsSumFactorizationIntegration<PhysSpaceTest,PhysSpaceTrial>::
evaluate_moments_op_u_v(
    const array<ValueTable<Function<1>::ValueType>,dim> &phi_1D_test,
    const array<ValueTable<Function<1>::ValueType>,dim> &phi_1D_trial,
    const array<Real,dim> &length_element_edge) const -> array<DynamicMultiArray<Real,3>,dim>
{
    array<DynamicMultiArray<Real,3>,dim> moments;

    for (int dir = 0 ; dir < dim ; ++dir)
    {
    	const auto & w_times_bernst = w_times_B_proj_1D_[dir];
    	const auto & phi_test  = phi_1D_test [dir];
    	const auto & phi_trial = phi_1D_trial[dir];

        const Size n_basis_projection = w_times_bernst.get_num_functions();
        const Size n_basis_test  = phi_test.get_num_functions();
        const Size n_basis_trial = phi_trial.get_num_functions();

        TensorSize<3> moments1D_tensor_size;
        moments1D_tensor_size[0] = n_basis_projection;
        moments1D_tensor_size[1] = n_basis_test;
        moments1D_tensor_size[2] = n_basis_trial;

        auto & moments1D = moments[dir];
        moments1D.resize(moments1D_tensor_size);

        const Size n_pts = w_times_bernst.get_num_points();
        Assert(phi_test.get_num_points() == n_pts,
               ExcDimensionMismatch(phi_test.get_num_points(),n_pts));
        Assert(phi_trial.get_num_points() == n_pts,
               ExcDimensionMismatch(phi_trial.get_num_points(),n_pts));

        vector<Real> phi_mu1_mu2(n_pts);

        const Real edge_length = length_element_edge[dir];

        Index flat_id_I = 0 ;
        for (int mu2 = 0 ; mu2 < n_basis_test ; ++mu2)
        {
            const auto phi_1D_mu2 = phi_test.get_function_view(mu2);

            for (int mu1 = 0 ; mu1 < n_basis_trial ; ++mu1)
            {
                const auto phi_1D_mu1 = phi_trial.get_function_view(mu1);

                for (int jpt = 0 ; jpt < n_pts ; ++jpt)
                    phi_mu1_mu2[jpt] = phi_1D_mu2[jpt](0) * phi_1D_mu1[jpt](0) * edge_length;

                for (int lambda = 0 ; lambda < n_basis_projection ; ++lambda)
                {
                    const auto w_B_lambda = w_times_bernst.get_function_view(lambda);
                    Real sum = 0.0;
                    for (int jpt = 0 ; jpt < n_pts ; ++jpt)
                        sum += w_B_lambda[jpt] * phi_mu1_mu2[jpt];
                    moments1D(flat_id_I++) = sum;
                } // end loop lambda
            } // end loop mu1
        } // end loop mu2
    } // end loop dir

    return moments;
}


template<class PhysSpaceTest, class PhysSpaceTrial>
inline
EllipticOperatorsSumFactorizationIntegration<PhysSpaceTest,PhysSpaceTrial>::
EllipticOperatorsSumFactorizationIntegration(
		const TensorSize<dim> &projection_degree,
		const Quadrature<dim> &projection_quadrature)
		:
		EllipticOperators<PhysSpaceTest,PhysSpaceTrial>(),
		proj_degree_(projection_degree),
		proj_quad_(projection_quadrature)
{
    //-----------------------------------------------------------------
    std::array<DenseMatrix,dim> B_proj_1D;
    for (int i = 0 ; i < dim ; ++i)
    {
    	n_bernst_1D_[i] = proj_degree_[i]+1;

        // values of the one-dimensional Bernstein basis at the projection pts
        B_proj_1D[i] = BernsteinBasis::evaluate(
        				   proj_degree_[i],
                           proj_quad_.get_points().get_data_direction(i));

        Assert(n_bernst_1D_[i] == B_proj_1D[i].size1(),
               ExcDimensionMismatch(n_bernst_1D_[i], B_proj_1D[i].size1()));

        Assert(proj_quad_.get_num_points_direction()[i] == B_proj_1D[i].size2(),
               ExcDimensionMismatch(proj_quad_.get_num_points_direction()[i], B_proj_1D[i].size2()));
    }
    n_basis_proj_ = n_bernst_1D_.flat_size();
    //-----------------------------------------------------------------


    //-----------------------------------------------------------------
    for (int i = 0 ; i < dim ; ++i)
    {
        const Size n_pts = proj_quad_.get_num_points_direction()[i];
        w_times_B_proj_1D_[i].resize(n_bernst_1D_[i],n_pts);

        const auto &B_i = B_proj_1D[i];

        auto &w_B_i = w_times_B_proj_1D_[i];

        const vector<Real> w_i = proj_quad_.get_weights().get_data_direction(i);

        for (Index ifn = 0 ; ifn < n_bernst_1D_[i] ; ++ifn)
        {
            auto w_B_ifn = w_B_i.get_function_view(ifn);

            for (Index jpt = 0 ; jpt < n_pts ; ++jpt)
                w_B_ifn[jpt] = w_i[jpt] * B_i(ifn,jpt);
        }
    }
    //-----------------------------------------------------------------



    //-----------------------------------------------------------------
    // evaluation of the Bernstein's mass matrix -- begin

    using boost::math::binomial_coefficient;

    B_proj_.resize(n_basis_proj_,n_basis_proj_);

    for (int dir = 0 ; dir < dim ; ++dir)
    {
    	B_proj_1D[dir].resize(n_bernst_1D_[dir],n_bernst_1D_[dir]);

    	const int q = proj_degree_[dir];
    	for (int i = 0 ; i < n_bernst_1D_[dir] ; ++i)
    	{
    		const Real tmp = binomial_coefficient<double>(q,i) / (2.0 * q + 1.0);

    		for (int j = 0 ; j < n_bernst_1D_[dir] ; ++j)
    		{
    			B_proj_1D[dir](i,j) =
    					tmp / binomial_coefficient<double>(2*q,i+j) * binomial_coefficient<double>(q,j) ;
    		}
    	}
    } // end loop dir


    using MAUtils = MultiArrayUtils<dim>;


    TensorIndex<dim> wgt_B_proj = MAUtils::compute_weight(n_bernst_1D_);
    for (Size row = 0 ; row < n_basis_proj_ ; ++row)
    {
    	const auto row_tensor_id = MAUtils::flat_to_tensor_index(row,wgt_B_proj);

    	for (Size col = 0 ; col < n_basis_proj_ ; ++col)
    	{
    		const auto col_tensor_id = MAUtils::flat_to_tensor_index(col,wgt_B_proj);

   			Real tmp = 1.0;
   			for (int k = 0 ; k < dim ; ++k)
   				tmp *= B_proj_1D[k](row_tensor_id[k],col_tensor_id[k]);

    		B_proj_(row,col) = tmp;
    	} //end loop col
    } //end loop row
    // evaluation of the Bernstein's mass matrix -- end
    //-----------------------------------------------------------------



    //-----------------------------------------------------------------
    // evaluation of the Bernstein's mass matrix inverse -- begin
    inv_B_proj_ = B_proj_.inverse();
    // evaluation of the Bernstein's mass matrix inverse -- end
    //-----------------------------------------------------------------
}


template<class PhysSpaceTest, class PhysSpaceTrial>
inline
void
EllipticOperatorsSumFactorizationIntegration<PhysSpaceTest,PhysSpaceTrial>::
eval_operator_u_v(
		const ElemTest &elem_test,
		const ElemTrial &elem_trial,
		const ValueVector<Real> &coeffs,
		DenseMatrix &operator_u_v) const
{
	//TODO: only the symmetric case is tested. In the non symmetric case, we need to check that
	// the physical space iterators have the same grid, map, reference space, index, etc.
	Assert(&elem_test == &elem_trial,ExcNotImplemented());



    //----------------------------------------------------
    // Assembly of the local mass matrix using sum-factorization -- begin
    const TimePoint start_assembly_mass_matrix = Clock::now();




    //--------------------------------------------------------------------------
    bool is_symmetric = this->test_same_space(elem_test,elem_trial);
    //--------------------------------------------------------------------------




    //--------------------------------------------------------------------------
    const auto start_initialization = Clock::now();

    //--------------------------------------------------------------------------
    // getting the number of basis along each coordinate direction for the projection space
    TensorSize<dim> n_basis_projection;
    for (int i = 0 ; i < dim ; ++i)
        n_basis_projection(i) = w_times_B_proj_1D_[i].get_num_functions();
    //--------------------------------------------------------------------------



    //--------------------------------------------------------------------------
    // getting the number of basis along each coordinate direction for the test and trial space

    const Index comp = 0; // only scalar spaces for the moment

    // test space -- begin
    TensorIndex<dim> degree_test = elem_test.get_physical_space()->get_reference_space()->get_degree()(comp);
    TensorSize<dim> n_basis_elem_test;
    for (int i = 0 ; i < dim ; ++i)
        n_basis_elem_test(i) = degree_test(i) + 1;

    const Size n_basis_test_flat = n_basis_elem_test.flat_size();
    Assert(n_basis_elem_test.flat_size()==elem_test.get_num_basis(),
           ExcDimensionMismatch(n_basis_elem_test.flat_size(),elem_test.get_num_basis()));

    const auto weight_basis_test = MultiArrayUtils<dim>::compute_weight(n_basis_elem_test);

    const TensorSize<dim> n_basis_test = n_basis_elem_test;
    // test space -- end


    // trial space -- begin
    TensorIndex<dim> degree_trial = elem_trial.get_physical_space()->get_reference_space()->get_degree()(comp);
    TensorSize<dim> n_basis_elem_trial;
    for (int i = 0 ; i < dim ; ++i)
        n_basis_elem_trial(i) = degree_trial(i) + 1;

    const Size n_basis_trial_flat = n_basis_elem_trial.flat_size();
    Assert(n_basis_elem_trial.flat_size()==elem_trial.get_num_basis(),
           ExcDimensionMismatch(n_basis_elem_trial.flat_size(),elem_trial.get_num_basis()));

    const auto weight_basis_trial = MultiArrayUtils<dim>::compute_weight(n_basis_elem_trial);

    const TensorSize<dim> n_basis_trial = n_basis_elem_trial;
    // trial space -- end

//    const Size n_basis = n_basis_elem.flat_size();
    //--------------------------------------------------------------------------


    using ValueType1D = Function<1>::ValueType;

    //--------------------------------------------------------------------------
    // getting the 1D values for the test space -- begin
    array< ValueTable<ValueType1D>,dim> phi_1D_test;
    {
    	const auto &ref_elem_accessor = elem_test.get_ref_space_accessor();

    	const auto &quad_points = ref_elem_accessor.get_quad_points();
    	const auto n_quad_points = quad_points.get_num_points_direction();

    	const auto &bspline_scalar_evaluators = ref_elem_accessor.get_scalar_evaluators()(comp);

    	for (int i = 0 ; i < dim ; ++i)
    		phi_1D_test[i].resize(n_basis_elem_test[i],n_quad_points[i]);

    	for (Index flat_fn_id = 0 ; flat_fn_id < n_basis_test_flat ; ++flat_fn_id)
    	{
    		const TensorIndex<dim> tensor_fn_id = MultiArrayUtils<dim>::flat_to_tensor_index(flat_fn_id,weight_basis_test);

    		const auto &bspline1D_values =
    				bspline_scalar_evaluators(tensor_fn_id)->get_derivative_components_view(0);

    		for (int i = 0 ; i < dim ; ++i)
    		{
    			auto phi_1D_ifn = phi_1D_test[i].get_function_view(tensor_fn_id[i]);

    			const auto &bsp_val = bspline1D_values[i];

    			for (int jpt = 0 ; jpt < n_quad_points[i] ; ++jpt)
    				phi_1D_ifn[jpt] = bsp_val(jpt);
    		}
    	}
    }
    // getting the 1D values for the test space -- end
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // getting the 1D values for the trial space -- begin
    array< ValueTable<ValueType1D>,dim> phi_1D_trial;
    {
    	const auto &ref_elem_accessor = elem_trial.get_ref_space_accessor();

    	const auto &quad_points = ref_elem_accessor.get_quad_points();
    	const auto n_quad_points = quad_points.get_num_points_direction();

    	const auto &bspline_scalar_evaluators = ref_elem_accessor.get_scalar_evaluators()(comp);

    	for (int i = 0 ; i < dim ; ++i)
    		phi_1D_trial[i].resize(n_basis_elem_trial[i],n_quad_points[i]);

    	for (Index flat_fn_id = 0 ; flat_fn_id < n_basis_trial_flat ; ++flat_fn_id)
    	{
    		const TensorIndex<dim> tensor_fn_id = MultiArrayUtils<dim>::flat_to_tensor_index(flat_fn_id,weight_basis_trial);

    		const auto &bspline1D_values =
    				bspline_scalar_evaluators(tensor_fn_id)->get_derivative_components_view(0);

    		for (int i = 0 ; i < dim ; ++i)
    		{
    			auto phi_1D_ifn =  phi_1D_trial[i].get_function_view(tensor_fn_id[i]);

    			const auto &bsp_val = bspline1D_values[i];

    			for (int jpt = 0 ; jpt < n_quad_points[i] ; ++jpt)
    				phi_1D_ifn[jpt] = bsp_val(jpt);
    		}
    	}
    }
    // getting the 1D values for the trial space -- end
    //--------------------------------------------------------------------------


    const auto end_initialization = Clock::now();
    const Duration elapsed_time_initialization = end_initialization - start_initialization;
    std::cout << "Elapsed_seconds initialization = " << elapsed_time_initialization.count() << std::endl;
    //--------------------------------------------------------------------------



    //----------------------------------------------------
    // Projection phase -- begin
    const TimePoint start_projection = Clock::now();


    // checks that the mapping used in the test space and in the trial space is the same
    Assert(elem_test.get_physical_space()->get_push_forward()->get_mapping() ==
    		elem_trial.get_physical_space()->get_push_forward()->get_mapping(),
    		ExcMessage("Test and trial spaces must have the same mapping (and the same grid)!"));


    //--------------------------------------------------------------------------
    // checks that the elements on the grid are the same
    using GridElem = CartesianGridElementAccessor<dim>;
    Assert(static_cast<const GridElem &>(elem_test.get_ref_space_accessor()) ==
    	   static_cast<const GridElem &>(elem_trial.get_ref_space_accessor()),
    	   ExcMessage("Different elements for test space and trial space."));
    //--------------------------------------------------------------------------


    // performs the projection of the function coeffs*det(DF) using as basis the Bernstein's polynomials
    const auto det_DF = elem_test.get_measures() ;

    Assert(det_DF.size() == coeffs.size(),
    		ExcDimensionMismatch(det_DF.size(), coeffs.size()));
    Size n_points = det_DF.size();

    ValueVector<Real> c_times_detDF(n_points);
    for (Index ipt = 0 ; ipt < n_points ; ++ipt)
    	c_times_detDF[ipt] = coeffs[ipt] * det_DF[ipt];

    // here we project the determinant Jacobian on a Bernstein polynomials space
    const auto K = projection_l2_bernstein_basis(c_times_detDF);

    const TimePoint end_projection = Clock::now();
    const Duration elapsed_time_projection = end_projection - start_projection;
    std::cout << "Elapsed seconds projection = " << elapsed_time_projection.count() << std::endl;
    // Projection phase -- end
    //----------------------------------------------------





    //----------------------------------------------------
    // precalculation of the I[i](lambda,mu1,mu2) terms (i.e. the moments)
    const auto start_compute_moments = Clock::now();

    const array<Real,dim> length_element_edge =
    		elem_test.get_ref_space_accessor().get_coordinate_lengths();

    const auto moments = evaluate_moments_op_u_v(phi_1D_test,phi_1D_trial,length_element_edge);

    const auto end_compute_moments = Clock::now();
    Duration elapsed_time_compute_moments = end_compute_moments - start_compute_moments;
    std::cout << "Elapsed seconds moments = " << elapsed_time_compute_moments.count() << std::endl;
    //----------------------------------------------------




    //----------------------------------------------------
    // Assembly of the local mass matrix using sum-factorization -- begin
    const auto start_sum_factorization = Clock::now();
    TensorSize<3> tensor_size_C_0;
    tensor_size_C_0[0] = n_basis_projection.flat_size(); // theta size
    tensor_size_C_0[1] = 1; // alpha size
    tensor_size_C_0[2] = 1; // beta size

    DynamicMultiArray<Real,3> C_0(tensor_size_C_0);
    const Size n_entries = tensor_size_C_0.flat_size();

    Assert(n_entries == K.flat_size(),ExcDimensionMismatch(n_entries,K.flat_size()));
    for (int entry_id = 0 ; entry_id < n_entries ; ++entry_id)
        C_0(entry_id) = K(entry_id);


    MassMatrixIntegrator<dim> integrate_mass_matrix;
    integrate_mass_matrix(is_symmetric,
    					  n_basis_projection,
                          n_basis_trial,
                          n_basis_test,
                          moments,
                          C_0,
                          operator_u_v);


    const auto end_sum_factorization = Clock::now();
    Duration elapsed_time_sum_factorization = end_sum_factorization - start_sum_factorization;
    std::cout << "Elapsed seconds sum-factorization = " << elapsed_time_sum_factorization.count() << std::endl;
    // Assembly of the local mass matrix using sum-factorization -- end
    //----------------------------------------------------

    Duration elapsed_time_assemble = elapsed_time_sum_factorization +
    		elapsed_time_compute_moments + elapsed_time_projection + elapsed_time_initialization ;
    std::cout << "Elapsed seconds assemblying = " << elapsed_time_assemble.count() << std::endl;


    const TimePoint end_assembly_mass_matrix = Clock::now();

    const_cast<Duration &>(elapsed_time_operator_u_v_) += end_assembly_mass_matrix - start_assembly_mass_matrix;
    std::cout << "Elapsed seconds operator u_v = " << elapsed_time_operator_u_v_.count() << std::endl;
    // Assembly of the local mass matrix using sum-factorization -- end
    //----------------------------------------------------
}


template<int dim>
void
PoissonProblemSumFactorization<dim>::
assemble()
{

    LogStream out;


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
//                            ValueFlags::gradient |
                            ValueFlags::measure |
                            ValueFlags::w_measure |
                            ValueFlags::point;

    elem->init_values(fill_flags, this->elem_quad);



    // number of points along each direction for the quadrature scheme.
    TensorSize<dim> n_quad_points = this->elem_quad.get_num_points_direction();

    using Clock = chrono::high_resolution_clock;
    using TimePoint = chrono::time_point<Clock>;


    DenseMatrix loc_mass_matrix_sf(n_basis,n_basis);


    using SpaceTest = Space;
    using SpaceTrial = Space;
    EllipticOperatorsSumFactorizationIntegration<SpaceTest,SpaceTrial>
    	elliptic_operators_sf(TensorSize<dim>(proj_deg_),quad_proj_);

    EllipticOperatorsStandardIntegration<SpaceTest,SpaceTrial> elliptic_operators_std;
    for (; elem != elem_end; ++elem)
    {
        loc_rhs.clear();

        elem->fill_values();


        //----------------------------------------------------
        // multiplicative coefficients of the mass matrix term.
        ValueVector<Real> c_mass(n_quad_points.flat_size());
        for ( auto & c : c_mass)
        	c = 1.0;
        //----------------------------------------------------



        //----------------------------------------------------
        // Assembly of the local mass matrix using sum-factorization -- begin
        const TimePoint start_assembly_mass_matrix = Clock::now();

        elliptic_operators_sf.eval_operator_u_v(*elem,*elem,c_mass,loc_mass_matrix_sf);

        const TimePoint end_assembly_mass_matrix = Clock::now();

        elapsed_time_assembly_mass_matrix_ += end_assembly_mass_matrix - start_assembly_mass_matrix;
        // Assembly of the local mass matrix using sum-factorization -- end
        //----------------------------------------------------




        //----------------------------------------------------
        // Assembly of the local mass matrix using the standard approach -- begin
        const TimePoint start_assembly_mass_matrix_old = Clock::now();

        elliptic_operators_std.eval_operator_u_v(*elem,*elem,c_mass,loc_mat);

        const TimePoint end_assembly_mass_matrix_old = Clock::now();
        elapsed_time_assembly_mass_matrix_old_ += end_assembly_mass_matrix_old - start_assembly_mass_matrix_old;
        // Assembly of the local mass matrix using the standard approach -- end
        //----------------------------------------------------


        //*/
        /*
        loc_dofs = elem->get_local_to_global();
        this->matrix->add_block(loc_dofs, loc_dofs, loc_mat);
        this->rhs->add_block(loc_dofs, loc_rhs);
        //*/

//        out<< "Local mass matrix sum-factorization=" << loc_mass_matrix_sf << endl << endl;
//        out<< "Local mass matrix original=" << loc_mat << endl << endl;
//        out<< "mass matrix difference=" << loc_mat - loc_mass_matrix_sf << endl << endl;
        const DenseMatrix m_diff = loc_mat - loc_mass_matrix_sf;
        out << "Frobenius norm of the difference="
            << m_diff.norm_frobenius() << endl;
//            << m_diff.norm_frobenius()/(m_diff.get_num_rows()*m_diff.get_num_cols()) << endl;
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
    do_test<1>();

    do_test<2>();

    do_test<3>();
//*/
    return  0;
}
