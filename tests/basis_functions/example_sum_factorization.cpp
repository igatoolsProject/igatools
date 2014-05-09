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
#include <igatools/contrib/table_handler.h>

#include <igatools/operators/elliptic_operators_std_integration.h>
#include <igatools/operators/elliptic_operators_sf_integration.h>

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





// [Problem class]
template<int dim,class DerivedClass>
class PoissonProblem
{
public:
    PoissonProblem(const TensorSize<dim> &n_knots, const int deg);
    void run();


    Real get_elapsed_time_eval_basis() const;
    Real get_elapsed_time_eval_mass_matrix() const;
    Real get_elapsed_time_eval_stiffness_matrix() const;
    Real get_elapsed_time_eval_rhs() const;
    Real get_elapsed_time_assemble_stiffness_matrix() const;
    Real get_elapsed_time_solve_linear_system() const;
    Real get_elapsed_time_fill_complete() const;

private:
    void assemble();
    void solve();
    void output();
    // [Problem class]

    // [type aliases]
protected:
    using RefSpace  = BSplineSpace<dim>;
    using PushFw    = PushForward<Transformation::h_grad, dim>;
    using Space     = PhysicalSpace<RefSpace, PushFw>;
    using SpaceTest = Space;
    using SpaceTrial = Space;

    using Duration = chrono::duration<Real>;
    using Clock = chrono::high_resolution_clock;
    using TimePoint = chrono::time_point<Clock>;

//    using Value     = typename Function<dim>::Value;
    // [type aliases]


    const int deg_;

    shared_ptr<Mapping<dim>> map;
    shared_ptr<Space>        space;

    const Quadrature<dim>   elem_quad;
    const Quadrature<dim-1> face_quad;

    const boundary_id dir_id = 0;

#if defined(USE_TRILINOS)
    const static LinearAlgebraPackage linear_algebra_package = LinearAlgebraPackage::trilinos;
#elif defined(USE_PETSC)
    const static LinearAlgebraPackage linear_algebra_package = LinearAlgebraPackage::petsc;
#endif
    using MatrixType = Matrix<linear_algebra_package>;
    using VectorType = Vector<linear_algebra_package>;
    using LinearSolverType = LinearSolver<linear_algebra_package>;


    std::shared_ptr<MatrixType> matrix;
    std::shared_ptr<VectorType> rhs;
    std::shared_ptr<VectorType> solution;

    Duration elapsed_time_eval_basis_;

    Duration elapsed_time_eval_mass_matrix_;

    Duration elapsed_time_eval_stiffness_matrix_;

    Duration elapsed_time_eval_rhs_;

    Duration elapsed_time_assemble_stiffness_matrix_;

    Duration elapsed_time_fill_complete_;

    Duration elapsed_time_solve_linear_system_;

    std::string filename_;
};

template<int dim,class DerivedClass>
Real
PoissonProblem<dim,DerivedClass>::
get_elapsed_time_eval_basis() const
{
    return elapsed_time_eval_basis_.count();
}

template<int dim,class DerivedClass>
Real
PoissonProblem<dim,DerivedClass>::
get_elapsed_time_eval_mass_matrix() const
{
    return elapsed_time_eval_mass_matrix_.count();
}

template<int dim,class DerivedClass>
Real
PoissonProblem<dim,DerivedClass>::
get_elapsed_time_eval_stiffness_matrix() const
{
    return elapsed_time_eval_stiffness_matrix_.count();
}

template<int dim,class DerivedClass>
Real
PoissonProblem<dim,DerivedClass>::
get_elapsed_time_eval_rhs() const
{
    return elapsed_time_eval_rhs_.count();
}

template<int dim,class DerivedClass>
Real
PoissonProblem<dim,DerivedClass>::
get_elapsed_time_assemble_stiffness_matrix() const
{
    return elapsed_time_assemble_stiffness_matrix_.count();
}

template<int dim,class DerivedClass>
Real
PoissonProblem<dim,DerivedClass>::
get_elapsed_time_solve_linear_system() const
{
    return elapsed_time_solve_linear_system_.count();
}

template<int dim,class DerivedClass>
Real
PoissonProblem<dim,DerivedClass>::
get_elapsed_time_fill_complete() const
{
    return elapsed_time_fill_complete_.count();
}

template<int dim,class DerivedClass>
PoissonProblem<dim,DerivedClass>::
PoissonProblem(const TensorSize<dim> &n_knots, const int deg)
    :
    deg_(deg),
    elem_quad(QGauss<dim>(deg+1)),
    face_quad(QGauss<dim-1>(deg+1)),
    filename_("poisson_problem-" + to_string(dim) + "d")
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
    map       = BallMapping<dim>::create(grid);
//    map       = IdentityMapping<dim,0>::create(grid);
    space     = Space::create(ref_space, PushFw::create(map));

    const auto n_basis = space->get_num_basis();
    matrix   = MatrixType::create(get_sparsity_pattern(const_pointer_cast<const Space>(space)));
    rhs      = VectorType::create(n_basis);
    solution = VectorType::create(n_basis);
}





template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
assemble()
{

    LogStream out;


    const Size n_basis = this->space->get_num_basis_per_element();
    DenseVector loc_rhs(n_basis);
    vector<Index> loc_dofs(n_basis);

    const Size n_qp = this->elem_quad.get_num_points();
    ConstantFunction<dim> f({0.5});
    vector< typename Function<dim>::ValueType > f_values(n_qp);



    //-----------------------------------------------------------------
    /**
     * Initialization of the container for the values of the function
     * that must be projected.
     */
    using ValueType = typename Function<dim>::ValueType;
    vector<ValueType> f_values_proj(this->elem_quad.get_num_points());
    //-----------------------------------------------------------------




    auto elem = this->space->begin();
    const auto elem_end = this->space->end();
    ValueFlags fill_flags = ValueFlags::value |
                            ValueFlags::gradient |
                            ValueFlags::measure |
                            ValueFlags::w_measure |
                            ValueFlags::point;

    elem->init_values(fill_flags, this->elem_quad);



    // number of points along each direction for the quadrature scheme.
    TensorSize<dim> n_quad_points = this->elem_quad.get_num_points_direction();



    DenseMatrix loc_mass_matrix(n_basis, n_basis);
    DenseMatrix loc_stiffness_matrix(n_basis, n_basis);


    const auto &elliptic_operators = static_cast<const DerivedClass &>(*this).get_elliptic_operators();

//    EllipticOperatorsStdIntegration<SpaceTest,SpaceTrial> elliptic_operators_std;
    for (; elem != elem_end; ++elem)
    {
        loc_rhs.clear();

        //----------------------------------------------------
        const TimePoint start_eval_basis = Clock::now();
        elem->fill_values();
        const TimePoint end_eval_basis = Clock::now();
        this->elapsed_time_eval_basis_ += end_eval_basis - start_eval_basis;

        auto points  = elem->get_points();
        auto phi     = elem->get_basis_values();
        auto grd_phi = elem->get_basis_gradients();
        auto w_meas  = elem->get_w_measures();
        //----------------------------------------------------


        //----------------------------------------------------
        f.evaluate(points, f_values);
        //----------------------------------------------------


        //----------------------------------------------------
        // multiplicative coefficients of the mass matrix term.
        ValueVector<Real> c_mass(n_quad_points.flat_size());
        for (auto & c : c_mass)
            c = 1.0;
        //----------------------------------------------------


        //----------------------------------------------------
        // Assembly of the local mass matrix -- begin
        const TimePoint start_eval_mass_matrix = Clock::now();

        elliptic_operators.eval_operator_u_v(*elem,*elem,c_mass,loc_mass_matrix);

        const TimePoint end_eval_mass_matrix = Clock::now();

        this->elapsed_time_eval_mass_matrix_ +=
            end_eval_mass_matrix - start_eval_mass_matrix;
        // Assembly of the local mass matrix -- end
        //----------------------------------------------------


        //----------------------------------------------------
        // multiplicative coefficients of the stiffness matrix term.
        vector<TMatrix<dim,dim>> c_stiffness(n_quad_points.flat_size());
        for (auto & c : c_stiffness)
            for (Index i = 0 ; i < dim ; ++i)
                c[i][i] = 1.0;
        //----------------------------------------------------


        //----------------------------------------------------
        // Assembly of the local stiffness matrix -- begin
        const TimePoint start_eval_stiffness_matrix = Clock::now();

        elliptic_operators.eval_operator_gradu_gradv(
            *elem,*elem,c_stiffness,loc_stiffness_matrix);

        const TimePoint end_eval_stiffness_matrix = Clock::now();
        this->elapsed_time_eval_stiffness_matrix_ +=
            end_eval_stiffness_matrix - start_eval_stiffness_matrix;
        // Assembly of the local stiffness matrix -- end
        //----------------------------------------------------


        //----------------------------------------------------
        // Assemblying the right hand side -- begin
        const TimePoint start_eval_rhs = Clock::now();
        for (int i = 0; i < n_basis; ++i)
        {
            auto phi_i = phi.get_function_view(i);
            for (int qp = 0; qp < n_qp; ++qp)
                loc_rhs(i) += scalar_product(phi_i[qp], f_values[qp]) * w_meas[qp];
        }
        const TimePoint end_eval_rhs = Clock::now();
        this->elapsed_time_eval_rhs_ += end_eval_rhs - start_eval_rhs;
        // Assemblying the right hand side -- end
        //----------------------------------------------------


        loc_dofs = elem->get_local_to_global();


        const TimePoint start_assemblying_matrix = Clock::now();
        this->matrix->add_block(loc_dofs, loc_dofs, loc_stiffness_matrix);
        const TimePoint end_assemblying_matrix = Clock::now();
        this->elapsed_time_assemble_stiffness_matrix_ += end_assemblying_matrix - start_assemblying_matrix;


        this->rhs->add_block(loc_dofs, loc_rhs);
    }

    const TimePoint start_fill_complete = Clock::now();
    this->matrix->fill_complete();
    const TimePoint end_fill_complete = Clock::now();
    this->elapsed_time_fill_complete_ = end_fill_complete - start_fill_complete;


    ConstantFunction<dim> g({0.0});
    std::map<Index, Real> values;
    const int dir_id = 0 ;
    project_boundary_values<Space,linear_algebra_package>(g, this->space, this->face_quad, dir_id, values);
    apply_boundary_values(values, *this->matrix, *this->rhs, *this->solution);
//*/

    out << "Dim=" << dim << "         space_deg=" << this->deg_ << endl;
    out << "Elapsed seconds eval mass matrix = "
        << this->elapsed_time_eval_mass_matrix_.count() << endl;
    out << "Elapsed seconds eval stiffness matrix = "
        << this->elapsed_time_eval_stiffness_matrix_.count() << endl;
    out << endl;


    // AssertThrow(false,ExcNotImplemented());
}


template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
solve()
{
    const Real tol = 1.0e-10;
    const Size max_iters = 10000000;
    LinearSolverType solver(LinearSolverType::SolverType::CG,tol,max_iters);
//    LinearSolverType solver(LinearSolverType::SolverType::LU,tol,max_iters);

    const TimePoint start_solve_linear_system = Clock::now();

    solver.solve(*matrix, *rhs, *solution);

    const TimePoint end_solve_linear_system = Clock::now();
    this->elapsed_time_solve_linear_system_ =
        end_solve_linear_system - start_solve_linear_system;
}



template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
output()
{
    const int n_plot_points = deg_+1;
    Writer<dim> writer(map, n_plot_points);

    writer.add_field(space, *solution, "solution");
//    string filename = "poisson_problem-" + to_string(dim) + "d" ;
    writer.save(filename_,"appended");
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
    using base_t = PoissonProblem< dim, PoissonProblemStandardIntegration<dim> >;
    using typename base_t::Space;
    using typename base_t::SpaceTest;
    using typename base_t::SpaceTrial;

    using EllipticOperatorsType = EllipticOperatorsStdIntegration<SpaceTest,SpaceTrial>;

    PoissonProblemStandardIntegration(const TensorSize<dim> &n_knots, const int space_deg);

    const EllipticOperatorsType &get_elliptic_operators() const;

private:

    EllipticOperatorsType elliptic_operators_std_;

};

template<int dim>
PoissonProblemStandardIntegration<dim>::
PoissonProblemStandardIntegration(const TensorSize<dim> &n_knots,const int space_deg)
    :
    base_t(n_knots,space_deg)
{
    this->filename_ += "-std";
}



template<int dim>
auto
PoissonProblemStandardIntegration<dim>::
get_elliptic_operators() const -> const EllipticOperatorsType &
{
    return elliptic_operators_std_;
}



template<int dim>
class PoissonProblemSumFactorization :
    public PoissonProblem< dim, PoissonProblemSumFactorization<dim> >
{
public:
    using base_t = PoissonProblem< dim, PoissonProblemSumFactorization<dim> >;
    using typename base_t::Space;
    using typename base_t::SpaceTest;
    using typename base_t::SpaceTrial;

    using EllipticOperatorsType = EllipticOperatorsSFIntegration<SpaceTest,SpaceTrial>;

    /** @name Constructors & destructor. */
    ///@{
    PoissonProblemSumFactorization() = delete ;

    /**
     * Constructor.
     * @param[in] n_knots Number of knots along each direction.
     * @param[in] space_deg Polynomial degree for the solution space.
     */
    PoissonProblemSumFactorization(const TensorSize<dim> &n_knots, const int space_deg);

    /** Copy constructor. Not allowed to be used. */
    PoissonProblemSumFactorization(const PoissonProblemSumFactorization<dim> &in) = delete;

    /** Move constructor. Not allowed to be used. */
    PoissonProblemSumFactorization(PoissonProblemSumFactorization<dim> &&in) = delete;

    /** Destructor. */
    ~PoissonProblemSumFactorization() = default;
    ///@}


    const EllipticOperatorsType &get_elliptic_operators() const;

private:
    EllipticOperatorsSFIntegration<SpaceTest,SpaceTrial>
    elliptic_operators_sf_;
};

template<int dim>
PoissonProblemSumFactorization<dim>::
PoissonProblemSumFactorization(const TensorSize<dim> &n_knots,const int space_deg)
    :
    base_t(n_knots,space_deg)
{
    this->filename_ += "-sf";
}



template<int dim>
auto
PoissonProblemSumFactorization<dim>::
get_elliptic_operators() const -> const EllipticOperatorsType &
{
    return elliptic_operators_sf_;
}


template <int dim>
void
do_test()
{
    using std::cout;
    using std::endl;

    const int n_knots = 11;

    TableHandler elapsed_time_table;

    string time_eval_basis = "Eval basis";

    string time_eval_rhs = "Eval rhs";

    string time_mass_sum_fac = "Eval mass sum_fac";
    string time_mass_orig = "Eval mass orig";

    string time_stiff_sum_fac = "Eval stiffness sum_fac";
    string time_stiff_orig = "Eval stiffness orig";

    string time_assemble = "Time assemble";

    string time_fill_complete = "Time fill_complete";

    string time_solve_lin_system = "Time solve lin.system";

    int degree_min = 1;
    int degree_max = 1;
    for (int degree = degree_min ; degree <= degree_max ; ++degree)
    {
        cout << "-----------------------------------" << endl;
        cout << "Sum-Factorization -- begin" << endl;
        PoissonProblemSumFactorization<dim> poisson_sf(TensorSize<dim>(n_knots),degree);
        poisson_sf.run();
        cout << "Sum-Factorization -- end" << endl;
        cout << "-----------------------------------" << endl;

        cout << endl;

        cout << "-----------------------------------" << endl;
        cout << "Standard Quadrature -- begin" << endl;
        PoissonProblemStandardIntegration<dim> poisson_std(TensorSize<dim>(n_knots),degree);
        poisson_std.run();
        cout << "Standard Quadrature -- end" << endl;
        cout << "-----------------------------------" << endl;

        cout << endl;

        //*/
        elapsed_time_table.add_value("Degree",degree);
        elapsed_time_table.add_value(time_eval_basis,poisson_sf.get_elapsed_time_eval_basis());

        elapsed_time_table.add_value(time_eval_rhs,poisson_sf.get_elapsed_time_eval_rhs());

        elapsed_time_table.add_value(time_mass_sum_fac,poisson_sf.get_elapsed_time_eval_mass_matrix());
        elapsed_time_table.add_value(time_mass_orig,poisson_std.get_elapsed_time_eval_mass_matrix());
        elapsed_time_table.add_value(time_stiff_sum_fac,poisson_sf.get_elapsed_time_eval_stiffness_matrix());
        elapsed_time_table.add_value(time_stiff_orig,poisson_std.get_elapsed_time_eval_stiffness_matrix());

        elapsed_time_table.add_value(time_assemble,poisson_sf.get_elapsed_time_assemble_stiffness_matrix());
        elapsed_time_table.add_value(time_fill_complete,poisson_sf.get_elapsed_time_fill_complete());
        elapsed_time_table.add_value(time_solve_lin_system,poisson_sf.get_elapsed_time_solve_linear_system());

    }
    elapsed_time_table.set_precision(time_eval_basis,10);
    elapsed_time_table.set_scientific(time_eval_basis,true);

    elapsed_time_table.set_precision(time_eval_rhs,10);
    elapsed_time_table.set_scientific(time_eval_rhs,true);

    elapsed_time_table.set_precision(time_mass_sum_fac,10);
    elapsed_time_table.set_scientific(time_mass_sum_fac,true);

    elapsed_time_table.set_precision(time_mass_orig,10);
    elapsed_time_table.set_scientific(time_mass_orig,true);

    elapsed_time_table.set_precision(time_stiff_sum_fac,10);
    elapsed_time_table.set_scientific(time_stiff_sum_fac,true);

    elapsed_time_table.set_precision(time_stiff_orig,10);
    elapsed_time_table.set_scientific(time_stiff_orig,true);

    elapsed_time_table.set_precision(time_assemble,10);
    elapsed_time_table.set_scientific(time_assemble,true);

    elapsed_time_table.set_precision(time_fill_complete,10);
    elapsed_time_table.set_scientific(time_fill_complete,true);

    elapsed_time_table.set_precision(time_solve_lin_system,10);
    elapsed_time_table.set_scientific(time_solve_lin_system,true);

    ofstream elapsed_time_file("sum_factorization_time_"+to_string(dim)+"D.txt");
    elapsed_time_table.write_text(elapsed_time_file);

}


int main(int argc,char **args)

{
#if defined(USE_TRILINOS)
#elif defined(USE_PETSC)
    PetscInitialize(&argc,&args,(char *)0,"Sum factorization example");
#endif

    do_test<1>();

    do_test<2>();

    do_test<3>();
//*/
#if defined(USE_TRILINOS)
#elif defined(USE_PETSC)
    auto ierr = PetscFinalize();
#endif

    return  0;
}
