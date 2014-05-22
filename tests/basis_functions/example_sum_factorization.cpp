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
    PoissonProblem(const int n_knots, const int deg);


    void run();


    Real get_elapsed_time_eval_basis() const;
    Real get_elapsed_time_eval_mass_matrix() const;
    Real get_elapsed_time_eval_stiffness_matrix() const;
    Real get_elapsed_time_eval_rhs() const;
    Real get_elapsed_time_assemble_stiffness_matrix() const;
    Real get_elapsed_time_solve_linear_system() const;
    Real get_elapsed_time_fill_complete() const;
    Real get_elapsed_time_total() const;

    int get_num_dofs() const;
    int get_num_iters() const;
    Real get_achieved_tol() const;

    static std::string get_filename();

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

    TimePoint start_poisson_;
    TimePoint   end_poisson_;

    const int deg_;

    shared_ptr<Mapping<dim>> map;
    shared_ptr<Space>        space;

    const Quadrature<dim>   elem_quad;
    const Quadrature<dim-1> face_quad;

    const boundary_id dir_id = 0;

#if defined(USE_PETSC)
    const static LinearAlgebraPackage linear_algebra_package = LinearAlgebraPackage::petsc;
#elif defined(USE_TRILINOS)
    const static LinearAlgebraPackage linear_algebra_package = LinearAlgebraPackage::trilinos;
#endif
    using MatrixType = Matrix<linear_algebra_package>;
    using VectorType = Vector<linear_algebra_package>;
    using LinearSolverType = LinearSolver<linear_algebra_package>;


    std::shared_ptr<MatrixType> matrix;
    std::shared_ptr<VectorType> rhs;
    std::shared_ptr<VectorType> solution;


    Duration elapsed_time_total_;

    Duration elapsed_time_eval_basis_;

    Duration elapsed_time_eval_mass_matrix_;

    Duration elapsed_time_eval_stiffness_matrix_;

    Duration elapsed_time_eval_rhs_;

    Duration elapsed_time_assemble_stiffness_matrix_;

    Duration elapsed_time_fill_complete_;

    Duration elapsed_time_solve_linear_system_;

    int num_dofs_;

    int num_iters_;

    Real achieved_tol_;
};



template<int dim,class DerivedClass>
Real
PoissonProblem<dim,DerivedClass>::
get_achieved_tol() const
{
    return achieved_tol_;
}

template<int dim,class DerivedClass>
Real
PoissonProblem<dim,DerivedClass>::
get_elapsed_time_total() const
{
    return elapsed_time_total_.count();
}

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
int
PoissonProblem<dim,DerivedClass>::
get_num_dofs() const
{
    return num_dofs_;
}

template<int dim,class DerivedClass>
int
PoissonProblem<dim,DerivedClass>::
get_num_iters() const
{
    return num_iters_;
}

template<int dim,class DerivedClass>
std::string
PoissonProblem<dim,DerivedClass>::
get_filename()
{
    return "poisson_problem-" + to_string(dim) + "D";
}


template<int dim,class DerivedClass>
PoissonProblem<dim,DerivedClass>::
PoissonProblem(const int n_knots, const int deg)
    :
    start_poisson_(Clock::now()),
    deg_(deg),
    elem_quad(QGauss<dim>(deg+1)),
    face_quad(QGauss<dim-1>(deg+1))
{
    LogStream out;
    out << "PoissonProblem constructor -- begin" << endl;

    BBox<dim> box;
    for (int i=0 ; i < dim ; ++i)
    {
        box[i][0] = 0.0;
        box[i][1] = 1.0;
    }
    box[0][1] = 0.5;

    auto grid = CartesianGrid<dim>::create(box,TensorSize<dim>(n_knots));
    auto ref_space = RefSpace::create(grid, deg);
    map       = BallMapping<dim>::create(grid);
//    map       = IdentityMapping<dim,0>::create(grid);
    space     = Space::create(ref_space, PushFw::create(map));

    num_dofs_ = space->get_num_basis();
    matrix   = MatrixType::create(get_sparsity_pattern(const_pointer_cast<const Space>(space)));
    rhs      = VectorType::create(num_dofs_);
    solution = VectorType::create(num_dofs_);
    out << "PoissonProblem constructor -- end" << endl;

}




template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
assemble()
{

    LogStream out;

    Duration elapsed_time_assemble;
    TimePoint start_assemble = Clock::now();


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


    TimePoint start_boundary_conditions = Clock::now();
    ConstantFunction<dim> g({0.0});
    std::map<Index, Real> values;
    const int dir_id = 0 ;
    project_boundary_values<Space,linear_algebra_package>(g, this->space, this->face_quad, dir_id, values);


    TimePoint start_apply_bc = Clock::now();
    apply_boundary_values(values, *this->matrix, *this->rhs, *this->solution);
    TimePoint end_apply_bc = Clock::now();
    Duration elapsed_time_apply_bc = end_apply_bc - start_apply_bc;



    TimePoint end_boundary_conditions = Clock::now();
    Duration elapsed_time_boundary_conditions = end_boundary_conditions - start_boundary_conditions;

//*/

    TimePoint end_assemble = Clock::now();

    elapsed_time_assemble += end_assemble - start_assemble;

    out << "Dim=" << dim << "         space_deg=" << this->deg_ << endl;
    out << "Elapsed seconds eval mass matrix = "
        << this->elapsed_time_eval_mass_matrix_.count() << endl;
    out << "Elapsed seconds eval stiffness matrix = "
        << this->elapsed_time_eval_stiffness_matrix_.count() << endl;
    out << "Elapsed seconds apply bc = "
        << elapsed_time_apply_bc.count() << endl;
    out << "Elapsed seconds boundary conditions = "
        << elapsed_time_boundary_conditions.count() << endl;
    out << "Elapsed seconds assemble() function = "
        << elapsed_time_assemble.count() << endl;
    out << endl;


    // AssertThrow(false,ExcNotImplemented());
}


template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
solve()
{
    const Real tol = 1.0e-6;
    const Size max_iters = 10000000;
    LinearSolverType solver(
        LinearSolverType::SolverType::CG,
        LinearSolverType::PreconditionerType::ILU,
        tol,
        max_iters);
//    LinearSolverType solver(LinearSolverType::SolverType::LU,tol,max_iters);

    const TimePoint start_solve_linear_system = Clock::now();

    solver.solve(*matrix, *rhs, *solution);

    num_iters_ = solver.get_num_iterations();
    achieved_tol_ = solver.get_achieved_tolerance();

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
    writer.save(PoissonProblem<dim,DerivedClass>::get_filename(),"appended");
}



template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
run()
{
    static_cast<DerivedClass &>(*this).assemble();
    solve();
    end_poisson_ = Clock::now();

//    output();


    elapsed_time_total_ = end_poisson_ - start_poisson_;
}



template<int dim>
class PoissonProblemStandardIntegration :
    public PoissonProblem< dim, PoissonProblemStandardIntegration<dim> >
{
public:
    using base_t = PoissonProblem< dim, PoissonProblemStandardIntegration<dim> >;
    using base_t::base_t;
    using typename base_t::Space;
    using typename base_t::SpaceTest;
    using typename base_t::SpaceTrial;

    using EllipticOperatorsType = EllipticOperatorsStdIntegration<SpaceTest,SpaceTrial>;

    PoissonProblemStandardIntegration(const TensorSize<dim> &n_knots, const int space_deg);

    const EllipticOperatorsType &get_elliptic_operators() const;

    static std::string get_filename();

private:

    EllipticOperatorsType elliptic_operators_std_;

};
/*
template<int dim>
PoissonProblemStandardIntegration<dim>::
PoissonProblemStandardIntegration(const TensorSize<dim> &n_knots,const int space_deg)
    :
    base_t(n_knots,space_deg)
{}
//*/
template<int dim>
std::string
PoissonProblemStandardIntegration<dim>::
get_filename()
{
    return base_t::get_filename() + "-std";
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
    using base_t::base_t;
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

    static std::string get_filename();

private:
    EllipticOperatorsSFIntegration<SpaceTest,SpaceTrial>
    elliptic_operators_sf_;
};
/*
template<int dim>
PoissonProblemSumFactorization<dim>::
PoissonProblemSumFactorization(const TensorSize<dim> &n_knots,const int space_deg)
    :
    base_t(n_knots,space_deg)
{}
//*/
template<int dim>
std::string
PoissonProblemSumFactorization<dim>::
get_filename()
{
    return base_t::get_filename() + "-sf";
}



template<int dim>
auto
PoissonProblemSumFactorization<dim>::
get_elliptic_operators() const -> const EllipticOperatorsType &
{
    return elliptic_operators_sf_;
}


template <class PoissonProblemSolver >
void
do_test(const int degree_min, const int degree_max,const int n_elems_per_direction)
{
    using std::cout;
    using std::endl;

    const int n_knots = n_elems_per_direction+1;

    TableHandler time_table;

    string time_eval_basis = "Basis";

    string time_eval_rhs = "RHS";

    string time_mass = "Mass";

    string time_stiff = "Stiffness";

    string time_assemble = "Loc-to-glob assemble";

    string time_fill_complete = "Fill-complete";

    string time_solve_lin_system = "Solve lin.system";

    string time_total = "Total";

    string achieved_tol = "Tol";

    for (int degree = degree_min ; degree <= degree_max ; ++degree)
    {
        PoissonProblemSolver poisson(n_knots,degree);
        poisson.run();


        //*/
        time_table.add_value("Degree",degree);
        time_table.add_value(time_eval_basis,poisson.get_elapsed_time_eval_basis());

        time_table.add_value(time_eval_rhs,poisson.get_elapsed_time_eval_rhs());

        time_table.add_value(time_mass,poisson.get_elapsed_time_eval_mass_matrix());
        time_table.add_value(time_stiff,poisson.get_elapsed_time_eval_stiffness_matrix());

        time_table.add_value(time_assemble,poisson.get_elapsed_time_assemble_stiffness_matrix());
        time_table.add_value(time_fill_complete,poisson.get_elapsed_time_fill_complete());
        time_table.add_value(time_solve_lin_system,poisson.get_elapsed_time_solve_linear_system());

        time_table.add_value(time_total,poisson.get_elapsed_time_total());

        time_table.add_value("Num dofs",poisson.get_num_dofs());
        time_table.add_value("Num iters",poisson.get_num_iters());
        time_table.add_value(achieved_tol,poisson.get_achieved_tol());

    }
    time_table.set_precision(time_eval_basis,10);
    time_table.set_scientific(time_eval_basis,true);

    time_table.set_precision(time_eval_rhs,10);
    time_table.set_scientific(time_eval_rhs,true);

    time_table.set_precision(time_mass,10);
    time_table.set_scientific(time_mass,true);

    time_table.set_precision(time_stiff,10);
    time_table.set_scientific(time_stiff,true);

    time_table.set_precision(time_assemble,10);
    time_table.set_scientific(time_assemble,true);

    time_table.set_precision(time_fill_complete,10);
    time_table.set_scientific(time_fill_complete,true);

    time_table.set_precision(time_solve_lin_system,10);
    time_table.set_scientific(time_solve_lin_system,true);

    time_table.set_precision(time_total,10);
    time_table.set_scientific(time_total,true);

    time_table.set_precision(achieved_tol,10);
    time_table.set_scientific(achieved_tol,true);




    ofstream elapsed_time_file(
        "time_" + PoissonProblemSolver::get_filename() +
        "_" + std::to_string(degree_min) +
        "_" + std::to_string(degree_max) +
        "_" + std::to_string(n_elems_per_direction) +
        ".txt");
    time_table.write_text(elapsed_time_file);

}


int main(int argc,char **args)

{
#if defined(USE_PETSC)
    PetscInitialize(&argc,&args,(char *)0,"Sum factorization example");
#endif


    int degree_min = 1;
    int degree_max = 1;
    int n_elems_per_direction = 50;


    cout << "-----------------------------------" << endl;
    cout << "Sum-Factorization -- begin" << endl;
    do_test< PoissonProblemSumFactorization<1> >(degree_min,degree_max,n_elems_per_direction);

    do_test< PoissonProblemSumFactorization<2> >(degree_min,degree_max,n_elems_per_direction);

    do_test< PoissonProblemSumFactorization<3> >(degree_min,degree_max,n_elems_per_direction);
    cout << "Sum-Factorization -- end" << endl;
    cout << "-----------------------------------" << endl;

    cout << endl;



    cout << "-----------------------------------" << endl;
    cout << "Standard Quadrature -- begin" << endl;
    do_test< PoissonProblemStandardIntegration<1> >(degree_min,degree_max,n_elems_per_direction);

    do_test< PoissonProblemStandardIntegration<2> >(degree_min,degree_max,n_elems_per_direction);

    do_test< PoissonProblemStandardIntegration<3> >(degree_min,degree_max,n_elems_per_direction);
    cout << "Standard Quadrature -- end" << endl;
    cout << "-----------------------------------" << endl;

    cout << endl;

    //*/
#if defined(USE_PETSC)
    auto ierr = PetscFinalize();
#endif

    return  0;
}
