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
//#include <igatools/utils/value_table.h>
//#include <igatools/utils/multi_array_utils.h>
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
     */
    PoissonProblemSumFactorization(const TensorSize<dim> &n_knots, const int space_deg);

    /** Copy constructor. Not allowed to be used. */
    PoissonProblemSumFactorization(const PoissonProblemSumFactorization<dim> &in) = delete;

    /** Move constructor. Not allowed to be used. */
    PoissonProblemSumFactorization(PoissonProblemSumFactorization<dim> &&in) = delete;

    /** Destructor. */
    ~PoissonProblemSumFactorization() = default;
    ///@}



    Real get_elapsed_time_assembly_mass_matrix_sf() const
    {
        return elapsed_time_assembly_mass_matrix_sf_.count();
    }

    Real get_elapsed_time_assembly_mass_matrix_std() const
    {
        return elapsed_time_assembly_mass_matrix_std_.count();
    }

    Real get_elapsed_time_assembly_stiffness_matrix_sf() const
    {
        return elapsed_time_assembly_stiffness_matrix_sf_.count();
    }

    Real get_elapsed_time_assembly_stiffness_matrix_std() const
    {
        return elapsed_time_assembly_stiffness_matrix_std_.count();
    }

    void assemble();

private:
    using base_t = PoissonProblem< dim, PoissonProblemSumFactorization<dim> >;

    int space_deg_;

    chrono::duration<Real> elapsed_time_assembly_mass_matrix_sf_;

    chrono::duration<Real> elapsed_time_assembly_mass_matrix_std_;

    chrono::duration<Real> elapsed_time_assembly_stiffness_matrix_std_;

    chrono::duration<Real> elapsed_time_assembly_stiffness_matrix_sf_;

//    chrono::duration<Real> elapsed_time_compute_I_;

    DenseMatrix eval_matrix_proj() const;

};

template<int dim>
PoissonProblemSumFactorization<dim>::
PoissonProblemSumFactorization(const TensorSize<dim> &n_knots,const int space_deg)
    :
    base_t(n_knots,space_deg),
    space_deg_(space_deg)
{}



template<int dim>
void
PoissonProblemSumFactorization<dim>::
assemble()
{

    LogStream out;


    const int n_basis = this->space->get_num_basis_per_element();
    DenseVector loc_rhs(n_basis);
    vector<Index> loc_dofs(n_basis);

    ConstantFunction<dim> f({1.0});



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

    using Clock = chrono::high_resolution_clock;
    using TimePoint = chrono::time_point<Clock>;


    DenseMatrix loc_mass_matrix_std(n_basis, n_basis);
    DenseMatrix loc_mass_matrix_sf(n_basis,n_basis);

    DenseMatrix loc_stiffness_matrix_std(n_basis, n_basis);
    DenseMatrix loc_stiffness_matrix_sf(n_basis, n_basis);


    using SpaceTest = Space;
    using SpaceTrial = Space;
    EllipticOperatorsSFIntegration<SpaceTest,SpaceTrial>
    elliptic_operators_sf;

    EllipticOperatorsStdIntegration<SpaceTest,SpaceTrial> elliptic_operators_std;
    for (; elem != elem_end; ++elem)
    {
        loc_rhs.clear();

        elem->fill_values();


        //----------------------------------------------------
        // multiplicative coefficients of the mass matrix term.
        ValueVector<Real> c_mass(n_quad_points.flat_size());
        for (auto & c : c_mass)
            c = 1.0;
        //----------------------------------------------------



        //----------------------------------------------------
        // Assembly of the local mass matrix using sum-factorization -- begin
        const TimePoint start_assembly_mass_matrix_sf = Clock::now();

        elliptic_operators_sf.eval_operator_u_v(*elem,*elem,c_mass,loc_mass_matrix_sf);

        const TimePoint end_assembly_mass_matrix_sf = Clock::now();

        elapsed_time_assembly_mass_matrix_sf_ =
            end_assembly_mass_matrix_sf - start_assembly_mass_matrix_sf;
        // Assembly of the local mass matrix using sum-factorization -- end
        //----------------------------------------------------




        //----------------------------------------------------
        // Assembly of the local mass matrix using the standard approach -- begin
        const TimePoint start_assembly_mass_matrix_std = Clock::now();

        elliptic_operators_std.eval_operator_u_v(*elem,*elem,c_mass,loc_mass_matrix_std);

        const TimePoint end_assembly_mass_matrix_std = Clock::now();
        elapsed_time_assembly_mass_matrix_std_ =
            end_assembly_mass_matrix_std - start_assembly_mass_matrix_std;
        // Assembly of the local mass matrix using the standard approach -- end
        //----------------------------------------------------



        //----------------------------------------------------
        // multiplicative coefficients of the stiffness matrix term.
        vector<TMatrix<dim,dim>> c_stiffness(n_quad_points.flat_size());
        for (auto & c : c_stiffness)
            for (Index i = 0 ; i < dim ; ++i)
                c[i][i] = 1.0;
        //----------------------------------------------------


        //----------------------------------------------------
        // Assembly of the local stiffness matrix using the sum-factorization approach -- begin
        const TimePoint start_assembly_stiffness_matrix_sf = Clock::now();

        elliptic_operators_sf.eval_operator_gradu_gradv(
            *elem,*elem,c_stiffness,loc_stiffness_matrix_sf);

        const TimePoint end_assembly_stiffness_matrix_sf = Clock::now();
        elapsed_time_assembly_stiffness_matrix_sf_ =
            end_assembly_stiffness_matrix_sf - start_assembly_stiffness_matrix_sf;
        // Assembly of the local stiffness matrix using the sum-factorization approach -- end
        //----------------------------------------------------



        //----------------------------------------------------
        // Assembly of the local stiffness matrix using the standard approach -- begin
        const TimePoint start_assembly_stiffness_matrix_std = Clock::now();

        elliptic_operators_std.eval_operator_gradu_gradv(
            *elem,*elem,c_stiffness,loc_stiffness_matrix_std);

        const TimePoint end_assembly_stiffness_matrix_std = Clock::now();
        elapsed_time_assembly_stiffness_matrix_std_ =
            end_assembly_stiffness_matrix_std - start_assembly_stiffness_matrix_std;
        // Assembly of the local stiffness matrix using the standard approach -- end
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
        const DenseMatrix mass_matrix_diff = loc_mass_matrix_std - loc_mass_matrix_sf;
        out << "Mass-matrix: maximum norm of the difference="
            << mass_matrix_diff.norm_max() << endl;


//        out<< "Local stiffness matrix standard=" << loc_stiffness_matrix_std << endl << endl;
        const DenseMatrix stiffness_matrix_diff = loc_stiffness_matrix_std - loc_stiffness_matrix_sf;
        out << "Siffness-matrix: maximum norm of the difference="
            << stiffness_matrix_diff.norm_max() << endl;

    }

    this->matrix->fill_complete();
//*/

    out << "Dim=" << dim << "         space_deg=" << space_deg_ << endl;
    out << "Elapsed seconds assembly mass matrix sum-factorization = "
        << elapsed_time_assembly_mass_matrix_sf_.count() << endl;
    out << "Elapsed seconds assembly mass matrix standard quadrature = "
        << elapsed_time_assembly_mass_matrix_std_.count() << endl;
    out << endl;


    // AssertThrow(false,ExcNotImplemented());
}


template <int dim>
void
do_test()
{
    const int n_knots = 2;

    TableHandler elapsed_time_table;

    string time_mass_sum_fac = "Time mass-matrix sum_fac";
    string time_mass_orig = "Time mass-matrix orig";

    string time_stiff_sum_fac = "Time stiff-matrix sum_fac";
    string time_stiff_orig = "Time stiff-matrix orig";

    int degree_min = 3;
    int degree_max = 3;
    for (int degree = degree_min ; degree <= degree_max ; ++degree)
    {
        PoissonProblemSumFactorization<dim> poisson_sf(TensorSize<dim>(n_knots),degree);
        poisson_sf.run();
        //*/
        elapsed_time_table.add_value("Degree",degree);
        elapsed_time_table.add_value(time_mass_sum_fac,poisson_sf.get_elapsed_time_assembly_mass_matrix_sf());
        elapsed_time_table.add_value(time_mass_orig,poisson_sf.get_elapsed_time_assembly_mass_matrix_std());
        elapsed_time_table.add_value(time_stiff_sum_fac,poisson_sf.get_elapsed_time_assembly_stiffness_matrix_sf());
        elapsed_time_table.add_value(time_stiff_orig,poisson_sf.get_elapsed_time_assembly_stiffness_matrix_std());

    }

    elapsed_time_table.set_precision(time_mass_sum_fac,10);
    elapsed_time_table.set_scientific(time_mass_sum_fac,true);

    elapsed_time_table.set_precision(time_mass_orig,10);
    elapsed_time_table.set_scientific(time_mass_orig,true);

    elapsed_time_table.set_precision(time_stiff_sum_fac,10);
    elapsed_time_table.set_scientific(time_stiff_sum_fac,true);

    elapsed_time_table.set_precision(time_stiff_orig,10);
    elapsed_time_table.set_scientific(time_stiff_orig,true);

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
