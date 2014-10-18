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
#include "../tests.h"

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
#include <igatools/linear_algebra/dof_tools.h>

#include <igatools/operators/elliptic_operators_std_integration.h>
#include <igatools/operators/elliptic_operators_sf_integration.h>

#include <boost/math/special_functions/binomial.hpp>

#include <numeric>

#include <chrono>

// [old includes]

// [unqualified names]
using namespace iga;

using functions::ConstantFunction;
using space_tools::project_boundary_values;
using dof_tools::apply_boundary_values;
using numbers::PI;

using std::cout;
// [unqualified names]





// [Problem class]
template<int dim,class DerivedClass>
class PoissonProblem
{
public:
    PoissonProblem(const int n_knots, const int deg);


    ValueFlags get_fill_flags() const;

    const DerivedClass &as_derived_class() const ;


    Real get_elapsed_time_eval_basis() const;
    Real get_elapsed_time_eval_mass_matrix() const;
    Real get_elapsed_time_eval_stiffness_matrix() const;
    Real get_elapsed_time_eval_rhs() const;
    Real get_elapsed_time_assemble_stiffness_matrix() const;
    Real get_elapsed_time_fill_complete() const;

    int get_num_dofs() const;

    void assemble();
    // [Problem class]

    // [type aliases]
protected:
    using RefSpace  = BSplineSpace<dim>;
    using PushFw    = PushForward<Transformation::h_grad, dim>;
    using Space     = PhysicalSpace<RefSpace, PushFw>;
    using SpaceTest = Space;
    using SpaceTrial = Space;

    using Duration = std::chrono::duration<Real>;
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    // [type aliases]

    TimePoint start_poisson_;
    TimePoint   end_poisson_;

    const int deg_;

    shared_ptr<Mapping<dim>> map;
    shared_ptr<Space>        space;

    const Quadrature<dim>   elem_quad;
    const Quadrature<dim-1> face_quad;

    const boundary_id dir_id = 0;

#if defined(USE_TRILINOS)
    static const auto la_pack = LAPack::trilinos;
#elif defined(USE_PETSC)
    static const auto la_pack = LAPack::petsc;
#endif

    using MatrixType = Matrix<la_pack>;
    using VectorType = Vector<la_pack>;


    std::shared_ptr<MatrixType> matrix;
    std::shared_ptr<VectorType> rhs;
    std::shared_ptr<VectorType> solution;


    Duration elapsed_time_eval_basis_;

    Duration elapsed_time_eval_mass_matrix_;

    Duration elapsed_time_eval_stiffness_matrix_;

    Duration elapsed_time_eval_rhs_;

    Duration elapsed_time_assemble_stiffness_matrix_;

    Duration elapsed_time_fill_complete_;

    int num_dofs_;
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
PoissonProblem<dim,DerivedClass>::
PoissonProblem(const int n_knots, const int deg)
    :
    start_poisson_(Clock::now()),
    deg_(deg),
    elem_quad(QGauss<dim>(deg+1)),
    face_quad(QGauss<dim-1>(deg+1)),
    elapsed_time_eval_basis_(0),
    elapsed_time_eval_mass_matrix_(0),
    elapsed_time_eval_stiffness_matrix_(0),
    elapsed_time_eval_rhs_(0),
    elapsed_time_assemble_stiffness_matrix_(0),
    elapsed_time_fill_complete_(0)
{


    LogStream out_screen;
    out_screen << "PoissonProblem constructor -- begin" << endl;

    BBox<dim> box;
    for (int i=0 ; i < dim ; ++i)
    {
        box[i][0] = 0.0;
        box[i][1] = 1.0;
    }
    box[0][1] = 0.5;

    auto grid = CartesianGrid<dim>::create(box,TensorSize<dim>(n_knots));
    auto ref_space = RefSpace::create(deg,grid);
    map       = BallMapping<dim>::create(grid);
//    map       = IdentityMapping<dim,0>::create(grid);
    space     = Space::create(ref_space, PushFw::create(map));

    num_dofs_ = space->get_num_basis();
    matrix   = MatrixType::create(*space->get_space_manager());
    rhs      = VectorType::create(num_dofs_);
    solution = VectorType::create(num_dofs_);
    out_screen << "PoissonProblem constructor -- end" << endl;

}


template<int dim,class DerivedClass>
const DerivedClass &
PoissonProblem<dim,DerivedClass>::
as_derived_class() const
{
    return static_cast<const DerivedClass &>(*this);
}

template<int dim,class DerivedClass>
ValueFlags
PoissonProblem<dim,DerivedClass>::
get_fill_flags() const
{
    return this->as_derived_class().get_fill_flags();
}

template<int dim,class DerivedClass>
void
PoissonProblem<dim,DerivedClass>::
assemble()
{
    LogStream out_screen;

    Duration elapsed_time_assemble;
    TimePoint start_assemble = Clock::now();



    const Size n_qp = this->elem_quad.get_num_points();
    ConstantFunction<dim> f( {0.5});
    using Value = typename Function<dim>::Value;
    ValueVector<Value> f_values(n_qp);



    //-----------------------------------------------------------------
    /**
     * Initialization of the container for the values of the function
     * that must be projected.
     */
    ValueVector<Value> f_values_proj(this->elem_quad.get_num_points());
    //-----------------------------------------------------------------




    auto elem = this->space->begin();
    const auto elem_end = this->space->end();

    ValueFlags fill_flags = this->get_fill_flags();

    elem->init_cache(fill_flags, this->elem_quad);



    // number of points along each direction for the quadrature scheme.
    TensorSize<dim> n_quad_points = this->elem_quad.get_num_points_direction();


    const auto &elliptic_operators = static_cast<const DerivedClass &>(*this).get_elliptic_operators();

    for (; elem != elem_end; ++elem)
    {
        //----------------------------------------------------
        const Size n_basis = elem->get_num_basis();
        DenseVector loc_rhs(n_basis);
        loc_rhs = 0.0;

        DenseMatrix loc_mass_matrix(n_basis, n_basis);
        loc_mass_matrix = 0.0;

        DenseMatrix loc_stiffness_matrix(n_basis, n_basis);
        loc_stiffness_matrix = 0.0;
        //----------------------------------------------------


        //----------------------------------------------------
        const TimePoint start_eval_basis = Clock::now();
        elem->fill_cache();
        const TimePoint end_eval_basis = Clock::now();
        this->elapsed_time_eval_basis_ += end_eval_basis - start_eval_basis;

        auto points  = elem->get_points();
        auto phi     = elem->get_basis_values();
//        auto grd_phi = elem->get_basis_gradients();
        auto w_meas  = elem->get_w_measures();
        //----------------------------------------------------


        //----------------------------------------------------
        f.evaluate(points, f_values);
        //----------------------------------------------------


        //----------------------------------------------------
        // multiplicative coefficients of the mass matrix term.
        ValueVector<Real> c_mass(n_quad_points.flat_size());
        for (auto &c : c_mass)
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
        iga::vector<TMatrix<dim,dim>> c_stiffness(n_quad_points.flat_size());
        for (auto &c : c_stiffness)
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
        /*
        for (int i = 0; i < n_basis; ++i)
        {
            auto phi_i = phi.get_function_view(i);
            for (int qp = 0; qp < n_qp; ++qp)
                loc_rhs(i) += scalar_product(phi_i[qp], f_values[qp]) * w_meas[qp];
        }
        //*/
        elliptic_operators.eval_operator_rhs_v(*elem,f_values,loc_rhs);

        const TimePoint end_eval_rhs = Clock::now();
        this->elapsed_time_eval_rhs_ += end_eval_rhs - start_eval_rhs;
        // Assemblying the right hand side -- end
        //----------------------------------------------------


        auto loc_dofs = elem->get_local_to_global();


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
    ConstantFunction<dim> g( {0.0});
    std::map<Index, Real> values;
    const int dir_id = 0 ;
    project_boundary_values<Space,la_pack>(g, this->space, this->face_quad, dir_id, values);


    TimePoint start_apply_bc = Clock::now();
    apply_boundary_values(values, *this->matrix, *this->rhs, *this->solution);
    TimePoint end_apply_bc = Clock::now();
    Duration elapsed_time_apply_bc = end_apply_bc - start_apply_bc;



    TimePoint end_boundary_conditions = Clock::now();
    Duration elapsed_time_boundary_conditions = end_boundary_conditions - start_boundary_conditions;

//*/

    TimePoint end_assemble = Clock::now();

    elapsed_time_assemble += end_assemble - start_assemble;

    out << "==========================================================" << endl;
    out << "Dim=" << dim << "         space_deg=" << this->deg_ << endl;

    out_screen << "Dim=" << dim << "         space_deg=" << this->deg_ << endl;
    out_screen << "Elapsed seconds eval mass matrix = "
               << this->elapsed_time_eval_mass_matrix_.count() << endl;
    out_screen << "Elapsed seconds eval stiffness matrix = "
               << this->elapsed_time_eval_stiffness_matrix_.count() << endl;
    out_screen << "Elapsed seconds apply bc = "
               << elapsed_time_apply_bc.count() << endl;
    out_screen << "Elapsed seconds boundary conditions = "
               << elapsed_time_boundary_conditions.count() << endl;
    out_screen << "Elapsed seconds assemble() function = "
               << elapsed_time_assemble.count() << endl;
    out_screen << endl;


    this->matrix->print(out);

    out << "==========================================================" << endl;
    out << endl;
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

    ValueFlags get_fill_flags() const;

private:

    EllipticOperatorsType elliptic_operators_std_;

};


template<int dim>
auto
PoissonProblemStandardIntegration<dim>::
get_elliptic_operators() const -> const EllipticOperatorsType &
{
    return elliptic_operators_std_;
}

template<int dim>
ValueFlags
PoissonProblemStandardIntegration<dim>::
get_fill_flags() const
{
    ValueFlags fill_flags = ValueFlags::value |
                            ValueFlags::gradient |
                            ValueFlags::measure |
                            ValueFlags::w_measure |
                            ValueFlags::point;

    return fill_flags;
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

    ValueFlags get_fill_flags() const;

private:
    EllipticOperatorsSFIntegration<SpaceTest,SpaceTrial>
    elliptic_operators_sf_;
};

template<int dim>
auto
PoissonProblemSumFactorization<dim>::
get_elliptic_operators() const -> const EllipticOperatorsType &
{
    return elliptic_operators_sf_;
}


template<int dim>
ValueFlags
PoissonProblemSumFactorization<dim>::
get_fill_flags() const
{
    ValueFlags fill_flags = ValueFlags::value |
//                            ValueFlags::gradient |
                            ValueFlags::map_inv_gradient |
                            ValueFlags::measure |
                            ValueFlags::w_measure |
                            ValueFlags::point;

    return fill_flags;
}


template <class PoissonProblemSolver >
void
do_test(const int degree_min, const int degree_max,const int n_elems_per_direction)
{
    const int n_knots = n_elems_per_direction+1;

    for (int degree = degree_min ; degree <= degree_max ; ++degree)
    {
        PoissonProblemSolver poisson(n_knots,degree);
        poisson.assemble();
    }
}


int main(int argc,char **args)

{
    int degree_min = 3;
    int degree_max = 3;
    int n_elems_per_direction = 2;

    int num_ref_lvls = 1;

    bool do_sum_factorization = true;
    bool do_std_quadrature = true;

    for (int ref_lvl=0 ; ref_lvl < num_ref_lvls ; ++ref_lvl, n_elems_per_direction *= 2)
    {
        if (do_sum_factorization)
        {
            cout << "-----------------------------------" << endl;
            cout << "Sum-Factorization -- begin" << endl;
            do_test< PoissonProblemSumFactorization<1> >(degree_min,degree_max,n_elems_per_direction);

            do_test< PoissonProblemSumFactorization<2> >(degree_min,degree_max,n_elems_per_direction);

            do_test< PoissonProblemSumFactorization<3> >(degree_min,degree_max,n_elems_per_direction);
            cout << "Sum-Factorization -- end" << endl;
            cout << "-----------------------------------" << endl << endl;
        }

        if (do_std_quadrature)
        {
            cout << "-----------------------------------" << endl;
            cout << "Standard Quadrature -- begin" << endl;
            do_test< PoissonProblemStandardIntegration<1> >(degree_min,degree_max,n_elems_per_direction);

            do_test< PoissonProblemStandardIntegration<2> >(degree_min,degree_max,n_elems_per_direction);

            do_test< PoissonProblemStandardIntegration<3> >(degree_min,degree_max,n_elems_per_direction);
            cout << "Standard Quadrature -- end" << endl;
            cout << "-----------------------------------" << endl << endl;
        }
    }

    return  0;
}
