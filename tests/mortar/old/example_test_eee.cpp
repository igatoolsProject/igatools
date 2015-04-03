//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
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

#include <igatools/base/function_lib.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/basis_functions/multi_patch_space.h>
#include <igatools/io/writer.h>
#include <iostream>

#include <igatools/base/quadrature_lib.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/linear_solver.h>
#include <igatools/linear_algebra/dof_tools.h>

#include <math.h>


using namespace iga;
using namespace std;
using functions::ConstantFunction;
using space_tools::project_boundary_values;
using dof_tools::apply_boundary_values;
using dof_tools::get_sparsity_pattern;



template <int dim>
shared_ptr<Mapping<dim,0>> create_map(const vector<Real> &control_pts)
{
    auto knots  = CartesianGrid<dim>::create(2);

    using RefSpaceMap = BSplineSpace<dim,dim,1>;

    const int p = 1;
    auto ref_space_map = RefSpaceMap::create(p,knots);

    Assert(control_pts.size() == ref_space_map->get_num_basis(),
           ExcDimensionMismatch(control_pts.size(),ref_space_map->get_num_basis()));

    return IgMapping<RefSpaceMap>::create(ref_space_map, control_pts);
}

template <int dim>
shared_ptr<Mapping<dim,0>>
                        create_geo_0()
{
    const int n_basis = dim * std::pow(2,dim);
    vector<Real> control_pts(n_basis);

    if (dim == 1)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;
    }
    else if (dim == 2)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
    }
    else if (dim == 3)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

    }
    return create_map<dim>(control_pts);
}



template <int dim>
shared_ptr<Mapping<dim,0>>
                        create_geo_1()
{
    const int n_basis = dim * std::pow(2,dim);
    vector<Real> control_pts(n_basis);

    if (dim == 1)
    {
        control_pts[0] = 1.0 ;
        control_pts[1] = 2.0 ;
    }
    else if (dim == 2)
    {
        control_pts[0] = 1.0 ;
        control_pts[1] = 2.0 ;

        control_pts[2] = 1.0 ;
        control_pts[3] = 2.0 ;

        control_pts[4] = 0.0 ;
        control_pts[5] = 0.0 ;

        control_pts[6] = 1.0 ;
        control_pts[7] = 1.0 ;
    }
    else if (dim == 3)
    {
        control_pts[0] = 1.0 ;
        control_pts[1] = 2.0 ;
        control_pts[2] = 1.0 ;

        control_pts[3] = 2.0 ;
        control_pts[4] = 1.0 ;
        control_pts[5] = 2.0 ;

        control_pts[6] = 1.0 ;
        control_pts[7] = 2.0 ;
        control_pts[8] = 0.0 ;

        control_pts[9] = 0.0 ;
        control_pts[10] = 1.0 ;
        control_pts[11] = 1.0 ;

        control_pts[12] = 0.0 ;
        control_pts[13] = 0.0 ;
        control_pts[14] = 1.0 ;

        control_pts[15] = 1.0 ;
        control_pts[16] = 0.0 ;
        control_pts[17] = 0.0 ;

        control_pts[18] = 0.0 ;
        control_pts[19] = 0.0 ;
        control_pts[20] = 1.0 ;

        control_pts[21] = 1.0 ;
        control_pts[22] = 1.0 ;
        control_pts[23] = 1.0 ;
    }

    return create_map<dim>(control_pts);
}



// Dirichlet_function
template <int dim, int dim_field, int rank=1>
class D_function : public Function<dim,dim_field,rank>
{
public:
    D_function() {};

    void evaluate(const std::vector< typename Function<dim,dim_field,rank>::Point> &points,
                  std::vector< typename Function<dim, dim_field,rank>::Value> &values) const

    {
        const int num_points = points.size() ;
        Assert(num_points == values.size(), ExcDimensionMismatch(num_points, values.size())) ;

        for (int i =0; i<num_points; ++i)
        {
            //auto pt_i = points[i];
            values[i]= 1.0;//sin(5*pt_i[0])*cos(6*pt_i[1]);
        }

    }
};



int main()
{

    LogStream out;
    const int dim(2);
    const int dim_field(1);
    auto map0=create_geo_0<dim>();
    auto map1=create_geo_1<dim>();

    using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
    //using RefSpaceField         = NURBSSpace<dim,dim_field,1>;
    shared_ptr<RefSpaceField>   ref_space_field0;
    shared_ptr<RefSpaceField>   ref_space_field1;
    const int deg0 = 1;
    const int deg1 = 1;
    ref_space_field0  = RefSpaceField::create(deg0, map0->get_grid());
    ref_space_field1  = RefSpaceField::create(deg1, map1->get_grid());

    using PushFw      = PushForward<Transformation::h_grad, dim,0>;
    using PhySpace    = PhysicalSpace<RefSpaceField, PushFw>;
    using PhySpace_ptr= shared_ptr<PhySpace>;
    vector<PhySpace_ptr> spaces;
    spaces.push_back(PhySpace::create(ref_space_field0, PushFw::create(map0),0));
    spaces.push_back(PhySpace::create(ref_space_field1, PushFw::create(map1),1));
    //auto space0       = PhySpace::create(ref_space_field0, PushFw::create(map0),0);
    //auto space1       = PhySpace::create(ref_space_field1, PushFw::create(map1),1);


    MultiPatchSpace<PhySpace> domain_multip;
    domain_multip.arrangement_open();
    for (auto crr_space: spaces)
    {
        domain_multip.add_patch(crr_space);
    }
    //domain_multip.add_patch(space0);

    //domain_multip.add_patch(space1);
    domain_multip.add_interface(InterfaceType::C0_strong,
                                space[0],1,
                                space[1],0);
    domain_multip.arrangement_close();
    domain_multip.print_info(out);
    auto dof_manager=domain_multip.get_dofs_manager();
    dof_manager->print_info(out);

    //dofs_manager->linear_constraints_open();
    //
    // To DO
    //
    //dofs_manager->linear_constraints_close();



    /////const boundary_id nh_id = 0;
    ////const boundary_id n_id = 1;
    ////const boundary_id dh_id = 2;
    const boundary_id d_id = 3;
    /////set<boundary_id> nh_faces = {1};
    set<boundary_id> d_faces    = {0};
    //// Change the face id 0->NH // 1->NNH // 2->DH // 3->DNH
    /////auto crr_nh=nh_faces.begin(), end_nh=nh_faces.end();
    /////for (; crr_nh!=end_nh; ++crr_nh)
    /////{
    /////   space0->get_grid()->set_boundary_id(*crr_nh,nh_id);
    /////}

    auto crr_d=d_faces.begin(), end_d=d_faces.end();
    for (; crr_d!=end_d; ++crr_d)
    {
        space0->get_grid()->set_boundary_id(*crr_d,d_id);
    }


    const Quadrature<dim>   elem_quad0(QGauss<dim>(deg0+1));
    const Quadrature<dim-1> face_quad0(QGauss<dim-1>(deg0+1));
    const Quadrature<dim>   elem_quad1(QGauss<dim>(deg1+1));
    const Quadrature<dim-1> face_quad1(QGauss<dim-1>(deg1+1));


    std::shared_ptr<Matrix<LAPack::trilinos>> matrix0;
    std::shared_ptr<Vector<LAPack::trilinos>> rhs0;
    std::shared_ptr<Vector<LAPack::trilinos>> solution0;

    const auto n_basis0 = space0->get_num_basis();
    matrix0   = Matrix<LAPack::trilinos>::create(get_sparsity_pattern<PhySpace>(space0));
    rhs0      = Vector<LAPack::trilinos>::create(n_basis0);
    solution0 = Vector<LAPack::trilinos>::create(n_basis0);

    DenseMatrix loc_mat0(n_basis0, n_basis0);
    DenseVector loc_rhs0(n_basis0);
    vector<Index> loc_dofs0(n_basis0);

    const int n_qp = elem_quad0.get_num_points();
    ConstantFunction<dim> f({0.0});
    using Value = typename Function<dim>::Value;
    vector<Value> f_values(n_qp);

    auto elem = space0->begin();
    const auto elem_end = space0->end();
    ValueFlags fill_flags = ValueFlags::value | ValueFlags::gradient |
                            ValueFlags::w_measure | ValueFlags::point;
    elem->init_values(fill_flags, elem_quad0);

    for (; elem != elem_end; ++elem)
    {
        loc_mat0.clear();
        loc_rhs0.clear();
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
        matrix->add_block(loc_dofs, loc_dofs, loc_mat);
        rhs->add_block(loc_dofs, loc_rhs);
    }

    matrix->FillComplete();

    //ConstantFunction<dim> g({0.0});
    D_function<dim, dim_field, 1> dfunction;
    //vector< Value > dfunc_at_points(points.size());
    //function.evaluate(points, d_func_at_points);
    std::map<Index, Real> values;
    //project_boundary_values<PhySpace,LinearAlgebraPackage::trilinos>(g,space0,face_quad0,dir_id,values);
    project_boundary_values<PhySpace,LAPack::trilinos>(dfunction,space0,face_quad0,d_id,values);
    apply_boundary_values(values, *matrix, *rhs, *solution);
    solution->print(out);


    using LinSolver = LinearSolver<LAPack::trilinos>;
    LinSolver solver(LinSolver::SolverType::CG);
    solver.solve(*matrix, *rhs, *solution);
    solution->print(out);

    const int n_plot_points = 2;
    Writer<dim> writer(map0, n_plot_points);

    writer.add_field(space0, *solution, "solution");
    string filename = "poisson_problem_P0" ;
    writer.save(filename);



    return 0;




}