
#include <igatools/base/function_lib.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/geometry/topology.h>
//#include <igatools/io/reader.h>
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
    auto knots = CartesianGrid<dim>::create(2);

    using RefSpaceMap = BSplineSpace<dim,dim,1>;

    const int p = 1;
    auto ref_space_map = RefSpaceMap::create(p,knots);

    Assert(control_pts.size() == ref_space_map->get_num_basis(),
           ExcDimensionMismatch(control_pts.size(),ref_space_map->get_num_basis()));

    return IgMapping<RefSpaceMap>::create(ref_space_map, control_pts);
}




template <int dim>
shared_ptr<Mapping<dim,0>> create_geo_0()
{
    const int n_basis = dim * std::pow(2,dim);
    vector<Real> control_pts(n_basis);

    if (dim == 1)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;
    }
    else if (dim == 2) // x coordinates of all the points, then y
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
shared_ptr<Mapping<dim,0>> create_geo_1()
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



//real multiply_vect(real i, real val) {i=i*val; return i}



int main()
{


    LogStream out;
    const int dim(2);
    const int dim_field(1);
    vector<int> degrees(1);
    const int deg(1);


    // Mapping $\R^2 \rightarrow \R^2$, scalar field
    vector<shared_ptr<Mapping<dim,0>>> maps_d21;
    maps_d21.push_back(create_geo_0<dim>());
    maps_d21.push_back(create_geo_1<dim>());

    using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
    using RefSpaceField_ptr     = shared_ptr<RefSpaceField>;
    using PushFw                = PushForward<Transformation::h_grad, dim,0>;
    using PhySpace              = PhysicalSpace<RefSpaceField, PushFw>;
    using PhySpace_ptr          = shared_ptr<PhySpace>;

    using Func = Function<dim, dim>;
    using Gradient = typename Func::Gradient;

    vector<RefSpaceField_ptr>   ref_spaces_field;
    vector<PhySpace_ptr>        spaces;

    MultiPatchSpace<PhySpace> domain_multip;
    domain_multip.arrangement_open();
    const boundary_id mortar_id = 7;
    vector<set<boundary_id>> mortar_faces;
    mortar_faces.push_back({1});
    mortar_faces.push_back({0});


    for (uint i=0; i!=maps_d21.size(); ++i)
    {
        ref_spaces_field.push_back(RefSpaceField::create(degrees[i], maps_d21[i]->get_grid()));
        spaces.push_back(PhySpace::create(ref_spaces_field[i], PushFw::create(maps_d21[i]),i));
        //      spaces[i]->print_info(out);
        if (i==1)
        {
            spaces[i]->refine_h(2);
        }// Problem here with a value different than 2
        //      spaces[i]->print_info(out);

        auto tmp_mortar_faces=mortar_faces[i];
        auto crr_mortar=tmp_mortar_faces.begin(), end_mortar=tmp_mortar_faces.end();
        for (; crr_mortar!=end_mortar; ++crr_mortar)
        {
            spaces[i]->get_grid()->set_boundary_id(*crr_mortar,mortar_id);
        }


        domain_multip.add_patch(spaces[i]);
    }

    domain_multip.arrangement_close();


    // multiplier space
    const int degree_multiplier=1;
    auto ref_space_multiplier = RefSpaceField::create(degree_multiplier,spaces[0]->get_grid());
    auto multiplier_space=PhySpace::create(ref_space_multiplier, PushFw::create(maps_d21[0]),1);





    vector< Index > joint_to_0;
    vector< Index > joint_to_1;
    auto grid0=spaces[0]->get_grid();
    auto grid1=spaces[1]->get_grid();

    vector<shared_ptr<CartesianGrid<dim>>> grid_J;
    grid_J.push_back(make_shared<CartesianGrid<dim>>(build_cartesian_grid_union(*grid0, *grid1, joint_to_0, joint_to_1)));  //pointeur for the iterator
    grid_J.push_back(make_shared<CartesianGrid<dim>>(build_cartesian_grid_union(*grid0, *grid1, joint_to_0, joint_to_1)));
    //map0->print_info(out);
    //map1->print_info(out);
    grid_J[0]->print_info(out);
    grid_J[1]->print_info(out);
    //auto bd_id0_0=grid_Js->get_boundary_id(0);
    //auto tmp0_bd_id0_0=grid_J[1]->get_boundary_id(0);
    //auto tmp0_bd_id0_1=grid_J[1]->get_boundary_id(1);
    //auto tmp0_bd_id0_2=grid_J[1]->get_boundary_id(2);
    //auto tmp0_bd_id0_3=grid_J[1]->get_boundary_id(3);
    //out<< to_string(tmp0_bd_id0_0)<< to_string(tmp0_bd_id0_1)<< to_string(tmp0_bd_id0_2)<< to_string(tmp0_bd_id0_3)<<endl;


    for (uint i=0; i!=maps_d21.size(); ++i)
    {
        auto tmp_mortar_faces=mortar_faces[i];
        auto crr_mortar=tmp_mortar_faces.begin(), end_mortar=tmp_mortar_faces.end();
        for (; crr_mortar!=end_mortar; ++crr_mortar)
        {
            grid_J[i]->set_boundary_id(*crr_mortar,mortar_id);
        }
    }



    //auto tmp_bd_id0_0=grid_J[1]->get_boundary_id(0);
    //auto tmp_bd_id0_1=grid_J[1]->get_boundary_id(1);
    //auto tmp_bd_id0_2=grid_J[1]->get_boundary_id(2);
    //auto tmp_bd_id0_3=grid_J[1]->get_boundary_id(3);
    //out<< to_string(tmp_bd_id0_0)<< to_string(tmp_bd_id0_1)<< to_string(tmp_bd_id0_2)<< to_string(tmp_bd_id0_3)<<endl;


    const Quadrature<dim>   elem_quad(QGauss<dim>(deg+1));
    const Quadrature<dim-1> face_quad(QGauss<dim-1>(deg+1));
    auto  face_quad_gbl  = extend_face_quad(face_quad,1);
    auto pts_unit_domain = face_quad_gbl.get_points().get_flat_cartesian_product();
    auto w_unit_domain   = face_quad_gbl.get_weights().get_flat_tensor_product();
    out<<to_string(w_unit_domain[0])<<endl;


    //ATTENTION INTEGRATION SUR FACE
    //  auto elem_Js = *grid_J[0]->begin();
    //  const auto elem_Js_end = *grid_J[0]->end();
    //  ValueFlags fill_flags = ValueFlags::w_measure;
    //  elem_Js.init_values(fill_flags, elem_quad);
    for (auto elem_Js: *grid_J[0])
    {
        //   for (; elem_Js != elem_Js_end; ++elem_Js){
        //      elem_Js.fill_values();
        if (elem_Js.is_boundary())
        {
            for (Index face_id = 0; face_id < UnitElement<dim>::n_faces; ++face_id)
            {
                if (elem_Js.is_boundary(face_id))
                {
                    auto    grid_tmp=elem_Js.get_grid();
                    out<<to_string(grid_tmp->get_boundary_id(face_id));
                    if (grid_tmp->get_boundary_id(face_id)==mortar_id)
                    {
                        auto elem_Js_id = elem_Js.get_flat_index();
                        auto quad_ref_domain = elem_Js.transform_points_unit_to_reference(pts_unit_domain);

                        //elem_Js.init_face_values(face_id, fill_flags, face_quad);
                        //elem_Js.fill_face_values(face_id);

                        // Problem should be filled without any cache
                        //const auto w_meas=elem_Js.get_measure();
                        const auto w_meas=static_cast<CartesianGridElement<dim>&>(elem_Js).get_measure(FaceTopology<dim>(1));
                        out<<w_meas<<endl;
                        //


                        auto elem_0_id = joint_to_0[elem_Js_id];
                        auto elem_1_id = joint_to_1[elem_Js_id];
                        auto elem0=spaces[0]->get_element(elem_0_id);
                        auto elem_m=multiplier_space->get_element(elem_0_id);
                        auto elem1=spaces[1]->get_element(elem_1_id);

                        vector<Gradient> grad_0;
                        maps_d21[0]->evaluate_gradients_at_points(pts_unit_domain, grad_0);

                        auto det=determinant<dim,dim>(grad_0[0]);
                        out<<det<<"ICI"<<endl;


                        auto basis0=elem0.evaluate_basis_values_at_points(pts_unit_domain);
                        auto basis1=elem1.evaluate_basis_values_at_points(pts_unit_domain);

                        auto phi_0 = basis0.get_function_view(0);
                        out<<phi_0[0]<<phi_0[1];
                        auto phi_1 = basis1.get_function_view(0);
                        out<<phi_1[0]<<phi_1[1];

                        //      out<<"elem0"  <<to_string(elem_0_id) << endl;
                        out<<"HERE";
                    }
                }
            }
        }
    }
    auto elem1=spaces[1]->get_element(1);
    //out << elem1.get_flat_index();

    return  0;
}
