
#include <igatools/base/function_lib.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/io/reader.h>
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




template<const int dim>
shared_ptr<Mapping<dim>> nurb_geometry_from_file0()
{
    string input_file = "nurb_geometry0.xml";
    auto map = get_mapping_from_file<dim>(input_file);

    const int n_plot_points = 10;
    Writer<dim> writer(map, n_plot_points);
    string filename = "view_nurb_geometry0";
    writer.save(filename);
    return map;
}


template<const int dim>
shared_ptr<Mapping<dim>> nurb_geometry_from_file1()
{
    string input_file = "nurb_geometry1.xml";
    auto map = get_mapping_from_file<dim>(input_file);

    const int n_plot_points = 10;
    Writer<dim> writer(map, n_plot_points);
    string filename = "view_nurb_geometry1";
    writer.save(filename);
    return map;
}




template <int dim, int dim_field>
void do_test(shared_ptr<Mapping<dim,0>> map0, shared_ptr<Mapping<dim,0>> map1)
{

    using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
    //using RefSpaceField         = NURBSSpace<dim,dim_field,1>;
    shared_ptr<RefSpaceField>   ref_space_field0;
    shared_ptr<RefSpaceField>   ref_space_field1;
    const int deg0 = 2;
    const int deg1 = 2;
    ref_space_field0  = RefSpaceField::create(deg0, map0->get_grid());
    ref_space_field1  = RefSpaceField::create(deg1, map1->get_grid());

    using PushFw      = PushForward<Transformation::h_grad, dim,0>;
    using PhySpace    = PhysicalSpace<RefSpaceField, PushFw>;
    auto space0       = PhySpace::create(ref_space_field0, PushFw::create(map0),0);
    auto space1       = PhySpace::create(ref_space_field1, PushFw::create(map1),1);

    //space0->refine_h(2);
    //space1->refine_h(2);

    space0->print_info(out);
    space1->print_info(out);
}


int main()
{
    LogStream out;
    const int dim(2);
    auto map0=nurb_geometry_from_file0<dim>();
    auto map1=nurb_geometry_from_file1<dim>();
    do_test<dim,1>(map0, map1, out);
    do_test<dim,2>(map0, map1, out);

}
