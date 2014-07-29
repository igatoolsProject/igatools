
#include "../tests.h"

#include <igatools/base/function_lib.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/io/writer.h>


using namespace iga;
using namespace std;




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




template <int dim, int dim_field>
void do_test(vector<shared_ptr<Mapping<dim,0>>> &maps, vector<int> degrees, LogStream &out)
{

    using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
    using RefSpaceField_ptr     = shared_ptr<RefSpaceField>;
    using PushFw                = PushForward<Transformation::h_grad, dim,0>;
    using PhySpace              = PhysicalSpace<RefSpaceField, PushFw>;
    using PhySpace_ptr          = shared_ptr<PhySpace>;

    vector<RefSpaceField_ptr>   ref_spaces_field;
    vector<PhySpace_ptr>        spaces;


    for (int i=0; i!=maps.size(); ++i)
    {
        ref_spaces_field.push_back(RefSpaceField::create(degrees[i], maps[i]->get_grid()));
        spaces.push_back(PhySpace::create(ref_spaces_field[i], PushFw::create(maps[i]),i));
        spaces[i]->get_grid()->print_info(out);
        spaces[i]->get_grid()->refine(3);
        spaces[i]->get_grid()->print_info(out);
        //spaces[i]->get_reference_space()->refine_h(3);
        spaces[i]->print_info(out);
        //spaces[i]->refine_h(3); // Problem here with a value different than 2
        spaces[i]->print_info(out);
    }
}




int main()
{
    const int dim1(1);
    const int dim2(2);
    const int dim3(3);
    vector<int> degrees(1);


    // Mapping $\R \rightarrow \R$, scalar field
    vector<shared_ptr<Mapping<dim1,0>>> maps_d11;
    maps_d11.push_back(create_geo_0<dim1>());
    maps_d11.push_back(create_geo_1<dim1>());
    do_test<dim1,1>(maps_d11, degrees, out);


    // Mapping $\R^2 \rightarrow \R^2$, vector field in $\R^2$
    vector<shared_ptr<Mapping<dim2,0>>> maps_d22;
    maps_d22.push_back(create_geo_0<dim2>());
    maps_d22.push_back(create_geo_1<dim2>());
    do_test<dim2,2>(maps_d22, degrees, out);
    // Mapping $\R^2 \rightarrow \R^2$, scalar field
    vector<shared_ptr<Mapping<dim2,0>>> maps_d21;
    maps_d21.push_back(create_geo_0<dim2>());
    maps_d21.push_back(create_geo_1<dim2>());
    do_test<dim2,1>(maps_d21, degrees, out);


    // Mapping $\R^3 \rightarrow \R^3$, vector field in $\R^3$
    vector<shared_ptr<Mapping<dim3,0>>> maps_d33;
    maps_d33.push_back(create_geo_0<dim3>());
    maps_d33.push_back(create_geo_1<dim3>());
    do_test<dim3,3>(maps_d33, degrees, out);
    // Mapping $\R^3 \rightarrow \R^3$, scalar field
    vector<shared_ptr<Mapping<dim3,0>>> maps_d31;
    maps_d31.push_back(create_geo_0<dim3>());
    maps_d31.push_back(create_geo_1<dim3>());
    do_test<dim3,1>(maps_d31, degrees, out);

    return 0;

}