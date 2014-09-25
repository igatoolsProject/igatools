
#include "../tests.h"

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


using namespace iga;
using namespace std;




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
//void do_test(shared_ptr<Mapping<dim,0>> map0, shared_ptr<Mapping<dim,0>> map1, LogStream& out){
void do_test(vector<shared_ptr<Mapping<dim,0>>> &maps, vector<int> degrees, LogStream &out)
{


    //Assert(maps.size() !=0 ,
    //     ExcDimensionMismatch(maps.size(),0);

    using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
    using RefSpaceField_ptr     = shared_ptr<RefSpaceField>;
    using PushFw                = PushForward<Transformation::h_grad, dim,0>;
    using PhySpace              = PhysicalSpace<RefSpaceField, PushFw>;
    using PhySpace_ptr          = shared_ptr<PhySpace>;

    vector<RefSpaceField_ptr>   ref_spaces_field;
    vector<PhySpace_ptr>        spaces;

    //shared_ptr<RefSpaceField>   ref_space_field0;
    //shared_ptr<RefSpaceField>   ref_space_field1;
    //Assert(maps.size() !=degrees.size(),
    //            ExcDimensionMismatch(maps.size(),degrees.size());
    for (int i=0; i!=maps.size(); ++i)
    {
        ref_spaces_field.push_back(RefSpaceField::create(degrees[i], maps[i]->get_grid()));
        spaces.push_back(PhySpace::create(ref_spaces_field[i], PushFw::create(maps[i]),i));
        spaces[i]->print_info(out);
        spaces[i]->refine_h(4);
        spaces[i]->print_info(out);
    }


    //ref_space_field0  = RefSpaceField::create(deg0, map0->get_grid());
    //ref_space_field1  = RefSpaceField::create(deg1, map1->get_grid());

    //using PushFw      = PushForward<Transformation::h_grad, dim,0>;
    //using PhySpace    = PhysicalSpace<RefSpaceField, PushFw>;

    //auto space0       = PhySpace::create(ref_space_field0, PushFw::create(map0),0);
    //auto space1       = PhySpace::create(ref_space_field1, PushFw::create(map1),1);

    //space0->refine_h(2); Pb of refinement
    //space1->refine_h(2);

    //space0->print_info(out);
    //space1->print_info(out);
}




int main()
{
    const int dim1(1);
    const int dim2(2);
    const int dim3(3);
    vector<int> degrees(1);
    //using RefSpaceMap_d1     = BSplineSpace<dim1,dim1,1>;
    //using IgMapping_d1_ptr   = shared_ptr<IgMapping<RefSpaceMap_d1>>;
    //vector<IgMapping_d1_ptr> maps_d11;


    vector<shared_ptr<Mapping<dim1,0>>> maps_d11;
    maps_d11.push_back(create_geo_0<dim1>());
    maps_d11.push_back(create_geo_1<dim1>());
    do_test<dim1,1>(maps_d11, degrees,out);
    //auto map0_d1=create_geo_0<dim1>();
    //auto map1_d1=create_geo_1<dim1>();
    //do_test<dim1,1>(map0_d1, map1_d1, out);
    //do_test<dim1,1>(maps_d11, out);


//  vector<shared_ptr<Mapping<dim2,0>>> maps_d22;
//  maps_d22.push_back(create_geo_0<dim2>());
//  maps_d22.push_back(create_geo_1<dim2>());
//  do_test<dim2,2>(maps_d22, degrees,out);
//  vector<shared_ptr<Mapping<dim2,0>>> maps_d21;
//  maps_d21.push_back(create_geo_0<dim2>());
//  maps_d21.push_back(create_geo_1<dim2>());
//  do_test<dim2,1>(maps_d21, degrees,out);


//  vector<shared_ptr<Mapping<dim3,0>>> maps_d33;
//  maps_d33.push_back(create_geo_0<dim3>());
//  maps_d33.push_back(create_geo_1<dim3>());
//  do_test<dim3,3>(maps_d33, degrees,out);
//  vector<shared_ptr<Mapping<dim3,0>>> maps_d31;
//  maps_d31.push_back(create_geo_0<dim3>());
//  maps_d31.push_back(create_geo_1<dim3>());
//  do_test<dim3,1>(maps_d31, degrees,out);

//  const int dim2(2);
    //using RefSpaceMap_d2     = BSplineSpace<dim2,dim2,1>;
    //using IgMapping_d2_ptr   = shared_ptr<IgMapping<RefSpaceMap_d2>>;
//  vector<shared_ptr<Mapping<dim2,0>>> maps_d22;
//  maps_d22.push_back(create_geo_0<dim2>);
//  maps_d22.push_back(create_geo_1<dim2>);
//  do_test<dim2,2>(maps_d22, out);

//  vector<shared_ptr<Mapping<dim2,0>>> maps_d21;
//  maps_d21.push_back(create_geo_0<dim2>);
//  maps_d21.push_back(create_geo_1<dim2>);
//  do_test<dim2,1>(maps_d21, out);

    //auto map0_d2=create_geo_0<dim2>();
    //auto map1_d2=create_geo_1<dim2>();
    //do_test<dim2,2>(map0_d2, map1_d2, out);
    //do_test<dim2,1>(map0_d2, map1_d2, out);


    //  const int dim3(3);
    //  auto map0_d3=create_geo_0<dim3>();
    //  auto map1_d3=create_geo_1<dim3>();
    //  do_test<dim3,3>(map0_d3, map1_d3, out);
    /////PB do_test<dim3,2>(map0_d3, map1_d3);
    //  do_test<dim3,1>(map0_d3, map1_d3, out);

}