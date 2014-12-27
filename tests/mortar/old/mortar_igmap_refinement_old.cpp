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
#include "../tests.h"

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


using namespace iga;
using namespace std;


//template<const int dim>
//shared_ptr<Mapping<dim>> nurb_geometry_from_file0()
//{
//    string input_file = "nurb_geometry0.xml";
//    auto map = get_mapping_from_file<dim>(input_file);
//
//    const int n_plot_points = 10;
//    Writer<dim> writer(map, n_plot_points);
//    string filename = "view_nurb_geometry0";
//    writer.save(filename);
//  return map;
//}


//template<const int dim>
//shared_ptr<Mapping<dim>> nurb_geometry_from_file1()
//{
//    string input_file = "nurb_geometry1.xml";
//    auto map = get_mapping_from_file<dim>(input_file);
//
//   const int n_plot_points = 10;
//    Writer<dim> writer(map, n_plot_points);
//    string filename = "view_nurb_geometry1"+ to_string(dim) + "d" ;;
//    writer.save(filename);
//  return map;
//}

template<const int dim>
shared_ptr<Mapping<dim>> geometry_from_file(int patch_nb)
{
    string input_file = "geometry"+ to_string(patch_nb)+".xml";
    auto map = get_mapping_from_file<dim>(input_file);

    const int n_plot_points = 10;
    Writer<dim> writer(map, n_plot_points);
    string filename = "paraview_geometry"+ to_string(patch_nb);
    writer.save(filename);
    return map;
}




//template <int dim, int dim_field>
//void do_test(shared_ptr<Mapping<dim,0>> map0, shared_ptr<Mapping<dim,0>> map1){
//
//  using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
//  //using RefSpaceField         = NURBSSpace<dim,dim_field,1>;
//  shared_ptr<RefSpaceField>   ref_space_field0;
//  shared_ptr<RefSpaceField>   ref_space_field1;
//  const int deg0 = 2;
//  const int deg1 = 2;
//  ref_space_field0  = RefSpaceField::create(deg0, map0->get_grid());
//  ref_space_field1  = RefSpaceField::create(deg1, map1->get_grid());
//
//    using PushFw      = PushForward<Transformation::h_grad, dim,0>;
//    using PhySpace    = PhysicalSpace<RefSpaceField, PushFw>;
//    auto space0       = PhySpace::create(ref_space_field0, PushFw::create(map0),0);
//  auto space1       = PhySpace::create(ref_space_field1, PushFw::create(map1),1);
//
//  //space0->refine_h(2);
//  //space1->refine_h(2);
//
//  //space0->print_info(out);
//  //space1->print_info(out);
//}



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
        spaces[i]->print_info(out);
        //spaces[i]->refine_h();
        spaces[i]->print_info(out);
    }
}



int main()
{
    const int patch_nb(2);
    const int dim2(2);
    vector<int> degrees(2);

    vector<shared_ptr<Mapping<dim2,0>>> maps_d22;
    for (int i=0; i!=patch_nb; ++i)
    {
        maps_d22.push_back(geometry_from_file<dim2>(i));
    }
    do_test<dim2,2>(maps_d22, degrees,out);

    vector<shared_ptr<Mapping<dim2,0>>> maps_d21;
    for (int i=0; i!=patch_nb; ++i)
    {
        maps_d21.push_back(geometry_from_file<dim2>(i));
    }
    do_test<dim2,1>(maps_d21, degrees,out);


    //LogStream out;
    //const int dim(2);
    //auto map0=nurb_geometry_from_file0<dim>();
    //auto map1=nurb_geometry_from_file1<dim>();
    //do_test<dim,1>(map0, map1);
    ////do_test<dim,2>(map0, map1, out);

}
