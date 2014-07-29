#include "../tests.h"

#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/io/reader.h>
#include <igatools/basis_functions/space_tools.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/io/writer.h>


using namespace iga;
using namespace std;




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
        spaces[i]->refine_h(3); //Problem for the refinement
        spaces[i]->print_info(out);
    }
}



int main()
{
    const int patch_nb(2);
    const int dim2(2);
    vector<int> degrees(2);


    // Mapping $\R^2 \rightarrow \R^2$, vector field in $\R^2$
    vector<shared_ptr<Mapping<dim2,0>>> maps_d22;
    for (int i=0; i!=patch_nb; ++i)
    {
        maps_d22.push_back(geometry_from_file<dim2>(i));
    }
    do_test<dim2,2>(maps_d22, degrees,out);

    // Mapping $\R^2 \rightarrow \R^2$, scalar field
    vector<shared_ptr<Mapping<dim2,0>>> maps_d21;
    for (int i=0; i!=patch_nb; ++i)
    {
        maps_d21.push_back(geometry_from_file<dim2>(i));
    }
    do_test<dim2,1>(maps_d21, degrees,out);

}
