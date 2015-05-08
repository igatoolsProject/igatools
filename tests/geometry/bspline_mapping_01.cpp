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
/*
 *  Test the IgMapping class on Bspline space
 *  The map is the identity of degree one.
 *
 *  author: pauletti
 *  date: 2013-10-04
 *
 */

#include "../tests.h"

#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>
#include <igatools/base/ig_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_space.h>
//#include <igatools/basis_functions/nurbs_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/base/function_element.h>

using namespace EpetraTools;





template <int dim,int codim>
void serialize_deserialize(std::shared_ptr<IgFunction<ReferenceSpace<dim, dim+codim>>> F)
{
    out.begin_item("Original IgFunction:");
    F->print_info(out);
    out.end_item();

    std::string template_strdetails = "_dim" + std::to_string(dim) + "_codim" + std::to_string(codim);
    std::string filename = "ig_function" + template_strdetails + ".xml";
    std::string tag_name = "IgFunction" + template_strdetails;
    {
        // serialize the IgFunction object to an xml file
        std::ofstream xml_ostream(filename);
        boost::archive::xml_oarchive xml_out(xml_ostream);
//        xml_out.template register_type<BSplineSpace<dim,dim+codim>>();

        xml_out << boost::serialization::make_nvp(tag_name.c_str(),F);
        xml_ostream.close();
    }

    F.reset();
    {
        // de-serialize the IgFunction object from an xml file
        std::ifstream xml_istream(filename);
        boost::archive::xml_iarchive xml_in(xml_istream);
//        xml_in.template register_type<BSplineSpace<dim,dim+codim>>();

        xml_in >> BOOST_SERIALIZATION_NVP(F);
        xml_istream.close();
    }
    out.begin_item("IgFunction after serialize-deserialize:");
    F->print_info(out);
    out.end_item();
//*/
}



template <int dim, int codim=0>
void bspline_map(const int deg = 1)
{
    OUTSTART

    const int sub_dim = dim;
    using Space = BSplineSpace<dim, dim+codim>;
    using RefSpace = ReferenceSpace<dim, dim+codim>;
    using Function = IgFunction<RefSpace>;
    using Mapping   = Mapping<dim, codim>;

    auto grid = CartesianGrid<dim>::create(2);
    auto space = Space::create(deg, grid);

    auto c_p = EpetraTools::create_vector(space, "active");
    auto &control_pts = *c_p;

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

    auto F = Function::create(space, c_p);


//    serialize_deserialize<dim,codim>(F);


    auto map = Mapping::create(F);

    auto quad = QGauss<dim>(3);
    auto flag =  ValueFlags::value| ValueFlags::gradient
                 | ValueFlags::hessian;

    map->template reset<sub_dim>(flag, quad);

    auto elem = map->begin();
    auto end  = map->end();
    const int s_id = 0;

    map->template init_cache<sub_dim>(elem);
    for (; elem != end; ++elem)
    {
        map->template fill_cache<sub_dim>(elem, s_id);
        out << "Values (x1,x2,...):" << endl;
        elem->template get_values<_Value,sub_dim>(s_id).print_info(out);
        out << endl;
        out << "Gradients:" << endl;
        elem->template get_values<_Gradient,sub_dim>(s_id).print_info(out);
        out << endl;
//        elem->template get_values<2,sub_dim>(s_id).print_info(out);
//        out << endl;
    }

    OUTEND
}

int main()
{
    out.depth_console(10);

    bspline_map<1>();
    bspline_map<2>();
    bspline_map<3>();

}
