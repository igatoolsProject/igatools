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
 *  Test for SubFunction class
 *  author: pauletti
 *  date: Oct 12, 2014
 */

#include "../tests.h"

#include <igatools/functions/functions_container.h>
#include <igatools/base/quadrature_lib.h>


using std::shared_ptr;
using std::static_pointer_cast;


void serialize_deserialize(std::shared_ptr<FunctionsContainer> funcs_container)
{
    OUTSTART

    out.begin_item("Original FunctionsContainer:");
    funcs_container->print_info(out);
    out.end_item();

    std::string filename = "functions_container.xml";
    std::string tag_name = "FunctionsContainer";
    {
        // serialize the PhysicalSpace object to an xml file
        std::ofstream xml_ostream(filename);
        OArchive xml_out(xml_ostream);

        xml_out << boost::serialization::make_nvp(tag_name.c_str(),funcs_container);
        xml_ostream.close();
    }

    funcs_container.reset();
    {
        // de-serialize the PhysicalSpace object from an xml file
        std::ifstream xml_istream(filename);
        IArchive xml_in(xml_istream);

        xml_in >> BOOST_SERIALIZATION_NVP(funcs_container);
        xml_istream.close();
    }
    out.begin_item("FunctionsContainer after serialize-deserialize:");
    funcs_container->print_info(out);
    out.end_item();
//*/

    OUTEND
}



template <int dim,int codim,int range>
using Func = Function<dim,codim,range,1>;

void print_container(std::shared_ptr<FunctionsContainer> funcs_container)
{
    const auto mappings_dim_2_codim_0 = funcs_container->template get_mappings_dim_codim<2,0>();
    out.begin_item("Mappings with dimension 2 and codimension 0:");
    int m_counter = 0;
    for (const auto &m : mappings_dim_2_codim_0)
    {
        auto mapping = m.first;
        auto name    = m.second; // this is the string we associated to the mapping object when we used insert_mapping()

        auto flag = ValueFlags::point | ValueFlags::value;
        QUniform<2> quad(2);
        mapping->reset(flag, quad);

        out.begin_item("Mapping before iterators");
        mapping->print_info(out);
        out.end_item();

        auto elem = mapping->begin();
        auto end  = mapping->end();

#if 0
        try
        {
            const shared_ptr<IgFunction<2, 0, 2, 1>> ig_fun_ptr = std::dynamic_pointer_cast<IgFunction<2, 0, 2, 1>> (mapping);
            if (ig_fun_ptr == nullptr)
            {
                out << "Identitity func start" << endl;
                const shared_ptr<IdentityFunction<2, 2>> id_fun_ptr = std::dynamic_pointer_cast<IdentityFunction<2, 2>> (mapping);
                const IdentityFunction<2, 2> &id_fun = *id_fun_ptr;
                out << "Identitity func stop" << endl;
            }
            else
            {
                out << "Igfunc start" << endl;
                const IgFunction<2, 0, 2, 1> &ig_fun = *ig_fun_ptr;
                out << "Igfunc stop" << endl;
            }
        }
        catch (const std::bad_weak_ptr &e)
        {
            out << "Catching " << e.what() << endl;
        }
#endif
        out << "Mapping[" << m_counter++ << "]   name= " << name << std::endl;

        const auto &funcs_dim_2_codim_0_range_1_rank_1 =
            funcs_container->template get_functions_associated_to_mapping<2,0,1,1>(mapping);
        out.begin_item("Functions<2,0,1,1>:");
        int f_counter = 0;
        for (const auto &f : funcs_dim_2_codim_0_range_1_rank_1)
        {
            out << "Function[" << f_counter++ << "]"
                << "   ID:" << f.first
                << "   name= " << f.second.name_ << std::endl;
        }
        out.end_item();

    }
}

void deserialize_only()
{
    std::string filename = "functions_container.xml";

    std::ifstream xml_istream(filename);
    IArchive xml_in(xml_istream);

    auto funcs_container = std::make_shared<FunctionsContainer>();
    xml_in >> BOOST_SERIALIZATION_NVP(funcs_container);
    xml_istream.close();

    out.begin_item("Inside deserialize_only()");
    print_container(funcs_container);
    out.end_item();
}


void do_test()
{
    int n_elem_per_side = 2;
    auto grid_1 = CartesianGrid<1>::create(n_elem_per_side+1);
    auto grid_2 = CartesianGrid<2>::create(n_elem_per_side+1);
    auto grid_3 = CartesianGrid<3>::create(n_elem_per_side+1);
//    create_fun<2, 0, 2>();

    auto func_identity_1_1 = IdentityFunction<1,1>::create(grid_1);
    auto func_identity_2_2 = IdentityFunction<2,2>::create(grid_2);
    auto func_identity_3_3 = IdentityFunction<3,3>::create(grid_3);



    const int deg = 3;
    auto bsp_space_1_1 = BSplineSpace<1,1,1>::create(deg, grid_1);
    auto bsp_space_2_1 = BSplineSpace<2,1,1>::create(deg, grid_2);
    auto bsp_space_3_1 = BSplineSpace<3,1,1>::create(deg, grid_3);
    auto bsp_space_2_2 = BSplineSpace<2,2,1>::create(deg, grid_2);
    auto bsp_space_3_3 = BSplineSpace<3,3,1>::create(deg, grid_3);
    auto bsp_space_2_3 = BSplineSpace<2,3,1>::create(deg, grid_2);


    Epetra_SerialComm comm;
    auto bsp_coeff_1_1 = EpetraTools::create_vector(
                             EpetraTools::create_map(*bsp_space_1_1, "active", comm));
    (*bsp_coeff_1_1)[0] = 1.;


    auto bsp_coeff_2_1 = EpetraTools::create_vector(
                             EpetraTools::create_map(*bsp_space_2_1, "active", comm));
    (*bsp_coeff_2_1)[0] = 1.;


    auto bsp_coeff_3_1 = EpetraTools::create_vector(
                             EpetraTools::create_map(*bsp_space_3_1, "active", comm));
    (*bsp_coeff_3_1)[0] = 1.;


    auto bsp_coeff_2_2 = EpetraTools::create_vector(
                             EpetraTools::create_map(*bsp_space_2_2, "active", comm));
    (*bsp_coeff_2_2)[0] = 1.;


    auto bsp_coeff_3_3 = EpetraTools::create_vector(
                             EpetraTools::create_map(*bsp_space_3_3, "active", comm));
    (*bsp_coeff_3_3)[0] = 1.;


    auto bsp_coeff_2_3 = EpetraTools::create_vector(
                             EpetraTools::create_map(*bsp_space_2_3, "active", comm));
    (*bsp_coeff_2_3)[0] = 1.;


    auto bsp_func_1_1 = IgFunction<1,0,1,1>::create(bsp_space_1_1, bsp_coeff_1_1);
    auto bsp_func_2_1 = IgFunction<2,0,1,1>::create(bsp_space_2_1, bsp_coeff_2_1);
    auto bsp_func_3_1 = IgFunction<3,0,1,1>::create(bsp_space_3_1, bsp_coeff_3_1);
    auto bsp_func_2_2 = IgFunction<2,0,2,1>::create(bsp_space_2_2, bsp_coeff_2_2);
    auto bsp_func_3_3 = IgFunction<3,0,3,1>::create(bsp_space_3_3, bsp_coeff_3_3);
    auto bsp_func_2_3 = IgFunction<2,0,3,1>::create(bsp_space_2_3, bsp_coeff_2_3);



    auto phys_space_1_1_1_0 =
        PhysicalSpace<1,1,1,0,Transformation::h_grad>::create(
            bsp_space_1_1,
            bsp_func_1_1->clone());

    auto phys_space_2_1_1_0 =
        PhysicalSpace<2,1,1,0,Transformation::h_grad>::create(
            bsp_space_2_1,
            bsp_func_2_2->clone());

    auto phys_space_3_1_1_0 =
        PhysicalSpace<3,1,1,0,Transformation::h_grad>::create(
            bsp_space_3_1,
            bsp_func_3_3->clone());

    auto phys_space_2_2_1_0 =
        PhysicalSpace<2,2,1,0,Transformation::h_grad>::create(
            bsp_space_2_2,
            bsp_func_2_2->clone());

    auto phys_space_3_3_1_0 =
        PhysicalSpace<3,3,1,0,Transformation::h_grad>::create(
            bsp_space_3_3,
            bsp_func_3_3->clone());

    auto phys_space_2_1_1_1 =
        PhysicalSpace<2,1,1,1,Transformation::h_grad>::create(
            bsp_space_2_1,
            bsp_func_2_3->clone());

    auto phys_space_2_3_1_1 =
        PhysicalSpace<2,3,1,1,Transformation::h_grad>::create(
            bsp_space_2_3,
            bsp_func_2_3->clone());

    auto phys_coeff_1_1_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(*phys_space_1_1_1_0, "active", comm));
    (*phys_coeff_1_1_1_0)[0] = 1.;


    auto phys_coeff_2_1_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(*phys_space_2_1_1_0, "active", comm));
    (*phys_coeff_2_1_1_0)[0] = 1.;


    auto phys_coeff_3_1_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(*phys_space_3_1_1_0, "active", comm));
    (*phys_coeff_3_1_1_0)[0] = 1.;


    auto phys_coeff_2_2_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(*phys_space_2_2_1_0, "active", comm));
    (*phys_coeff_2_2_1_0)[0] = 1.;


    auto phys_coeff_3_3_1_0 = EpetraTools::create_vector(
                                  EpetraTools::create_map(*phys_space_3_3_1_0, "active", comm));
    (*phys_coeff_3_3_1_0)[0] = 1.;


    auto phys_coeff_2_1_1_1 = EpetraTools::create_vector(
                                  EpetraTools::create_map(*phys_space_2_1_1_1, "active", comm));
    (*phys_coeff_2_1_1_1)[0] = 1.;


    auto phys_coeff_2_3_1_1 = EpetraTools::create_vector(
                                  EpetraTools::create_map(*phys_space_2_3_1_1, "active", comm));
    (*phys_coeff_2_3_1_1)[0] = 1.;

    auto phys_func_1_1_1_0 = IgFunction<1,0,1,1>::create(phys_space_1_1_1_0,phys_coeff_1_1_1_0);
    auto phys_func_2_1_1_0 = IgFunction<2,0,1,1>::create(phys_space_2_1_1_0,phys_coeff_2_1_1_0);
    auto phys_func_3_1_1_0 = IgFunction<3,0,1,1>::create(phys_space_3_1_1_0,phys_coeff_3_1_1_0);
    auto phys_func_2_2_1_0 = IgFunction<2,0,2,1>::create(phys_space_2_2_1_0,phys_coeff_2_2_1_0);
    auto phys_func_3_3_1_0 = IgFunction<3,0,3,1>::create(phys_space_3_3_1_0,phys_coeff_3_3_1_0);
    auto phys_func_2_1_1_1 = IgFunction<2,1,1,1>::create(phys_space_2_1_1_1,phys_coeff_2_1_1_1);
    auto phys_func_2_3_1_1 = IgFunction<2,1,3,1>::create(phys_space_2_3_1_1,phys_coeff_2_3_1_1);


    auto funcs_container = std::make_shared<FunctionsContainer>();

    funcs_container->insert_mapping(
        std::const_pointer_cast<MapFunction<1,1>>(
            phys_func_1_1_1_0->get_ig_space()->get_ptr_const_map_func()),
        "map_1_1_1_0");

    funcs_container->insert_mapping(
        std::const_pointer_cast<MapFunction<2,2>>(
            phys_func_2_1_1_0->get_ig_space()->get_ptr_const_map_func()),
        "map_2_1_1_0");

    funcs_container->insert_mapping(
        std::const_pointer_cast<MapFunction<3,3>>(
            phys_func_3_1_1_0->get_ig_space()->get_ptr_const_map_func()),
        "map_3_1_1_0");

    funcs_container->insert_mapping(
        std::const_pointer_cast<MapFunction<2,2>>(
            phys_func_2_2_1_0->get_ig_space()->get_ptr_const_map_func()),
        "map_2_2_1_0");

    funcs_container->insert_mapping(
        std::const_pointer_cast<MapFunction<3,3>>(
            phys_func_3_3_1_0->get_ig_space()->get_ptr_const_map_func()),
        "map_3_3_1_0");

    funcs_container->insert_mapping(
        std::const_pointer_cast<MapFunction<2,3>>(
            phys_func_2_1_1_1->get_ig_space()->get_ptr_const_map_func()),
        "map_2_1_1_1");

    funcs_container->insert_mapping(
        std::const_pointer_cast<MapFunction<2,3>>(
            phys_func_2_3_1_1->get_ig_space()->get_ptr_const_map_func()),
        "map_2_3_1_1");

    funcs_container->insert_mapping(func_identity_1_1,"map_identity_1_1");
    funcs_container->insert_mapping(func_identity_2_2,"map_identity_2_2");
    funcs_container->insert_mapping(func_identity_3_3,"map_identity_3_3");

    funcs_container->insert_function(
        std::const_pointer_cast<MapFunction<1,1>>(
            phys_func_1_1_1_0->get_ig_space()->get_ptr_const_map_func()),
        static_pointer_cast<Func<1,0,1>>(phys_func_1_1_1_0),
        "phys_func_1_1_1_0");

    funcs_container->insert_function(
        std::const_pointer_cast<MapFunction<2,2>>(
            phys_func_2_1_1_0->get_ig_space()->get_ptr_const_map_func()),
        static_pointer_cast<Func<2,0,1>>(phys_func_2_1_1_0),
        "phys_func_2_1_1_0");

    funcs_container->insert_function(
        std::const_pointer_cast<MapFunction<3,3>>(
            phys_func_3_1_1_0->get_ig_space()->get_ptr_const_map_func()),
        static_pointer_cast<Func<3,0,1>>(phys_func_3_1_1_0),
        "phys_func_3_1_1_0");

    funcs_container->insert_function(
        std::const_pointer_cast<MapFunction<2,2>>(
            phys_func_2_2_1_0->get_ig_space()->get_ptr_const_map_func()),
        static_pointer_cast<Func<2,0,2>>(phys_func_2_2_1_0),
        "phys_func_2_2_1_0");

    funcs_container->insert_function(
        std::const_pointer_cast<MapFunction<3,3>>(
            phys_func_3_3_1_0->get_ig_space()->get_ptr_const_map_func()),
        static_pointer_cast<Func<3,0,3>>(phys_func_3_3_1_0),
        "phys_func_3_3_1_0");

    funcs_container->insert_function(
        std::const_pointer_cast<MapFunction<2,3>>(
            phys_func_2_1_1_1->get_ig_space()->get_ptr_const_map_func()),
        static_pointer_cast<Func<2,1,1>>(phys_func_2_1_1_1),
        "phys_func_2_1_1_1");

    funcs_container->insert_function(
        std::const_pointer_cast<MapFunction<2,3>>(
            phys_func_2_3_1_1->get_ig_space()->get_ptr_const_map_func()),
        static_pointer_cast<Func<2,1,3>>(phys_func_2_3_1_1),
        "phys_func_2_3_1_1");

    serialize_deserialize(funcs_container);
    print_container(funcs_container);

    deserialize_only();

}



int main()
{
    do_test();

    return 0;
}
