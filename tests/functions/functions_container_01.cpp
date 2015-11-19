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

#ifdef SERIALIZATION
void serialize_deserialize(std::shared_ptr<FunctionsContainer> funcs_container)
{
  OUTSTART

  out.begin_item("Original FunctionsContainer:");
  funcs_container->print_info(out);
  out.end_item();

  std::string filename = "functions_container.xml";
  std::string tag_name = "FunctionsContainer";
  {
    // serialize the PhysicalSpaceBasis object to an xml file
    std::ofstream xml_ostream(filename);
    OArchive xml_out(xml_ostream);

    xml_out << funcs_container;
  }

  funcs_container.reset();
  {
    // de-serialize the PhysicalSpaceBasis object from an xml file
    std::ifstream xml_istream(filename);
    IArchive xml_in(xml_istream);

    xml_in >> funcs_container;
  }
  out.begin_item("FunctionsContainer after serialize-deserialize:");
  funcs_container->print_info(out);
  out.end_item();
//*/

  OUTEND
}



void deserialize_only()
{
  std::string filename = "functions_container.xml";

  std::ifstream xml_istream(filename);
  IArchive xml_in(xml_istream);

  auto funcs_container = std::make_shared<FunctionsContainer>();
  xml_in >> funcs_container;

  out.begin_item("Inside deserialize_only()");
  print_container(funcs_container);
  out.end_item();
}
#endif // SERIALIZATION


template <int dim,int codim,int range>
using Func = Function<dim,codim,range,1>;

void print_container(const std::shared_ptr<FunctionsContainer> &funcs_container)
{
  const auto domains_dim_2_codim_0 = funcs_container->template get_domains_dim_codim<2,0>();
  out.begin_item("Mappings with dimension 2 and codimension 0:");
  int m_counter = 0;
  for (const auto &domain : domains_dim_2_codim_0)
  {
#if 0
    auto flag = ValueFlags::point | ValueFlags::value;
    QUniform<2> quad(2);
    domain->reset(flag, quad);

    out.begin_item("Mapping before iterators");
    domain->print_info(out);
    out.end_item();

    auto elem = domain->begin();
    auto end  = domain->end();
#endif

#if 0
    try
    {
      const shared_ptr<IgFunction<2, 0, 2, 1>> ig_fun_ptr = std::dynamic_pointer_cast<IgFunction<2, 0, 2, 1>> (domain);
      if (ig_fun_ptr == nullptr)
      {
        out << "Identitity func start" << endl;
        const shared_ptr<IdentityFunction<2, 2>> id_fun_ptr = std::dynamic_pointer_cast<IdentityFunction<2, 2>> (domain);
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
    out << "Mapping[" << m_counter++ << "]   ID= " << domain->get_object_id()
        << "   name= " << domain->get_name() << std::endl;

    const auto &funcs_dim_2_codim_0_range_1_rank_1 =
      funcs_container->template get_functions_with_same_domain<2,0,1,1>(*domain);
    out.begin_item("Functions<2,0,1,1>:");
    int f_counter = 0;
    for (const auto &f : funcs_dim_2_codim_0_range_1_rank_1)
    {
      out << "Function[" << f_counter++ << "]"
          << "   ID:" << f.first
          << "   name= " << f.second->get_name() << std::endl;
    }
    out.end_item();

  }
}


void do_test()
{
  int n_elem_per_side = 2;
  auto grid_1 = Grid<1>::const_create(n_elem_per_side+1);
  auto grid_2 = Grid<2>::const_create(n_elem_per_side+1);
  auto grid_3 = Grid<3>::const_create(n_elem_per_side+1);
//    create_fun<2, 0, 2>();

#if 0
  auto func_identity_1_1 = IdentityFunction<1,1>::const_create(grid_1);
  auto func_identity_2_2 = IdentityFunction<2,2>::const_create(grid_2);
  auto func_identity_3_3 = IdentityFunction<3,3>::const_create(grid_3);
#endif


  const int deg = 3;
  auto bsp_space_1_1 = BSpline<1,1,1>::const_create(SplineSpace<1,1,1>::const_create(deg,grid_1));
  auto bsp_space_2_1 = BSpline<2,1,1>::const_create(SplineSpace<2,1,1>::const_create(deg,grid_2));
  auto bsp_space_3_1 = BSpline<3,1,1>::const_create(SplineSpace<3,1,1>::const_create(deg,grid_3));
  auto bsp_space_2_2 = BSpline<2,2,1>::const_create(SplineSpace<2,2,1>::const_create(deg,grid_2));
  auto bsp_space_3_3 = BSpline<3,3,1>::const_create(SplineSpace<3,3,1>::const_create(deg,grid_3));
  auto bsp_space_2_3 = BSpline<2,3,1>::const_create(SplineSpace<2,3,1>::const_create(deg,grid_2));

  Epetra_SerialComm comm;
  auto bsp_coeff_1_1 = EpetraTools::create_vector(*bsp_space_1_1,DofProperties::active,comm);
  (*bsp_coeff_1_1)[0] = 1.;


  auto bsp_coeff_2_1 = EpetraTools::create_vector(*bsp_space_2_1,DofProperties::active,comm);
  (*bsp_coeff_2_1)[0] = 1.;


  auto bsp_coeff_3_1 = EpetraTools::create_vector(*bsp_space_3_1,DofProperties::active,comm);
  (*bsp_coeff_3_1)[0] = 1.;


  auto bsp_coeff_2_2 = EpetraTools::create_vector(*bsp_space_2_2,DofProperties::active,comm);
  (*bsp_coeff_2_2)[0] = 1.;


  auto bsp_coeff_3_3 = EpetraTools::create_vector(*bsp_space_3_3,DofProperties::active,comm);
  (*bsp_coeff_3_3)[0] = 1.;


  auto bsp_coeff_2_3 = EpetraTools::create_vector(*bsp_space_2_3,DofProperties::active,comm);
  (*bsp_coeff_2_3)[0] = 1.;


  auto bsp_func_1_1 = IgGridFunction<1,1>::const_create(bsp_space_1_1, *bsp_coeff_1_1);
  auto bsp_func_2_1 = IgGridFunction<2,1>::const_create(bsp_space_2_1, *bsp_coeff_2_1);
  auto bsp_func_3_1 = IgGridFunction<3,1>::const_create(bsp_space_3_1, *bsp_coeff_3_1);
  auto bsp_func_2_2 = IgGridFunction<2,2>::const_create(bsp_space_2_2, *bsp_coeff_2_2);
  auto bsp_func_3_3 = IgGridFunction<3,3>::const_create(bsp_space_3_3, *bsp_coeff_3_3);
  auto bsp_func_2_3 = IgGridFunction<2,3>::const_create(bsp_space_2_3, *bsp_coeff_2_3);



  auto phys_space_1_1_1_0 =
    PhysicalSpaceBasis<1,1,1,0>::const_create(
      bsp_space_1_1,
      Domain<1,0>::const_create(bsp_func_1_1,"map_1_1_1_0"));

  auto phys_space_2_1_1_0 =
    PhysicalSpaceBasis<2,1,1,0>::const_create(
      bsp_space_2_1,
      Domain<2,0>::const_create(bsp_func_2_2,"map_2_1_1_0"));

  auto phys_space_3_1_1_0 =
    PhysicalSpaceBasis<3,1,1,0>::const_create(
      bsp_space_3_1,
      Domain<3,0>::const_create(bsp_func_3_3));

  auto phys_space_2_2_1_0 =
    PhysicalSpaceBasis<2,2,1,0>::const_create(
      bsp_space_2_2,
      Domain<2,0>::const_create(bsp_func_2_2,"map_2_2_1_0"));

  auto phys_space_3_3_1_0 =
    PhysicalSpaceBasis<3,3,1,0>::const_create(
      bsp_space_3_3,
      Domain<3,0>::const_create(bsp_func_3_3,"map_3_3_1_0"));

  auto phys_space_2_1_1_1 =
    PhysicalSpaceBasis<2,1,1,1>::const_create(
      bsp_space_2_1,
      Domain<2,1>::const_create(bsp_func_2_3,"map_2_1_1_1"));

  auto phys_space_2_3_1_1 =
    PhysicalSpaceBasis<2,3,1,1>::const_create(
      bsp_space_2_3,
      Domain<2,1>::const_create(bsp_func_2_3,"map_2_3_1_1"));

  auto phys_coeff_1_1_1_0 = EpetraTools::create_vector(*phys_space_1_1_1_0,DofProperties::active,comm);
  (*phys_coeff_1_1_1_0)[0] = 1.;


  auto phys_coeff_2_1_1_0 = EpetraTools::create_vector(*phys_space_2_1_1_0,DofProperties::active,comm);
  (*phys_coeff_2_1_1_0)[0] = 1.;


  auto phys_coeff_3_1_1_0 = EpetraTools::create_vector(*phys_space_3_1_1_0,DofProperties::active,comm);
  (*phys_coeff_3_1_1_0)[0] = 1.;


  auto phys_coeff_2_2_1_0 = EpetraTools::create_vector(*phys_space_2_2_1_0,DofProperties::active,comm);
  (*phys_coeff_2_2_1_0)[0] = 1.;


  auto phys_coeff_3_3_1_0 = EpetraTools::create_vector(*phys_space_3_3_1_0,DofProperties::active,comm);
  (*phys_coeff_3_3_1_0)[0] = 1.;


  auto phys_coeff_2_1_1_1 = EpetraTools::create_vector(*phys_space_2_1_1_1,DofProperties::active,comm);
  (*phys_coeff_2_1_1_1)[0] = 1.;


  auto phys_coeff_2_3_1_1 = EpetraTools::create_vector(*phys_space_2_3_1_1,DofProperties::active,comm);
  (*phys_coeff_2_3_1_1)[0] = 1.;

  auto phys_func_1_1_1_0 = IgFunction<1,0,1,1>::const_create(phys_space_1_1_1_0,*phys_coeff_1_1_1_0,DofProperties::active,"phys_func_1_1_1_0");
  auto phys_func_2_1_1_0 = IgFunction<2,0,1,1>::const_create(phys_space_2_1_1_0,*phys_coeff_2_1_1_0,DofProperties::active,"phys_func_2_1_1_0");
  auto phys_func_3_1_1_0 = IgFunction<3,0,1,1>::const_create(phys_space_3_1_1_0,*phys_coeff_3_1_1_0,DofProperties::active,"phys_func_3_1_1_0");
  auto phys_func_2_2_1_0 = IgFunction<2,0,2,1>::const_create(phys_space_2_2_1_0,*phys_coeff_2_2_1_0,DofProperties::active,"phys_func_2_2_1_0");
  auto phys_func_3_3_1_0 = IgFunction<3,0,3,1>::const_create(phys_space_3_3_1_0,*phys_coeff_3_3_1_0,DofProperties::active,"phys_func_3_3_1_0");
  auto phys_func_2_1_1_1 = IgFunction<2,1,1,1>::const_create(phys_space_2_1_1_1,*phys_coeff_2_1_1_1,DofProperties::active,"phys_func_1_1_1_1");
  auto phys_func_2_3_1_1 = IgFunction<2,1,3,1>::const_create(phys_space_2_3_1_1,*phys_coeff_2_3_1_1,DofProperties::active,"phys_func_2_3_1_1");


  auto funcs_container = std::make_shared<FunctionsContainer>();

  funcs_container->insert_domain(phys_func_1_1_1_0->get_domain());

  funcs_container->insert_domain(phys_func_2_1_1_0->get_domain());

  funcs_container->insert_domain(phys_func_3_1_1_0->get_domain());

  funcs_container->insert_domain(phys_func_2_2_1_0->get_domain());

  funcs_container->insert_domain(phys_func_3_3_1_0->get_domain());

  funcs_container->insert_domain(phys_func_2_1_1_1->get_domain());

  funcs_container->insert_domain(phys_func_2_3_1_1->get_domain());

#if 0
  funcs_container->insert_domain(func_identity_1_1,"map_identity_1_1");
  funcs_container->insert_domain(func_identity_2_2,"map_identity_2_2");
  funcs_container->insert_domain(func_identity_3_3,"map_identity_3_3");
#endif

  funcs_container->insert_function(
    phys_func_1_1_1_0);

  funcs_container->insert_function(
    phys_func_2_1_1_0);

  funcs_container->insert_function(
    phys_func_3_1_1_0);

  funcs_container->insert_function(
    phys_func_2_2_1_0);

  funcs_container->insert_function(
    phys_func_3_3_1_0);

  funcs_container->insert_function(
    phys_func_2_1_1_1);

  funcs_container->insert_function(
    phys_func_2_3_1_1);

  out.begin_item("Functions container:");
  funcs_container->print_info(out);
  out.end_item();


#ifdef SERIALIZATION
  serialize_deserialize(funcs_container);
//    funcs_container->print_info(out);
//    print_container(funcs_container);

  deserialize_only();
#endif // SERIALIZATION
}



int main()
{
  do_test();

  return 0;
}
