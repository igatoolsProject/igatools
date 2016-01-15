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
 *
 *  Test for the evaluation of physical basis basis functions
 *  values and gradients with an ig function for the map
 *  author:
 *  date: 2014-11-08
 *
 */

#include "../tests.h"

#include <igatools/functions/ig_grid_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/nurbs.h>
#include <igatools/basis_functions/physical_basis.h>
#include <igatools/basis_functions/physical_basis_element.h>
#include <igatools/basis_functions/physical_basis_handler.h>



template <class Basis>
void serialize_deserialize(std::shared_ptr<Basis> basis)
{
  OUTSTART

  out.begin_item("Original PhysicalBasis:");
  basis->print_info(out);
  out.end_item();


  std::string template_args =
    "_dim" + std::to_string(Basis::dim) +
    "_range" + std::to_string(Basis::range) +
    "_rank" + std::to_string(Basis::rank) +
    "_codim" + std::to_string(Basis::codim);
  std::string filename = "phys_basis" + template_args + ".xml";
  std::string tag_name = "PhysicalBasis" + template_args;
  {
    // serialize the PhysicalBasis object to an xml file
    std::ofstream xml_ostream(filename);
    OArchive xml_out(xml_ostream);

    xml_out << boost::serialization::make_nvp(tag_name.c_str(),basis);
    xml_ostream.close();
  }

  basis.reset();
  {
    // de-serialize the PhysicalBasis object from an xml file
    std::ifstream xml_istream(filename);
    IArchive xml_in(xml_istream);

    xml_in >> BOOST_SERIALIZATION_NVP(basis);
    xml_istream.close();
  }
  out.begin_item("PhysicalBasis after serialize-deserialize:");
  basis->print_info(out);
  out.end_item();
//*/

  OUTEND
}


template<int dim, int codim=0>
auto
create_function(shared_ptr<const BSpline<dim, dim + codim>> basis)
{
  IgCoefficients control_pts;

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

  using Function = IgGridFunction<dim,dim+codim>;
  return Function::const_create(basis, control_pts);
}


template<int dim,int range=dim,int rank=1,int codim=0>
auto
create_phys_basis(shared_ptr<const BSpline<dim,range,rank>> ref_basis)
{
  using Basis = PhysicalBasis<dim,range,rank,codim>;

  return Basis::const_create(
           ref_basis,
           Domain<dim,codim>::const_create(create_function(ref_basis)),
           Transformation::h_grad);
}


template <int dim, int order = 0, int range=dim, int rank=1, int codim = 0>
void elem_values(const int n_knots = 2, const int deg=1)
{
  OUTSTART
  const int k = dim;
  using BspBasis = BSpline<dim, range, rank>;


  auto grid  = Grid<dim>::const_create(n_knots);

  auto ref_basis = BspBasis::const_create(
                     SplineSpace<dim,range,rank>::const_create(deg,grid));


  auto basis = create_phys_basis(ref_basis);

//  serialize_deserialize(basis);


  const int n_qp = 3;
  auto quad = QGauss<k>::create(n_qp);

  using basis_element::Flags;
  auto flag = Flags::value |
              Flags::gradient |
              Flags::hessian |
              Flags::divergence |
              Flags::w_measure;

  auto elem_handler = basis->create_cache_handler();
  elem_handler->template set_flags<k>(flag);

  using basis_element::_Value;
  using basis_element::_Gradient;
  using basis_element::_Hessian;
  using basis_element::_Divergence;

  auto elem = basis->begin();
  auto end = basis->end();
  elem_handler->init_element_cache(elem,quad);
  int elem_id = 0;
  for (; elem != end; ++elem)
  {
    elem_handler->fill_element_cache(elem);
    out.begin_item("Element " + std::to_string(elem_id));
    elem->print_info(out);
//    elem->print_cache_info(out);

    out.begin_item("Values: ");
    elem->template get_basis_data<_Value, k>(0,DofProperties::active).print_info(out);
    out.end_item();

    out.begin_item("Gradients: ");
    elem->template get_basis_data<_Gradient, k>(0,DofProperties::active).print_info(out);
    out.end_item();

    out.begin_item("Hessians: ");
    elem->template get_basis_data<_Hessian, k>(0,DofProperties::active).print_info(out);
    out.end_item();

    out.begin_item("Divergences: ");
    elem->template get_basis_data<_Divergence,k>(0,DofProperties::active).print_info(out);
    out.end_item();

    out.begin_item("W * Measures: ");
    elem->template get_w_measures<k>(0).print_info(out);
    out.end_item();

    out.end_item();

    ++elem_id;
  }


  OUTEND
}

int main()
{
  out.depth_console(10);

  elem_values<1>();
  elem_values<2>();
  elem_values<3>();

}
