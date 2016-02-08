//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

/**
 *  Common code to use in the phys_basis_iterators_* tests
 *
 *   author: pauletti
 *   date: 2015-04-17
 *
 */

#include "../tests.h"

#include <igatools/base/quadrature_lib.h>
//#include <igatools/functions/function.h>
#include <igatools/functions/grid_function_lib.h>

#include <igatools/basis_functions/bspline.h>

#include <igatools/basis_functions/physical_basis.h>
#include <igatools/basis_functions/physical_basis_element.h>
#include <igatools/basis_functions/physical_basis_handler.h>

#include <igatools/basis_functions/space_tools.h>

using space_tools::get_boundary_dofs;

template <int dim, int range=1, int rank=1, int codim = 0>
shared_ptr<PhysicalBasis<dim,range,rank,codim> >
create_phys_basis(const shared_ptr<Grid<dim>> &grid,
                  const shared_ptr<GridFunction<dim,dim+codim>> &grid_func,
                  const int deg=1)
{
  using BSpBasis = BSpline<dim, range, rank>;
  using PhysBasis = PhysicalBasis<dim,range,rank,codim>;
  auto bsp_basis = BSpBasis::create(SplineSpace<dim,range,rank>::create(deg,grid));
  return PhysBasis::create(bsp_basis, Domain<dim,codim>::create(grid_func), Transformation::h_grad);
}


struct DofProp
{
  static const std::string interior;
  static const std::string dirichlet;
  static const std::string neumman;
};

const std::string DofProp::interior = "interior";
const std::string DofProp::dirichlet  = "dirichlet";
const std::string DofProp::neumman  = "neumman";

enum  bc : boundary_id
{
  dir=0, neu
};

template <int dim, int range=1, int rank=1, int codim = 0>
shared_ptr<PhysicalBasis<dim,range,rank,codim> >
create_basis_prop(const shared_ptr<Grid<dim>> &grid,
                  const shared_ptr<GridFunction<dim,dim+codim>> &grid_func,
                  const int deg=1)
{

  using BSpBasis = BSpline<dim, range, rank>;
  using PhysBasis = PhysicalBasis<dim,range,rank,codim>;
  auto space = SplineSpace<dim,range,rank>::create(deg,grid);
  auto bsp_basis = BSpBasis::create(space);
  auto phys_basis = PhysBasis::create(bsp_basis, Domain<dim,codim>::create(grid_func));

//  std::set<boundary_id>  dir_ids = {bc::dir};
  std::set<int> dir_faces;
  for (int face_id = 1 ; face_id < UnitElement<dim>::n_faces ; ++face_id)
    dir_faces.insert(face_id);
  auto dir_dofs = get_boundary_dofs<PhysBasis>(phys_basis, dir_faces);


  const std::set<int> neu_faces{0};
  auto neu_dofs = get_boundary_dofs<PhysBasis>(phys_basis, neu_faces);


  auto dof_dist = space->get_dof_distribution();
  auto int_dofs = dof_dist->get_interior_dofs();


  SafeSTLVector<Index> common(dim*range);
  auto end1 =
    std::set_intersection(neu_dofs.begin(), neu_dofs.end(),
                          dir_dofs.begin(), dir_dofs.end(), common.begin());
  common.resize(end1-common.begin());
  for (auto &id : common)
    neu_dofs.erase(id);

  dof_dist->add_dofs_property(DofProp::interior);
  dof_dist->add_dofs_property(DofProp::dirichlet);
  dof_dist->add_dofs_property(DofProp::neumman);


  dof_dist->set_dof_property_status(DofProp::interior, int_dofs, true);
  dof_dist->set_dof_property_status(DofProp::dirichlet, dir_dofs, true);
  dof_dist->set_dof_property_status(DofProp::neumman, neu_dofs, true);

  return phys_basis;
}



template <int dim, int k, int range=1, int rank=1, int codim = 0>
void elem_values(shared_ptr<PhysicalBasis<dim,range,rank,codim>> phys_basis,
                 const int n_qp = 1,
                 const string &prop = DofProperties::active,
                 const bool no_boundary=true)
{

  using std::to_string;
  out << "elem_values<" << dim << "," << k << "," << range << "," << rank << "," << codim << ">("
      << n_qp << "," << prop <<  "," << no_boundary << ")" << std::endl;


//    using PhysBasis = PhysicalBasis<dim,range,rank,codim>;
//    using Handler = typename PhysBasis::Handler;

  auto quad = QGauss<k>::create(n_qp);
  using basis_element::Flags;
  auto flag = Flags::value|
              Flags::gradient |
              Flags::hessian |
              Flags::w_measure;
  //|
  //      Flags::point |
  //    ValueFlags::w_measure;

  auto elem_filler = phys_basis->create_cache_handler();
  elem_filler->template set_flags<k>(flag);

  auto elem = phys_basis->begin();
  auto end  = phys_basis->end();
  elem_filler->template init_cache<k>(elem,quad);

  using basis_element::_Value;
  using basis_element::_Gradient;
  using basis_element::_Hessian;

  int elem_id = 0;
  for (; elem != end; ++elem)
  {
    const auto &grid_elem = elem->get_grid_element();
    if ((no_boundary) || (grid_elem.is_boundary()))
    {
      out.begin_item("Element " + std::to_string(elem_id));
      out << "Element index: " << grid_elem.get_index() << endl;
      for (auto &s_id : UnitElement<dim>::template elems_ids<k>())
      {
        if ((no_boundary) || (grid_elem.is_boundary(s_id)))
        {
          out.begin_item("Sub element id:  " + std::to_string(s_id));
          elem_filler->template fill_cache<k>(*elem, s_id);

          out.begin_item("Values: ");
          elem->template get_basis_data<_Value, k>(s_id,prop).print_info(out);
          out.end_item();

          out.begin_item("Gradients: ");
          elem->template get_basis_data<_Gradient, k>(s_id,prop).print_info(out);
          out.end_item();

          out.begin_item("Hessians: ");
          elem->template get_basis_data<_Hessian, k>(s_id,prop).print_info(out);
          out.end_item();

          out.begin_item("W * Measures: ");
          elem->template get_w_measures<k>(s_id).print_info(out);
          out.end_item();
          out.end_item();
        }
      }
      out.end_item();
    }
    ++elem_id;
  }


}
