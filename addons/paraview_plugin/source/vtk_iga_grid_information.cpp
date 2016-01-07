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

#include <paraview_plugin/vtk_iga_grid_information.h>


IGA_NAMESPACE_OPEN

VtkGridInformation::
VtkGridInformation(const NumCellsContainer_ &num_cells,
                   const VtkGridType &grid_type)
  :
  grid_type_(grid_type),
  cells_per_element_(num_cells)
{
#ifndef NDEBUG
  for (const auto &it : num_cells)
    Assert(it > 0, ExcMessage("The number of cells per visualization "
                              "element must be > 0 in every direction."));
#endif
};



auto
VtkGridInformation::
create(const NumCellsContainer_ &num_cells,
       const VtkGridType &grid_type) -> SelfPtr_
{
  return SelfPtr_(new Self_(num_cells, grid_type));
}



bool
VtkGridInformation::
update(SelfPtr_ grid_info)
{
  if (grid_type_ != grid_info->get_grid_type())
  {
    grid_type_ = grid_info->get_grid_type();
    return true;
  }

  const auto &n_cells = grid_info->cells_per_element_;
  for (int dir = 0; dir < 3; ++dir)
  {
    if (n_cells[dir] != cells_per_element_[dir])
    {
      cells_per_element_ = n_cells;
      return true;
    }
  }

  return false;
}




const VtkGridType &
VtkGridInformation::
get_grid_type() const
{
  return grid_type_;
};



bool
VtkGridInformation::
is_structured() const
{
  return grid_type_ == VtkGridType::Structured;
};



bool
VtkGridInformation::
is_quadratic() const
{
  return grid_type_ == VtkGridType::UnstructuredQuadratic;
};



template <int dim>
SafeSTLArray <Size, dim>
VtkGridInformation::
get_num_cells_per_element() const
{
  SafeSTLArray <Size, dim> n_vis_elements;
  for (int dir = 0; dir < dim; ++dir)
    n_vis_elements[dir] = cells_per_element_[dir];
  return n_vis_elements;
}



void
VtkGridInformation::
print_info (LogStream &out) const
{
    out.begin_item("VtkGridInformation");
    out << "VTK grid type: ";
    switch (grid_type_)
    {
        case VtkGridType::Structured:
            out << "Structured Grid";
            break;
        case VtkGridType::UnstructuredLinear:
            out << "Unstructured Linear Grid";
            break;
        case VtkGridType::UnstructuredQuadratic:
            out << "Unstructured Quadratic Grid";
            break;
        default:
            break;
    }
    out << std::endl;
    out << "Number of VTK cells per Bezier element: " << cells_per_element_ << std::endl;
    out.end_item();
};



VtkControlGridInformation::
VtkControlGridInformation(const bool structured)
  :
  grid_type_(structured ? VtkGridType::Structured :
            VtkGridType::UnstructuredLinear)
{}



auto
VtkControlGridInformation::
create(const bool structured) -> SelfPtr_
{
  return SelfPtr_(new Self_(structured));
}



bool
VtkControlGridInformation::
update(SelfPtr_ grid_info)
{
  if (grid_type_ != grid_info->get_grid_type())
  {
    grid_type_ = grid_info->get_grid_type();
    return true;
  }

  return false;
}



bool
VtkControlGridInformation::
is_structured() const
{
  return grid_type_ == VtkGridType::Structured;
}




const VtkGridType &
VtkControlGridInformation::
get_grid_type() const
{
  return grid_type_;
};



void
VtkControlGridInformation::
print_info (LogStream &out) const
{
    out.begin_item("VtkControlGridInformation");
    out << "VTK grid type: ";
    switch (grid_type_)
    {
        case VtkGridType::Structured:
            out << "Structured Grid";
            break;
        case VtkGridType::UnstructuredLinear:
            out << "Unstructured Linear Grid";
            break;
        case VtkGridType::UnstructuredQuadratic:
            out << "Unstructured Quadratic Grid";
            break;
        default:
            break;
    }
    out << std::endl;
    out.end_item();
};

template SafeSTLArray <Size, 1> VtkGridInformation::get_num_cells_per_element<1>() const;
template SafeSTLArray <Size, 2> VtkGridInformation::get_num_cells_per_element<2>() const;
template SafeSTLArray <Size, 3> VtkGridInformation::get_num_cells_per_element<3>() const;

IGA_NAMESPACE_CLOSE
