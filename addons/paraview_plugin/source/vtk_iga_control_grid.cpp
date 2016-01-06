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

#include "../include/paraview_plugin/vtk_iga_control_grid.h"

#include <igatools/functions/ig_grid_function.h>
#include <igatools/geometry/domain.h>
#include <igatools/basis_functions/dof_distribution.h>

#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyVertex.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include "../include/paraview_plugin/vtk_iga_grid_information.h"


IGA_NAMESPACE_OPEN

template <class Domain>
auto
VtkIgaControlGrid<Domain>::
create(const DomainPtr_ domain,
       const ControlGridInfoPtr_ grid_info) -> VtkGridPtr_
{
  Assert(domain != nullptr, ExcNullPtr());
  Assert(grid_info != nullptr, ExcNullPtr());

  const auto ig_grid_fun = std::dynamic_pointer_cast<const IgGridFun_>
      (domain->get_grid_function());
  Assert(ig_grid_fun != nullptr, ExcNullPtr());

  const auto &coefs = ig_grid_fun->get_coefficients();

  const auto &dofs = ig_grid_fun->get_basis()->get_dof_distribution();
  const auto &dofs_table = dofs->get_num_dofs_table();

  const Size n_total_pts = dofs_table.get_component_size(0);

  auto points = vtkSmartPointer <vtkPoints>::New();
  points->SetNumberOfPoints(n_total_pts);

  // Filling the control point coordinates by retrieving the
  // values of the control point coefficients.

  // The points are initialized to zero.
  const double point_void[3] =
  { 0.0, 0.0, 0.0 };
  for (int i_pt = 0; i_pt < n_total_pts; ++i_pt)
    points->SetPoint(i_pt, point_void);

  Index comp;
  Index local_id;
  double point_tmp[3];
  for (const auto &it : coefs)
  {
    dofs->global_to_comp_local(it.first, comp, local_id);
    points->GetPoint(local_id, point_tmp);
    point_tmp[comp] = it.second;
    points->SetPoint(local_id, point_tmp);
  }

  if (!grid_info->is_structured() || dim == 1)  // Unstructured grid
    return Self_::create_grid_vtu(ig_grid_fun, points);
  else // Structured grid
    return Self_::create_grid_vts(ig_grid_fun, points);
}



template <class Domain>
auto
VtkIgaControlGrid<Domain>::
create_grid_vts(const IgGridFunPtr_ ig_grid_fun,
                const vtkSmartPointer<vtkPoints> points) -> VtkGridPtr_
{
  const auto &dofs = ig_grid_fun->get_basis()->get_dof_distribution();
  const auto &dofs_table = dofs->get_num_dofs_table();
  const auto n_pts_dir = dofs_table[0];

  auto grid_dim = TensorSize <3> (1);
  for (int dir = 0; dir < dim; ++dir)
    grid_dim[dir] = n_pts_dir[dir];

  auto grid = vtkSmartPointer <vtkStructuredGrid>::New();
  grid->SetDimensions(grid_dim[0], grid_dim[1], grid_dim[2]);
  grid->SetPoints(points);

  return grid;
}



template <class Domain>
auto
VtkIgaControlGrid<Domain>::
create_grid_vtu(const IgGridFunPtr_ ig_grid_fun,
                const vtkSmartPointer<vtkPoints> points) -> VtkGridPtr_
{
  const auto &dofs = ig_grid_fun->get_basis()->get_dof_distribution();
  const auto &dofs_table = dofs->get_num_dofs_table();
  const auto n_pts_dir = dofs_table[0];

  auto grid = vtkSmartPointer <vtkUnstructuredGrid>::New();
  grid->SetPoints(points);

  // Creating poly vertices.
  vtkSmartPointer <vtkPolyVertex> vertices = vtkSmartPointer <
                                             vtkPolyVertex>::New();

  const Size n_points = points->GetNumberOfPoints();
  auto vertex_points = vertices->GetPointIds();
  vertex_points->SetNumberOfIds(n_points);
  for (int i_pt = 0; i_pt < n_points; ++i_pt)
    vertex_points->SetId(i_pt, i_pt);

  // Create a cell array to store the lines and vertices,
  vtkSmartPointer <vtkCellArray> cells = vtkSmartPointer <
                                         vtkCellArray>::New();
  cells->InsertNextCell(vertices);

  const TensorSizedContainer <dim> ids(n_pts_dir);
  TensorIndex <dim> tid;

  UnitElement <dim> elem;

  for (const auto &dir : elem.active_directions)
  {
    // Building the point id container for the face.

    const Index face_id = dir * 2;

    // Creating the cartesian product container for the face.
    TensorSize <dim - 1> n_pts_face_dir;

    const auto face_elem = elem.template get_elem <dim - 1> (face_id);

    Index i = 0;
    for (const auto &face_dir : face_elem.active_directions)
      n_pts_face_dir[i++] = n_pts_dir[face_dir];

    const TensorSizedContainer <dim - 1> ids_face(
      n_pts_face_dir);

    const Size &n_pts_line = n_pts_dir[dir];

    // Creating the offset of the line points.
    Index offset = 1;
    for (int dir2 = 0; dir2 < dir; ++dir2)
      offset *= n_pts_dir[dir2];

    tid[dir] = 0;
    for (int l_id = 0; l_id < ids_face.flat_size(); ++l_id)
    {
      // Getting the flat index of the initial point of the line.
      const auto fid = ids_face.flat_to_tensor(l_id);
      auto fid_it = fid.cbegin();
      for (const auto &face_dir : face_elem.active_directions)
        tid[face_dir] = *fid_it++;
      Index pt_id = ids.tensor_to_flat(tid);

      vtkSmartPointer <vtkPolyLine> line = vtkSmartPointer <
                                           vtkPolyLine>::New();
      auto line_points = line->GetPointIds();
      line_points->SetNumberOfIds(n_pts_line);

      for (int i_pt = 0; i_pt < n_pts_line; ++i_pt, pt_id += offset)
        line_points->SetId(i_pt, pt_id);

      cells->InsertNextCell(line);

    }
  }

  const Size n_cells = cells->GetNumberOfCells();

  int cell_types[n_cells];
  cell_types[0] = VTK_POLY_VERTEX;

  for (int i = 1; i < n_cells; ++i)
    cell_types[i] = VTK_POLY_LINE;

  grid->Allocate(n_cells, 0);
  grid->SetCells(cell_types, cells);

  return grid;
}

template class VtkIgaControlGrid<Domain<1, 0>>;
template class VtkIgaControlGrid<Domain<1, 1>>;
template class VtkIgaControlGrid<Domain<1, 2>>;
template class VtkIgaControlGrid<Domain<2, 0>>;
template class VtkIgaControlGrid<Domain<2, 1>>;
template class VtkIgaControlGrid<Domain<3, 0>>;

IGA_NAMESPACE_CLOSE
