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

#include <paraview_plugin/grid_information.h>

#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include <vtkPolyVertex.h>

#include <boost/range/irange.hpp>

#include <igatools/functions/ig_function.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <paraview_plugin/control_grid_generator.h>


IGA_NAMESPACE_OPEN

template <int dim, int codim>
VtkIgaControlGridGenerator<dim, codim>::
VtkIgaControlGridGenerator (const MapFunPtr_ mapping,
                         const ControlGridInfoPtr_ grid_info)
    :
                         ig_fun_ (std::dynamic_pointer_cast <IgFun_> (mapping)),
                         grid_info_ (grid_info)
{
    Assert (mapping != nullptr, ExcNullPtr());
    Assert (ig_fun_ != nullptr, ExcNullPtr());
    Assert (grid_info != nullptr, ExcNullPtr());

    // TODO: include assert here to verify that all the components share the
    // same space.
}
template <int dim, int codim>
auto
VtkIgaControlGridGenerator<dim, codim>::
get_grid (const MapFunPtr_ mapping,
             const ControlGridInfoPtr_ grid_info) -> VtkGridPtr_
{
    VtkIgaControlGridGenerator generator (mapping, grid_info);
    return generator.create_grid();
}



template <int dim, int codim>
auto
VtkIgaControlGridGenerator<dim, codim>::
create_grid () const -> VtkGridPtr_
{
    const auto space = ig_fun_->get_ig_space ();
    const auto &coefs = ig_fun_->get_coefficients ();
    const auto &dofs = space->get_ptr_const_dof_distribution ();

    const auto &dofs_table = dofs->get_num_dofs_table ();

    const Size n_total_pts = dofs_table.get_component_size (0);

    auto points = vtkSmartPointer <vtkPoints>::New ();
    points->SetNumberOfPoints (n_total_pts);

    // Filling the control point coordinates by retrieving the
    // values of the control point coefficients.

    // The points are initialized to zero.
    const double point_void[3] =
            { 0.0, 0.0, 0.0 };
    for (int i_pt = 0; i_pt < n_total_pts; ++i_pt)
        points->SetPoint (i_pt, point_void);

    Index comp;
    Index local_id;
    double point_tmp[3];
    for (const auto &it : coefs)
    {
        dofs->global_to_comp_local (it.first, comp, local_id);
        points->GetPoint (local_id, point_tmp);
        point_tmp[comp] = it.second;
        points->SetPoint (local_id, point_tmp);
    }

    if (!grid_info_->is_structured () || dim == 1) // Unstructured grid
        return this->create_grid_vtu(points);
    else // Structured grid
        return this->create_grid_vts(points);
}






template <int dim, int codim>
auto
VtkIgaControlGridGenerator<dim, codim>::
create_grid_vts(const vtkSmartPointer<vtkPoints> points) const -> VtkGridPtr_
{
    const auto space = ig_fun_->get_ig_space ();
    const auto &dofs = space->get_ptr_const_dof_distribution ();

    const auto &dofs_table = dofs->get_num_dofs_table ();
    const auto n_pts_dir = dofs_table[0];

    auto grid_dim = TensorSize <3> (1);
    for (const auto &dir : boost::irange (0, dim))
        grid_dim[dir] = n_pts_dir[dir];

    auto grid = vtkSmartPointer <vtkStructuredGrid>::New ();
    grid->SetDimensions (grid_dim[0], grid_dim[1], grid_dim[2]);
    grid->SetPoints (points);

    return grid;
}



template <int dim, int codim>
auto
VtkIgaControlGridGenerator<dim, codim>::
create_grid_vtu(const vtkSmartPointer<vtkPoints> points) const -> VtkGridPtr_
{
    const auto space = ig_fun_->get_ig_space ();
    const auto &dofs = space->get_ptr_const_dof_distribution ();

    const auto &dofs_table = dofs->get_num_dofs_table ();
    const auto n_pts_dir = dofs_table[0];

    auto grid = vtkSmartPointer <vtkUnstructuredGrid>::New ();
    grid->SetPoints (points);

    // Creating poly vertices.
    vtkSmartPointer <vtkPolyVertex> vertices = vtkSmartPointer <
            vtkPolyVertex>::New ();

    const Size n_points = points->GetNumberOfPoints ();
    auto vertex_points = vertices->GetPointIds ();
    vertex_points->SetNumberOfIds (n_points);
    for (const auto &i_pt : boost::irange (0, n_points))
        vertex_points->SetId (i_pt, i_pt);

    // Create a cell array to store the lines and vertices,
    vtkSmartPointer <vtkCellArray> cells = vtkSmartPointer <
            vtkCellArray>::New ();
    cells->InsertNextCell (vertices);

    const TensorSizedContainer <dim> ids (n_pts_dir);
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

        const TensorSizedContainer <dim - 1> ids_face (
                n_pts_face_dir);

        const Size &n_pts_line = n_pts_dir[dir];

        // Creating the offset of the line points.
        Index offset = 1;
        for (int dir2 = 0; dir2 < dir; ++dir2)
            offset *= n_pts_dir[dir2];

        tid[dir] = 0;
        for (const auto &l_id : boost::irange (0, ids_face.flat_size ()))
        {
            // Getting the flat index of the initial point of the line.
            const auto fid = ids_face.flat_to_tensor (l_id);
            auto fid_it = fid.cbegin ();
            for (const auto &face_dir : face_elem.active_directions)
                tid[face_dir] = *fid_it++;
            Index pt_id = ids.tensor_to_flat (tid);

            vtkSmartPointer <vtkPolyLine> line = vtkSmartPointer <
                    vtkPolyLine>::New ();
            auto line_points = line->GetPointIds ();
            line_points->SetNumberOfIds (n_pts_line);

            for (int i_pt = 0; i_pt < n_pts_line; ++i_pt, pt_id += offset)
                line_points->SetId (i_pt, pt_id);

            cells->InsertNextCell (line);

        }
    }

    const Size n_cells = cells->GetNumberOfCells ();

    int cell_types[n_cells];
    cell_types[0] = VTK_POLY_VERTEX;

    for (int i = 1; i < n_cells; ++i)
        cell_types[i] = VTK_POLY_LINE;

    grid->Allocate (n_cells, 0);
    grid->SetCells (cell_types, cells);

    return grid;
}

template class VtkIgaControlGridGenerator<1, 0>;
template class VtkIgaControlGridGenerator<1, 1>;
template class VtkIgaControlGridGenerator<1, 2>;
template class VtkIgaControlGridGenerator<2, 0>;
template class VtkIgaControlGridGenerator<2, 1>;
template class VtkIgaControlGridGenerator<3, 0>;



IGA_NAMESPACE_CLOSE
