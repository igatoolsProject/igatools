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

#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <paraview_plugin/vtk_iga_solid_grid.h>

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/domain_element.h>
#include <igatools/utils/multi_array_utils.h>
#include <igatools/functions/function_handler.h>
#include <igatools/functions/function_element.h>
#include <paraview_plugin/vtk_iga_grid_information.h>

using std::shared_ptr;
using namespace boost::fusion;

using std::remove_reference;

IGA_NAMESPACE_OPEN

namespace paraview_plugin
{

template <class Domain>
VtkIgaSolidGrid<Domain>::
PointsTopology::
PointsTopology(const std::shared_ptr<const Grid<dim>> cartesian_grid,
           const GridInfoPtr_ grid_info)
           :
           connectivity_(create_element_connectivity(grid_info)),
           quad_ (create_visualization_quadrature(grid_info)),
           n_vis_elements_ (grid_info->get_num_cells_per_element<dim>())
{
  Assert (cartesian_grid != nullptr, ExcNullPtr());
  Assert (grid_info != nullptr, ExcNullPtr());

  this->fill_points_map_mask(cartesian_grid, grid_info);
};



template <class Domain>
void
VtkIgaSolidGrid<Domain>::
PointsTopology::
fill_points_map_mask(const std::shared_ptr<const Grid<dim>> cartesian_grid,
                     const GridInfoPtr_ grid_info)
{
  Assert (cartesian_grid != nullptr, ExcNullPtr());
  Assert (grid_info != nullptr, ExcNullPtr());

  // Points map and mask
  map_.clear();
  mask_.clear();

  const Size n_bezier_elements = cartesian_grid->get_num_all_elems();

  if (n_bezier_elements == 0)
      VtkIgaWarningMacro("Grid \"" + cartesian_grid->get_name() +
                         "\" with ObjectId=\"" +
                         std::to_string(cartesian_grid->get_object_id()) +
                         "\" has 0 Bezier elements");

  map_.resize(n_bezier_elements);

  const int n_quad_points = quad_->get_num_points();

  if (!grid_info->is_structured() || dim == 1) // VTK unstructured grid
  {
    if (grid_info->is_quadratic()) // VTK quadratic elements
    {
      // Number of visualization points per Bezier element in each direction.
      // const auto n_pts_dir_per_bezier_elem = quad_->get_num_coords_direction();

      // Iteration along all the quadrature point.
      // Only the points in an edge are added to the mask.
      const auto &quad_points_1d = quad_->get_points_1d();
      for (int i_pt = 0; i_pt < n_quad_points; ++i_pt)
      {
        const auto tensor_id = quad_points_1d.flat_to_tensor(i_pt);

        // If the number of odd values is > 1. This means that the point is
        // not in an edge.
        // For the case 1D, this conditions is never fulfilled, no
        // points will be removed (this is what is desired).
        Size n_odd_values = 0;
        for (int dir = 0; dir < dim; ++dir)
          n_odd_values += tensor_id[dir] % 2;

        if (n_odd_values < 2) // It's on an edge.
          mask_.push_back(i_pt);
      } // end loop i_pt

      // Total number of visualization points per Bezier element.
      const Size n_pts_per_bezier_elem = mask_.size();

      Index point_id = 0;
      for (auto &pm_el : map_)
      {
        pm_el.resize(n_pts_per_bezier_elem);
        for (auto &pm : pm_el)
          pm = point_id++;
      } // end loop pm_el

      n_total_points_ = point_id;
    } // end VTK quadratic elements
    else // VTK linear elements
    {
      // Total number of visualization points per Bezier element.
      const Size n_pts_per_bezier_elem = n_quad_points;

      mask_.resize(n_pts_per_bezier_elem);
      Index id = 0;
      for (auto &pm : mask_)
        pm = id++;

      Index point_id = 0;
      for (auto &pm_el : map_)
      {
        pm_el.resize(n_pts_per_bezier_elem);
        for (auto &pm : pm_el)
          pm = point_id++;
      } // map_
      n_total_points_ = point_id;
    } // end VTK linear elements
  }
  else // end VTK structured grid
  {
    // Total number of visualization points per Bezier element.
    const Size n_pts_per_bezier_elem = n_quad_points;
    // Number of visualization points per Bezier element in each direction.
    const auto n_pts_dir_per_bezier_elem = quad_->get_num_coords_direction();

    mask_.resize(n_pts_per_bezier_elem);
    Index id = 0;
    for (auto &pm : mask_)
      pm = id++;

    n_total_points_ = n_pts_per_bezier_elem * n_bezier_elements;

    TensorSize <dim> n_pts_per_mesh; // Number of points per direction of
    // VTK structured grid.
    const auto n_intervals = cartesian_grid->get_num_intervals();
    for (int dir = 0; dir < dim; ++dir)
      n_pts_per_mesh[dir] = n_intervals[dir]
                            * n_pts_dir_per_bezier_elem[dir];

    // Tensorial index of the first point in Bezier element.
    TensorIndex <dim> pt_mesh_t_offset;

    // Tensorial index of the point.
    TensorIndex <dim> pt_mesh_t_id;

    // Tensorial index of the point referred to the number of points in a single element.
    TensorIndex <dim> pt_elem_t_id;

    const auto w_elem_pts = MultiArrayUtils <dim>::compute_weight(
                              n_pts_dir_per_bezier_elem);
    const auto w_mesh_pts = MultiArrayUtils <dim>::compute_weight(
                              n_pts_per_mesh);

    int i_el = 0;
    for (const auto &elem : *cartesian_grid)
    {
      auto &pmi = map_[i_el];
      pmi.resize(n_pts_per_bezier_elem);

      const auto &elem_t_id = elem.get_index().get_tensor_index();

      // Computing the tensor index of the first point of the element.
      for (int dir = 0; dir < dim; ++dir)
        pt_mesh_t_offset[dir] = elem_t_id[dir]
                                * n_pts_dir_per_bezier_elem[dir];

      Index point_id = 0;
      for (auto &pm : pmi)
      {
        // Computing the tensor index of the point referred to the number
        // of points into a single Bezier element.
        pt_elem_t_id = MultiArrayUtils <dim>::flat_to_tensor_index(
                         point_id++, w_elem_pts);

        // Computing the tensor index of the point referred to the number
        // of points into the whole mesh.
        for (int dir = 0; dir < dim; ++dir)
          pt_mesh_t_id[dir] = pt_mesh_t_offset[dir]
                              + pt_elem_t_id[dir];

        pm = MultiArrayUtils <dim>::tensor_to_flat_index(
               pt_mesh_t_id, w_mesh_pts);
      } // end loop pm
      ++i_el;
    } //end loop elem
  }
};



template <class Domain>
auto
VtkIgaSolidGrid<Domain>::
PointsTopology::
create_element_connectivity(const GridInfoPtr_ grid_info) -> Connectivity_
{
  Assert (grid_info != nullptr, ExcNullPtr());

  if (!grid_info->is_structured() || dim == 1) // VTK unstructured grid
  {
    if (grid_info->is_quadratic()) // VTK quadratic elements
      return Self_::create_quadratic_element_connectivity<dim>(grid_info);
    else // VTK linear elements
      return Self_::create_linear_element_connectivity(grid_info);
  }
  else // end VTK structured grid
      return Connectivity_();
}



template <class Domain>
Size
VtkIgaSolidGrid<Domain>::
PointsTopology::
get_num_pts_per_bezier_elem() const
{
    return mask_.size();
}



template <class Domain>
Size
VtkIgaSolidGrid<Domain>::
PointsTopology::
get_num_bezier_elems() const
{
    return map_.size();
}



template <class Domain>
Size
VtkIgaSolidGrid<Domain>::
PointsTopology::
get_flat_num_cells_per_bezier_elem() const
{
    return n_vis_elements_.flat_size();
}



template <class Domain>
Size
VtkIgaSolidGrid<Domain>::
PointsTopology::
get_num_total_pts() const
{
    return n_total_points_;
}



template <class Domain>
Size
VtkIgaSolidGrid<Domain>::
PointsTopology::
get_num_pts_per_single_vtk_cell() const
{
    return connectivity_[0].size();
}



template <class Domain>
Size
VtkIgaSolidGrid<Domain>::
PointsTopology::
get_num_vtk_cells_per_bezier_elem (const Index &dir) const
{
    Assert (dir >= 0 && dir < dim, ExcIndexRange(dir, 0, dim));
    return n_vis_elements_[dir];
}



template <class Domain>
const SafeSTLVector<Index> &
VtkIgaSolidGrid<Domain>::
PointsTopology::
get_mask() const
{
    return mask_;
}



template <class Domain>
auto
VtkIgaSolidGrid<Domain>::
PointsTopology::
get_connectivity() const -> const Connectivity_ &
{
    return connectivity_;
}



template <class Domain>
auto
VtkIgaSolidGrid<Domain>::
PointsTopology::
map_cbegin() const -> Map_::const_iterator
{
    return map_.cbegin();
}



template <class Domain>
auto
VtkIgaSolidGrid<Domain>::
PointsTopology::
get_quadrature() const ->QuadPtr_
{
    return quad_;
}




template <class Domain>
auto
VtkIgaSolidGrid<Domain>::
PointsTopology::
create_linear_element_connectivity(const GridInfoPtr_ grid_info) -> Connectivity_
{
  Assert (grid_info != nullptr, ExcNullPtr());
  Connectivity_ conn;

  // Number of vertices in a dim-dimensional square.
  static constexpr int n_vertices =
    UnitElement <dim>::template num_elem <0> ();

  const auto n_vis_elems = grid_info->get_num_cells_per_element<dim>();
  TensorSize <dim> n_elem_bound_per_dir;
  Size n_cells_per_bezier = 1;
  for (int dir = 0; dir < dim; ++dir)
  {
    n_elem_bound_per_dir[dir] = n_vis_elems[dir] + 1;
    n_cells_per_bezier *= n_vis_elems[dir];
  }

  // This grid is going to help in building the connectivity.
  // Every element of the grid refers to a cell.
  const auto cells_grid = Grid<dim>::const_create(n_elem_bound_per_dir);

  conn.resize(n_cells_per_bezier);

  const Size n_points_per_single_cell = n_vertices;

  // Creating the connectivity ---------------------------------------------//

  // Building the offsets container. According to the vtk elements connectivity.
  using T_ = SafeSTLArray < SafeSTLArray<int, dim>, n_vertices>;
  const T_ delta_idx =
  dim == 1 ? T_({{ 0 },
    { 1 }
  }) :   // dim = 1

    dim == 2 ? T_({{ 0, 0 },
    { 1, 0 },
    { 1, 1 },
    { 0, 1 }
  }) :   // dim = 2

    T_(
  {
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 1, 1, 0 },
    { 0, 1, 0 },
    { 0, 0, 1 },
    { 1, 0, 1 },
    { 1, 1, 1 },
    { 0, 1, 1 }
  }); // dim = 3

  const TensorIndex <dim> weight_points =
    MultiArrayUtils <dim>::compute_weight(n_elem_bound_per_dir);

  TensorIndex <dim> vtk_vertex_tensor_idx;

  auto conn_el = conn.begin();
  auto cell = cells_grid->begin();
  const auto end = cells_grid->end();
  for (; cell != end; ++cell, ++conn_el)
  {
    conn_el->resize(n_points_per_single_cell);

    const auto &vtk_elem_tensor_idx = cell->get_index().get_tensor_index();

    auto conn = conn_el->begin();
    for (int iVertex = 0; iVertex < n_points_per_single_cell;
         ++iVertex, ++conn)
    {
      for (int i = 0; i < dim; ++i)
        vtk_vertex_tensor_idx[i] = vtk_elem_tensor_idx[i]
                                   + delta_idx[iVertex][i];

      *conn = MultiArrayUtils <dim>::tensor_to_flat_index
              (vtk_vertex_tensor_idx, weight_points);
    } // end loop iVertex
  } // end loop cell
  //--------------------------------------------------------------------------//

  return conn;
}



template <class Domain>
template <int aux_dim>
auto
VtkIgaSolidGrid<Domain>::
PointsTopology::
create_quadratic_element_connectivity(const GridInfoPtr_ grid_info,
  typename std::enable_if_t<aux_dim == 1> *) -> Connectivity_
{
  Assert (grid_info != nullptr, ExcNullPtr());

  Connectivity_ conn;

  //  This is the connectivity pattern of the 1D VTK quadratic element.
  //
  //         0 -- 2 -- 1

  // Number of cells per Bezier element.
  Size n_cells_per_bezier = 1;
  const auto n_vis_elems = grid_info->get_num_cells_per_element<dim>();
  for (int dir = 0; dir < dim; ++dir)
      n_cells_per_bezier *= n_vis_elems[dir];

  conn.clear();
  conn.resize(n_cells_per_bezier);

  Index point_id = 0;
  for (auto &conn_el : conn)
  {
    conn_el = {point_id, point_id+2, point_id+1};
    point_id += 2;
  }

  return conn;
}



template <class Domain>
template <int aux_dim>
auto
VtkIgaSolidGrid<Domain>::
PointsTopology::
create_quadratic_element_connectivity(const GridInfoPtr_ grid_info,
  typename std::enable_if_t<aux_dim == 2> *) -> Connectivity_
{
  Assert (grid_info != nullptr, ExcNullPtr());

  Connectivity_ conn;

  //  This is the connectivity pattern of the 2D VTK quadratic element.
  //
  //         3 -- 6 -- 2
  //         |         |
  //         7         5
  //         |         |
  //         0 -- 4 -- 1

  // Number of cells per Bezier element.
  const auto n_vis_elems = grid_info->get_num_cells_per_element<dim>();
  Size n_cells_per_bezier = 1;
  for (int dir = 0; dir < dim; ++dir)
      n_cells_per_bezier *= n_vis_elems[dir];

  conn.clear();
  conn.resize(n_cells_per_bezier);

  TensorSize <aux_dim> n_full_points_per_dir;
  for (int dir = 0; dir < aux_dim; ++dir)
    n_full_points_per_dir[dir] = n_vis_elems[dir] * 2 + 1;

  static const int n_points_per_single_cell = 8;

  TensorSize <aux_dim> n_elem_bound_per_dir;
  for (int dir = 0; dir < aux_dim; ++dir)
    n_elem_bound_per_dir[dir] = n_vis_elems[dir] + 1;

  // This grid is going to help in building the connectivity.
  // Every element of the grid refers to a cell.
  const auto cells_grid = Grid <aux_dim>::const_create(n_elem_bound_per_dir);

  // This array constaints the offsets of the points along the
  // first (u) direction.
  //   Offset for the line:  3 -- 6 -- 2 -> offset = 2
  //   Offset for the line:  7         5 -> offset = 1
  //   Offset for the line:  0 -- 4 -- 1 -> offset = 2

  const SafeSTLArray <Index, 8> offsets_u =
  {
    { 2, 2, 2, 2, 2, 1, 2, 1 }
  };

  // This container if for the offsets of the points along the
  // second (v) direction.
  const Index offsets_v = n_full_points_per_dir[0] * 2
                          - n_vis_elems[0];

  SafeSTLVector <Index> vtk_vertex_id_0(n_points_per_single_cell);
  vtk_vertex_id_0[0] = 0;
  vtk_vertex_id_0[4] = 1;
  vtk_vertex_id_0[1] = 2;

  vtk_vertex_id_0[7] = n_full_points_per_dir[0];
  vtk_vertex_id_0[5] = vtk_vertex_id_0[7] + 1;

  vtk_vertex_id_0[3] = vtk_vertex_id_0[7] + n_full_points_per_dir[0]
                       - n_vis_elems[0];
  vtk_vertex_id_0[6] = vtk_vertex_id_0[3] + 1;
  vtk_vertex_id_0[2] = vtk_vertex_id_0[3] + 2;

  auto conn_el = conn.begin();
  for (const auto &cell : *cells_grid)
  {
    conn_el->resize(n_points_per_single_cell);

    const auto &vtk_elem_tensor_idx = cell.get_index().get_tensor_index();

    auto conn = conn_el->begin();
    for (int i_pt = 0; i_pt < n_points_per_single_cell; ++i_pt, ++conn)
    {
      *conn = vtk_vertex_id_0[i_pt]
              + offsets_u[i_pt] * vtk_elem_tensor_idx[0]
              + offsets_v * vtk_elem_tensor_idx[1];
    } // end loop i_pt
    ++conn_el;
  } // end loop cell

  return conn;
}



template <class Domain>
template <int aux_dim>
auto
VtkIgaSolidGrid<Domain>::
PointsTopology::
create_quadratic_element_connectivity(const GridInfoPtr_ grid_info,
  typename std::enable_if_t<aux_dim == 3> *) -> Connectivity_
{
  Assert (grid_info != nullptr, ExcNullPtr());

  Connectivity_ conn;

  //  This is the connectivity pattern of the 3D VTK quadratic element.
  //
  //         w = 0             w = 1             w = 2
  //
  //     3 -- 10 --  2    19 -------- 18     7 -- 14 --  6
  //     |           |     |           |     |           |
  //    11           9     |           |    15          13
  //     |           |     |           |     |           |
  //     0 --  8 --  1    16 -------- 17     4 -- 12 --  5
  //

  // Number of cells per Bezier element.
  const auto n_vis_elems = grid_info->get_num_cells_per_element<dim>();
  Size n_cells_per_bezier = 1;
  for (int dir = 0; dir < dim; ++dir)
      n_cells_per_bezier *= n_vis_elems[dir];

  conn.clear();
  conn.resize(n_cells_per_bezier);

  TensorSize <aux_dim> n_full_points_per_dir;
  for (int dir = 0; dir < aux_dim; ++dir)
    n_full_points_per_dir[dir] = n_vis_elems[dir] * 2 + 1;

  static const int n_points_per_single_cell = 20;

  TensorSize <aux_dim> n_elem_bound_per_dir;
  for (int dir = 0; dir < aux_dim; ++dir)
    n_elem_bound_per_dir[dir] = n_vis_elems[dir] + 1;

  // This grid is going to help in building the connectivity.
  // Every element of the grid refers to a cell.
  const auto cells_grid = Grid <dim>::const_create(n_elem_bound_per_dir);

  // This array constaints the offsets of the points along the
  // first (u) direction.
  //   Offset for the line:  3 -- 10 --  2   -> offset = 2
  //   Offset for the line: 11 --------  9   -> offset = 1
  //   Offset for the line:  0 --  8 --  1   -> offset = 2

  //   Offset for the line: 19 -------- 18   -> offset = 1
  //   Offset for the line: 16 -------- 17   -> offset = 1

  //   Offset for the line:  7 -- 14 --  6   -> offset = 2
  //   Offset for the line: 15 -------- 13   -> offset = 1
  //   Offset for the line:  4 -- 12 --  5   -> offset = 2

  const SafeSTLArray <Index, 20> offsets_u =
  {
    {
      2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 2,
      1, 1, 1,
      1, 1
    }
  };

  TensorSize <aux_dim> n_full_points_per_dir_red;
  for (int dir = 0; dir < aux_dim; ++dir)
    n_full_points_per_dir_red[dir] = n_full_points_per_dir[dir]
                                     - n_vis_elems[dir];

  // This array constaints the offsets of the points along the
  // second (v) direction.
  //   Offset for the line:  3 -- 10 --  2   -> offset = ov0
  //   Offset for the line: 11 --------  9   -> offset = ov0
  //   Offset for the line:  0 --  8 --  1   -> offset = ov0

  //   Offset for the line: 19 -------- 18   -> offset = ov1
  //   Offset for the line: 16 -------- 17   -> offset = ov1

  //   Offset for the line:  7 -- 14 --  6   -> offset = ov0
  //   Offset for the line: 15 -------- 13   -> offset = ov0
  //   Offset for the line:  4 -- 12 --  5   -> offset = ov0

  const Index ov0 = n_full_points_per_dir[0]
                    + n_full_points_per_dir_red[0];
  const Index ov1 = n_full_points_per_dir_red[0];
  const SafeSTLArray <Index, 20> offsets_v =
  {
    {
      ov0, ov0, ov0, ov0, ov0, ov0, ov0, ov0,
      ov0, ov0, ov0, ov0, ov0, ov0,
      ov0, ov0, ov1, ov1, ov1, ov1
    }
  };

  // Offset along the third (w) direction.
  const Size ow0 = n_full_points_per_dir[0] * n_full_points_per_dir_red[1]
                   + n_full_points_per_dir_red[0] * n_vis_elems[1];
  const Size ow1 = n_full_points_per_dir_red[0]
                   * n_full_points_per_dir_red[1];

  const Index offsets_w = ow0 + ow1;

  SafeSTLVector <Index> vtk_vertex_id_0(n_points_per_single_cell);
  vtk_vertex_id_0[0] = 0;
  vtk_vertex_id_0[8] = 1;
  vtk_vertex_id_0[1] = 2;
  vtk_vertex_id_0[11] = vtk_vertex_id_0[0] + n_full_points_per_dir[0];
  vtk_vertex_id_0[9] = vtk_vertex_id_0[11] + 1;
  vtk_vertex_id_0[3] = vtk_vertex_id_0[11] + n_full_points_per_dir_red[0];
  vtk_vertex_id_0[10] = vtk_vertex_id_0[3] + 1;
  vtk_vertex_id_0[2] = vtk_vertex_id_0[3] + 2;

  vtk_vertex_id_0[16] = vtk_vertex_id_0[0] + ow0;
  vtk_vertex_id_0[17] = vtk_vertex_id_0[16] + 1;
  vtk_vertex_id_0[19] = vtk_vertex_id_0[16]
                        + n_full_points_per_dir_red[0];
  vtk_vertex_id_0[18] = vtk_vertex_id_0[19] + 1;

  vtk_vertex_id_0[4] = vtk_vertex_id_0[16] + ow1;
  vtk_vertex_id_0[12] = vtk_vertex_id_0[4] + 1;
  vtk_vertex_id_0[5] = vtk_vertex_id_0[4] + 2;
  vtk_vertex_id_0[15] = vtk_vertex_id_0[4] + n_full_points_per_dir[0];
  vtk_vertex_id_0[13] = vtk_vertex_id_0[15] + 1;
  vtk_vertex_id_0[7] = vtk_vertex_id_0[15] + n_full_points_per_dir_red[0];
  vtk_vertex_id_0[14] = vtk_vertex_id_0[7] + 1;
  vtk_vertex_id_0[6] = vtk_vertex_id_0[7] + 2;

  auto conn_el = conn.begin();
  for (const auto &cell : *cells_grid)
  {
    conn_el->resize(n_points_per_single_cell);

    const auto &vtk_elem_tensor_idx = cell.get_index().get_tensor_index();

    auto conn = conn_el->begin();
    for (int i_pt = 0; i_pt < n_points_per_single_cell; ++i_pt, ++conn)
    {
      *conn = vtk_vertex_id_0[i_pt]
              + offsets_u[i_pt] * vtk_elem_tensor_idx[0]
              + offsets_v[i_pt] * vtk_elem_tensor_idx[1]
              + offsets_w * vtk_elem_tensor_idx[2];
    } // end loop i_pt
    ++conn_el;
  } // end loop cell

  return conn;
}



template <class Domain>
auto
VtkIgaSolidGrid<Domain>::
PointsTopology::
create_visualization_quadrature(const GridInfoPtr_ grid_info) ->
QuadPtr_
{
  Assert (grid_info != nullptr, ExcNullPtr());

  const auto n_vis_elems = grid_info->get_num_cells_per_element<dim>();
  const Size n_pts_per_cell_dir =  grid_info->is_quadratic() ? 3 : 2;
  TensorSize <dim> n_points;
  for (int dir = 0; dir < dim; ++dir)
    n_points[dir] = n_vis_elems[dir] * (n_pts_per_cell_dir - 1) + 1;

  return QUniform <dim>::create(n_points);

#if 0
  // TODO: Right now is not possible to use non tensor product quadrature
  // for cache elements in igatools. It will be in the future.
  // The code below defines a non tensor product quadrature for VTK quadratic
  // cells.

  const auto &uniform_array_points = uniform_quad.get_points_1d();
  const Size n_total_points = uniform_quad.get_num_points();

  // The points that are not in an edge are removed from the quadrature points.
  // The condition to be in an edge is that just one or zero directions are
  // active. (This also works for the 1D case, where no point is removed.)
  SafeSTLVector<Index> point_indices;
  for (int i_pt = 0; i_pt < n_total_points; ++i_pt)
  {
    const auto t_id = uniform_array_points.flat_to_tensor(i_pt);

    int active_directions = 0;
    for (int dir = 0; dir < dim; ++dir)
      if ((t_id[dir] > 0) && (t_id[dir] < (n_points[dir] - 1)))
        ++active_directions;

    if (active_directions < 2)
      point_indices.push_back(i_pt);
  }

  ValueVector<Points<dim>> quad_points(point_indices.size());
  auto qp = quad_points.begin();
  for (const auto &id : point_indices)
    *qp++ = uniform_quad.get_point(id);

  BBox<dim> box;

  return shared_ptr<Quadrature<dim>> (new Quadrature<dim> (quad_points,
                                                           uniform_quad.get_weights_1d(), box));
#endif
}



template <class Domain>
auto
VtkIgaSolidGrid<Domain>::
create(const DomainPtr_ domain,
       const ObjContPtr_t_ objs_container,
       const GridInfoPtr_ grid_info,
       const bool is_physical) -> VtkGridPtr_
{
  Assert (domain != nullptr, ExcNullPtr());
  Assert (grid_info != nullptr, ExcNullPtr());
  Assert (objs_container != nullptr, ExcNullPtr());

  // 1D geometries must be unstructured.
  if (!grid_info->is_structured() || dim == 1)  // VTK unstructured grid.
    return Self_::create_grid_vtu(domain, objs_container, grid_info,
                                  is_physical);
  else // VTK structured grid.
    return Self_::create_grid_vts(domain, objs_container, grid_info,
                                  is_physical);
}



template <class Domain>
auto
VtkIgaSolidGrid<Domain>::
create_grid_vts(const DomainPtr_ domain,
                const ObjContPtr_t_ objs_container,
                const GridInfoPtr_ grid_info,
                const bool is_physical) ->
VtkGridPtr_
{
  Assert (domain != nullptr, ExcNullPtr());
  Assert (objs_container != nullptr, ExcNullPtr());
  Assert (grid_info != nullptr, ExcNullPtr());

  PointsTopology topology (domain->get_grid_function()->get_grid(), grid_info);

  const auto points = Self_::create_points(domain, topology);

  if (dim == 1)
      VtkIgaWarningMacro("It is not possible to create vtk structured "
                         "grids for 1D geometries.");

  auto vts_grid = vtkSmartPointer <vtkStructuredGrid>::New();

  const auto n_intervals = domain->get_grid_function()->get_grid()->get_num_intervals();
  SafeSTLArray<int, 3> grid_dim({{1, 1, 1}});
  for (int dir = 0; dir < dim; ++dir)
    grid_dim[dir] = (topology.get_num_vtk_cells_per_bezier_elem(dir) + 1) * n_intervals[dir];
  vts_grid->SetDimensions(grid_dim[0], grid_dim[1], grid_dim[2]);

  vts_grid->SetPoints(points);

  if (is_physical)
      Self_::create_point_data_physical(domain, objs_container,
                                        topology, vts_grid);
  else
      Self_::create_point_data_parametric(domain, objs_container,
                                          topology, vts_grid);

  return vts_grid;
}



template <class Domain>
auto
VtkIgaSolidGrid<Domain>::
create_grid_vtu(const DomainPtr_ domain,
                const ObjContPtr_t_ objs_container,
                const GridInfoPtr_ grid_info,
                const bool is_physical) ->
VtkGridPtr_
{
  Assert (domain != nullptr, ExcNullPtr());
  Assert (objs_container != nullptr, ExcNullPtr());
  Assert (grid_info != nullptr, ExcNullPtr());

  PointsTopology topology (domain->get_grid_function()->get_grid(), grid_info);

  const auto points = Self_::create_points(domain, topology);

  auto vtu_grid = vtkSmartPointer <vtkUnstructuredGrid>::New();

  vtu_grid->SetPoints(points);

  const Size n_bezier_elements = topology.get_num_bezier_elems();

  // Total number of cells.
  const int n_total_cells = n_bezier_elements * topology.get_flat_num_cells_per_bezier_elem();

  auto cell_ids = vtkSmartPointer <vtkIdTypeArray>::New();

  const Size tuple_size = topology.get_num_pts_per_single_vtk_cell() + 1;
  cell_ids->SetNumberOfComponents(tuple_size);
  cell_ids->SetNumberOfTuples(n_total_cells);

  std::vector<vtkIdType> tuple(tuple_size);
  tuple[0] = tuple_size - 1;

  Index cell_id = 0;
  for (int be = 0; be < n_bezier_elements; ++be)
  {
    const int vtk_vertex_id_offset = topology.get_num_pts_per_bezier_elem() * be;

    // Iterating along the cells of one Bezier element.
    for (const auto &cc : topology.get_connectivity())
    {
      int point_id = 1;
      for (const auto &c : cc)
        tuple[point_id++] = c + vtk_vertex_id_offset;

      cell_ids->SetTupleValue(cell_id++, tuple.data());
    }
  }

  const auto cells = vtkSmartPointer <vtkCellArray>::New();
  cells->SetCells(cell_ids->GetNumberOfTuples(), cell_ids);

  vtu_grid->Allocate(cells->GetNumberOfCells(), 0);
  const int vtk_enum_type =
    grid_info->is_quadratic() ?
    dim == 3 ? VTK_QUADRATIC_HEXAHEDRON :
    dim == 2 ? VTK_QUADRATIC_QUAD :
    VTK_QUADRATIC_EDGE
    :
    dim == 3 ? VTK_HEXAHEDRON :
    dim == 2 ? VTK_QUAD :
    VTK_LINE;

  vtu_grid->SetCells(vtk_enum_type, cells);


  if (is_physical)
      Self_::create_point_data_physical(domain, objs_container,
                                        topology, vtu_grid);
  else
      Self_::create_point_data_parametric(domain, objs_container,
                                          topology, vtu_grid);

  return vtu_grid;
}



template <class Domain>
vtkSmartPointer<vtkPoints>
VtkIgaSolidGrid<Domain>::
create_points(const DomainPtr_ domain,
              const PointsTopology &points_top)
{
  Assert(domain != nullptr, ExcNullPtr());

  const auto points = vtkSmartPointer <vtkPoints>::New();
  points->SetNumberOfPoints(points_top.get_num_total_pts());

  SafeSTLArray<Real, 3> point_tmp(0.0);

  auto domain_cache_handler = domain->create_cache_handler();
  domain_cache_handler->template set_flags<dim>(domain_element::Flags::point);


  auto elem = domain->cbegin();
  const auto end = domain->cend();

  domain_cache_handler->init_cache(elem, points_top.get_quadrature());

  for (auto pm_el = points_top.map_cbegin(); elem != end; ++elem, ++pm_el)
  {
    domain_cache_handler->template fill_cache<dim>(elem,0);

    const auto element_vertices_tmp = elem->template get_points<dim>(0);
    auto pm = pm_el->cbegin();
    for (const auto &mask : points_top.get_mask())
    {
      const auto &point = element_vertices_tmp[mask];
      for (int dir = 0; dir < space_dim; ++dir)
        point_tmp[dir] = point[dir];
      points->SetPoint(*pm++, point_tmp.data());
    }
  }

  return points;
}



template <class Domain>
void
VtkIgaSolidGrid<Domain>::
create_point_data_physical(const DomainPtr_ domain,
                           const ObjContPtr_t_ objs_container,
                           const PointsTopology &points_top,
                           const VtkGridPtr_ vtk_grid)
{
  Assert (domain != nullptr, ExcNullPtr());
  Assert (objs_container != nullptr, ExcNullPtr());
  Assert (vtk_grid != nullptr, ExcNullPtr());

  vtkPointData *point_data = vtk_grid->GetPointData();

  ValidFuncs_ valid_f_ptr_types;
  for_each(valid_f_ptr_types, [&](const auto &ptr_type)
  {
    using FunctionType = typename remove_reference<decltype(ptr_type)>::type::element_type;


    // Constant functions from the objects container.
    for (const auto &id : objs_container->template get_const_object_ids<FunctionType>())
    {
        const auto func = objs_container->template get_const_object<FunctionType>(id);

        // Checking if func corresponds to the domain of the generator.
        if (func->get_domain() == domain)
        {
            using Value = typename FunctionType::Value;
            using _D0 = typename function_element::template _D<0>;
            using Flags = function_element::Flags;

            const vtkIdType n_tuples = points_top.get_num_bezier_elems() * points_top.get_num_pts_per_bezier_elem();

            auto arr = vtkSmartPointer <vtkDoubleArray>::New();

            arr->SetName(func->get_name().c_str());
            arr->SetNumberOfComponents(Value::size);
            arr->SetNumberOfTuples(n_tuples);

            auto func_handler = func->create_cache_handler();
            func_handler->set_element_flags(Flags::D0);

            auto elem = func->cbegin();
            const auto end = func->cend();

            func_handler->init_cache(*elem, points_top.get_quadrature());

            auto pnm_it = points_top.map_cbegin();
            for (; elem != end; ++elem, ++pnm_it)
            {
                func_handler->fill_element_cache(*elem);

                auto pnm = pnm_it->cbegin();
                auto values = elem->template get_values_from_cache<_D0, dim>(0);
                for (const auto &pm : points_top.get_mask())
                    arr->SetTuple(*pnm++, values[pm].get_flat_values().data());
            }

            point_data->AddArray(arr.Get());

        }
    }
  });

}



template <class Domain>
void
VtkIgaSolidGrid<Domain>::
create_point_data_parametric(const DomainPtr_ domain,
                             const ObjContPtr_t_ objs_container,
                             const PointsTopology &points_top,
                             const VtkGridPtr_ vtk_grid)
{
  Assert (domain != nullptr, ExcNullPtr());
  Assert (objs_container != nullptr, ExcNullPtr());
  Assert (vtk_grid != nullptr, ExcNullPtr());

  vtkPointData *point_data = vtk_grid->GetPointData();

  ValidGridFuncs_ valid_gf_ptr_types;
  for_each(valid_gf_ptr_types, [&](const auto &ptr_type)
  {
    using GridFunctionType = typename remove_reference<decltype(ptr_type)>::type::element_type;

    // Constant functions from the objects container.
    for (const auto &id : objs_container->template get_const_object_ids<GridFunctionType>())
    {
        const auto grid_func = objs_container->template get_const_object<GridFunctionType>(id);
        const auto grid = domain->get_grid_function()->get_grid();

        // Checking if func corresponds to the domain of the generator.
        if (grid_func->get_grid() == grid)
        {
            using Value = typename GridFunctionType::Value;
            using _D0 = typename grid_function_element::template _D<0>;
            using Flags = grid_function_element::Flags;

            const vtkIdType n_tuples = points_top.get_num_bezier_elems() * points_top.get_num_pts_per_bezier_elem();

            auto arr = vtkSmartPointer <vtkDoubleArray>::New();

            arr->SetName(grid_func->get_name().c_str());
            arr->SetNumberOfComponents(Value::size);
            arr->SetNumberOfTuples(n_tuples);

            auto grid_func_handler = grid_func->create_cache_handler();
            grid_func_handler->set_element_flags(Flags::D0);

            auto elem = grid_func->cbegin();
            const auto end = grid_func->cend();

            grid_func_handler->init_cache(*elem, points_top.get_quadrature());

            auto pnm_it = points_top.map_cbegin();
            for (; elem != end; ++elem, ++pnm_it)
            {
                grid_func_handler->fill_element_cache(*elem);

                auto pnm = pnm_it->cbegin();
                auto values = elem->template get_values_from_cache<_D0, dim>(0);
                for (const auto &pm : points_top.get_mask())
                    arr->SetTuple(*pnm++, values[pm].get_flat_values().data());
            }

            point_data->AddArray(arr.Get());
        }
    }
  });

}

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#include <paraview_plugin/vtk_iga_solid_grid.inst>
