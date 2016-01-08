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

#include <paraview_plugin/vtk_iga_knot_grid.h>

#include <boost/range/irange.hpp>

#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/domain_element.h>

#include <paraview_plugin/vtk_iga_grid_information.h>

using std::shared_ptr;

IGA_NAMESPACE_OPEN

template <class Domain>
auto
VtkIgaKnotGrid<Domain>::
create(const DomainPtr_ domain,
       const GridInfoPtr_ grid_info) -> VtkGridPtr_
{
  return Self_::create_grid<Domain::dim>(domain, grid_info);
}



template <class Domain>
template<int aux_dim>
auto
VtkIgaKnotGrid<Domain>::
create_grid(const DomainPtr_ domain,
            const GridInfoPtr_ grid_info) ->
EnableIf<aux_dim == 1, VtkGridPtr_>
{
  // Implementation for 1D case.
  // In this case the grid consists on a set of points corresponding
  // to the knots.

  Assert(domain != nullptr, ExcNullPtr());
  Assert(grid_info != nullptr, ExcNullPtr());

  const auto cartesian_grid = domain->get_grid_function()->get_grid();
  const auto &n_intervals = cartesian_grid->get_num_intervals();

  const auto quad = QUniform<dim>::create(2);
  const Size n_vtk_points = n_intervals[0] + 1;
  const Size n_vtk_cells = n_vtk_points;

  const auto vtk_points = vtkSmartPointer <vtkPoints>::New();
  vtk_points->SetNumberOfPoints(n_vtk_points);

  const auto vtk_cell_ids = vtkSmartPointer <vtkIdTypeArray>::New();

  const Size tuple_size = 2;
  vtk_cell_ids->SetNumberOfComponents(tuple_size);
  vtk_cell_ids->SetNumberOfTuples(n_vtk_cells);

  auto elem = domain->cbegin();
  const auto end = domain->cend();

  SafeSTLArray<Real, 3> point_tmp(0.0);
  std::array<vtkIdType, 2> tuple({{ 1, 0 }});

  auto domain_cache_handler = domain->create_cache_handler();
  domain_cache_handler->template set_flags<dim>(domain_element::Flags::point);
  domain_cache_handler->init_cache(elem, quad);

  Index pt_id = 0;
  // Adding first point of every element.
  for (; elem != end; ++elem, ++pt_id)
  {
    domain_cache_handler->template fill_cache<dim>(elem, 0);
    const auto points = elem->template get_points<dim>(0);
    const auto &pp = points[0];
    for (const auto &dir : boost::irange(0, space_dim))
      point_tmp[dir] = pp[dir];
    vtk_points->SetPoint(pt_id, point_tmp.data());
    tuple[1] = pt_id;
    vtk_cell_ids->SetTupleValue(pt_id, tuple.data());
  }

  // Adding last point of the last element.
  const auto points = elem->template get_points<dim>(0);
  const auto &pp = points[1];
  for (const auto &dir : boost::irange(0, space_dim))
    point_tmp[dir] = pp[dir];
  vtk_points->SetPoint(pt_id, point_tmp.data());
  tuple[1] = pt_id;
  vtk_cell_ids->SetTupleValue(pt_id, tuple.data());

  // Creating grid.
  auto grid = vtkSmartPointer <vtkUnstructuredGrid>::New();
  grid->SetPoints(vtk_points);

  // Creating cells.
  const auto vtk_cells = vtkSmartPointer <vtkCellArray>::New();
  vtk_cells->SetCells(n_vtk_cells, vtk_cell_ids);
  grid->Allocate(vtk_cells->GetNumberOfCells(), 0);
  const int vtk_enum_type = VTK_VERTEX;
  grid->SetCells(vtk_enum_type, vtk_cells);

  return grid;
}



template <class Domain>
template<int aux_dim>
auto
VtkIgaKnotGrid<Domain>::
create_grid(const DomainPtr_ domain,
            const GridInfoPtr_ grid_info) ->
EnableIf<aux_dim == 2 || aux_dim == 3, VtkGridPtr_>
{
  // Implementation for 2D and 3D cases.

  const auto num_visualization_elements =
  grid_info->get_num_cells_per_element<dim>();

  TensorSize <dim> n_vis_elements;
  for (int dir = 0; dir < dim; ++dir)
  {
    n_vis_elements[dir] = num_visualization_elements[dir];
    Assert(n_vis_elements[dir] > 0,
    ExcMessage("The number of visualization elements must be > 0."))
  }

  const auto cartesian_grid = domain->get_grid_function()->get_grid();
  const auto &n_intervals = cartesian_grid->get_num_intervals();

  const Size n_points_per_single_cell = grid_info->is_quadratic() ? 3 : 2;
  const SafeSTLVector <int> vtk_line_connectivity =
  grid_info->is_quadratic() ?
  SafeSTLVector <int> ({{ 0, 2, 1 } }) :
    SafeSTLVector <int> ({{ 0, 1 } });

  // A 1D quadrature is built in every direction.
  SafeSTLArray <shared_ptr<Quadrature <1>>,dim> quadratures;
  for (const auto &dir : boost::irange(0, dim))
    quadratures[dir] = QUniform<1>::create(
      n_vis_elements[dir] * (n_points_per_single_cell - 1) + 1);

  const Size n_pts_per_cell_dir = grid_info->is_quadratic() ? 3 : 2;
  TensorSize <aux_dim> n_points;
  for (int dir = 0; dir < aux_dim; ++dir)
    n_points[dir] = n_vis_elements[dir] * (n_pts_per_cell_dir - 1) + 1;

  SafeSTLVector <Real> zero_vec(1, 0.0);
  SafeSTLVector <Real> one_vec(1, 1.0);

  // Computing total number of vtk points and cells.
  Size n_vtk_cells = 0;
  Size n_vtk_points = 0;

  for (const auto &dir : boost::irange(0, dim))
  {
    const Index face_id = dir * 2;

    Size n_knot_lines = 1;
    const auto &face_elem = UnitElement <dim>::template get_elem <dim - 1>
    (face_id);
    for (const auto &face_dir : face_elem.active_directions)
      n_knot_lines *= (n_intervals[face_dir] + 1);
    const Size n_cells_in_knot_line = n_intervals[dir]
    * n_vis_elements[dir];
    n_vtk_cells += n_knot_lines * n_cells_in_knot_line;
    n_vtk_points +=
    n_knot_lines
    * (n_cells_in_knot_line
    * (n_points_per_single_cell - 1) + 1);
  }

  // Creating points.
  const auto vtk_points = vtkSmartPointer <vtkPoints>::New();
  vtk_points->SetNumberOfPoints(n_vtk_points);

  SafeSTLArray<Real, 3> point_tmp(0.0);

  // Creating cells.
  const auto vtk_cell_ids = vtkSmartPointer <vtkIdTypeArray>::New();
  const Size tuple_size = n_points_per_single_cell + 1;
  vtk_cell_ids->SetNumberOfComponents(tuple_size);
  vtk_cell_ids->SetNumberOfTuples(n_vtk_cells);

  std::vector<vtkIdType> tuple(tuple_size);
  tuple[0] = n_points_per_single_cell;

  // Global counters for points and cell tuples.
  Index vtk_pt_id = 0;
  Index vtk_tuple_id = 0;

  // Creating knot lines along each direction.
  for (const auto &dir : boost::irange(0, dim))
  {
    const Index face_id = dir * 2;
    const auto &face_elem = UnitElement <dim>::template get_elem <dim - 1>
    (face_id);


    // Looping along all the knot coordinates of the face.
    typename Grid<dim>::template SubGridMap<dim-1> elem_map;
    const auto sub_grid =
    cartesian_grid->template get_sub_grid <dim-1>(face_id,elem_map);
    const auto &face_coords_t_size = sub_grid->get_num_knots_dim();
    const Size n_pts_face = face_coords_t_size.flat_size();
    const auto face_coords_t_w = MultiArrayUtils<dim-1>::compute_weight(face_coords_t_size);


    // Points and weights needed for creating quadratures with points
    // only on one edge (along the direction dir).
    const auto &quad_dir = quadratures[dir];
    TensorSize <dim> n_quad_points(1);
    n_quad_points[dir] = quad_dir->get_num_points();
    TensorProductArray <dim> quad_points_1d(n_quad_points);
    TensorProductArray <dim> quad_weights_1d(n_quad_points);
    quad_points_1d.copy_data_direction(
      dir, quad_dir->get_coords_direction(0));

    // Iterating along every knot coordinates of the face
    auto grid = domain->get_grid_function()->get_grid();
    auto elem = domain->cbegin();


    TensorIndex<dim> elem_t_id(0);
    for (const auto &i_pt : boost::irange(0, n_pts_face))
    {
      // Creating a quadrature along a certain edge of the element.
      const auto t_id_face = MultiArrayUtils<dim-1>::flat_to_tensor_index(i_pt,face_coords_t_w);
      auto ad = face_elem.active_directions.cbegin();
      for (const auto &t : t_id_face)
      {
        if (t == n_intervals[*ad]) // Last point in this face direction.
        {
          elem_t_id[*ad] = t - 1;
          quad_points_1d.copy_data_direction(*ad, one_vec);
        }
        else
        {
          elem_t_id[*ad] = t;
          quad_points_1d.copy_data_direction(*ad, zero_vec);
        }
        ++ad;
      }
      const auto quad = Quadrature<dim>::create(quad_points_1d,quad_weights_1d,BBox<dim>());

      auto domain_cache_handler = domain->create_cache_handler();
      domain_cache_handler->template set_flags<dim>(domain_element::Flags::point);
      domain_cache_handler->init_cache(elem, quad);

      elem_t_id[dir] = 0;

      auto &elem_t_id_dir = elem_t_id[dir];

      // Iterating along a single knot line in the direction dir.
      for (int itv = 0; itv < n_intervals[dir]; ++itv, ++elem_t_id_dir)
      {
        const int elem_f_id = grid->tensor_to_flat_element_id(elem_t_id);
        elem->move_to(ElementIndex<dim>(elem_f_id,elem_t_id));

        domain_cache_handler->template fill_cache<dim>(*elem, 0);
        const auto &physical_points = elem->template get_points<dim>(0);

        // Filling points.
        Index vtk_pt_id_0 = vtk_pt_id;
        for (const auto &pp : physical_points)
        {
          for (const auto &dir2 : boost::irange(0, space_dim))
            point_tmp[dir2] = pp[dir2];
          vtk_points->SetPoint(vtk_pt_id++, point_tmp.data());
        } // pp
        --vtk_pt_id; // The first point of the next element will be filled
        // into the position of the last point of the current
        // element.

        // Filling connectivity.
        for (int c_id = 0; c_id < n_vis_elements[dir];
             ++c_id,
             vtk_pt_id_0 +=
               (n_points_per_single_cell - 1))
        {
          Index id = 1;
          for (const auto &c : vtk_line_connectivity)
            tuple[id++] = vtk_pt_id_0 + c;

          vtk_cell_ids->SetTupleValue(vtk_tuple_id++, tuple.data());
        } // c_id

      } // itv
      ++vtk_pt_id; // Advancing point position for the next knot line.
    } // i_pt

  } // dir

  // Creating cells.
  const auto vtk_cells = vtkSmartPointer <vtkCellArray>::New();
  vtk_cells->SetCells(n_vtk_cells, vtk_cell_ids);

  // Creating grid.
  auto grid = vtkSmartPointer <vtkUnstructuredGrid>::New();
  grid->SetPoints(vtk_points);
  grid->Allocate(vtk_cells->GetNumberOfCells(), 0);
  const int vtk_enum_type =
    grid_info->is_quadratic() ? VTK_QUADRATIC_EDGE : VTK_LINE;
  grid->SetCells(vtk_enum_type, vtk_cells);

  return grid;
}

IGA_NAMESPACE_CLOSE

#include <paraview_plugin/vtk_iga_knot_grid.inst>
