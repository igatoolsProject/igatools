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

#include <paraview_plugin/iga_vtk.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkHexahedron.h>
#include <vtkQuad.h>
#include <vtkLine.h>
#include <vtkPolyLine.h>

#include <igatools/functions/functions_container.h>
#include <igatools/base/quadrature_lib.h>

using std::string;
using std::shared_ptr;
using std::pair;


IGA_NAMESPACE_OPEN


IGAVTK::
IGAVTK ()
  :
  file_name_(),
  file_path_(),
  num_visualization_elements_(TensorSize<3>(1)),
  quadratic_cells_(false),
  unstructured_grid_(true),
  funcs_container_(std::make_shared<FunctionsContainer>())
{};



void
IGAVTK::
set_file (const string& file_name, const string& file_path)
{
  file_name_ = file_name;
  file_path_ = file_path;
  // TODO: Check here that the file exists.
};



void
IGAVTK::
set_visualization_element_properties (const int* const num_visualization_elements,
                                      const int& grid_type)
{
  for (int dir = 0; dir < num_visualization_elements_.size(); ++dir)
    num_visualization_elements_[dir] = *(num_visualization_elements + dir);

  // Grid type 0 : Unstructured grid : quadratic elements.
  // Grid type 1 : Unstructured grid : linear elements.
  // Grid type 2 : Structured grid.
  Assert (grid_type >= 0 && grid_type <= 2, ExcIndexRange(grid_type, 0, 3));

  if (grid_type == 2) // Structured grid.
  {
    unstructured_grid_ = false;
    quadratic_cells_ = false;
  }
  else  // Unstructured grid.
  {
    unstructured_grid_ = true;
    if (grid_type == 0)
      quadratic_cells_ = true;
    else
      quadratic_cells_ = false;
  }
};



void
IGAVTK::
clear ()
{
  AssertThrow (false, ExcNotImplemented ());
};



void
IGAVTK::
generate_vtk_grids(const bool& create_control_mesh,
                   const bool& create_knot_mesh,
                   const bool& create_parametric_mesh,
                   const bool& create_physical_mesh,
                   vtkMultiBlockDataSet* const mb) const
{
  Assert (file_name_ != "", ExcMessage ("Not specified file name."));
  Assert (file_path_ != "", ExcMessage ("Not specified file path."));

  Assert (num_visualization_elements_[0] > 0,
          ExcMessage ("Number of visualization elements must be > 0 in every "
                      "direction."));
  Assert (num_visualization_elements_[1] > 0,
          ExcMessage ("Number of visualization elements must be > 0 in every "
                      "direction."));
  Assert (num_visualization_elements_[2] > 0,
          ExcMessage ("Number of visualization elements must be > 0 in every "
                      "direction."));

  const auto n_functions = this->get_number_functions();
  const Size &n_identity_funcs = n_functions[0];
  const Size &n_not_identity_funcs = n_functions[1];
  const Size &n_ig_funcs = n_functions[2];

  unsigned int block_index = 0;
  if (create_parametric_mesh)
  {
    vtkMultiBlockDataSet* block =
        vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index));
    this->fill_vtk_grids(block, n_identity_funcs, 0,
                         true, false, create_knot_mesh);
    ++block_index;
  }

  if (create_physical_mesh)
  {
    vtkMultiBlockDataSet* block =
        vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index));
    this->fill_vtk_grids(block, n_not_identity_funcs, n_ig_funcs,
                         false, create_control_mesh, create_knot_mesh);
    ++block_index;
  }
};



void
IGAVTK::
fill_vtk_grids(vtkMultiBlockDataSet *const mb,
               const Size &n_solid_mesh,
               const Size &n_control_mesh,
               const bool is_identity,
               const bool create_control_mesh,
               const bool create_knot_mesh) const
{

  if (n_solid_mesh == 0)
    return;

  Size n_blocks = 1;
  if (create_control_mesh)
    ++n_blocks;
  if (create_knot_mesh)
    ++n_blocks;
  mb->SetNumberOfBlocks(n_blocks);

  unsigned int block_index = 0;
  auto solid_block = vtkSmartPointer<vtkMultiBlockDataSet>::New ();
  solid_block->SetNumberOfBlocks(n_solid_mesh);
  mb->SetBlock(block_index, solid_block);
  mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), "Solid mesh");
  ++block_index;

  vtkSmartPointer<vtkMultiBlockDataSet> control_block;
  if (create_control_mesh)
  {
    control_block = vtkSmartPointer<vtkMultiBlockDataSet>::New ();
    control_block->SetNumberOfBlocks(n_control_mesh);
    mb->SetBlock(block_index, control_block);
    mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), "Control mesh");
    ++block_index;
  }

  vtkSmartPointer<vtkMultiBlockDataSet> knot_block;
  if (create_knot_mesh)
  {
    knot_block = vtkSmartPointer<vtkMultiBlockDataSet>::New ();
    knot_block->SetNumberOfBlocks(n_solid_mesh);
    mb->SetBlock(block_index, knot_block);
    mb->GetMetaData(block_index)->Set(vtkCompositeDataSet::NAME(), "Knot mesh");
    ++block_index;
  }


  Index solid_mesh_index = 0;
  Index control_mesh_index = 0;

  this->template fill_vtk_grid<1, 0>(mb, is_identity,
                                     create_control_mesh, create_knot_mesh,
                                     solid_mesh_index, control_mesh_index);
  this->template fill_vtk_grid<2, 0>(mb, is_identity,
                                     create_control_mesh, create_knot_mesh,
                                     solid_mesh_index, control_mesh_index);
  this->template fill_vtk_grid<1, 1>(mb, is_identity,
                                     create_control_mesh, create_knot_mesh,
                                     solid_mesh_index, control_mesh_index);
  this->template fill_vtk_grid<3, 0>(mb, is_identity,
                                     create_control_mesh, create_knot_mesh,
                                     solid_mesh_index, control_mesh_index);
  this->template fill_vtk_grid<2, 1>(mb, is_identity,
                                     create_control_mesh, create_knot_mesh,
                                     solid_mesh_index, control_mesh_index);
  this->template fill_vtk_grid<1, 2>(mb, is_identity,
                                     create_control_mesh, create_knot_mesh,
                                     solid_mesh_index, control_mesh_index);
};



template <int dim, int codim>
void
IGAVTK::
fill_vtk_grid(vtkMultiBlockDataSet *const mb,
              const bool is_identity,
              const bool create_control_mesh,
              const bool create_knot_mesh,
              Index& solid_mesh_index,
              Index& control_mesh_index) const
{

  const auto mappings = funcs_container_->template get_all_mappings<dim, codim>();

  if (mappings.size () == 0)
    return;

//   using Map = Mapping<dim, codim>;
  using IdFun = IdentityFunction<dim, dim + codim>;
  using MapFunPtr = MapFunPtr_<dim, codim>;

  unsigned int block_index = 0;
  vtkMultiBlockDataSet* solid_block =
    vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index));
  ++block_index;

  vtkMultiBlockDataSet* control_block;
  if (create_control_mesh)
  {
    control_block = vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index));
    Assert (control_block != nullptr, ExcNullPtr());
    ++block_index;
  }

  vtkMultiBlockDataSet* knot_block;
  if (create_knot_mesh)
  {
    knot_block = vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index));
    Assert (knot_block != nullptr, ExcNullPtr());
    ++block_index;
  }

  for (const auto &m : mappings)
  {
    const MapFunPtr mapping = m.first;
    const string &name = m.second;

    if ((std::dynamic_pointer_cast<IdFun>(mapping) != nullptr) != is_identity)
      continue;

    solid_block->GetMetaData(solid_mesh_index)->Set(vtkCompositeDataSet::NAME(),
                                                    name.c_str());
    this->generate_solid_mesh_grids<dim, codim>
      (mapping, unstructured_grid_, solid_mesh_index, solid_block);

    if (create_control_mesh)
    {
      control_block->GetMetaData(control_mesh_index)->Set(vtkCompositeDataSet::NAME(),
                                                        name.c_str());
      this->generate_control_mesh_grids<dim, codim>(mapping, control_mesh_index,
                                                    control_block);
      ++control_mesh_index;
    }

    if (create_knot_mesh)
    {
      knot_block->GetMetaData(solid_mesh_index)->Set(vtkCompositeDataSet::NAME(),
                                                     name.c_str());
      this->generate_knot_mesh_grids<dim, codim>(mapping, solid_mesh_index,
                                                 knot_block);
    }

    ++solid_mesh_index;
  }
};



template <int dim, int codim>
void
IGAVTK::
generate_solid_mesh_grids(const MapFunPtr_<dim, codim> mapping,
                          const bool unstructured_grid,
                          const Index& vtk_block_id,
                          vtkMultiBlockDataSet* const vtk_block) const
{
  TensorSize<dim> n_vis_elements;
  for (int dir = 0; dir < dim; ++dir)
    n_vis_elements[dir] = num_visualization_elements_[dir];

  // An implicit cast it is done from
  // vtkSmartPointer<vtkUnstructuredGrid> or vtkSmartPointer<vtkStructuredGrid>
  // to vtkSmartPointer<vtkPointSet>.
  vtkSmartPointer<vtkPointSet> grid;

  if (unstructured_grid) // VTK unstructured grid.
    grid = this->create_solid_vtu_grid<dim, codim>(mapping, n_vis_elements,
                                                   quadratic_cells_);
  else // VTK structured grid.
    grid = this->create_solid_vts_grid<dim, codim>(mapping, n_vis_elements);

  Assert (grid != nullptr, ExcNullPtr());
  vtk_block->SetBlock (vtk_block_id, grid);
};



template <int dim, int codim>
vtkSmartPointer<vtkUnstructuredGrid>
IGAVTK::
create_solid_vtu_grid(const MapFunPtr_<dim, codim> mapping,
                      const TensorSize<dim> &n_vis_elements,
                      const bool quadratic_cells) const
{
  Assert (mapping != nullptr, ExcNullPtr());

#ifndef NDEBUG
  for (int dir = 0; dir < dim; ++dir)
    Assert (n_vis_elements[dir] > 0,
            ExcMessage ("Number of visualization elements must be > 0 in every "
                        "direction."));
#endif

  vtkSmartPointer<vtkUnstructuredGrid> vtu_grid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkPointData* const point_data = vtu_grid->GetPointData();
  const auto points = this->create_points_solid_vtk_grid<dim, codim>
      (mapping, n_vis_elements, false, quadratic_cells, point_data);

  vtu_grid->SetPoints(points);


  const Size n_bezier_elements = mapping->get_grid()->get_num_all_elems ();
  const auto cells = this->template create_cells_solid_vtu_grid<dim>
    (n_vis_elements, n_bezier_elements, quadratic_cells);

  vtu_grid->Allocate(cells->GetNumberOfCells (), 0);
  const int vtk_enum_type = quadratic_cells ?
                              dim == 3 ? VTK_QUADRATIC_HEXAHEDRON :
                              dim == 2 ? VTK_QUADRATIC_QUAD :
                                        VTK_QUADRATIC_EDGE
                            :
                              dim == 3 ? VTK_HEXAHEDRON :
                              dim == 2 ? VTK_QUAD :
                                        VTK_LINE;

  vtu_grid->SetCells(vtk_enum_type, cells);

  return vtu_grid;
};



template <int dim, int codim>
vtkSmartPointer<vtkPointSet>
IGAVTK::
create_solid_vts_grid(const MapFunPtr_<dim, codim> mapping,
                      const TensorSize<dim> &n_vis_elements) const
{
  Assert (mapping != nullptr, ExcNullPtr());

#ifndef NDEBUG
  for (int dir = 0; dir < dim; ++dir)
    Assert (n_vis_elements[dir] > 0,
            ExcMessage ("Number of visualization elements must be > 0 in every "
                        "direction."));
#endif


  if (dim != 1)
  {

    vtkSmartPointer<vtkStructuredGrid> vts_grid =
        vtkSmartPointer<vtkStructuredGrid>::New();

    vtkPointData* const point_data = vts_grid->GetPointData();
    const auto points = this->create_points_solid_vtk_grid<dim, codim>
      (mapping, n_vis_elements, true, false, point_data);

    const auto n_intervals = mapping->get_grid()->get_num_intervals();
    int grid_dim[3] = {1, 1, 1};
    for (int dir = 0; dir < dim; ++dir)
      grid_dim[dir] = (n_vis_elements[dir] + 1) * n_intervals[dir];
    vts_grid->SetDimensions(grid_dim[0], grid_dim[1], grid_dim[2]);

    vts_grid->SetPoints(points);

    return vts_grid;
  }
  else // 1D
  {
    // The 1D structure grid are not visualized in VTK.
    // Therefore, an unstructured grid is created instead.

    auto vtu_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    vtkPointData* const point_data = vtu_grid->GetPointData();
    const auto points = this->create_points_solid_vtk_grid<dim, codim>
      (mapping, n_vis_elements, true, false, point_data);

    vtu_grid->SetPoints(points);

    vtkSmartPointer<vtkPolyLine> polyLine =  vtkSmartPointer<vtkPolyLine>::New();

    const Size n_points = points->GetNumberOfPoints();
    polyLine->GetPointIds()->SetNumberOfIds(n_points);
    for(int i = 0; i < n_points; ++i)
      polyLine->GetPointIds()->SetId(i, i);

    // Create a cell array to store the lines in and add the lines to it
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->InsertNextCell(polyLine);

    vtu_grid->Allocate(cells->GetNumberOfCells (), 0);
    const int vtk_enum_type = VTK_POLY_LINE;
    vtu_grid->SetCells(vtk_enum_type, cells);

    return vtu_grid;
  }
};



template <int dim, int codim>
vtkSmartPointer<vtkPoints>
IGAVTK::
create_points_solid_vtk_grid(const MapFunPtr_<dim, codim> mapping,
                             const TensorSize<dim> &n_vis_elements,
                             const bool structured_grid,
                             const bool quadratic_cells,
                             vtkPointData* const point_data) const
{
  Assert (mapping != nullptr, ExcNullPtr());

  Assert (!(structured_grid && quadratic_cells),
          ExcMessage("Not possible to create quadratic cells in VTK structured "
                     "grid."));

#ifndef NDEBUG
  for (int dir = 0; dir < dim; ++dir)
    Assert (n_vis_elements[dir] > 0,
            ExcMessage ("Number of visualization elements must be > 0 in every "
                        "direction."));
#endif

  static const int space_dim = dim + codim;

  const auto quad = IGAVTK::create_visualization_quadrature<dim>
      (n_vis_elements, quadratic_cells);
  const auto cartesian_grid = mapping->get_grid();

  SafeSTLVector<SafeSTLVector<Index>> points_map;
  SafeSTLVector<Index> points_mask;

  Size n_total_points{0};
  IGAVTK::create_points_numbering_map<dim>(cartesian_grid, quad,
      structured_grid, quadratic_cells, points_map, points_mask, n_total_points);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints (n_total_points);

  double point_tmp[3] = {0.0, 0.0, 0.0};

  auto flag = ValueFlags::point | ValueFlags::value;
  mapping->reset(flag, *quad);

  auto elem = mapping->begin();
  const auto end = mapping->end();

  const auto topology = Topology<dim>();
  mapping->init_cache(elem, topology);

  auto pm_el = points_map.cbegin();
  for (; elem != end; ++elem, ++pm_el)
  {
    mapping->fill_cache(elem, topology, 0);

    auto element_vertices_tmp = elem->template get_values<_Value, dim>(0);
    auto pm = pm_el->cbegin();
    for (const auto &mask : points_mask)
    {
      const auto &point = element_vertices_tmp[mask];
      for (int dir = 0; dir < space_dim ; ++dir)
        point_tmp[dir] = point[dir];
      points->SetPoint (*pm++, point_tmp);
    }
  }

  this->create_point_data<dim, codim>
    (mapping, *quad, points_map, points_mask, point_data);

  return points;
};



template <int dim, int codim>
void
IGAVTK::
generate_knot_mesh_grids(const MapFunPtr_<dim, codim> mapping,
                         const Index& vtk_block_id,
                         vtkMultiBlockDataSet* const vtk_block) const
{
  static const int space_dim = dim + codim;

  const auto cartesian_grid = mapping->get_grid();
  const auto& n_intervals = cartesian_grid->get_num_intervals();

  if (dim == 1)
  {
    const QUniform<dim> quad (2);
    const Size n_vtk_points = n_intervals[0] + 1;
    const Size n_vtk_cells = n_vtk_points;

    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
    vtk_points->SetNumberOfPoints (n_vtk_points);

    vtkSmartPointer<vtkIdTypeArray> vtk_cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();
    const Size tuple_size = 2;
    vtk_cell_ids->SetNumberOfComponents (tuple_size);
    vtk_cell_ids->SetNumberOfTuples (n_vtk_cells);

    auto flag = ValueFlags::point | ValueFlags::value;
    auto elem = mapping->begin();
    auto end  = mapping->end();
    const auto topology = Topology<dim>();

    double point_tmp[3] = {0.0, 0.0, 0.0};
    vtkIdType tuple[2] = {1, 0};

    mapping->reset(flag, quad);
    mapping->init_cache(elem, topology);

    Index pt_id = 0;
    for (; elem != end; ++elem)
    {
      mapping->fill_cache(elem, topology, 0);
      auto points = elem->template get_values<_Value, dim>(0);
      const auto &pp = points[0];
      for (int dir = 0; dir < space_dim; ++dir)
        point_tmp[dir] = pp[dir];
      vtk_points->SetPoint (pt_id, point_tmp);
      tuple[1] = pt_id;
      vtk_cell_ids->SetTupleValue(pt_id, tuple);
      ++pt_id;
    }
    auto points = elem->template get_values<_Value, dim>(0);
    const auto &pp = points[1];
    for (int dir = 0; dir < space_dim; ++dir)
      point_tmp[dir] = pp[dir];
    tuple[1] = pt_id;
    vtk_cell_ids->SetTupleValue(pt_id, tuple);
    vtk_points->SetPoint (pt_id, point_tmp);

    // Creating grid.
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(vtk_points);

    // Creating cells.
    vtkSmartPointer<vtkCellArray> vtk_cells = vtkSmartPointer<vtkCellArray>::New();
    vtk_cells->SetCells (n_vtk_cells, vtk_cell_ids);
    grid->Allocate(vtk_cells->GetNumberOfCells (), 0);
    const int vtk_enum_type = VTK_VERTEX;
    grid->SetCells(vtk_enum_type, vtk_cells);

    vtk_block->SetBlock (vtk_block_id, grid);

    return;
  }

  // TODO: improve performance.

  AssertThrow (dim != 1, ExcNotImplemented());

  const Size n_points_per_single_cell = quadratic_cells_ ? 3 : 2;
  const SafeSTLVector<int> vtk_line_connectivity = quadratic_cells_ ?
    SafeSTLVector<int>({{0, 2, 1}}) : SafeSTLVector<int>({{0, 1}});


  SafeSTLArray<shared_ptr<Quadrature<1>>, dim> quadratures;
  TensorSize<1> n_vis_elems;
  for (int dir = 0; dir < dim; ++dir)
  {
    n_vis_elems[0] = num_visualization_elements_[dir];
    quadratures[dir] = IGAVTK::create_visualization_quadrature<1>(n_vis_elems, quadratic_cells_);
  }

  SafeSTLVector<Real> zero_vec{{0.0}};
  SafeSTLVector<Real> one_vec{{1.0}};

  // Computing total number of vtk points and cells.
  Size n_vtk_cells = 0;
  Size n_vtk_points = 0;
  for (int dir = 0; dir < dim; ++dir)
  {
    Size n_knot_lines = 1;
    for (int dir2 = 0; dir2 < dim; ++dir2)
    {
      if (dir2 != dir)
      {
        n_knot_lines *= (n_intervals[dir2] + 1);
      }
    }
    const Size n_cells_in_knot_line = n_intervals[dir] * num_visualization_elements_[dir];
    n_vtk_cells += n_knot_lines * n_cells_in_knot_line;
    n_vtk_points += n_knot_lines * (n_cells_in_knot_line * (n_points_per_single_cell - 1) + 1);
  }

  vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
  vtk_points->SetNumberOfPoints (n_vtk_points);

  vtkSmartPointer<vtkIdTypeArray> vtk_cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();
  const Size tuple_size = n_points_per_single_cell + 1;
  vtk_cell_ids->SetNumberOfComponents (tuple_size);
  vtk_cell_ids->SetNumberOfTuples (n_vtk_cells);

  Index vtk_pt_id = 0;
  double point_tmp[3] = {0.0, 0.0, 0.0};

  Index vtk_tuple_id = 0;
  vtkIdType* tuple = new vtkIdType[tuple_size];
  tuple[0] = n_points_per_single_cell;

  for (int dir = 0; dir < dim; ++dir) // Creating knot lines along each direction.
  {
    const Index face_id = dir * 2;
    auto &k_elem = UnitElement<dim>::template get_elem<dim-1>(face_id);
    const auto &active_directions = k_elem.active_directions;

    auto flag = ValueFlags::point | ValueFlags::value;
    auto elem = mapping->begin();
    const auto topology = Topology<dim>();

    // Looping along all the knot coordinates of the face.
    using InterGridMap = typename CartesianGrid<dim>::template InterGridMap<dim-1>;
    InterGridMap elem_map;
    const auto &sub_grid = cartesian_grid->template get_sub_grid<dim-1>(face_id, elem_map);
    const auto &face_coords_tensor = sub_grid->get_knot_coordinates();
    const Size n_pts_face = face_coords_tensor.flat_size();

    const auto &quad_dir = quadratures[dir];
    TensorSize<dim> n_quad_points(1);
    n_quad_points[dir] = quad_dir->get_num_points();
    TensorProductArray<dim> quad_points_1d(n_quad_points);
    TensorProductArray<dim> quad_weights_1d(n_quad_points);
    quad_points_1d.copy_data_direction(dir, quad_dir->get_coords_direction(0));

    // Iterating along every knot coordinates of the face
    TensorIndex<dim> elem_t_id(0);
    for (int i_pt = 0; i_pt < n_pts_face; ++i_pt)
    {
      const auto t_id_face = face_coords_tensor.flat_to_tensor(i_pt);
      auto ad = active_directions.cbegin();
      for (const auto &t : t_id_face)
      {
        if (t == n_intervals[*ad])
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
      Quadrature<dim> quad (quad_points_1d, quad_weights_1d);

      mapping->reset(flag, quad);
      mapping->init_cache(elem, topology);
      elem_t_id[dir] = 0;
      for (int itv = 0; itv < n_intervals[dir]; ++itv)
      {
        elem.move_to(cartesian_grid->tensor_to_flat(elem_t_id));
        elem_t_id[dir] += 1;

        mapping->fill_cache(elem, topology, 0);
        auto physical_points = elem->template get_values<_Value, dim>(0);

        Index vtk_pt_id_0 = vtk_pt_id;
        for (const auto &pp : physical_points)
        {
          for (int dir2 = 0; dir2 < space_dim; ++dir2)
            point_tmp[dir2] = pp[dir2];
          vtk_points->SetPoint (vtk_pt_id++, point_tmp);
        } // pp
        --vtk_pt_id;

        for (int c_id = 0; c_id < num_visualization_elements_[dir];  ++c_id,
              vtk_pt_id_0 += (n_points_per_single_cell-1))
        {
          Index id = 1;
          for (const auto &c : vtk_line_connectivity)
            tuple[id++] = vtk_pt_id_0 + c;

          vtk_cell_ids->SetTupleValue(vtk_tuple_id++, tuple);
        } // c_id
      } // itv
      ++vtk_pt_id;
    } // i_pt

  } // dir

  // Creating grid.
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(vtk_points);

  // Creating cells.
  vtkSmartPointer<vtkCellArray> vtk_cells = vtkSmartPointer<vtkCellArray>::New();
  vtk_cells->SetCells (n_vtk_cells, vtk_cell_ids);
  grid->Allocate(vtk_cells->GetNumberOfCells (), 0);
  const int vtk_enum_type = quadratic_cells_ ? VTK_QUADRATIC_EDGE : VTK_LINE;
  grid->SetCells(vtk_enum_type, vtk_cells);

  vtk_block->SetBlock (vtk_block_id, grid);
};



template <int dim, int codim>
void
IGAVTK::
generate_control_mesh_grids(const MapFunPtr_<dim, codim> mapping,
                            const Index& vtk_block_id,
                            vtkMultiBlockDataSet* const vtk_block) const
{
  static const int space_dim = dim + codim;
  using IgFun_ = IgFunction<dim, 0, space_dim, 1>;
  const auto ig_func = std::dynamic_pointer_cast<IgFun_>(mapping);
  Assert (ig_func != nullptr, ExcNullPtr());

  const auto space = ig_func->get_ig_space();
  const auto& coefs = ig_func->get_coefficients();
  const auto& dofs = space->get_dof_distribution();

  // TODO: include assert here to verify that all the components share the
  // same space.
  const auto &dofs_table = dofs->get_num_dofs_table();
  const auto n_pts_dir = dofs_table[0];

  auto points = vtkSmartPointer<vtkPoints>::New();

  const Size n_total_pts = dofs_table.get_component_size(0);
  points->SetNumberOfPoints (n_total_pts);

  // Setting all the points void.
  const double point_void[3] = {0.0, 0.0, 0.0};
  for (int i_pt = 0; i_pt < n_total_pts; ++i_pt)
      points->SetPoint (i_pt, point_void);

  Index comp;
  Index local_id;
  double point_tmp[3];
  for (const auto &it : coefs)
  {
    dofs->global_to_comp_local(it.first, comp, local_id);
    points->GetPoint (local_id, point_tmp);
    point_tmp[comp] = it.second;
    points->SetPoint (local_id, point_tmp);
  }

  if (dim == 1)
  {
    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    grid->SetPoints(points);

    vtkSmartPointer<vtkPolyLine> polyLine =  vtkSmartPointer<vtkPolyLine>::New();

    const Size n_points = points->GetNumberOfPoints();
    polyLine->GetPointIds()->SetNumberOfIds(n_points);
    for(int i = 0; i < n_points; ++i)
      polyLine->GetPointIds()->SetId(i, i);

    // Create a cell array to store the lines in and add the lines to it
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->InsertNextCell(polyLine);

    grid->Allocate(cells->GetNumberOfCells (), 0);
    const int vtk_enum_type = VTK_POLY_LINE;
    grid->SetCells(vtk_enum_type, cells);

    vtk_block->SetBlock (vtk_block_id, grid);
  }
  else
  {
    double grid_dim[3] = {1, 1, 1};
    for (int dir = 0; dir < dim; ++dir)
      grid_dim[dir] = n_pts_dir[dir];;

    auto grid = vtkSmartPointer<vtkStructuredGrid>::New();
    grid->SetDimensions(grid_dim[0], grid_dim[1], grid_dim[2]);
    grid->SetPoints(points);
    vtk_block->SetBlock (vtk_block_id, grid);
  }

};



void
IGAVTK::
parse_file ()
{
  // TODO: this is a temporary solution for avoiding runtime error in deserialization.
//   ifstream xml_istream(file_name_);
//   IArchive xml_in(xml_istream);
//   xml_in >> BOOST_SERIALIZATION_NVP (funcs_container_);
//   xml_istream.close();
//   Assert here if funcs_container_ is void.

  funcs_container_ = std::make_shared<FunctionsContainer>();
  this->template create_geometries<1>();
  this->template create_geometries<2>();
  this->template create_geometries<3>();
};



SafeSTLArray<Size, 3>
IGAVTK::
get_number_functions () const
{
  SafeSTLArray<Size, 3> number_functions(0);
  Size &n_identity = number_functions[0];
  Size &n_not_identity = number_functions[1];
  Size &n_ig_functions = number_functions[2];

  for (const auto& m : funcs_container_->template get_all_mappings<1, 0>())
    if (std::dynamic_pointer_cast<IdentityFunction<1, 1>>(m.first) != nullptr)
      ++n_identity;
    else
    {
      ++n_not_identity;
      if (std::dynamic_pointer_cast<IgFunction<1, 0, 1, 1>>(m.first) != nullptr)
        ++n_ig_functions;
    }

  for (const auto& m : funcs_container_->template get_all_mappings<2, 0>())
    if (std::dynamic_pointer_cast<IdentityFunction<2, 2>>(m.first) != nullptr)
      ++n_identity;
    else
    {
      ++n_not_identity;
      if (std::dynamic_pointer_cast<IgFunction<2, 0, 2, 1>>(m.first) != nullptr)
        ++n_ig_functions;
    }

  for (const auto& m : funcs_container_->template get_all_mappings<1, 1>())
    if (std::dynamic_pointer_cast<IdentityFunction<1, 2>>(m.first) != nullptr)
      ++n_identity;
    else
    {
      ++n_not_identity;
      if (std::dynamic_pointer_cast<IgFunction<1, 0, 2, 1>>(m.first) != nullptr)
        ++n_ig_functions;
    }

  for (const auto& m : funcs_container_->template get_all_mappings<3, 0>())
    if (std::dynamic_pointer_cast<IdentityFunction<3, 3>>(m.first) != nullptr)
      ++n_identity;
    else
    {
      ++n_not_identity;
      if (std::dynamic_pointer_cast<IgFunction<3, 0, 3, 1>>(m.first) != nullptr)
        ++n_ig_functions;
    }

  for (const auto& m : funcs_container_->template get_all_mappings<2, 1>())
    if (std::dynamic_pointer_cast<IdentityFunction<2, 3>>(m.first) != nullptr)
      ++n_identity;
    else
    {
      ++n_not_identity;
      if (std::dynamic_pointer_cast<IgFunction<2, 0, 3, 1>>(m.first) != nullptr)
        ++n_ig_functions;
    }

  for (const auto& m : funcs_container_->template get_all_mappings<1, 2>())
    if (std::dynamic_pointer_cast<IdentityFunction<1, 3>>(m.first) != nullptr)
      ++n_identity;
    else
    {
      ++n_not_identity;
      if (std::dynamic_pointer_cast<IgFunction<1, 0, 3, 1>>(m.first) != nullptr)
        ++n_ig_functions;
    }

  return number_functions;
};



template <>
void
IGAVTK::
create_VTK_quadratic_element_connectivity<1> (
     const TensorSize<1> &n_vis_elements,
     SafeSTLVector<SafeSTLVector<Index>> &connectivity,
     Size &n_points_per_bezier_element)
{
  //  This is the connectivity pattern of the 1D VTK quadratic element.
  //
  //         0 -- 2 -- 1

  // Number of cells per Bezier element.
  const Size n_cells_per_bezier = n_vis_elements.flat_size();

  connectivity.clear();
  connectivity.resize(n_cells_per_bezier);

  Index point_id = 0;
  for (auto &conn_el : connectivity)
  {
    conn_el = {point_id, point_id+2, point_id+1};
    point_id += 2;
  }

  n_points_per_bezier_element = n_vis_elements[0] * 2 + 1;
};



template <>
void
IGAVTK::
create_VTK_quadratic_element_connectivity<2> (
     const TensorSize<2> &n_vis_elements,
     SafeSTLVector<SafeSTLVector<Index>> &connectivity,
     Size &n_points_per_bezier_element)
{
  //  This is the connectivity pattern of the 2D VTK quadratic element.
  //
  //         3 -- 6 -- 2
  //         |         |
  //         7         5
  //         |         |
  //         0 -- 4 -- 1

  static const int dim = 2;

  // Number of cells per Bezier element.
  const Size n_cells_per_bezier = n_vis_elements.flat_size();

  connectivity.clear();
  connectivity.resize(n_cells_per_bezier);

    TensorSize<dim> n_full_points_per_dir;
    for (int dir = 0; dir < dim; ++dir)
      n_full_points_per_dir[dir] = n_vis_elements[dir] * 2 + 1;

  static const int n_points_per_single_cell = 8;

  TensorSize<dim> n_elem_bound_per_dir;
  for (int dir = 0; dir < dim; ++dir)
    n_elem_bound_per_dir[dir] = n_vis_elements[dir] + 1;

  // This grid is going to help in building the connectivity.
  // Every element of the grid refers to a cell.
  const auto cells_grid = CartesianGrid<dim>::create (n_elem_bound_per_dir);


  // This array constaints the offsets of the points along the
  // first (u) direction.
  //   Offset for the line:  3 -- 6 -- 2 -> offset = 2
  //   Offset for the line:  7         5 -> offset = 1
  //   Offset for the line:  0 -- 4 -- 1 -> offset = 2

  const SafeSTLArray<Index, 8> offsets_u = {{2, 2, 2, 2, 2, 1, 2, 1}};

  // This container if for the offsets of the points along the
  // second (v) direction.
  const Index offsets_v = n_full_points_per_dir[0] * 2 - n_vis_elements[0];

  SafeSTLVector<Index> vtk_vertex_id_0(n_points_per_single_cell);
  vtk_vertex_id_0[0] = 0;
  vtk_vertex_id_0[4] = 1;
  vtk_vertex_id_0[1] = 2;

  vtk_vertex_id_0[7] = n_full_points_per_dir[0];
  vtk_vertex_id_0[5] = vtk_vertex_id_0[7] + 1;

  vtk_vertex_id_0[3] = vtk_vertex_id_0[7] + n_full_points_per_dir[0] - n_vis_elements[0];
  vtk_vertex_id_0[6] = vtk_vertex_id_0[3] + 1;
  vtk_vertex_id_0[2] = vtk_vertex_id_0[3] + 2;

  auto conn_el = connectivity.begin ();
  auto cell = cells_grid->begin();
  const auto end = cells_grid->end();
  for (; cell != end; ++cell, ++conn_el)
  {
    conn_el->resize(n_points_per_single_cell);

    SafeSTLArray<Index,dim> vtk_elem_tensor_idx = cell->get_tensor_index();
    auto conn = conn_el->begin();
    for (int i_pt = 0; i_pt < n_points_per_single_cell; ++i_pt, ++conn)
      *conn = vtk_vertex_id_0[i_pt]
            + offsets_u[i_pt] * vtk_elem_tensor_idx[0]
            + offsets_v       * vtk_elem_tensor_idx[1];
  }

  n_points_per_bezier_element =
    (n_vis_elements[0] * 2 + 1) * (n_vis_elements[1] * 2 + 1)  - n_cells_per_bezier;
};



template <>
void
IGAVTK::
create_VTK_quadratic_element_connectivity<3> (
     const TensorSize<3> &n_vis_elements,
     SafeSTLVector<SafeSTLVector<Index>> &connectivity,
     Size &n_points_per_bezier_element)
{
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

  static const int dim = 3;

  // Number of cells per Bezier element.
  const Size n_cells_per_bezier = n_vis_elements.flat_size();

  connectivity.clear();
  connectivity.resize(n_cells_per_bezier);

    TensorSize<dim> n_full_points_per_dir;
    for (int dir = 0; dir < dim; ++dir)
      n_full_points_per_dir[dir] = n_vis_elements[dir] * 2 + 1;

  static const int n_points_per_single_cell = 20;

  TensorSize<dim> n_elem_bound_per_dir;
  for (int dir = 0; dir < dim; ++dir)
    n_elem_bound_per_dir[dir] = n_vis_elements[dir] + 1;

  // This grid is going to help in building the connectivity.
  // Every element of the grid refers to a cell.
  const auto cells_grid = CartesianGrid<dim>::create (n_elem_bound_per_dir);

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

  const SafeSTLArray<Index, 20> offsets_u =
    {{2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1}};

  TensorSize <dim> n_full_points_per_dir_red;
  for (int dir = 0; dir < dim; ++dir)
    n_full_points_per_dir_red[dir] = n_full_points_per_dir[dir] - n_vis_elements[dir];

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

  const Index ov0 = n_full_points_per_dir[0] + n_full_points_per_dir_red[0];
  const Index ov1 = n_full_points_per_dir_red[0];
  const SafeSTLArray<Index, 20> offsets_v =
    {{ov0, ov0, ov0, ov0, ov0, ov0, ov0, ov0, ov0, ov0, ov0, ov0, ov0, ov0,
      ov0, ov0, ov1, ov1, ov1, ov1}};

  // Offset along the third (w) direction.
  const Size ow0 = n_full_points_per_dir[0] * n_full_points_per_dir_red[1]
                  + n_full_points_per_dir_red[0] * n_vis_elements[1];
  const Size ow1 = n_full_points_per_dir_red[0] * n_full_points_per_dir_red[1];

  const Index offsets_w = ow0 + ow1;

  SafeSTLVector<Index> vtk_vertex_id_0(n_points_per_single_cell);
  vtk_vertex_id_0[0]  = 0;
  vtk_vertex_id_0[8]  = 1;
  vtk_vertex_id_0[1]  = 2;
  vtk_vertex_id_0[11] = vtk_vertex_id_0[0]  + n_full_points_per_dir[0];
  vtk_vertex_id_0[9]  = vtk_vertex_id_0[11] + 1;
  vtk_vertex_id_0[3]  = vtk_vertex_id_0[11] + n_full_points_per_dir_red[0];
  vtk_vertex_id_0[10] = vtk_vertex_id_0[3]  + 1;
  vtk_vertex_id_0[2]  = vtk_vertex_id_0[3]  + 2;

  vtk_vertex_id_0[16]  = vtk_vertex_id_0[0]  + ow0;
  vtk_vertex_id_0[17]  = vtk_vertex_id_0[16] + 1;
  vtk_vertex_id_0[19]  = vtk_vertex_id_0[16] + n_full_points_per_dir_red[0];
  vtk_vertex_id_0[18]  = vtk_vertex_id_0[19] + 1;

  vtk_vertex_id_0[4]  = vtk_vertex_id_0[16] + ow1;
  vtk_vertex_id_0[12] = vtk_vertex_id_0[4]  + 1;
  vtk_vertex_id_0[5]  = vtk_vertex_id_0[4]  + 2;
  vtk_vertex_id_0[15] = vtk_vertex_id_0[4]  + n_full_points_per_dir[0];
  vtk_vertex_id_0[13] = vtk_vertex_id_0[15] + 1;
  vtk_vertex_id_0[7]  = vtk_vertex_id_0[15] + n_full_points_per_dir_red[0];
  vtk_vertex_id_0[14] = vtk_vertex_id_0[7]  + 1;
  vtk_vertex_id_0[6]  = vtk_vertex_id_0[7]  + 2;

  auto conn_el = connectivity.begin ();
  auto cell = cells_grid->begin();
  const auto end = cells_grid->end();
  for (; cell != end; ++cell, ++conn_el)
  {
    conn_el->resize(n_points_per_single_cell);

    SafeSTLArray<Index,dim> vtk_elem_tensor_idx = cell->get_tensor_index();
    auto conn = conn_el->begin();
    for (int i_pt = 0; i_pt < n_points_per_single_cell; ++i_pt, ++conn)
      *conn = vtk_vertex_id_0[i_pt]
            + offsets_u[i_pt] * vtk_elem_tensor_idx[0]
            + offsets_v[i_pt] * vtk_elem_tensor_idx[1]
            + offsets_w       * vtk_elem_tensor_idx[2];
  }

  n_points_per_bezier_element =
    (n_vis_elements[0] * 2 + 1) * (n_vis_elements[1] * 2 + 1) 
    * (n_vis_elements[2] * 2 + 1) - n_cells_per_bezier * 4
    - n_vis_elements[0] * n_vis_elements[1]
    - n_vis_elements[0] * n_vis_elements[2]
    - n_vis_elements[1] * n_vis_elements[2];
};



template <int dim>
vtkSmartPointer<vtkCellArray>
IGAVTK::
create_cells_solid_vtu_grid (const TensorSize<dim> &n_vis_elements,
                             const Size &n_bezier_elements,
                             const bool quadratic_cells) const
{
  AssertThrow (n_bezier_elements > 0, ExcMessage("0 Bezier elements found."));

#ifndef NDEBUG
  for (int dir = 0; dir < dim; ++dir)
    Assert (n_vis_elements[dir] > 0,
            ExcMessage ("Number of visualization elements must be > 0 in every "
                        "direction."));
#endif

  // Number of cells per Bezier element.
  const Size n_cells_per_bezier = n_vis_elements.flat_size();

  // These variables are filled inside the if-else blocks --------------------//
  // Total number of VTK points into a single Bezier element.
  Size n_points_per_bezier_element = 0;
  // Connectivity of the VTK cells referrered to the VTK points numbering.
  SafeSTLVector<SafeSTLVector<Index>> connectivity_base;
  //--------------------------------------------------------------------------//


  if (quadratic_cells) // VTK quadratic cells
  {
    IGAVTK::create_VTK_quadratic_element_connectivity<dim>
      (n_vis_elements, connectivity_base, n_points_per_bezier_element);
  } // VTK quadratic cells

  else // VTK linear cells
  {
    // Number of vertices in a dim-dimensional square.
    static constexpr int n_vertices = UnitElement<dim>::template num_elem<0>();

    TensorSize<dim> n_elem_bound_per_dir;
    for (int dir = 0; dir < dim; ++dir)
      n_elem_bound_per_dir[dir] = n_vis_elements[dir] + 1;

    // This grid is going to help in building the connectivity.
    // Every element of the grid refers to a cell.
    const auto cells_grid = CartesianGrid<dim>::create (n_elem_bound_per_dir);

    connectivity_base.resize(n_cells_per_bezier);

    n_points_per_bezier_element = n_elem_bound_per_dir.flat_size ();
    const Size n_points_per_single_cell = n_vertices;

    // Creating the connectivity base ----------------------------------------//

    // Building the offsets container. According to the vtk elements connectivity.
    using T_ = SafeSTLArray < SafeSTLArray<int, dim>, n_vertices>;
    const T_ delta_idx =
      dim == 1 ? T_({{0},       {1}})                          :   // dim = 1
      dim == 2 ? T_({{0, 0},    {1, 0},    {1, 1},    {0, 1}}) :   // dim = 2
                T_({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}); // dim = 3

    const TensorIndex<dim> weight_points =
      MultiArrayUtils< dim >::compute_weight(n_elem_bound_per_dir);

    TensorIndex<dim> vtk_vertex_tensor_idx;

    auto conn_el = connectivity_base.begin ();
    auto cell = cells_grid->begin();
    const auto end = cells_grid->end();
    for (; cell != end; ++cell, ++conn_el)
    {
      conn_el->resize(n_points_per_single_cell);

      SafeSTLArray<Index,dim> vtk_elem_tensor_idx = cell->get_tensor_index();

      auto conn = conn_el->begin ();
      for (int iVertex = 0; iVertex < n_points_per_single_cell; ++iVertex, ++conn)
      {
          for (int i = 0; i < dim; ++i)
            vtk_vertex_tensor_idx[i] = vtk_elem_tensor_idx[i]
                                      + delta_idx[iVertex][i];

          *conn = MultiArrayUtils<dim>::tensor_to_flat_index
            (vtk_vertex_tensor_idx, weight_points);
      }
    }
    //--------------------------------------------------------------------------//

  } // VTK linear cells


  // Total number of cells.
  const int n_total_cells = n_bezier_elements * n_cells_per_bezier;

  auto cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();

  const Size n_points_per_single_cell = connectivity_base[0].size();
  const Size tuple_size = n_points_per_single_cell + 1;
  cell_ids->SetNumberOfComponents (tuple_size);
  cell_ids->SetNumberOfTuples (n_total_cells);

  vtkIdType* tuple = new vtkIdType[tuple_size];
  tuple[0] = tuple_size - 1;

  Index cell_id = 0;
  for (int be = 0; be < n_bezier_elements; ++be)
  {
    const int vtk_vertex_id_offset = n_points_per_bezier_element * be;

    // Iterating along the cells of one Bezier element.
    for (const auto& cc : connectivity_base)
    {
      int point_id = 1;
      for (const auto& c : cc)
        tuple[point_id++] = c + vtk_vertex_id_offset;

      cell_ids->SetTupleValue(cell_id++, tuple);
    }
  }

  auto cells = vtkSmartPointer<vtkCellArray>::New();
  cells->SetCells (cell_ids->GetNumberOfTuples() , cell_ids);

  return cells;
};



template <>
void
IGAVTK::
create_point_data<1, 0> (const shared_ptr<Function<1, 0, 1, 1>> map,
                         const Quadrature<1> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         const SafeSTLVector<Index>& points_mask,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<1, 0, 1, 1>(map, quad, points_map, points_mask, point_data);
//   this->template create_point_data<1, 0, 2, 1>(map, quad, points_map, points_mask, point_data);
//   this->template create_point_data<1, 0, 3, 1>(map, quad, points_map, points_mask, point_data);
};



template <>
void
IGAVTK::
create_point_data<2, 0> (const shared_ptr<Function<2, 0, 2, 1>> map,
                         const Quadrature<2> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         const SafeSTLVector<Index>& points_mask,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<2, 0, 1, 1>(map, quad, points_map, points_mask, point_data);
  this->template create_point_data<2, 0, 2, 1>(map, quad, points_map, points_mask, point_data);
};



template <>
void
IGAVTK::
create_point_data<1, 1> (const shared_ptr<Function<1, 0, 2, 1>> map,
                         const Quadrature<1> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         const SafeSTLVector<Index>& points_mask,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<1, 1, 1, 1>(map, quad, points_map, points_mask, point_data);
  this->template create_point_data<1, 1, 2, 1>(map, quad, points_map, points_mask, point_data);
};



template <>
void
IGAVTK::
create_point_data<3, 0> (const shared_ptr<Function<3, 0, 3, 1>> map,
                         const Quadrature<3> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         const SafeSTLVector<Index>& points_mask,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<3, 0, 1, 1>(map, quad, points_map, points_mask, point_data);
  this->template create_point_data<3, 0, 3, 1>(map, quad, points_map, points_mask, point_data);
};



template <>
void
IGAVTK::
create_point_data<2, 1> (const shared_ptr<Function<2, 0, 3, 1>> map,
                         const Quadrature<2> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         const SafeSTLVector<Index>& points_mask,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<2, 1, 1, 1>(map, quad, points_map, points_mask, point_data);
  this->template create_point_data<2, 1, 3, 1>(map, quad, points_map, points_mask, point_data);
};



template <>
void
IGAVTK::
create_point_data<1, 2> (const shared_ptr<Function<1, 0, 3, 1>> map,
                         const Quadrature<1> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         const SafeSTLVector<Index>& points_mask,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<1, 2, 1, 1>(map, quad, points_map, points_mask, point_data);
  this->template create_point_data<1, 2, 3, 1>(map, quad, points_map, points_mask, point_data);
};



template <int dim, int codim, int range, int rank>
void
IGAVTK::
create_point_data (const shared_ptr<Function<dim, 0, dim + codim, 1>> mapping,
                   const Quadrature<dim> &quad,
                   const SafeSTLVector<SafeSTLVector<Index>>& point_num_map,
                   const SafeSTLVector<Index> &points_mask,
                   vtkPointData* const point_data) const
{
  const auto &funcs_map = funcs_container_->template
    get_functions_associated_to_mapping<dim, codim, range, rank>(mapping);

  using Value = typename Function<dim, codim, range, rank>::Value;

  static constexpr Size n_comp = Value::size;
  const Size n_pts = quad.get_num_points();
  Real tuple[n_comp];

  auto flag = ValueFlags::value | ValueFlags::point;
  for (const auto& it : funcs_map)
  {
    const auto &fun = it.first;
    const auto &name = it.second;

    const int n_bezier_elements = point_num_map.size ();
    const vtkIdType n_tuples = n_bezier_elements * n_pts;

    vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();

    arr->SetName(name.c_str());
    arr->SetNumberOfComponents(n_comp);
    arr->SetNumberOfTuples(n_tuples);

    fun->reset(flag, quad);

    auto elem = fun->begin();
    const auto end = fun->end();

    const auto topology = Topology<dim>();
    fun->init_cache(elem, topology);

    auto pnm_it = point_num_map.cbegin();
    for (; elem != end; ++elem, ++pnm_it)
    {
      fun->fill_cache(elem, topology, 0);

      auto pnm = pnm_it->cbegin();
      auto values = elem->template get_values<_Value, dim>(0);
      for (const auto &pm : points_mask)
      {
        this->template tensor_to_tuple<Value>(values[pm], tuple);
        arr->SetTuple(*pnm++, tuple);
      }
    }

    if (point_data->HasArray(name.c_str()));
    {
      // TODO: Throw warning.
    }
    point_data->AddArray(arr.Get());
  }
};



template <int dim>
shared_ptr<Quadrature<dim>>
IGAVTK::
create_visualization_quadrature (const TensorSize<dim>& n_elements_per_direction,
                                 const bool is_quadratic)
{
  const Size n_pts_per_cell_dir = is_quadratic ? 3 : 2;
  TensorSize<dim> n_points;
  for (int dir = 0; dir < dim; ++dir)
  {
    Assert (n_elements_per_direction[dir] > 0, ExcMessage ("Wrong number of visualization elements."));
    n_points[dir] = n_elements_per_direction[dir] * (n_pts_per_cell_dir - 1) + 1;
  }

  QUniform<dim> uniform_quad (n_points);

  if (!is_quadratic)
  {
    return std::make_shared<Quadrature<dim>> (uniform_quad);
  }

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
      if ( (t_id[dir] > 0) && (t_id[dir] < (n_points[dir] - 1)) )
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

  return std::make_shared<Quadrature<dim>> (uniform_quad);
};



template <int dim>
void
IGAVTK::
create_points_numbering_map (const shared_ptr<const CartesianGrid<dim>> cartesian_grid,
                             const shared_ptr<Quadrature<dim>> quad,
                             const bool is_structured,
                             const bool is_quadratic,
                             SafeSTLVector<SafeSTLVector<Index>> &points_map,
                             SafeSTLVector<Index> &points_mask,
                             Size &n_total_points)
{
  points_map.clear();
  points_mask.clear();

  const Size n_bezier_elements = cartesian_grid->get_num_all_elems();

  points_map.resize(n_bezier_elements);

  if (is_structured) // VTK structured grid
  {
    // Total number of visualization points per Bezier element.
    const Size n_pts_per_bezier_elem = quad->get_num_points();
    // Number of visualization points per Bezier element in each direction.
    const auto n_pts_dir_per_bezier_elem = quad->get_num_coords_direction();

    points_mask.resize(n_pts_per_bezier_elem);
    Index id = 0;
    for (auto &pm : points_mask)
      pm = id++;

    n_total_points = n_pts_per_bezier_elem * n_bezier_elements;

    TensorSize<dim>  n_pts_per_mesh;   // Number of points per direction of
                                        // VTK structured grid.
    const auto n_intervals = cartesian_grid->get_num_intervals();
    for (int dir = 0 ; dir < dim ; ++dir)
      n_pts_per_mesh[dir] = n_intervals[dir] * n_pts_dir_per_bezier_elem[dir];

    TensorIndex<dim> elem_t_id;        // Tensorial index of the Bezier element.
    TensorIndex<dim> pt_mesh_t_offset; // Tensorial index of the first point in
                                        // Bezier element.
    TensorIndex<dim> pt_mesh_t_id;     // Tensorial index of the point.
    TensorIndex<dim> pt_elem_t_id;     // Tensorial index of the point referred
                                        // to the number of points in a
                                        // single element.

    const auto w_elem_pts = MultiArrayUtils<dim>::compute_weight(n_pts_dir_per_bezier_elem);
    const auto w_mesh_pts = MultiArrayUtils<dim>::compute_weight(n_pts_per_mesh);

    for (int i_el = 0; i_el < n_bezier_elements; ++i_el)
    {
      auto &pmi = points_map[i_el];
      pmi.resize(n_pts_per_bezier_elem);

      // Computing the tensor index of the first point of the element.
      elem_t_id = cartesian_grid->flat_to_tensor(i_el);
      for (int dir = 0 ; dir < dim ; ++dir)
        pt_mesh_t_offset[dir] = elem_t_id[dir] * n_pts_dir_per_bezier_elem[dir];

      Index point_id = 0;
      for (auto &pm : pmi)
      {
        // Computing the tensor index of the point referred to the number
        // of points into a single Bezier element.
        pt_elem_t_id = MultiArrayUtils<dim>::flat_to_tensor_index(point_id++, w_elem_pts);

        // Computing the tensor index of the point referred to the number
        // of points into the whole mesh.
        for (int dir = 0 ; dir < dim ; ++dir)
          pt_mesh_t_id[dir] = pt_mesh_t_offset[dir] + pt_elem_t_id[dir];

        pm = MultiArrayUtils<dim>::tensor_to_flat_index(pt_mesh_t_id, w_mesh_pts);
      }
    }
  }
  else // VTK unstructured grid
  {
    if (is_quadratic) // VTK quadratic elements
    {
      // Number of visualization points per Bezier element in each direction.
      const auto n_pts_dir_per_bezier_elem = quad->get_num_coords_direction();

      // Iteration along all the quadrature point.
      // Only the points in an edge are added to the mask.
      const auto &quad_points_1d = quad->get_points_1d();
      const Size n_quad_points = quad->get_num_points();
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
          points_mask.push_back(i_pt);
      }

      // Total number of visualization points per Bezier element.
      const Size n_pts_per_bezier_elem = points_mask.size();

      Index point_id = 0;
      for (auto &pm_el : points_map)
      {
        pm_el.resize(n_pts_per_bezier_elem);
        for (auto &pm : pm_el)
          pm = point_id++;
      } // points_map
      n_total_points = point_id;

    } // VTK quadratic elements

    else // VTK linear elements
    {
      // Total number of visualization points per Bezier element.
      const Size n_pts_per_bezier_elem = quad->get_num_points();

      points_mask.resize(n_pts_per_bezier_elem);
      Index id = 0;
      for (auto &pm : points_mask)
        pm = id++;

      Index point_id = 0;
      for (auto &pm_el : points_map)
      {
        pm_el.resize(n_pts_per_bezier_elem);
        for (auto &pm : pm_el)
          pm = point_id++;
      } // points_map
      n_total_points = point_id;

    } // VTK linear elements
  }
};


IGA_NAMESPACE_CLOSE

#include <igatools/io/reader.h>

IGA_NAMESPACE_OPEN

template <int dim>
void
IGAVTK::
create_geometries ()
{
  using std::dynamic_pointer_cast;
  using std::const_pointer_cast;

  static const int codim = dim == 1 ? 1 : 0;
  static const int space_dim = dim + codim;
  static const int range = space_dim;
  static const int rank = 1;
  static const Transformation transf = Transformation::h_grad;

  using Fun_ = Function<dim, 0, range, rank>;
  using FunPhys_ = Function<dim, codim, range, rank>;
  using IgFun_ = IgFunction<dim, 0, range, rank>;
  using IgFunPhys_ = IgFunction<dim, codim, range, rank>;
  using PhysSpace_ = PhysicalSpace<dim, range, rank, codim, transf>;
  using RefSpace_ = ReferenceSpace<dim, range, rank>;
  using IdFun_ = IdentityFunction<dim, dim>;
  using MapFunction_ = Function<dim, 0, space_dim, 1>;
  using Space_ = Space<dim, 0, space_dim, 1>;
  using Grid_ = CartesianGrid<dim>;

  // File names;
  const string fname_0 = "patch_0_" + std::to_string (dim) + "D.xml";
  const string fname_1 = "patch_1_" + std::to_string (dim) + "D.xml";
  const string fname_2 = "patch_2_" + std::to_string (dim) + "D.xml";
  const string fname_3 = "patch_3_" + std::to_string (dim) + "D.xml";

  // Reading maps
  const shared_ptr<MapFunction_> map_0 = get_mapping_from_file<dim, codim> (fname_0);
  const shared_ptr<MapFunction_> map_1 = get_mapping_from_file<dim, codim> (fname_1);
  const shared_ptr<MapFunction_> map_2 = get_mapping_from_file<dim, codim> (fname_2);
  const shared_ptr<MapFunction_> map_3 = get_mapping_from_file<dim, codim> (fname_3);

  // Transforming mapping to ig mapping.
  const shared_ptr<IgFun_> ig_func_0 = dynamic_pointer_cast<IgFun_>(map_0);
  const shared_ptr<IgFun_> ig_func_1 = dynamic_pointer_cast<IgFun_>(map_1);
  const shared_ptr<IgFun_> ig_func_2 = dynamic_pointer_cast<IgFun_>(map_2);
  const shared_ptr<IgFun_> ig_func_3 = dynamic_pointer_cast<IgFun_>(map_3);

  // Getting ig spaces.
  const shared_ptr<const Space_> space_0 = ig_func_0->get_ig_space ();
  const shared_ptr<const Space_> space_1 = ig_func_1->get_ig_space ();
  const shared_ptr<const Space_> space_2 = ig_func_2->get_ig_space ();
  const shared_ptr<const Space_> space_3 = ig_func_3->get_ig_space ();

  // Getting reference spaces.
  const shared_ptr<RefSpace_> ref_space_0 = const_pointer_cast<RefSpace_>(dynamic_pointer_cast<const RefSpace_>(space_0));
  const shared_ptr<RefSpace_> ref_space_1 = const_pointer_cast<RefSpace_>(dynamic_pointer_cast<const RefSpace_>(space_1));
  const shared_ptr<RefSpace_> ref_space_2 = const_pointer_cast<RefSpace_>(dynamic_pointer_cast<const RefSpace_>(space_2));
  const shared_ptr<RefSpace_> ref_space_3 = const_pointer_cast<RefSpace_>(dynamic_pointer_cast<const RefSpace_>(space_3));

  // Building physical spaces.
  const shared_ptr<PhysSpace_> phys_space_0 = PhysSpace_::create (ref_space_0, map_0);
  const shared_ptr<PhysSpace_> phys_space_1 = PhysSpace_::create (ref_space_1, map_1);
  const shared_ptr<PhysSpace_> phys_space_2 = PhysSpace_::create (ref_space_2, map_2);
  const shared_ptr<PhysSpace_> phys_space_3 = PhysSpace_::create (ref_space_3, map_3);

  // Getting grids;
  const shared_ptr<Grid_> grid_0 = ref_space_0->get_grid ();
  const shared_ptr<Grid_> grid_1 = ref_space_1->get_grid ();
  const shared_ptr<Grid_> grid_2 = ref_space_2->get_grid ();
  const shared_ptr<Grid_> grid_3 = ref_space_3->get_grid ();

  // Creating identity functions.
  const auto id_map_0 = IdFun_::create (grid_0);
  const auto id_map_1 = IdFun_::create (grid_1);
  const auto id_map_2 = IdFun_::create (grid_2);
  const auto id_map_3 = IdFun_::create (grid_3);

  Epetra_SerialComm comm;

  // Creating coefficients for ig physical space functions.
  auto phys_coeff_0 = EpetraTools::create_vector(
                        EpetraTools::create_map(phys_space_0, "active", comm));
  const auto dofs_dist_0 = space_0->get_dof_distribution();
  const Real val0 = 10.0;
  Index counter = 0;
  for (const auto& d : dofs_dist_0->get_dofs_view ())
  {
    double val = (counter++) * 1.5 + val0;
    phys_coeff_0->ReplaceGlobalValues (1, &val, &d);
//     phys_coeff_0->ReplaceGlobalValues (1, &val0, &d);
  }

  auto phys_coeff_1 = EpetraTools::create_vector(
                        EpetraTools::create_map(phys_space_1, "active", comm));
  const auto dofs_dist_1 = space_1->get_dof_distribution();
  const Real val1 = 11.0;
  counter = 0;
  for (const auto& d : dofs_dist_1->get_dofs_view ())
  {
    double val = (counter++) * 1.5 + val1;
    phys_coeff_1->ReplaceGlobalValues (1, &val, &d);
//     phys_coeff_1->ReplaceGlobalValues (1, &val1, &d);
  }

  auto phys_coeff_2 = EpetraTools::create_vector(
                        EpetraTools::create_map(phys_space_2, "active", comm));
  const auto dofs_dist_2 = space_2->get_dof_distribution();
  const Real val2 = 12.0;
  counter = 0;
  for (const auto& d : dofs_dist_2->get_dofs_view ())
  {
    double val = (counter++) * 1.5 + val2;
    phys_coeff_2->ReplaceGlobalValues (1, &val, &d);
//     phys_coeff_2->ReplaceGlobalValues (1, &val2, &d);
  }

  auto phys_coeff_3 = EpetraTools::create_vector(
                        EpetraTools::create_map(phys_space_3, "active", comm));
  const auto dofs_dist_3 = space_3->get_dof_distribution();
  const Real val3 = 13.0;
  counter = 0;
  for (const auto& d : dofs_dist_3->get_dofs_view ())
  {
    double val = (counter++) * 1.5 + val3;
    phys_coeff_3->ReplaceGlobalValues (1, &val, &d);
//     phys_coeff_3->ReplaceGlobalValues (1, &val3, &d);
  }


  // Creating coefficients for ig reference space functions.
  auto ref_coeff_0 = EpetraTools::create_vector(
                        EpetraTools::create_map(ref_space_0, "active", comm));
  const Real ref_val0 = 20.0;
  counter = 0;
  for (const auto& d : dofs_dist_0->get_dofs_view ())
  {
    double val = (counter++) * 1.5 + ref_val0;
    ref_coeff_0->ReplaceGlobalValues (1, &val, &d);
//     ref_coeff_0->ReplaceGlobalValues (1, &ref_val0, &d);
  }

  auto ref_coeff_1 = EpetraTools::create_vector(
                        EpetraTools::create_map(ref_space_1, "active", comm));
  const Real ref_val1 = 21.0;
  counter = 0;
  for (const auto& d : dofs_dist_1->get_dofs_view ())
  {
    double val = (counter++) * 1.5 + ref_val1;
    ref_coeff_1->ReplaceGlobalValues (1, &val, &d);
//     ref_coeff_1->ReplaceGlobalValues (1, &ref_val1, &d);
  }

  auto ref_coeff_2 = EpetraTools::create_vector(
                        EpetraTools::create_map(ref_space_2, "active", comm));
  const Real ref_val2 = 22.0;
  counter = 0;
  for (const auto& d : dofs_dist_2->get_dofs_view ())
  {
    double val = (counter++) * 1.5 + ref_val2;
    ref_coeff_2->ReplaceGlobalValues (1, &val, &d);
//     ref_coeff_2->ReplaceGlobalValues (1, &ref_val2, &d);
  }

  auto ref_coeff_3 = EpetraTools::create_vector(
                        EpetraTools::create_map(ref_space_3, "active", comm));
  const Real ref_val3 = 23.0;
  counter = 0;
  for (const auto& d : dofs_dist_3->get_dofs_view ())
  {
    double val = (counter++) * 1.5 + ref_val3;
    ref_coeff_3->ReplaceGlobalValues (1, &val, &d);
//     ref_coeff_3->ReplaceGlobalValues (1, &ref_val3, &d);
  }

  // Creating ig functions for physical spaces.
  auto ps_func_0 = dynamic_pointer_cast<FunPhys_>(IgFunPhys_::create (phys_space_0, phys_coeff_0));
  auto ps_func_1 = dynamic_pointer_cast<FunPhys_>(IgFunPhys_::create (phys_space_1, phys_coeff_1));
  auto ps_func_2 = dynamic_pointer_cast<FunPhys_>(IgFunPhys_::create (phys_space_2, phys_coeff_2));
  auto ps_func_3 = dynamic_pointer_cast<FunPhys_>(IgFunPhys_::create (phys_space_3, phys_coeff_3));

  // Creating ig functions for reference spaces.
  auto rf_func_0 = dynamic_pointer_cast<Fun_>(IgFun_::create (ref_space_0, ref_coeff_0));
  auto rf_func_1 = dynamic_pointer_cast<Fun_>(IgFun_::create (ref_space_1, ref_coeff_1));
  auto rf_func_2 = dynamic_pointer_cast<Fun_>(IgFun_::create (ref_space_2, ref_coeff_2));
  auto rf_func_3 = dynamic_pointer_cast<Fun_>(IgFun_::create (ref_space_3, ref_coeff_3));

  // Adding all the stuff to the functions container.

  // Inserting geometries.
  const string map_name_0 = "map_0_" + std::to_string (dim) + "D";
  const string map_name_1 = "map_1_" + std::to_string (dim) + "D";
  const string map_name_2 = "map_2_" + std::to_string (dim) + "D";
  const string map_name_3 = "map_3_" + std::to_string (dim) + "D";
  funcs_container_->insert_mapping(map_0, map_name_0);
  funcs_container_->insert_mapping(map_1, map_name_1);
  funcs_container_->insert_mapping(map_2, map_name_2);
  funcs_container_->insert_mapping(map_3, map_name_3);

  const string id_map_name_0 = "id_map_0_" + std::to_string (dim) + "D";
  const string id_map_name_1 = "id_map_1_" + std::to_string (dim) + "D";
  const string id_map_name_2 = "id_map_2_" + std::to_string (dim) + "D";
  const string id_map_name_3 = "id_map_3_" + std::to_string (dim) + "D";
  funcs_container_->insert_mapping(id_map_0, id_map_name_0);
  funcs_container_->insert_mapping(id_map_1, id_map_name_1);
  funcs_container_->insert_mapping(id_map_2, id_map_name_2);
  funcs_container_->insert_mapping(id_map_3, id_map_name_3);

  // Inserting associated functions.
  const string fun_map_name_0 = "phys_func_0_" + std::to_string (dim) + "D";
  const string fun_map_name_1 = "phys_func_1_" + std::to_string (dim) + "D";
  const string fun_map_name_2 = "phys_func_2_" + std::to_string (dim) + "D";
  const string fun_map_name_3 = "phys_func_3_" + std::to_string (dim) + "D";
  funcs_container_->insert_function(map_0, ps_func_0, fun_map_name_0);
  funcs_container_->insert_function(map_1, ps_func_1, fun_map_name_1);
  funcs_container_->insert_function(map_2, ps_func_2, fun_map_name_2);
  funcs_container_->insert_function(map_3, ps_func_3, fun_map_name_3);

#if 0 // The combination dim=1, codim=1, range=2, rank=1 is not instantiated.
  const string id_fun_map_name_0 = "ref_func_0_" + std::to_string (dim) + "D";
  const string id_fun_map_name_1 = "ref_func_1_" + std::to_string (dim) + "D";
  const string id_fun_map_name_2 = "ref_func_2_" + std::to_string (dim) + "D";
  const string id_fun_map_name_3 = "ref_func_3_" + std::to_string (dim) + "D";
  funcs_container_->insert_function(id_map_0, rf_func_0, id_fun_map_name_0);
  funcs_container_->insert_function(id_map_1, rf_func_1, id_fun_map_name_1);
  funcs_container_->insert_function(id_map_2, rf_func_2, id_fun_map_name_2);
  funcs_container_->insert_function(id_map_3, rf_func_3, id_fun_map_name_3);
#endif
};

IGA_NAMESPACE_CLOSE