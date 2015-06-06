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

#include <igatools/functions/functions_container.h>
#include <igatools/base/quadrature_lib.h>



using namespace iga;
using std::string;
using std::shared_ptr;
using std::pair;



IGAVTK::
IGAVTK ()
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
set_number_visualization_elements (const int* const num_visualization_elements)
{

  static const int dim = 3;
  for (int dir = 0; dir < dim; ++dir)
    num_visualization_elements_[dir] = *(num_visualization_elements + dir);
};



void
IGAVTK::
clear ()
{
  AssertThrow (false, ExcNotImplemented ());
};



void
IGAVTK::
generate_vtk_grids(const int& grid_type,
                   const bool& create_control_mesh,
                   const bool& create_knot_mesh,
                   const bool& create_parametric_mesh,
                   const bool& create_physical_mesh,
                   vtkMultiBlockDataSet* const mb) const
{
  Assert (file_name_ != "", ExcMessage ("Not specified file name."));
  Assert (file_path_ != "", ExcMessage ("Not specified file path."));
  Assert (num_visualization_elements_[0] > 0,
          ExcMessage ("Number of visualization elements must be > 0 in every direction."));
  Assert (num_visualization_elements_[1] > 0,
          ExcMessage ("Number of visualization elements must be > 0 in every direction."));
  Assert (num_visualization_elements_[2] > 0,
          ExcMessage ("Number of visualization elements must be > 0 in every direction."));

  const bool unstructured = grid_type == 0;

  const auto names = this->get_map_names ();

  unsigned int block_index = 0;
  if (create_control_mesh)
  {
    const auto& n = names[0];
    vtkMultiBlockDataSet* internal_mb =
      vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index++));
    Assert (internal_mb != nullptr, ExcNullPtr ());
    internal_mb->SetNumberOfBlocks(n.size ());

    unsigned int internal_block_index = 0;
    for (const auto& name : n)
      internal_mb->GetMetaData(internal_block_index++)->Set(vtkCompositeDataSet::NAME(), name.c_str());

    internal_block_index = 0;
    this->generate_control_mesh_grids<2, 0>(internal_mb, internal_block_index);
    this->generate_control_mesh_grids<1, 1>(internal_mb, internal_block_index);
    this->generate_control_mesh_grids<3, 0>(internal_mb, internal_block_index);
    this->generate_control_mesh_grids<2, 1>(internal_mb, internal_block_index);
    this->generate_control_mesh_grids<1, 2>(internal_mb, internal_block_index);
  }

  if (create_knot_mesh)
  {
    vtkMultiBlockDataSet* internal_mb =
      vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index++));
    Assert (internal_mb != nullptr, ExcNullPtr ());
    int number_blocks = 0;
    if (create_parametric_mesh)
      number_blocks += names[1].size();
    if (create_physical_mesh)
      number_blocks += names[2].size();

    internal_mb->SetNumberOfBlocks(number_blocks);

    unsigned int internal_block_index = 0;
    if (create_parametric_mesh)
      for (const auto& name : names[1])
        internal_mb->GetMetaData(internal_block_index++)->Set(vtkCompositeDataSet::NAME(), name.c_str());

    if (create_physical_mesh)
      for (const auto& name : names[2])
        internal_mb->GetMetaData(internal_block_index++)->Set(vtkCompositeDataSet::NAME(), name.c_str());

    internal_block_index = 0;
    if (create_parametric_mesh)
    {
      this->generate_knot_mesh_grids<2, 0>(internal_mb, internal_block_index, true);
  //     this->generate_knot_mesh_grids<1, 1>(internal_mb, internal_block_index, true);
      this->generate_knot_mesh_grids<3, 0>(internal_mb, internal_block_index, true);
      this->generate_knot_mesh_grids<2, 1>(internal_mb, internal_block_index, true);
  //     this->generate_knot_mesh_grids<1, 2>(internal_mb, internal_block_index, true);
    }

    if (create_physical_mesh)
    {
      this->generate_knot_mesh_grids<2, 0>(internal_mb, internal_block_index, false);
  //     this->generate_knot_mesh_grids<1, 1>(internal_mb, internal_block_index, false);
      this->generate_knot_mesh_grids<3, 0>(internal_mb, internal_block_index, false);
      this->generate_knot_mesh_grids<2, 1>(internal_mb, internal_block_index, false);
  //     this->generate_knot_mesh_grids<1, 2>(internal_mb, internal_block_index, false);
    }
  }

  if (create_parametric_mesh)
  {
    const auto& n = names[1];
    vtkMultiBlockDataSet* internal_mb =
      vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index++));
    Assert (internal_mb != nullptr, ExcNullPtr ());
    internal_mb->SetNumberOfBlocks(n.size ());

    unsigned int internal_block_index = 0;
    for (const auto& name : n)
      internal_mb->GetMetaData(internal_block_index++)->Set(vtkCompositeDataSet::NAME(), name.c_str());

    internal_block_index = 0;
    const bool is_parametric = true;
    this->generate_solid_mesh_grids<2, 0>(internal_mb, internal_block_index, unstructured, is_parametric);
    this->generate_solid_mesh_grids<1, 1>(internal_mb, internal_block_index, unstructured, is_parametric);
    this->generate_solid_mesh_grids<3, 0>(internal_mb, internal_block_index, unstructured, is_parametric);
    this->generate_solid_mesh_grids<2, 1>(internal_mb, internal_block_index, unstructured, is_parametric);
    this->generate_solid_mesh_grids<1, 2>(internal_mb, internal_block_index, unstructured, is_parametric);
  }

  if (create_physical_mesh)
  {
    const auto& n = names[2];
    vtkMultiBlockDataSet* internal_mb =
      vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index++));
    Assert (internal_mb != nullptr, ExcNullPtr ());
    internal_mb->SetNumberOfBlocks(n.size ());
    const bool is_parametric = false;

    unsigned int internal_block_index = 0;
    for (const auto& name : n)
      internal_mb->GetMetaData(internal_block_index++)->Set(vtkCompositeDataSet::NAME(), name.c_str());

    internal_block_index = 0;
    this->generate_solid_mesh_grids<2, 0>(internal_mb, internal_block_index, unstructured, is_parametric);
    this->generate_solid_mesh_grids<1, 1>(internal_mb, internal_block_index, unstructured, is_parametric);
    this->generate_solid_mesh_grids<3, 0>(internal_mb, internal_block_index, unstructured, is_parametric);
    this->generate_solid_mesh_grids<2, 1>(internal_mb, internal_block_index, unstructured, is_parametric);
    this->generate_solid_mesh_grids<1, 2>(internal_mb, internal_block_index, unstructured, is_parametric);
  }
};




template <int dim, int codim>
void
IGAVTK::
generate_solid_mesh_grids(vtkMultiBlockDataSet* const mb,
                          unsigned int& id,
                          const bool unstructured,
                          const bool is_parametric) const
{
  static const int space_dim = dim + codim;

  const auto mappings = funcs_container_->template get_all_mappings<dim, codim>();

  for (const auto &m : mappings)
  {
      using Fun_ = Function<dim, 0, dim + codim, 1>;
      shared_ptr<Fun_> mapping = m.first;

      if ((std::dynamic_pointer_cast<IdentityFunction<dim, space_dim>>(mapping) != nullptr) != is_parametric)
        continue;

      vtkSmartPointer<vtkPoints>   points = vtkSmartPointer<vtkPoints>::New();

      // Generating visualization quadrature.
      // TODO: quadratures for 1D, 2D and 3D should be moved outside this
      // function.

      TensorSize<dim> n_points;
      for (int dir = 0; dir < dim; ++dir)
      {
        n_points[dir] = num_visualization_elements_[dir] + 1;
        Assert (n_points[dir] > 1, ExcMessage ("Wrong number of visualization points."));
      }
      QUniform<dim> quad (n_points);

      const int n_bezier_elements = mapping->get_grid()->get_num_all_elems();
      const int n_points_per_bezier_element = quad.get_num_points();
      const int total_num_points = n_points_per_bezier_element * n_bezier_elements;

      // Setting the points --------------------------------------------------//
      points->SetNumberOfPoints (total_num_points);

      auto flag = ValueFlags::point | ValueFlags::value;
      mapping->reset(flag, quad);

      auto elem = mapping->begin();
      const auto end = mapping->end();

      const auto topology = Topology<dim>();
      mapping->init_cache(elem, topology);

      const auto point_num_map =
        create_points_numbering_map(mapping->get_grid(),
                                    n_points, unstructured);

      Assert (point_num_map.size() == n_bezier_elements,
              ExcDimensionMismatch(point_num_map.size(), n_bezier_elements));
#ifndef NDEBUG
      for (const auto& it : point_num_map)
        Assert (it.size() == n_points_per_bezier_element,
              ExcDimensionMismatch(it.size(), n_points_per_bezier_element));
#endif

      double point_tmp[3] = {0.0, 0.0, 0.0};
      auto pnm_it = point_num_map.cbegin();
      for (; elem != end; ++elem, ++pnm_it)
      {
        mapping->fill_cache(elem, topology, 0);

        auto element_vertices_tmp = elem->template get_values<_Value, dim>(0);
        auto pnm = pnm_it->cbegin();
        for (const auto& p : element_vertices_tmp)
        {
          for (int dir = 0; dir < dim ; ++dir)
            point_tmp[dir] = p[dir];
          points->SetPoint (*pnm++, point_tmp);
        }
      }
      //----------------------------------------------------------------------//

      if (unstructured) // Creating an unstuctured grid.
      {
        // Setting the cells -------------------------------------------------//
        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        int n_cells_per_bezier = 1;
        for (int dir = 0; dir < dim; ++dir)
          n_cells_per_bezier *= num_visualization_elements_[dir];
        const int n_total_cells = n_bezier_elements * n_cells_per_bezier;

        const auto cell_ids = IGAVTK::create_vtu_cell_ids<dim>
          (n_points, n_bezier_elements);
        cells->SetCells (n_total_cells , cell_ids);

        //--------------------------------------------------------------------//
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        grid->Allocate(cells->GetNumberOfCells (), 0);
        grid->SetPoints(points);
        const int vtk_enum_type = dim == 3 ? VTK_HEXAHEDRON :
                                  dim == 2 ? VTK_QUAD :
                                            VTK_LINE;
        grid->SetCells(vtk_enum_type, cells);
        mb->SetBlock (id, grid);

        auto point_data = grid->GetPointData();
        this->template create_point_data<dim, codim>
        (mapping, quad, point_num_map, point_data);
      }
      else  // Creating a stuctured grid.
      {
        vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();

        const auto grid_elem = mapping->get_grid()->get_num_intervals ();
        if (dim == 1)
          grid->SetDimensions((num_visualization_elements_[0] + 1) * grid_elem[0], 1, 1);
        else if (dim == 2)
          grid->SetDimensions((num_visualization_elements_[0] + 1) * grid_elem[0],
                              (num_visualization_elements_[1] + 1) * grid_elem[1], 1);
        else if (dim == 3)
          grid->SetDimensions((num_visualization_elements_[0] + 1) * grid_elem[0],
                              (num_visualization_elements_[1] + 1) * grid_elem[1],
                              (num_visualization_elements_[2] + 1) * grid_elem[2]);

        grid->SetPoints(points);
        mb->SetBlock (id, grid);

        auto point_data = grid->GetPointData();
        this->template create_point_data<dim, codim>
        (mapping, quad, point_num_map, point_data);
      }

      ++id;
  }
};



template <int dim, int codim>
void
IGAVTK::
generate_knot_mesh_grids(vtkMultiBlockDataSet* const mb,
                         unsigned int& id,
                         const bool is_parametric) const
{
  static const int space_dim = dim + codim;

  const auto mappings = funcs_container_->template get_all_mappings<dim, codim>();

  using Mapping = Function<dim, 0, dim+codim, 1>;
  using PhysPoint = typename Mapping::Value;
  using ParamPoint = Points<dim>;

  AssertThrow (dim != 1, ExcNotImplemented());

  SafeSTLArray<shared_ptr<Quadrature<1>>, dim> quadratures;
  for (int dir = 0; dir < dim; ++dir)
    quadratures[dir] = QUniform<1>::create (num_visualization_elements_[dir] + 1);

  SafeSTLVector<Real> zero_vec{{0.0}};
  SafeSTLVector<Real> one_vec{{1.0}};

  const Size n_points_per_single_cell = 2;

  for (const auto &m : mappings)
  {
    using Fun_ = Function<dim, 0, dim + codim, 1>;
    shared_ptr<Fun_> mapping = m.first;

    if ((std::dynamic_pointer_cast<IdentityFunction<dim, space_dim>>(mapping) != nullptr) != is_parametric)
      continue;

    const auto cartesian_grid = mapping->get_grid();
    const auto& n_intervals = cartesian_grid->get_num_intervals();

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
      const Size n_cells_in_knot_line = n_intervals[dir] * (num_visualization_elements_[dir]);
      n_vtk_cells += n_knot_lines * n_cells_in_knot_line;
      n_vtk_points += n_knot_lines * (n_cells_in_knot_line + 1);
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
      const auto &face_coords = face_coords_tensor.get_flat_cartesian_product();

      const auto &quad_dir = quadratures[dir];
      TensorSize<dim> n_quad_points(1);
      n_quad_points[dir] = quad_dir->get_num_points();
      TensorProductArray<dim> quad_points_1d(n_quad_points);
      TensorProductArray<dim> quad_weights_1d(n_quad_points);
      quad_points_1d.copy_data_direction(dir, quad_dir->get_coords_direction(0));

      // Iterating along every knot coordinates of the face
      Index flat_id = 0;
      TensorIndex<dim> elem_t_id(0);
      for (const auto &cpm : face_coords)
      {
        const auto t_id_face = face_coords_tensor.flat_to_tensor(flat_id++);
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
            for (int dir2 = 0; dir2 < dim; ++dir2)
              point_tmp[dir2] = pp[dir2];
            vtk_points->SetPoint (vtk_pt_id++, point_tmp);
          }
          --vtk_pt_id;

          for (int c_id = 0; c_id < physical_points.get_num_points()-1;  ++c_id,
                vtk_pt_id_0 += (n_points_per_single_cell-1))
          {
            for (int i = 0; i < n_points_per_single_cell; ++i)
              tuple[i+1] = vtk_pt_id_0 + i;
            vtk_cell_ids->SetTupleValue(vtk_tuple_id++, tuple);
          }
        }
        ++vtk_pt_id;
      } // face_coords

    } // dir

    // Creating grid.
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(vtk_points);

    // Creating cells.
    vtkSmartPointer<vtkCellArray> vtk_cells = vtkSmartPointer<vtkCellArray>::New();
    vtk_cells->SetCells (n_vtk_cells, vtk_cell_ids);
    grid->Allocate(vtk_cells->GetNumberOfCells (), 0);
    const int vtk_enum_type = VTK_LINE;
    grid->SetCells(vtk_enum_type, vtk_cells);

    mb->SetBlock (id, grid);
    ++id;

  } // mappings
};



template <int dim, int codim>
void
IGAVTK::
generate_control_mesh_grids(vtkMultiBlockDataSet* const mb,
                            unsigned int& id) const
{
  const auto mappings = funcs_container_->template get_all_mappings<dim, codim>();

  static const int space_dim = dim + codim;
  using Fun_   =   Function<dim, 0, space_dim, 1>;
  using IgFun_ = IgFunction<dim, 0, space_dim, 1>;

  for (const auto &m : mappings)
  {
      shared_ptr<Fun_> mapping = m.first;
      const auto name = m.second;

      const auto ig_func = std::dynamic_pointer_cast<IgFun_>(mapping);
      if (ig_func == nullptr)
        continue;

      const auto space = ig_func->get_ig_space();
      const auto& coefs = ig_func->get_coefficients();
      const auto& dofs = space->get_dof_distribution();

      // TODO: include assert here to verify that all the components share the 
      // same space.
      const auto &dofs_table = dofs->get_num_dofs_table();
      const auto n_pts_dir = dofs_table[0];

      auto grid = vtkSmartPointer<vtkStructuredGrid>::New();
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
        grid->SetDimensions(n_pts_dir[0], 1, 1);
      else if (dim == 2)
        grid->SetDimensions(n_pts_dir[0], n_pts_dir[1], 1);
      else if (dim == 3)
        grid->SetDimensions(n_pts_dir[0], n_pts_dir[1], n_pts_dir[2]);

      grid->SetPoints(points);
      mb->SetBlock (id, grid);

      ++id;
  }
}



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
//   this->template create_geometries<2>();
  this->template create_geometries<3>();
};



SafeSTLArray<SafeSTLVector<string>, 3>
IGAVTK::
get_map_names () const
{
  SafeSTLArray<SafeSTLVector<string>, 3> names;
  auto& all_names = names[0];
  auto& parametric_names = names[1];
  auto& mapped_names = names[2];

  for (const auto& m : funcs_container_->template get_all_mappings<2, 0>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (std::dynamic_pointer_cast<IdentityFunction<2, 2>>(mapping) != nullptr)
      parametric_names.push_back(name);
    else
      mapped_names.push_back(name);
    all_names.push_back(name);
  }

  for (const auto& m : funcs_container_->template get_all_mappings<1, 1>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (std::dynamic_pointer_cast<IdentityFunction<1, 2>>(mapping) != nullptr)
      parametric_names.push_back(name);
    else
      mapped_names.push_back(name);
    all_names.push_back(name);
  }

  for (const auto& m : funcs_container_->template get_all_mappings<3, 0>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (std::dynamic_pointer_cast<IdentityFunction<3, 3>>(mapping) != nullptr)
      parametric_names.push_back(name);
    else
      mapped_names.push_back(name);
    all_names.push_back(name);
  }

  for (const auto& m : funcs_container_->template get_all_mappings<2, 1>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (std::dynamic_pointer_cast<IdentityFunction<2, 3>>(mapping) != nullptr)
      parametric_names.push_back(name);
    else
      mapped_names.push_back(name);
    all_names.push_back(name);
  }

  for (const auto& m : funcs_container_->template get_all_mappings<1, 2>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (std::dynamic_pointer_cast<IdentityFunction<1, 3>>(mapping) != nullptr)
      parametric_names.push_back(name);
    else
      mapped_names.push_back(name);
    all_names.push_back(name);
  }

  return names;
};



template <int dim>
vtkSmartPointer<vtkIdTypeArray>
IGAVTK::
create_vtu_cell_ids (const TensorSize<dim>& n_points_per_direction,
                     const Size& n_bezier_elements)
{
  vtkSmartPointer<vtkIdTypeArray> cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();

  const auto grid = CartesianGrid<dim>::create (n_points_per_direction);

  static constexpr iga::Size n_points_per_single_cell = pow (2, dim);
  const int n_cells_per_bezier = grid->get_num_all_elems ();
  const int n_points_per_bezier_element = n_points_per_direction.flat_size ();
  const int n_total_cells = n_bezier_elements * n_cells_per_bezier;

  // Creating the connectivity base ------------------------------------------//
  Connectivity_t_<dim> connectivity_base (n_cells_per_bezier);

  // Building the offsets container. According to the vtk elements connectivity.
  using T_ = SafeSTLArray < SafeSTLArray<int, dim>, n_points_per_single_cell>;
  const  T_ delta_idx =
    dim == 1 ? T_({{0},       {1}})                          :   // dim = 1
    dim == 2 ? T_({{0, 0},    {1, 0},    {1, 1},    {0, 1}}) :   // dim = 2
               T_({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                   {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}); // dim = 3

  const TensorIndex<dim> weight_points =
    MultiArrayUtils< dim >::compute_weight(n_points_per_direction);

  TensorIndex<dim> vtk_vertex_tensor_idx;

  auto conn_el = connectivity_base.begin ();
  auto elem = grid->begin();
  const auto elem_end = grid->end();
  for (; elem != elem_end; ++elem, ++conn_el)
  {
    SafeSTLArray<Index,dim> vtk_elem_tensor_idx = elem->get_tensor_index();

    auto conn = conn_el->begin ();
    for (int iVertex = 0; iVertex < n_points_per_single_cell; ++iVertex, ++conn)
    {
        for (int i = 0; i < dim; ++i)
          vtk_vertex_tensor_idx[i] = vtk_elem_tensor_idx[i] + delta_idx[iVertex][i];

        *conn = MultiArrayUtils<dim>::tensor_to_flat_index(vtk_vertex_tensor_idx, weight_points);
    }
  }
  //--------------------------------------------------------------------------//

  cell_ids->SetNumberOfComponents (n_points_per_single_cell + 1);
  cell_ids->SetNumberOfTuples (n_total_cells);

  const int tuple_size = n_points_per_single_cell + 1;
  int cell_id = 0;
  for (int be = 0; be < n_bezier_elements; ++be)
  {
    const int vtk_vertex_id_offset = n_points_per_bezier_element * be;
    for (const auto& cc : connectivity_base)
    {
      vtkIdType* tuple = new vtkIdType[tuple_size];
      tuple[0] = n_points_per_single_cell;
      int point_id = 1;
      for (const auto& c : cc)
        tuple[point_id++] = c + vtk_vertex_id_offset;
      cell_ids->SetTupleValue(cell_id++, tuple);
    }
  }

  return cell_ids;
};


template <int dim>
SafeSTLVector<SafeSTLVector<Index>>
IGAVTK::
create_points_numbering_map (const shared_ptr<const CartesianGrid<dim>> grid,
                             const TensorSize<dim>& n_pts_per_vtk_elem,
                             const bool is_unstructured) const
{
  const Size num_elems = grid->get_num_all_elems();

  Size n_total_pts_per_element = 1;
  for (int dir = 0; dir < dim; ++dir)
    n_total_pts_per_element *= n_pts_per_vtk_elem[dir];

  SafeSTLVector<SafeSTLVector<Index>>
  points_map (num_elems, SafeSTLVector<Index> (n_total_pts_per_element));

  if (is_unstructured) // Unstructured grid.
  {
    Index point_id = 0;
    for (auto &npm_el : points_map)
      for (auto &npm : npm_el)
        npm = point_id++;
  }
  else // Structured grid.
  {
    TensorSize<dim>  n_pts_per_mesh;   // Number of points per direction of
                                       // vtk mesh.
    const auto grid_elems = grid->get_num_intervals();

    for (int dir = 0 ; dir < dim ; ++dir)
      n_pts_per_mesh[dir] = grid_elems[dir] * n_pts_per_vtk_elem[dir];

    TensorIndex<dim> elem_t_id;        // Tensorial index of the Bezier element.
    TensorIndex<dim> pt_mesh_t_offset; // Tensorial index of the first point in
                                       // Bezier element.
    TensorIndex<dim> pt_mesh_t_id;     // Tensorial index of the point.
    TensorIndex<dim> pt_elem_t_id;     // Tensorial index of the point referred
                                       // to the number of points in a
                                       // single element.

    const auto w_elem_pts = MultiArrayUtils<dim>::compute_weight(n_pts_per_vtk_elem);
    const auto w_mesh_pts = MultiArrayUtils<dim>::compute_weight(n_pts_per_mesh);

    Index i_el = 0;
    for (auto &nmp_el : points_map)
    {
      // Computing the tensor index of the first point of the element.
      elem_t_id = grid->flat_to_tensor(i_el++);
      for (int dir = 0 ; dir < dim ; ++dir)
        pt_mesh_t_offset[dir] = elem_t_id[dir] * n_pts_per_vtk_elem[dir];

      Index point_id = 0;
      for (auto &nmp : nmp_el) // Iterating along the points in a Bezier element.
      {
        // Computing the tensor index of the point referred to the number
        // of points into a single Bezier element.
        pt_elem_t_id = MultiArrayUtils<dim>::flat_to_tensor_index(point_id++, w_elem_pts);

        // Computing the tensor index of the point referred to the number
        // of points into the whole mesh.
        for (int dir = 0 ; dir < dim ; ++dir)
          pt_mesh_t_id[dir] = pt_mesh_t_offset[dir] + pt_elem_t_id[dir];

        nmp = MultiArrayUtils<dim>::tensor_to_flat_index(pt_mesh_t_id, w_mesh_pts);
      }
    }
  }

  return points_map;
};



template <>
void
IGAVTK::
create_point_data<2, 0> (const shared_ptr<Function<2, 0, 2, 1>> map,
                         const Quadrature<2> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<2, 0, 1, 1>(map, quad, points_map, point_data);
  this->template create_point_data<2, 0, 2, 1>(map, quad, points_map, point_data);
};



template <>
void
IGAVTK::
create_point_data<1, 1> (const shared_ptr<Function<1, 0, 2, 1>> map,
                         const Quadrature<1> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<1, 1, 1, 1>(map, quad, points_map, point_data);
  this->template create_point_data<1, 1, 2, 1>(map, quad, points_map, point_data);
};



template <>
void
IGAVTK::
create_point_data<3, 0> (const shared_ptr<Function<3, 0, 3, 1>> map,
                         const Quadrature<3> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<3, 0, 1, 1>(map, quad, points_map, point_data);
  this->template create_point_data<3, 0, 3, 1>(map, quad, points_map, point_data);
};



template <>
void
IGAVTK::
create_point_data<2, 1> (const shared_ptr<Function<2, 0, 3, 1>> map,
                         const Quadrature<2> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<2, 1, 1, 1>(map, quad, points_map, point_data);
  this->template create_point_data<2, 1, 3, 1>(map, quad, points_map, point_data);
};



template <>
void
IGAVTK::
create_point_data<1, 2> (const shared_ptr<Function<1, 0, 3, 1>> map,
                         const Quadrature<1> &quad,
                         const SafeSTLVector<SafeSTLVector<Index>>& points_map,
                         vtkPointData* const point_data) const
{
  this->template create_point_data<1, 2, 1, 1>(map, quad, points_map, point_data);
  this->template create_point_data<1, 2, 3, 1>(map, quad, points_map, point_data);
};



template <int dim, int codim, int range, int rank>
void
IGAVTK::
create_point_data (const shared_ptr<Function<dim, 0, dim + codim, 1>> mapping,
                   const Quadrature<dim> &quad,
                   const SafeSTLVector<SafeSTLVector<Index>>& point_num_map,
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
      Index pt_id = 0;
      for (const auto& v : values)
      {
//         for (int c = 0; c < n_comp; ++c)
//           *(tuple + c) = v[c];
        this->template tensor_to_tuple<Value>(v, tuple);
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

#include <igatools/io/reader.h>


template <int dim>
void
IGAVTK::
create_geometries ()
{
  using std::dynamic_pointer_cast;
  using std::const_pointer_cast;

  static const int codim = 0;
  static const int space_dim = dim + codim;
  static const int range = dim;
  static const int rank = 1;
  static const Transformation transf = Transformation::h_grad;

  using Fun_ = Function<dim, codim, range, rank>;
  using IgFun_ = IgFunction<dim, codim, range, rank>;
  using PhysSpace_ = PhysicalSpace<dim, range, rank, codim, transf>;
  using RefSpace_ = ReferenceSpace<dim, range, rank>;
  using IdFun_ = IdentityFunction<dim, space_dim>;

  // File names;
  const string fname_0 = "patch_0_" + std::to_string (dim) + "D.xml";
  const string fname_1 = "patch_1_" + std::to_string (dim) + "D.xml";
  const string fname_2 = "patch_2_" + std::to_string (dim) + "D.xml";
  const string fname_3 = "patch_3_" + std::to_string (dim) + "D.xml";

  // Reading maps
  const auto map_0 = get_mapping_from_file<dim, codim> (fname_0);
  const auto map_1 = get_mapping_from_file<dim, codim> (fname_1);
  const auto map_2 = get_mapping_from_file<dim, codim> (fname_2);
  const auto map_3 = get_mapping_from_file<dim, codim> (fname_3);

  // Transforming mapping to ig mapping.
  const auto ig_func_0 = dynamic_pointer_cast<IgFun_>(map_0);
  const auto ig_func_1 = dynamic_pointer_cast<IgFun_>(map_1);
  const auto ig_func_2 = dynamic_pointer_cast<IgFun_>(map_2);
  const auto ig_func_3 = dynamic_pointer_cast<IgFun_>(map_3);

  // Getting ig spaces.
  const auto space_0 = ig_func_0->get_ig_space ();
  const auto space_1 = ig_func_1->get_ig_space ();
  const auto space_2 = ig_func_2->get_ig_space ();
  const auto space_3 = ig_func_3->get_ig_space ();

  // Getting reference spaces.
  const auto ref_space_0 = const_pointer_cast<RefSpace_>(dynamic_pointer_cast<const RefSpace_>(space_0));
  const auto ref_space_1 = const_pointer_cast<RefSpace_>(dynamic_pointer_cast<const RefSpace_>(space_1));
  const auto ref_space_2 = const_pointer_cast<RefSpace_>(dynamic_pointer_cast<const RefSpace_>(space_2));
  const auto ref_space_3 = const_pointer_cast<RefSpace_>(dynamic_pointer_cast<const RefSpace_>(space_3));

  // Building physical spaces.
  const auto phys_space_0 = PhysSpace_::create (ref_space_0, map_0);
  const auto phys_space_1 = PhysSpace_::create (ref_space_1, map_1);
  const auto phys_space_2 = PhysSpace_::create (ref_space_2, map_2);
  const auto phys_space_3 = PhysSpace_::create (ref_space_3, map_3);

  // Getting grids;
  const auto grid_0 = ref_space_0->get_grid ();
  const auto grid_1 = ref_space_1->get_grid ();
  const auto grid_2 = ref_space_2->get_grid ();
  const auto grid_3 = ref_space_3->get_grid ();

  // Creating identity functions.
  auto id_map_0 = IdFun_::create (grid_0);
  auto id_map_1 = IdFun_::create (grid_1);
  auto id_map_2 = IdFun_::create (grid_2);
  auto id_map_3 = IdFun_::create (grid_3);

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
  auto ps_func_0 = dynamic_pointer_cast<Fun_>(IgFun_::create (phys_space_0, phys_coeff_0));
  auto ps_func_1 = dynamic_pointer_cast<Fun_>(IgFun_::create (phys_space_1, phys_coeff_1));
  auto ps_func_2 = dynamic_pointer_cast<Fun_>(IgFun_::create (phys_space_2, phys_coeff_2));
  auto ps_func_3 = dynamic_pointer_cast<Fun_>(IgFun_::create (phys_space_3, phys_coeff_3));

  // Creating ig functions for reference spaces.
  auto rf_func_0 = dynamic_pointer_cast<Fun_>(IgFun_::create (ref_space_0, ref_coeff_0));
  auto rf_func_1 = dynamic_pointer_cast<Fun_>(IgFun_::create (ref_space_1, ref_coeff_1));
  auto rf_func_2 = dynamic_pointer_cast<Fun_>(IgFun_::create (ref_space_2, ref_coeff_2));
  auto rf_func_3 = dynamic_pointer_cast<Fun_>(IgFun_::create (ref_space_3, ref_coeff_3));

  // Adding all the stuff to the functions container.
  funcs_container_ = std::make_shared<FunctionsContainer>();

  // Inserting geometries.
  funcs_container_->insert_mapping(map_0, "map_0");
  funcs_container_->insert_mapping(map_1, "map_1");
  funcs_container_->insert_mapping(map_2, "map_2");
  funcs_container_->insert_mapping(map_3, "map_3");

  funcs_container_->insert_mapping(id_map_0, "id_map_0");
  funcs_container_->insert_mapping(id_map_1, "id_map_1");
  funcs_container_->insert_mapping(id_map_2, "id_map_2");
  funcs_container_->insert_mapping(id_map_3, "id_map_3");

  // Inserting associated functions.
  funcs_container_->insert_function(map_0, ps_func_0, "phys_func_0");
  funcs_container_->insert_function(map_1, ps_func_1, "phys_func_0");
  funcs_container_->insert_function(map_2, ps_func_2, "phys_func_0");
  funcs_container_->insert_function(map_3, ps_func_3, "phys_func_0");

  funcs_container_->insert_function(id_map_0, rf_func_0, "ref_func_0");
  funcs_container_->insert_function(id_map_1, rf_func_1, "ref_func_1");
  funcs_container_->insert_function(id_map_2, rf_func_2, "ref_func_2");
  funcs_container_->insert_function(id_map_3, rf_func_3, "ref_func_3");
};