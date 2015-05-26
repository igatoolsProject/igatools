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
set_number_visualization_points (const int* const num_visualization_points)
{

  static const int dim = 3;
  for (int dir = 0; dir < dim; ++dir)
    num_visualization_points_[dir] = *(num_visualization_points + dir);
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
                   const bool& create_parametric_mesh,
                   const bool& create_physical_mesh,
                   vtkMultiBlockDataSet* const mb) const
{
  std::cout << "HOla" << std::endl;
  Assert (file_name_ != "", ExcMessage ("Not specified file name."));
  Assert (file_path_ != "", ExcMessage ("Not specified file path."));
  Assert (num_visualization_points_[0] > 1,
          ExcMessage ("Number of visualization points must be > 1 in every direction."));
  Assert (num_visualization_points_[1] > 1,
          ExcMessage ("Number of visualization points must be > 1 in every direction."));
  Assert (num_visualization_points_[2] > 1,
          ExcMessage ("Number of visualization points must be > 1 in every direction."));

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

//     this->generate_control_mesh_grid<2, 0>(internal_mb, internal_block_index);
//     this->generate_control_mesh_grid<1, 1>(internal_mb, internal_block_index);
//     this->generate_control_mesh_grid<3, 0>(internal_mb, internal_block_index);
//     this->generate_control_mesh_grid<2, 1>(internal_mb, internal_block_index);
//     this->generate_control_mesh_grid<1, 2>(internal_mb, internal_block_index);
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
//     this->generate_parametric_mesh_grid<2, 0>(internal_mb, internal_block_index, unstructured);
//     this->generate_parametric_mesh_grid<1, 1>(internal_mb, internal_block_index, unstructured);
//     this->generate_parametric_mesh_grid<3, 0>(internal_mb, internal_block_index, unstructured);
//     this->generate_parametric_mesh_grid<2, 1>(internal_mb, internal_block_index, unstructured);
//     this->generate_parametric_mesh_grid<1, 2>(internal_mb, internal_block_index, unstructured);
  }

  if (create_physical_mesh)
  {
    const auto& n = names[2];
    vtkMultiBlockDataSet* internal_mb =
      vtkMultiBlockDataSet::SafeDownCast(mb->GetBlock(block_index++));
    Assert (internal_mb != nullptr, ExcNullPtr ());
    internal_mb->SetNumberOfBlocks(n.size ());

    unsigned int internal_block_index = 0;
    for (const auto& name : n)
      internal_mb->GetMetaData(internal_block_index++)->Set(vtkCompositeDataSet::NAME(), name.c_str());
//     this->generate_physical_mesh_grid<2, 0>(internal_mb, internal_block_index, unstructured);
//     this->generate_physical_mesh_grid<1, 1>(internal_mb, internal_block_index, unstructured);
//     this->generate_physical_mesh_grid<3, 0>(internal_mb, internal_block_index, unstructured);
//     this->generate_physical_mesh_grid<2, 1>(internal_mb, internal_block_index, unstructured);
//     this->generate_physical_mesh_grid<1, 2>(internal_mb, internal_block_index, unstructured);
  }
};


template <int dim, int codim>
void
IGAVTK::
generate_physical_mesh_grids(vtkMultiBlockDataSet* const mb,
                             unsigned int& id,
                             const bool unstructured) const
{
#if 0
  const auto mappings = funcs_container_->template get_all_mappings<dim, codim>();

  for (const auto &m : mappings)
  {
      using Fun_ = Function<dim, 0, dim + codim, 1>;
      shared_ptr<Fun_> mapping = m.first;

      if (is_identity_mapping<dim, codim>(mapping) != identity_map)
        continue;

      auto name    = m.second;

      const int vtk_enum_type = dim == 3 ? VTK_HEXAHEDRON :
                                dim == 2 ? VTK_QUAD :
                                           VTK_LINE;

      vtkSmartPointer<vtkPoints>   points = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkCellArray> cellsArray = vtkSmartPointer<vtkCellArray>::New();

      // Generating visualization quadrature.
      // TODO: quadratures for 1D, 2D and 3D should be moved outside this
      // function.

      TensorSize<dim> n_points;
      for (int dir = 0; dir < dim; ++dir)
      {
        n_points[dir] = num_visualization_points_[dir];
        Assert (n_points[dir] > 1, ExcMessage ("Wrong number of visualization points."));
      }
      QUniform<dim> quad (n_points);

      const int n_bezier_elements = mapping->get_grid ()->get_num_all_elems ();
      const int n_points_per_bezier_element = quad.get_num_points ();
      const int total_num_points = n_points_per_bezier_element * n_bezier_elements;

      // Setting the points ------------------------------------------------------//
      points->SetNumberOfPoints (total_num_points);

      auto flag = ValueFlags::point | ValueFlags::value;
      mapping->reset(flag, quad);

      using IgFun = IgFunction<dim, 0, dim+codim, 1>;
      const auto ig_fun = std::dynamic_pointer_cast<IgFun>(mapping);
      AssertThrow (ig_fun != nullptr, ExcNullPtr ());
      LogStream out;
      ig_fun->print_info (out);

      try
      {
        const auto kk = *ig_fun;
      }
      catch (const std::bad_weak_ptr& e)
      {
        std::cout << "Capturing ig mapping exception " << e.what() << std::endl;
      }


      using ElementIterator = typename  Fun_::ElementIterator;
      ElementIterator elem;
      ElementIterator end;
      try
      {
        elem = mapping->begin();
        end = mapping->end();
      }
      catch (const std::bad_weak_ptr& e)
      {
        std::cout << "Capturing element exception " << e.what() << std::endl;
      }

      const auto topology = Topology<dim>();

      mapping->init_cache(elem, topology);

      double point_tmp[3]; point_tmp[0] = 0.0; point_tmp[1] = 0.0; point_tmp[2] = 0.0;
      int point_id = 0;
      for (; elem != end; ++elem)
      {
        mapping->fill_cache(elem, topology, 0);

        auto element_vertices_tmp = elem->template get_values<_Value, dim>(0);
        for (const auto& p : element_vertices_tmp)
        {
          for (int dir = 0; dir < dim ; ++dir)
            point_tmp[dir] = p[dir];
          points->SetPoint (point_id++, point_tmp);
        }
      }
      //----------------------------------------------------------------------//

      if (unstructured) // Creating an unstuctured grid.
      {
        // Setting the cells -------------------------------------------------//
        int n_cells_per_bezier = 1;
        for (int dir = 0; dir < dim; ++dir)
          n_cells_per_bezier *= num_visualization_points_[dir] - 1;
        const int n_total_cells = n_bezier_elements * n_cells_per_bezier;

        const auto cell_ids = IGAVTK::create_vtu_cell_ids<dim>
          (n_points, n_bezier_elements);
        cellsArray->SetCells (n_total_cells , cell_ids);

        //--------------------------------------------------------------------//
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        grid->Allocate(cellsArray->GetNumberOfCells (), 0);
        grid->SetPoints(points);
        grid->SetCells(vtk_enum_type, cellsArray);
        mb->SetBlock (id, grid);
      }
      else  // Creating a stuctured grid.
      {
        vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();

        // TODO: Revisit this.
        if (dim == 1)
          grid->SetDimensions(num_visualization_points_[0], 1, 1);
        else if (dim == 2)
          grid->SetDimensions(num_visualization_points_[0],
                              num_visualization_points_[1], 1);
        else if (dim == 3)
          grid->SetDimensions(num_visualization_points_[0],
                              num_visualization_points_[1],
                              num_visualization_points_[2]);

        grid->SetPoints(points);
        mb->SetBlock (id, grid);
      }

      ++id;
  }
#endif

};


vtkSmartPointer<vtkStructuredGrid>
IGAVTK::
make_grid (int i) const
{
  const auto grid = vtkSmartPointer<vtkStructuredGrid>::New();

  const auto points = vtkSmartPointer<vtkPoints>::New();
  double x, y, z;

  x = 0.0;
  y = 0.0;
  z = 0.0;

  for(unsigned int k = 0; k < 2; k++)
  {
    z += 2.0 * i;
    for(unsigned int j = 0; j < 3; j++)
    {
      y += 1.0 * i;
      for(unsigned int i = 0; i < 2; i++)
      {
        x += .5 * i;
        points->InsertNextPoint(x, y, z);
      }
    }
  }

  // Specify the dimensions of the grid
  grid->SetDimensions(2,3,2);
  grid->SetPoints(points);

  return grid;
};



void
IGAVTK::
parse_file ()
{
  std::cout << "parsing" << std::endl;
//   funcs_container_ = std::make_shared<FunctionsContainer> ();
  std::cout << "parsing" << std::endl;
  ifstream xml_istream(file_name_);
  std::cout << "parsing" << std::endl;
  IArchive xml_in(xml_istream);
  std::cout << "parsing " << file_name_ << std::endl;
  xml_in >> BOOST_SERIALIZATION_NVP (funcs_container_);
  std::cout << "parsing" << std::endl;
  xml_istream.close();
  std::cout << "parsing" << std::endl;

  // Assert here if funcs_container_ is void.
};


template <int dim, int codim>
bool
IGAVTK::
is_identity_mapping (shared_ptr<Function<dim, 0, dim+codim, 1>> map) const
{
  using IdFun_ = IdentityFunction<dim, dim+codim>;
  const auto id_func = std::dynamic_pointer_cast<IdFun_>(map);
  return (id_func != nullptr);
};




SafeSTLArray<SafeSTLVector<string>, 3>
IGAVTK::
get_map_names () const
{
  SafeSTLArray<SafeSTLVector<string>, 3> names;
  auto& all_names = names[0];
  auto& identity_names = names[1];
  auto& mapped_names = names[2];

  for (const auto& m : funcs_container_->template get_all_mappings<2, 0>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (this->is_identity_mapping<2, 0>(mapping))
      identity_names.push_back(name);
    else
      mapped_names.push_back(name);
    all_names.push_back(name);
  }

  for (const auto& m : funcs_container_->template get_all_mappings<1, 1>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (this->is_identity_mapping<1, 1>(mapping))
      identity_names.push_back(name);
    else
      mapped_names.push_back(name);
    all_names.push_back(name);
  }

  for (const auto& m : funcs_container_->template get_all_mappings<3, 0>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (this->is_identity_mapping<3, 0>(mapping))
      identity_names.push_back(name);
    else
      mapped_names.push_back(name);
    all_names.push_back(name);
  }

  for (const auto& m : funcs_container_->template get_all_mappings<2, 1>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (this->is_identity_mapping<2, 1>(mapping))
      identity_names.push_back(name);
    else
      mapped_names.push_back(name);
    all_names.push_back(name);
  }

  for (const auto& m : funcs_container_->template get_all_mappings<1, 2>())
  {
    auto mapping = m.first;
    auto name    = m.second;
    if (this->is_identity_mapping<1, 2>(mapping))
      identity_names.push_back(name);
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