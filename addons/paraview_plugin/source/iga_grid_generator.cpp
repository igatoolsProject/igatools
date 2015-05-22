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

#include <paraview_plugin/iga_grid_generator.h>

#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>

#include <vtkHexahedron.h>
#include <vtkQuad.h>
#include <vtkLine.h>

#include <vtkTypeInt32Array.h>
#include <vtkPointData.h>

#include <igatools/base/quadrature_lib.h>
#include <igatools/base/identity_function.h>
#include <igatools/io/reader.h>
#include <igatools/utils/multi_array_utils.h>
#include <igatools/functions/functions_container.h>

using namespace iga;
using std::string;
using std::shared_ptr;


IGAVTKGridGenerator::
IGAVTKGridGenerator (const string& file_name, const string& file_path,
                     const int* const num_visualization_points)
  :
  file_name_ (file_name),
  file_path_ (file_path),
  quad_1D_ (create_quadrature<1>(num_visualization_points)),
  quad_2D_ (create_quadrature<2>(num_visualization_points)),
  quad_3D_ (create_quadrature<3>(num_visualization_points))
{
  // TODO: to check here if the file exists.
  // TODO: to check here if the file is correct.

  ifstream xml_istream(file_name_);
  IArchive xml_in(xml_istream);
  xml_in >> BOOST_SERIALIZATION_NVP (funcs_container_);
  xml_istream.close();

  const auto& kk_2 = funcs_container_->get_data_dim_codim<2, 0>();
  LogStream out;
  kk_2.print_info (out);

  const auto kk_3 = funcs_container_->get_data_dim_codim<3, 0>();

  // TODO: set here max_dim;
  max_dim_ = 3;

};



auto
IGAVTKGridGenerator::
create (const string& file_name, const string& file_path,
        const int* num_visualization_points) ->
SelfPtr_t_
{
  return SelfPtr_t_ (new Self_t_ (file_name, file_path, num_visualization_points));
};




#if 0
std::pair<int, int>
IGAVTKGridGenerator::
get_dimensions_from_xml (const ptree& xml_tree) const
{
  AssertThrow (xml_element_is_present (xml_tree, "IgMapping"),
               ExcMessage ("IgMapping element is not present."));
  AssertThrow (xml_element_is_unique (xml_tree, "IgMapping"),
               ExcMessage ("IgMapping element is not unique."));

  const auto& igm_xml = get_xml_element (xml_tree, "IgMapping");

  const auto& atts = get_xml_element_attributes (igm_xml);

  std::pair<int, int> dims;
  auto& dim = dims.first;
  auto& codim = dims.second;

  Assert (xml_element_is_present (atts, "Dim"),
          ExcMessage ("Dim attribute is not present."));
  Assert (xml_element_is_unique (atts, "Dim"),
          ExcMessage ("Dim attribute is not unique."));
  dim = atts.get<int> ("Dim");
  AssertThrow (dim > 0 && dim < 4, ExcMessage ("Dim=" + std::to_string (dim) +
                                               " is not a valid value."));

  Assert (xml_element_is_present (atts, "Codim"),
          ExcMessage ("Codim attribute is not present."));
  Assert (xml_element_is_unique (atts, "Codim"),
          ExcMessage ("Codim attribute is not unique."));
  codim = atts.get<int> ("Codim");
  AssertThrow (codim >= 0 && codim < 3,
               ExcMessage ("Codim=" + std::to_string (codim) +
               " is not a valid value."));

  AssertThrow ((dim + codim) > 0 && (dim + codim) < 4,
               ExcMessage ("Invalid space dimension with Dim=" +
               std::to_string (codim) + " and Codim=" + std::to_string (codim) +
               "."));

  return dims;
};
#endif


int
IGAVTKGridGenerator::
fill_solid_output (vtkInformation* outInfo) const
{
#if 0

  // get the ouptut
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
                    outInfo->Get(vtkDataObject::DATA_OBJECT()));
//   output->Reset ();
//   output->Initialize ();

  vtkSmartPointer<vtkPoints>   points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> cellsArray = vtkSmartPointer<vtkCellArray>::New();

  IGAVTKGridGenerator::FillSolidGridVisitor visitor; 
  const bool success = apply_visitor(visitor, function_variant_, points, cellsArray, &quadrature_variant_);
  if (success)
  {
    const int vtk_enum_type = dim_ == 3 ? VTK_HEXAHEDRON :
                              dim_ == 2 ? VTK_QUAD :
                                          VTK_LINE;
    output->Allocate (cellsArray->GetNumberOfCells (), 0);
    output->SetPoints(points);
    output->SetCells(vtk_enum_type, cellsArray);
    
//     std::cout << "Hola" << std::endl;
//   const auto aterials = vtkTypeInt32Array::New();
//   materials->SetNumberOfComponents(3);
//   materials->SetNumberOfTuples(points->GetNumberOfPoints());
//   materials->SetName("material id");
//   for (int i = 0; i < points->GetNumberOfPoints(); ++i)
//   {
//     double* a = new double[3];
//     a[0] = 0;
//     a[1] = 1;
//     a[2] = i;
//     materials->InsertTuple (i, a);
//   }
//   output->GetPointData ()->AddArray(materials);
//   materials->FastDelete ();
//     std::cout << "Adios" << std::endl;
    
    
    return 1;
  }
  else
    return 0;
#endif
};



int
IGAVTKGridGenerator::
fill_identity_output (vtkInformation* outInfo) const
{
#if 0

  // get the ouptut
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
                    outInfo->Get(vtkDataObject::DATA_OBJECT()));
//   output->Reset ();
//   output->Initialize ();

  vtkSmartPointer<vtkPoints>   points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> cellsArray = vtkSmartPointer<vtkCellArray>::New();

  IGAVTKGridGenerator::FillIdentityGridVisitor visitor; 
  const bool success = apply_visitor(visitor, function_variant_, points, cellsArray, &quadrature_variant_);
  if (success)
  {
    const int vtk_enum_type = dim_ == 3 ? VTK_HEXAHEDRON :
                              dim_ == 2 ? VTK_QUAD :
                                          VTK_LINE;
    output->Allocate (cellsArray->GetNumberOfCells (), 0);
    output->SetPoints(points);
    output->SetCells(vtk_enum_type, cellsArray);
    return 1;
  }
  else
    return 0;
#endif
};



int
IGAVTKGridGenerator::
fill_control_mesh_output (vtkInformation* outInfo) const
{
#if 0
  // get the ouptut
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
                    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkPoints>   points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> cellsArray = vtkSmartPointer<vtkCellArray>::New();

  IGAVTKGridGenerator::FillControlMeshGridVisitor visitor; 
  const bool success = apply_visitor(visitor, function_variant_, points, cellsArray);
  if (success)
  {
    const int vtk_enum_type = dim_ == 3 ? VTK_HEXAHEDRON :
                              dim_ == 2 ? VTK_QUAD :
                                          VTK_LINE;
    output->Allocate (cellsArray->GetNumberOfCells (), 0);
    output->SetPoints(points);
    output->SetCells(vtk_enum_type, cellsArray);
    return 1;
  }
  else
    return 0;
#endif
}



template <int dim>
auto
IGAVTKGridGenerator::
create_quadrature (const int* const num_visualization_points) -> QuadPtr_<dim>
{
  TensorSize<dim> n_points;

  for (int dir = 0; dir < dim; ++dir)
  {
    Assert (*(num_visualization_points + dir) < 2,
      ExcMessage ("Number of visualization points per direction must be > 1."));
    n_points[dir] = *(num_visualization_points + dir);
  }

  return QUniform<dim>::create (n_points);
};



template <class FunctionPtr_t_>
auto
IGAVTKGridGenerator::FillSolidGridVisitor::operator()
  (const FunctionPtr_t_ func_ptr,
   const vtkSmartPointer<vtkPoints> points,
   const vtkSmartPointer<vtkCellArray> cellsArray,
   const QuadraturePtrVariant* const quad_var_ptr) ->
result_type
{
#if 0
  // TODO: to manager here the return value.
  // Return false in not sucessfull cases.

  typedef typename FunctionPtr_t_::element_type Function_t_;
  static const int st_dim = Function_t_::dim;

  typedef shared_ptr<QuadratureTensorProduct<st_dim>> QuadraturePtr_t_;
  const auto quad = get<QuadraturePtr_t_> (*quad_var_ptr);
  Assert (quad != nullptr, ExcNullPtr ());

  const auto tmp = quad->get_num_coords_direction();
  TensorSize <st_dim> n_points_per_direction;
  for (int dir = 0; dir < st_dim; ++dir)
    n_points_per_direction[dir] = tmp[dir]; 

  const int n_bezier_elements = func_ptr->get_grid ()->get_num_active_elems ();
  const int n_points_per_bezier_element = quad->get_num_points ();
  const int total_num_points = n_points_per_bezier_element * n_bezier_elements;

  // Setting the points ------------------------------------------------------//
  points->SetNumberOfPoints (total_num_points);

  func_ptr->reset(ValueFlags::value | ValueFlags::point, *quad);

  auto m_elem = func_ptr->begin();
  auto m_end  = func_ptr->end();

  const auto topology = Int<st_dim>();

  func_ptr->init_cache(m_elem, topology);

  double point_tmp[3]; point_tmp[0] = 0.0; point_tmp[1] = 0.0; point_tmp[2] = 0.0;
  int point_id = 0;
  for (; m_elem != m_end; ++m_elem)
  {
    func_ptr->fill_cache(m_elem, topology, 0);

    auto element_vertices_tmp = m_elem->template get_values<0,st_dim>(0);
    for (const auto& p : element_vertices_tmp)
    {
      for (int dir = 0; dir < st_dim ; ++dir)
        point_tmp[dir] = p[dir];
      points->SetPoint (point_id++, point_tmp);
    }
  }
  //--------------------------------------------------------------------------//

  // Setting the cells -------------------------------------------------------//
  int n_cells_per_bezier = 1;
  for (int dir = 0; dir < st_dim; ++dir)
    n_cells_per_bezier *= n_points_per_direction[dir] - 1;
  const int n_total_cells = n_bezier_elements * n_cells_per_bezier;

  const auto cell_ids = IGAVTKGridGenerator::create_cell_ids<st_dim>
    (n_points_per_direction, n_bezier_elements);
  cellsArray->SetCells (n_total_cells , cell_ids);
  //--------------------------------------------------------------------------//

  return true;
#endif
};



template <class FunctionPtr_t_>
auto
IGAVTKGridGenerator::FillIdentityGridVisitor::operator()
  (const FunctionPtr_t_ func_ptr,
   const vtkSmartPointer<vtkPoints> points,
   const vtkSmartPointer<vtkCellArray> cellsArray,
   const QuadraturePtrVariant* const quad_var_ptr) ->
result_type
{
#if 0
  // TODO: to manager here the return value.
  // Return false in not sucessfull cases.

  typedef typename FunctionPtr_t_::element_type Function_t_;
  static const int st_dim = Function_t_::dim;
  static const int st_range = Function_t_::range;

  typedef shared_ptr<QuadratureTensorProduct<st_dim>> QuadraturePtr_t_;
  const auto quad = get<QuadraturePtr_t_> (*quad_var_ptr);
  Assert (quad != nullptr, ExcNullPtr ());

  const auto tmp = quad->get_num_coords_direction();
  TensorSize <st_dim> n_points_per_direction;
  for (int dir = 0; dir < st_dim; ++dir)
    n_points_per_direction[dir] = tmp[dir]; 

  const auto grid = func_ptr->get_grid ();
  typedef IdentityFunction<st_dim, st_range> IdentityFunction_t_;
  const auto id_func_ptr = IdentityFunction_t_::create (grid); 

  const int n_bezier_elements = id_func_ptr->get_grid ()->get_num_active_elems ();
  const int n_points_per_bezier_element = quad->get_num_points ();
  const int total_num_points = n_points_per_bezier_element * n_bezier_elements;

  // Setting the points ------------------------------------------------------//
  points->SetNumberOfPoints (total_num_points);

  id_func_ptr->reset(ValueFlags::value | ValueFlags::point, *quad);

  auto m_elem = id_func_ptr->begin();
  auto m_end  = id_func_ptr->end();

  const auto topology = Int<st_dim>();

  id_func_ptr->init_cache(m_elem, topology);

  double point_tmp[3]; point_tmp[0] = 0.0; point_tmp[1] = 0.0; point_tmp[2] = 0.0;
  int point_id = 0;
  for (; m_elem != m_end; ++m_elem)
  {
    id_func_ptr->fill_cache(m_elem, topology, 0);

    auto element_vertices_tmp = m_elem->template get_values<0,st_dim>(0);
    for (const auto& p : element_vertices_tmp)
    {
      for (int dir = 0; dir < st_dim ; ++dir)
        point_tmp[dir] = p[dir];
      points->SetPoint (point_id++, point_tmp);
    }
  }
  //--------------------------------------------------------------------------//

  // Setting the cells -------------------------------------------------------//
  int n_cells_per_bezier = 1;
  for (int dir = 0; dir < st_dim; ++dir)
    n_cells_per_bezier *= n_points_per_direction[dir] - 1;
  const int n_total_cells = n_bezier_elements * n_cells_per_bezier;

  const auto cell_ids = IGAVTKGridGenerator::create_cell_ids<st_dim>
    (n_points_per_direction, n_bezier_elements);
  cellsArray->SetCells (n_total_cells , cell_ids);
  //--------------------------------------------------------------------------//

  return true;
#endif
};



template <class FunctionPtr_t_>
auto
IGAVTKGridGenerator::FillControlMeshGridVisitor::operator()
  (const FunctionPtr_t_ func_ptr,
   const vtkSmartPointer<vtkPoints> points,
   const vtkSmartPointer<vtkCellArray> cellsArray) ->
result_type
{
#if 0
  // TODO: to manager here the return value.
  // Return false in not sucessfull cases.

  typedef typename FunctionPtr_t_::element_type Function_t_;
  static const int st_dim  = Function_t_::dim;
  static const int st_range = Function_t_::range;
  static const int st_rank  = Function_t_::rank;

  typedef ReferenceSpace<st_dim, st_range, st_rank> RefSpace_t_;
  typedef IgFunction<RefSpace_t_> IgFunction_t_;
  const auto ig_func_ptr = std::dynamic_pointer_cast<IgFunction_t_>(func_ptr);
  const auto space = ig_func_ptr->get_iga_space ();
  const auto& ctrl_pts = ig_func_ptr->get_coefficients ();

  // Getting the grid dimensions ---------------------------------------------//
  TensorSize<st_dim> n_points_per_direction;
  // TODO: this must be changed. How to manage different components?
  const int comp_id = 0;
  int n_cells = 1;
  for (int dir = 0; dir < st_dim; ++dir)
  {
    n_points_per_direction[dir] = space->get_num_basis (comp_id, dir);
    n_cells *= (n_points_per_direction[dir] - 1);
  }
  const int n_points = n_points_per_direction.flat_size ();

  AssertThrow (n_points * st_dim == ctrl_pts.size (),
               ExcDimensionMismatch (n_points * st_dim, ctrl_pts.size ()));
  //--------------------------------------------------------------------------//

  // Setting the points ------------------------------------------------------//
  points->SetNumberOfPoints (n_points);

  double point_tmp[3]; point_tmp[0] = 0.0; point_tmp[1] = 0.0; point_tmp[2] = 0.0;
  for (int point_id = 0; point_id < n_points; ++point_id)
  {
    for (int dir = 0; dir < st_dim ; ++dir)
      point_tmp[dir] = ctrl_pts[point_id + dir * n_points];
    points->SetPoint (point_id, point_tmp);
  }
  //--------------------------------------------------------------------------//

  // Setting the cells -------------------------------------------------------//
  const auto cell_ids = IGAVTKGridGenerator::create_cell_ids<st_dim>
    (n_points_per_direction, 1);
  cellsArray->SetCells (n_cells , cell_ids);
  //--------------------------------------------------------------------------//

  return true;
#endif
};



template <int dim>
vtkSmartPointer<vtkIdTypeArray>
IGAVTKGridGenerator::
create_cell_ids (const TensorSize<dim>& n_points_per_direction,
                 const Size& n_bezier_elements)
{
#if 0
  vtkSmartPointer<vtkIdTypeArray> cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();

  const auto grid = CartesianGrid<dim>::create (n_points_per_direction);

  static constexpr iga::Size n_points_per_single_cell = pow (2, dim);
  const int n_cells_per_bezier = grid->get_num_active_elems ();
  const int n_points_per_bezier_element = n_points_per_direction.flat_size ();
  const int n_total_cells = n_bezier_elements * n_cells_per_bezier;

  const auto connectivity_base = IGAVTKGridGenerator::create_connectivity_base<dim>
    (n_points_per_direction);
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
#endif
};



template <int dim>
auto
IGAVTKGridGenerator::
create_connectivity_base_vtu (const TensorSize<dim>& n_points_per_direction) ->
Connectivity_t_<dim>
{
  vtkSmartPointer<vtkIdTypeArray> cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();

  const auto grid = CartesianGrid<dim>::create (n_points_per_direction);

  static constexpr iga::Size n_points_per_cell = pow (2, dim);
  const int n_cells_per_bezier = grid->get_num_active_elems ();

  Connectivity_t_<dim> connectivity_base (n_cells_per_bezier);

  // Building the offsets container.
  special_array < special_array<int, dim>, n_points_per_cell> delta_idx;

  if (dim == 1)
  {
      delta_idx[0][0] = 0;
      delta_idx[1][0] = 1;
  }
  else if (dim == 2)
  {
      delta_idx[0][0] = 0;
      delta_idx[0][1] = 0;

      delta_idx[1][0] = 1;
      delta_idx[1][1] = 0;

      delta_idx[2][0] = 1;
      delta_idx[2][1] = 1;

      delta_idx[3][0] = 0;
      delta_idx[3][1] = 1;
  }
  else if (dim == 3)
  {
      delta_idx[0][0] = 0;
      delta_idx[0][1] = 0;
      delta_idx[0][2] = 0;

      delta_idx[1][0] = 1;
      delta_idx[1][1] = 0;
      delta_idx[1][2] = 0;

      delta_idx[2][0] = 1;
      delta_idx[2][1] = 1;
      delta_idx[2][2] = 0;

      delta_idx[3][0] = 0;
      delta_idx[3][1] = 1;
      delta_idx[3][2] = 0;

      delta_idx[4][0] = 0;
      delta_idx[4][1] = 0;
      delta_idx[4][2] = 1;

      delta_idx[5][0] = 1;
      delta_idx[5][1] = 0;
      delta_idx[5][2] = 1;

      delta_idx[6][0] = 1;
      delta_idx[6][1] = 1;
      delta_idx[6][2] = 1;

      delta_idx[7][0] = 0;
      delta_idx[7][1] = 1;
      delta_idx[7][2] = 1;
  }

  TensorIndex<dim> weight_points =
    MultiArrayUtils< dim >::compute_weight(n_points_per_direction);

  TensorIndex<dim> vtk_vertex_tensor_idx;

  auto conn_el = connectivity_base.begin ();
  auto elem = grid->begin();
  const auto elem_end = grid->end();
  for (; elem != elem_end; ++elem, ++conn_el)
  {
    special_array<Index,dim> vtk_elem_tensor_idx = elem->get_tensor_index();

    auto conn = conn_el->begin ();
    for (int iVertex = 0; iVertex < n_points_per_cell; ++iVertex, ++conn)
    {
        for (int i = 0; i < dim; ++i)
          vtk_vertex_tensor_idx[i] = vtk_elem_tensor_idx[i] + delta_idx[iVertex][i];

        *conn = MultiArrayUtils<dim>::tensor_to_flat_index(vtk_vertex_tensor_idx, weight_points);
    }
  }

  return connectivity_base;
};
