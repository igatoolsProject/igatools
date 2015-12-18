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

#include <igatools/io/writer.h>

#ifdef XML_IO

#include <igatools/basis_functions/physical_space_element.h>
#include <igatools/utils/multi_array_utils.h>
//#include <igatools/functions/identity_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/grid_function_lib.h>

#include <igatools/io/xml_document.h>
#include <igatools/io/xml_element.h>


#include <fstream>

using std::shared_ptr;
using std::string;
using std::ofstream;
using std::ios;
using std::endl;
using std::to_string;

#include <boost/detail/endian.hpp>


IGA_NAMESPACE_OPEN

namespace
{

template<int dim,int codim>
inline
std::shared_ptr<const Domain<dim,codim> >
create_domain_from_grid(const shared_ptr<const Grid<dim>> &grid)
{
  const int space_dim = dim+codim;
  using F = grid_functions::LinearGridFunction<dim,space_dim>;

  using Grad = typename F::Gradient;
  using Val = typename F::Value;

  Grad A;
  for (int i = 0 ; i < dim ; ++i)
    A[i][i] = Tdouble(1.0);

  Val b;

  auto domain = Domain<dim,codim>::const_create(F::const_create(grid,A,b));
  return domain;
}

}



template<int dim, int codim, class T>
Writer<dim, codim, T>::
Writer(const shared_ptr<const Grid<dim>> &grid,
       const Index num_points_direction)
  :
  Writer(create_domain_from_grid<dim,codim>(grid),
        QUniform<dim>::create(num_points_direction))
{}


template<int dim, int codim, class T>
Writer<dim, codim, T>::
Writer(const std::shared_ptr<const GridFunction<dim,dim+codim>> &grid_function,
       const Index num_points_direction)
  :
  Writer(Domain<dim,codim>::const_create(grid_function),num_points_direction)
{}

template<int dim, int codim, class T>
Writer<dim, codim, T>::
Writer(const std::shared_ptr<const Domain<dim,codim>> &domain,
       const Index num_points_direction)
  :
  Writer(domain,QUniform<dim>::create(num_points_direction))
{}


template<int dim, int codim, class T>
Writer<dim, codim, T>::
Writer(const shared_ptr<const Domain<dim,codim> > &domain,
       const shared_ptr<const Quadrature<dim> > &quadrature)
  :
  domain_(domain),
  quad_plot_(quadrature),
  num_points_direction_(quad_plot_->get_num_coords_direction()),
  n_iga_elements_(domain->get_grid_function()->get_grid()->get_num_all_elems()),
  n_points_per_iga_element_(quad_plot_->get_num_points()),
  n_vtk_points_(n_iga_elements_*n_points_per_iga_element_),
  sizeof_Real_(sizeof(T)),
  sizeof_int_(sizeof(int)),
  sizeof_uchar_(sizeof(unsigned char)),
  offset_(0)
{
  Assert(domain_ != nullptr, ExcNullPtr());
  Assert(quad_plot_ != nullptr, ExcNullPtr());

#if defined( BOOST_LITTLE_ENDIAN )
  byte_order_ = "LittleEndian";
#elif defined( BOOST_BIG_ENDIAN )
  byte_order_ = "BigEndian";
#else
  AssertThrow(false, ExcMessage("Unsupported Endian-ness"));
#endif

  //--------------------------------------------------------------------------
  Assert(sizeof_Real_ == 8 || sizeof_Real_ == 4,
         ExcMessage("The size of the Real type can be only 8 o 4 bytes."));

  if (sizeof_Real_ == 8)
  {
    string_Real_ = "Float64";
    precision_ = 15;
  }
  else if (sizeof_Real_ == 4)
  {
    string_Real_ = "Float32";
    precision_ = 8;
  }
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  Assert(sizeof_int_ == 4,
         ExcMessage("The size of the int type can be only 4 bytes."));
  string_int_ = "UInt32";
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  Assert(sizeof_uchar_ == 1,
         ExcMessage("The size of the unsigned char type can be only 1 byte."));
  string_uchar_ = "UInt8";
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------
  appended_data_ << '_';
  //--------------------------------------------------------------------------


  if (dim == 1)
  {
    vtk_element_type_ = 3; // VTK_LINE
  }
  else if (dim == 2)
  {
    vtk_element_type_ = 9; // VTK_QUAD
  }
  else if (dim == 3)
  {
    vtk_element_type_ = 12; // VTK_HEXAHEDRON
  }

  //--------------------------------------------------------------------------
  n_vtk_elements_per_iga_element_ = 1;
  for (int i = 0; i < dim; i++)
  {
    Assert(num_points_direction_[i] >= 2, ExcLowerRange(num_points_direction_[i], 2));

    num_subelements_direction_[i] = num_points_direction_[i] - 1;

    n_vtk_elements_per_iga_element_ *= num_subelements_direction_[i];
  }
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  n_vtk_elements_ = n_vtk_elements_per_iga_element_ * n_iga_elements_;
  //--------------------------------------------------------------------------
}





template<int dim, int codim, class T>
void Writer<dim, codim, T>::
fill_points_and_connectivity(
  SafeSTLVector<SafeSTLVector<SafeSTLArray<T,3> > > &points_in_iga_elements,
  SafeSTLVector<SafeSTLVector<SafeSTLArray<int,n_vertices_per_vtk_element_> > >
  &vtk_elements_connectivity) const
{
  auto domain_cache_handler = domain_->create_cache_handler();
  domain_cache_handler->template set_flags<dim>(domain_element::Flags::point);


  auto elem = domain_->cbegin();
  auto end  = domain_->cend();

  domain_cache_handler->init_cache(*elem,quad_plot_);

  int elem_id = 0;
  for (; elem != end; ++elem, ++elem_id)
  {
    domain_cache_handler->template fill_cache<dim>(*elem,0);

    this->get_subelements(
      *elem,
      elem_id,
      vtk_elements_connectivity[elem_id],
      points_in_iga_elements[elem_id]);
  }
}



template<int dim, int codim, class T>
void Writer<dim, codim, T>::
get_subelements(
  const DomainElement<dim,codim> &elem,
  const int elem_flat_id,
  SafeSTLVector< SafeSTLArray<int,n_vertices_per_vtk_element_ > > &vtk_elements_connectivity,
  SafeSTLVector< SafeSTLArray<T,3> > &points_phys_iga_element) const
{
  Assert(Size(points_phys_iga_element.size()) == n_points_per_iga_element_,
         ExcDimensionMismatch(points_phys_iga_element.size(), n_points_per_iga_element_));

  Assert(Size(vtk_elements_connectivity.size())== n_vtk_elements_per_iga_element_,
         ExcDimensionMismatch(vtk_elements_connectivity.size(), n_vtk_elements_per_iga_element_));


  auto element_vertices_tmp = elem.template get_points<dim>(0);

  const T zero = T(0.0);

  const int space_dim = dim + codim;

  // here we evaluate the position of the evaluation points in the physical domain
  for (int ipt = 0; ipt < n_points_per_iga_element_; ++ipt)
  {
    for (int i = 0; i < space_dim; ++i)
      points_phys_iga_element[ipt][i] = element_vertices_tmp[ipt][i];

    for (int i = space_dim; i < 3; ++i)
      points_phys_iga_element[ipt][i] = zero;
  }


  const int iga_element_id = elem_flat_id;

  SafeSTLVector< SafeSTLArray<int,dim> > delta_idx(n_vertices_per_vtk_element_);

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
  //--------------------------------------------------------------------------



  TensorIndex<dim> weight_points =
    MultiArrayUtils< dim >::compute_weight(num_points_direction_);



  //--------------------------------------------------------------------------
  // grid defining the vtk elements inside the iga element

  const auto vtk_elements_grid = Grid<dim>::create(num_points_direction_);
  auto vtk_elem = vtk_elements_grid->begin();
  const auto vtk_elem_end = vtk_elements_grid->end();

  int vtk_vertex_id_offset = n_points_per_iga_element_ * iga_element_id;
  int vtk_elem_flat_id = 0;
  for (; vtk_elem != vtk_elem_end; ++vtk_elem, ++vtk_elem_flat_id)
  {
    auto vtk_elem_tensor_idx = vtk_elem->get_index().get_tensor_index();

    for (int iVertex = 0; iVertex < n_vertices_per_vtk_element_; ++iVertex)
    {
      TensorIndex<dim> vtk_vertex_tensor_idx;
      for (int i = 0; i < dim; ++i)
        vtk_vertex_tensor_idx[i] = vtk_elem_tensor_idx[i] + delta_idx[iVertex][i];

      const int vtk_vertex_local_id = MultiArrayUtils<dim>::tensor_to_flat_index(vtk_vertex_tensor_idx, weight_points);

      vtk_elements_connectivity[vtk_elem_flat_id][iVertex] = vtk_vertex_local_id + vtk_vertex_id_offset;
    }
  }
  //--------------------------------------------------------------------------
}



template<int dim, int codim, class T>
void
Writer<dim, codim, T>::
add_element_data(const SafeSTLVector<double> &element_data,
                 const std::string &name)
{
  cell_data_double_.emplace_back(CellData<double>(element_data, name));

  const string type = cell_data_double_.back().type_;
  if (type == "scalar")
  {
    names_cell_data_scalar_.emplace_back(name);
  }
  else if (type == "vector")
  {
    names_cell_data_vector_.emplace_back(name);
  }
  else if (type == "tensor")
  {
    names_cell_data_tensor_.emplace_back(name);
  }
}



template<int dim, int codim, class T>
void
Writer<dim, codim, T>::
add_element_data(const SafeSTLVector<int> &element_data,
                 const std::string &name)
{
  cell_data_int_.emplace_back(CellData<int>(element_data, name));

  const string type = cell_data_int_.back().type_;
  if (type == "scalar")
  {
    names_cell_data_scalar_.emplace_back(name);
  }
  else if (type == "vector")
  {
    names_cell_data_vector_.emplace_back(name);
  }
  else if (type == "tensor")
  {
    names_cell_data_tensor_.emplace_back(name);
  }
}



template<int dim, int codim, class T>
void
Writer<dim, codim, T>::
add_point_data(const int n_values_per_point,
               const std::string &type,
               const SafeSTLVector<SafeSTLVector<SafeSTLVector<T>>> &data_iga_elements,
               const std::string &name)
{
  Assert(data_iga_elements.size() == n_iga_elements_,
         ExcDimensionMismatch(data_iga_elements.size(), n_iga_elements_));
  Assert(type == "scalar" || type == "vector" || type == "tensor",
         ExcMessage("The point_data type can only be \"scalar\", \"vector\" or \"tensor\" (and not \"" + type + "\")"));

  shared_ptr<SafeSTLVector<T>> data_ptr(new SafeSTLVector<T>(n_iga_elements_ * n_points_per_iga_element_ * n_values_per_point));
  auto &data = *data_ptr;

  Index pos = 0;
  for (const auto &data_element : data_iga_elements)
  {
    Assert(data_element.size() == n_points_per_iga_element_,
           ExcDimensionMismatch(data_element.size(), n_points_per_iga_element_));

    for (const auto &data_point : data_element)
    {
      Assert(data_point.size() == n_values_per_point,
             ExcDimensionMismatch(data_point.size(), n_values_per_point));

      for (const double &value : data_point)
      {
        data[pos++] = value;
      }
    }
  }
  fields_.emplace_back(PointData(name,type,n_iga_elements_,
                                 n_points_per_iga_element_,
                                 n_values_per_point,
                                 data_ptr));

  if (type == "scalar")
  {
    names_point_data_scalar_.emplace_back(name);
  }
  else if (type == "vector")
  {
    names_point_data_vector_.emplace_back(name);
  }
  else if (type == "tensor")
  {
    names_point_data_tensor_.emplace_back(name);
  }

}





template<int dim, int codim, class T>
void Writer<dim, codim, T>::
save(const string &filename, const bool pretty_print) const
{

  SafeSTLVector< SafeSTLVector< SafeSTLArray<T,3> > >
  points_in_iga_elements(n_iga_elements_, SafeSTLVector< SafeSTLArray<T,3> >(n_points_per_iga_element_));

  SafeSTLVector< SafeSTLVector< SafeSTLArray< int, n_vertices_per_vtk_element_> > >
  vtk_elements_connectivity(n_iga_elements_);
  for (auto &iga_elem_connectivity : vtk_elements_connectivity)
    iga_elem_connectivity.resize(n_vtk_elements_per_iga_element_);

  this->fill_points_and_connectivity(points_in_iga_elements, vtk_elements_connectivity);

  const string vtu_filename = filename + ".vtu";
  const auto xml_doc = this->save_xml(points_in_iga_elements, vtk_elements_connectivity);
  xml_doc->write_to_file(vtu_filename, pretty_print);
}



template<int dim, int codim, class T>
shared_ptr<XMLDocument>
Writer<dim, codim, T>::
save_xml(const SafeSTLVector< SafeSTLVector< SafeSTLArray<T,3> > > &points_in_iga_elements,
         const SafeSTLVector< SafeSTLVector< SafeSTLArray<int,n_vertices_per_vtk_element_> > >
         &vtk_elements_connectivity) const
{
  const auto xml_doc = XMLDocument::create_void_document("VTKFile");

  const auto xml_elem = xml_doc->get_document_element();
  xml_elem->add_attribute("type", "UnstructuredGrid");
  // xml_elem->add_attribute("type_order", byte_order_);

  const auto vtu_elem = xml_doc->create_new_element("UnstructuredGrid");

  const auto piece_elem = xml_doc->create_new_element("Piece");
  piece_elem->add_attribute("NumberOfPoints", n_vtk_points_);
  piece_elem->add_attribute("NumberOfCells", n_vtk_elements_);
  const auto pts_elem = xml_doc->create_new_element("Points");

  SafeSTLVector<T> pt_coords(3 * n_iga_elements_ * n_points_per_iga_element_);

  Index i = 0;
  for (const auto &point_in_iga_element : points_in_iga_elements)
    for (const auto &point : point_in_iga_element)
    {
      pt_coords[i++] = point[0];
      pt_coords[i++] = point[1];
      pt_coords[i++] = point[2];
    }
  const auto data1 = xml_doc->create_vector_element("DataArray", pt_coords, precision_);
  data1->add_attribute("type", string_Real_);
  data1->add_attribute("NumberOfComponents", 3);
  data1->add_attribute("format", "ascii");
  pts_elem->append_child_element(data1);

  piece_elem->append_child_element(pts_elem);
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  // writing the <Cells> section
  const auto cells_elem = xml_doc->create_new_element("Cells");

  // Writing connectivity
  SafeSTLVector<int> conn_vec(n_iga_elements_ *
                              n_vtk_elements_per_iga_element_ *
                              n_vertices_per_vtk_element_);
  i = 0;
  for (const auto &iga_elem_connectivity : vtk_elements_connectivity)
    for (const auto &vtk_elem_connectivity : iga_elem_connectivity)
      for (const auto &point_id : vtk_elem_connectivity)
        conn_vec[i++] = point_id;

  const auto data2 = xml_doc->create_vector_element("DataArray", conn_vec, precision_);
  data2->add_attribute("Name", "connectivity");
  data2->add_attribute("type", string_int_);
  data2->add_attribute("format", "ascii");
  cells_elem->append_child_element(data2);


  // Writing offsets
  i = 0;
  SafeSTLVector<int> offsets_vec(n_vtk_elements_);
  for (int vtk_elem_id = 1; vtk_elem_id <= n_vtk_elements_; ++vtk_elem_id)
    offsets_vec[i++] = n_vertices_per_vtk_element_ * vtk_elem_id;
  const auto data3 = xml_doc->create_vector_element("DataArray", offsets_vec, precision_);
  data3->add_attribute("Name", "offsets");
  data3->add_attribute("type", string_int_);
  data3->add_attribute("format", "ascii");
  cells_elem->append_child_element(data3);


  // Writing types
  i = 0;
  SafeSTLVector<int> types_vec(n_vtk_elements_);
  for (int vtk_elem_id = 1; vtk_elem_id <= n_vtk_elements_; ++vtk_elem_id)
    types_vec[i++] = static_cast<int>(vtk_element_type_);
  const auto data4 = xml_doc->create_vector_element("DataArray", types_vec, precision_);
  data4->add_attribute("Name", "types");
  data4->add_attribute("type", string_uchar_);
  data4->add_attribute("format", "ascii");
  cells_elem->append_child_element(data4);

  piece_elem->append_child_element(cells_elem);
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  // writing the <PointData> section
  const auto pt_data_elem = xml_doc->create_new_element("PointData");

  // writing the <PointData> section
  string point_data_optional_attr;
  if (!names_point_data_scalar_.empty())
  {
    string attr = "";
    for (const string &name : names_point_data_scalar_)
      attr += name + " ";
    pt_data_elem->add_attribute("Scalars", attr);
  }
  if (!names_point_data_vector_.empty())
  {
    string attr = "";
    for (const string &name : names_point_data_vector_)
      attr += name + " ";
    pt_data_elem->add_attribute("Vectors", attr);
  }
  if (!names_point_data_tensor_.empty())
  {
    string attr = "";
    for (const string &name : names_point_data_tensor_)
      attr += name + " ";
    pt_data_elem->add_attribute("Tensors", attr);
  }

  for (const auto &point_data : fields_)
  {
    const auto data = xml_doc->create_vector_element("DataArray", *(point_data.values_), precision_);
    data->add_attribute("Name", point_data.name_);
    data->add_attribute("NumberOfComponents", point_data.num_components_);
    data->add_attribute("type", string_Real_);
    data->add_attribute("format", "ascii");
    pt_data_elem->append_child_element(data);
  }
  piece_elem->append_child_element(pt_data_elem);
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  // writing the <CellData> section
  const auto cell_data_elem = xml_doc->create_new_element("CellData");

  if (!names_cell_data_scalar_.empty())
  {
    string attr = "";
    for (const string &name : names_cell_data_scalar_)
      attr += name + " ";
    cell_data_elem->add_attribute("Scalars", attr);
  }
  if (!names_cell_data_vector_.empty())
  {
    string attr = "";
    for (const string &name : names_cell_data_vector_)
      attr += name + " ";
    cell_data_elem->add_attribute("Vectors", attr);
  }
  if (!names_cell_data_tensor_.empty())
  {
    string attr = "";
    for (const string &name : names_cell_data_tensor_)
      attr += name + " ";
    cell_data_elem->add_attribute("Tensors", attr);
  }

  for (const auto &cell_data : cell_data_double_)
  {
    const auto data = xml_doc->create_vector_element("DataArray", *(cell_data.values_), precision_);
    data->add_attribute("Name", cell_data.name_);
    data->add_attribute("NumberOfComponents", cell_data.num_components_);
    data->add_attribute("type", string_Real_);
    data->add_attribute("format", "ascii");
    pt_data_elem->append_child_element(data);
  }

  for (const auto &cell_data : cell_data_int_)
  {
    const auto data = xml_doc->create_vector_element("DataArray", *(cell_data.values_), precision_);
    data->add_attribute("Name", cell_data.name_);
    data->add_attribute("NumberOfComponents", cell_data.num_components_);
    data->add_attribute("type", string_int_);
    data->add_attribute("format", "ascii");
    pt_data_elem->append_child_element(data);
  }
  piece_elem->append_child_element(cell_data_elem);
  //--------------------------------------------------------------------------

  vtu_elem->append_child_element(piece_elem);
  xml_elem->append_child_element(vtu_elem);

  return xml_doc;
}



template<int dim, int codim, class T>
void
Writer<dim, codim, T>::
print_info(LogStream &out) const
{
  SafeSTLVector< SafeSTLVector< SafeSTLArray<T,3> > >
  points_in_iga_elements(n_iga_elements_, SafeSTLVector< SafeSTLArray<T,3> >(n_points_per_iga_element_));

  SafeSTLVector< SafeSTLVector< SafeSTLArray< int, n_vertices_per_vtk_element_> > >
  vtk_elements_connectivity(n_iga_elements_);
  for (auto &iga_elem_connectivity : vtk_elements_connectivity)
    iga_elem_connectivity.resize(n_vtk_elements_per_iga_element_);

  this->fill_points_and_connectivity(points_in_iga_elements, vtk_elements_connectivity);

  const auto xml_doc = save_xml(points_in_iga_elements, vtk_elements_connectivity);
  xml_doc->print_info(out);
}


template<int dim, int codim, class T>
int
Writer<dim, codim, T>::
get_num_points_per_iga_element() const
{
  return n_points_per_iga_element_;
}


template<int dim, int codim, class T>
int
Writer<dim, codim, T>::
get_num_vtk_elements_per_iga_element() const
{
  return n_vtk_elements_per_iga_element_;
}




template<int dim, int codim, class T>
int
Writer<dim, codim, T>::
get_num_iga_elements() const
{
  return n_iga_elements_;
}



template<int dim, int codim, class T>
int
Writer<dim, codim, T>::
get_num_vtk_elements() const
{
  return n_vtk_elements_;
}


IGA_NAMESPACE_CLOSE

#include <igatools/io/writer.inst>

#endif // XML_IO

