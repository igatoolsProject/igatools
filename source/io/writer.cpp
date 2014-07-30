//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
#include <igatools/base/exceptions.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/utils/multi_array_utils.h>

#include <vector>
#include <sstream>
#include <fstream>
#include <utility>

using std::vector;
using std::array;
using std::shared_ptr;
using std::make_shared;
using std::string;
using std::stringstream;
using std::ofstream;
using std::fstream;
using std::ios;
using std::pair;
using std::to_string;
using std::endl;

#include <boost/detail/endian.hpp>

IGA_NAMESPACE_OPEN

//TODO: Add patch id as a cell field

template<int dim, int codim, class T>
Writer<dim, codim, T>::
Writer(const shared_ptr<Grid> grid)
    :
    Writer(IdentityMapping<dim, codim>::create(grid),
           shared_ptr< QUniform<dim> >(new QUniform<dim>(2)))
{}



template<int dim, int codim, class T>
Writer<dim, codim, T>::
Writer(const shared_ptr<Grid> grid,
       const Index n_points_direction = 2)
    :
    Writer(IdentityMapping<dim, codim>::create(grid),
           shared_ptr< QUniform<dim> >(new QUniform<dim>(n_points_direction)))
{}



template<int dim, int codim, class T>
Writer<dim, codim, T>::
Writer(const shared_ptr<const Map> map,
       const Index n_points_direction = 2)
    :
    Writer(map,
           shared_ptr< QUniform<dim> >(new QUniform<dim>(n_points_direction)))
{}



template<int dim, int codim, class T>
Writer<dim, codim, T>::
Writer(const shared_ptr<const Mapping<dim,codim> > map,
       const shared_ptr<const Quadrature<dim> > quadrature)
    :
    grid_(map->get_grid()),
    map_(map),
    quad_plot_(*quadrature),
    num_points_direction_(quad_plot_.get_num_points_direction()),
    n_iga_elements_(grid_->get_num_elements()),
    n_points_per_iga_element_(quad_plot_.get_num_points()),
    n_vtk_points_(n_iga_elements_*n_points_per_iga_element_),
    sizeof_Real_(sizeof(T)),
    sizeof_int_(sizeof(int)),
    sizeof_uchar_(sizeof(unsigned char)),
    offset_(0)
{
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
int
Writer<dim, codim, T>::
get_num_iga_elements() const
{
    return n_iga_elements_;
}



template<int dim, int codim, class T>
int Writer<dim, codim, T>::get_num_vtk_elements() const
{
    return n_vtk_elements_;
}



template<int dim, int codim, class T>
int Writer<dim, codim, T>::get_num_points_per_iga_element() const
{
    return n_points_per_iga_element_;
}



template<int dim, int codim, class T>
int Writer<dim, codim, T>::get_num_vtk_elements_per_iga_element() const
{
    return n_vtk_elements_per_iga_element_;
}


template<int dim, int codim, class T>
void Writer<dim, codim, T>::
add_point_data(const int n_values_per_point,
               const std::string &type,
               const std::vector<std::vector<std::vector<T>>> &data_iga_elements,
               const std::string &name)
{
    Assert(data_iga_elements.size() == n_iga_elements_,
           ExcDimensionMismatch(data_iga_elements.size(), n_iga_elements_));
    Assert(type == "scalar" || type == "vector" || type == "tensor",
           ExcMessage("The point_data type can only be \"scalar\", \"vector\" or \"tensor\" (and not \"" + type + "\")"));

    shared_ptr<vector<T>> data_ptr(new vector<T>(n_iga_elements_ * n_points_per_iga_element_ * n_values_per_point));
    auto &data = *data_ptr;

    int pos = 0;
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
template<class Space, LAPack la_pack>
void Writer<dim, codim, T>::
add_field(shared_ptr<Space> space_,
          const Vector<la_pack> &coefs,
          const string &name)
{
    // Compromise to keep type safe but avoid the user for writing
    // pedantically correct but comprehensible undesirable casting
    shared_ptr<const Space> space = std::const_pointer_cast<const Space> (space_);

    //--------------------------------------------------------------------------
    Assert(space_dim <= 3,
           ExcMessage("The maximum allowed physical domain for VTK file is 3."));
    Assert(space->get_num_basis() == coefs.size(),
           ExcDimensionMismatch(space->get_num_basis(), coefs.size()));
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // get the fields to write and assign them to the vtkUnstructuredGrid object

    auto element     = space->begin();
    auto element_end = space->end();

    element->init_cache(ValueFlags::value, quad_plot_);



    const int n_elements = grid_->get_num_elements();
    const int n_pts_per_elem = quad_plot_.get_num_points();

    static const int dim_phys_range = Space::range;
    static const int rank = Space::rank;


    const int n_values_per_pt =
        dim_phys_range == 1 ? 1 : std::pow(dim_phys_range, rank);
    shared_ptr< vector<T> > data_ptr(new vector<T>(n_elements * n_pts_per_elem * n_values_per_pt));
    auto &data = *data_ptr;
    if (rank == 0)
    {
        int pos = 0;
        for (int iElement = 0; element != element_end; ++element, ++iElement)
        {
            element->fill_cache();
            const auto field_values = element->evaluate_field(
                                          coefs.get_local_coefs(element->get_local_to_global()));

            for (int iPt = 0; iPt < n_pts_per_elem; ++iPt)
                data[pos++] = field_values[iPt][0];
        }

        fields_.emplace_back(PointData(name,"scalar",n_elements,n_pts_per_elem, n_values_per_pt, data_ptr));
        names_point_data_scalar_.emplace_back(name);
    }
    else if (rank == 1)
    {
        int pos = 0;
        for (int iElement = 0; element != element_end; ++element, ++iElement)
        {
            element->fill_cache();

            const auto field_values = element->evaluate_field(
                                          coefs.get_local_coefs(element->get_local_to_global()));

            for (int iPt = 0; iPt < n_pts_per_elem; ++iPt)
            {
                const auto &field_value_ipt = field_values[ iPt ];
                for (int i = 0; i < dim_phys_range; ++i)
                    data[pos++] = field_value_ipt[i];
            }
        }

        fields_.emplace_back(PointData(name,"vector",n_elements,n_pts_per_elem, n_values_per_pt, data_ptr));
        names_point_data_vector_.emplace_back(name);
    }
    else if (rank == 2)
    {
        int pos = 0;
        for (int iElement = 0; element != element_end; ++element, ++iElement)
        {
            element->fill_cache();

            const auto field_values = element->evaluate_field(
                                          coefs.get_local_coefs(element->get_local_to_global()));

            for (int iPt = 0; iPt < n_pts_per_elem; ++iPt)
            {
                const auto &field_value_ipt = field_values[ iPt ];
                for (int i = 0; i < dim_phys_range; ++i)
                {
                    const auto &field_value_ipt_i = field_value_ipt[i];

                    for (int j = 0; j < dim_phys_range; ++j)
                        data[pos++] = field_value_ipt_i[j];
                }
            }
        }

        fields_.emplace_back(PointData(name,"tensor",n_elements,n_pts_per_elem,n_values_per_pt,data_ptr));
        names_point_data_tensor_.emplace_back(name);
    }

    //--------------------------------------------------------------------------
}



template<int dim, int codim, class T>
void Writer<dim, codim, T>::fill_points_and_connectivity(
    std::vector< std::vector< std::array<T,3> > > &points_in_iga_elements,
    std::vector< std::vector< std::array< int,n_vertices_per_vtk_element_> > >
    &vtk_elements_connectivity) const
{

    auto element = map_->begin();
    const auto element_end = map_->end();

    element->init_cache(ValueFlags::map_value, quad_plot_);

    for (; element != element_end; ++element)
    {
        const int iga_elem_id = element->get_flat_index();

        element->fill_cache();
        get_subelements(element,
                        vtk_elements_connectivity[iga_elem_id],
                        points_in_iga_elements[iga_elem_id]);
    }
}



template<int dim, int codim, class T>
void Writer<dim, codim, T>::
get_subelements(
    const typename Mapping< dim, codim>::ElementIterator elem,
    vector< array< int, n_vertices_per_vtk_element_ > > &vtk_elements_connectivity,
    vector< array<T,3> > &points_phys_iga_element) const
{

    Assert(Size(points_phys_iga_element.size()) == n_points_per_iga_element_,
           ExcDimensionMismatch(points_phys_iga_element.size(), n_points_per_iga_element_));

    Assert(Size(vtk_elements_connectivity.size())== n_vtk_elements_per_iga_element_,
           ExcDimensionMismatch(vtk_elements_connectivity.size(), n_vtk_elements_per_iga_element_));


    auto element_vertices_tmp = elem->get_map_values();

    const T zero = T(0.0);

    // here we evaluate the position of the evaluation points in the physical domain
    for (int ipt = 0; ipt < n_points_per_iga_element_; ++ipt)
    {
        for (int i = 0; i < space_dim; ++i)
            points_phys_iga_element[ipt][i] = element_vertices_tmp[ipt][i];

        for (int i = space_dim; i < 3; ++i)
            points_phys_iga_element[ipt][i] = zero;
    }


    const int iga_element_id = elem->get_flat_index();

    vector< array<int,dim> > delta_idx(n_vertices_per_vtk_element_);


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

    const auto  vtk_elements_grid = CartesianGrid<dim>::create(num_points_direction_);
    auto vtk_elem = vtk_elements_grid->begin();
    const auto vtk_elem_end = vtk_elements_grid->end();

    int vtk_vertex_id_offset = n_points_per_iga_element_ * iga_element_id;
    for (; vtk_elem != vtk_elem_end; ++vtk_elem)
    {
        int vtk_elem_flat_id = vtk_elem->get_flat_index();
        array<Index,dim> vtk_elem_tensor_idx = vtk_elem->get_tensor_index();

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
void Writer<dim, codim, T>::
add_element_data(const std::vector<double> &element_data,
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
void Writer<dim, codim, T>::
add_element_data(const std::vector<int> &element_data,
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
template<class Out>
void Writer<dim, codim, T>::save_ascii(Out &file,
                                       const std::vector< std::vector< std::array<T,3> > > &points_in_iga_elements,
                                       const std::vector< std::vector< std::array< int,n_vertices_per_vtk_element_> > >
                                       &vtk_elements_connectivity) const
{
    const string tab1("\t");
    const string tab2 = tab1 + tab1;
    const string tab3 = tab2 + tab1;
    const string tab4 = tab3 + tab1;
    const string tab5 = tab4 + tab1;

    file << "<?xml version=\"1.0\"?>" << endl;
    file << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"" << byte_order_ << "\">" << endl;

    file << tab1 << "<UnstructuredGrid>" << endl;

    file << tab2 << "<Piece NumberOfPoints=\"" << to_string(n_vtk_points_) << "\" NumberOfCells=\""<< to_string(n_vtk_elements_) << "\">" << endl;

    file << tab3 << "<Points>" << endl;
    file << tab4 << "<DataArray type=\"" << string_Real_ << "\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

    for (const auto &point_in_iga_element : points_in_iga_elements)
        for (const auto &point : point_in_iga_element)
            file << tab5 << point[0] << " " << point[1] << " " << point[2] << endl;

    file << tab4 << "</DataArray>" << endl;
    file << tab3 << "</Points>" << endl;

    file << tab3 << "<Cells>" << endl;
    file << tab4 << "<DataArray Name=\"connectivity\" type=\"" << string_int_ << "\" format=\"ascii\">" << endl;
    file << tab5;
    for (const auto &iga_elem_connectivity : vtk_elements_connectivity)
        for (const auto &vtk_elem_connectivity : iga_elem_connectivity)
            for (const auto &point_id : vtk_elem_connectivity)
                file << point_id << " ";
    file << endl;
    file << tab4 << "</DataArray>" << endl;

    file << tab4 << "<DataArray Name=\"offsets\" type=\"" << string_int_ << "\" format=\"ascii\">" << endl;
    file << tab5;
    for (int vtk_elem_id = 1; vtk_elem_id <= n_vtk_elements_; ++vtk_elem_id)
        file << n_vertices_per_vtk_element_ * vtk_elem_id << " ";
    file << endl;
    file << tab4 << "</DataArray>" << endl;

    file << tab4 << "<DataArray Name=\"types\" type=\"" << string_uchar_ << "\" format=\"ascii\">" << endl;
    file << tab5;
    for (int vtk_elem_id = 1; vtk_elem_id <= n_vtk_elements_; ++vtk_elem_id)
        file << static_cast<int>(vtk_element_type_) << " ";
    file << endl;
    file << tab4 << "</DataArray>" << endl;
    file << tab3 << "</Cells>" << endl;


    //--------------------------------------------------------------------------
    // writing the <PointData> section
    string point_data_optional_attr;
    if (!names_point_data_scalar_.empty())
    {
        point_data_optional_attr += " Scalars=\"";
        for (const string &name : names_point_data_scalar_)
            point_data_optional_attr += name + " ";
        point_data_optional_attr+= "\"";
    }
    if (!names_point_data_vector_.empty())
    {
        point_data_optional_attr += " Vectors=\"";
        for (const string &name : names_point_data_vector_)
            point_data_optional_attr += name + " ";
        point_data_optional_attr+= "\"";
    }
    if (!names_point_data_tensor_.empty())
    {
        point_data_optional_attr += " Tensors=\"";
        for (const string &name : names_point_data_tensor_)
            point_data_optional_attr += name + " ";
        point_data_optional_attr+= "\"";
    }

    file << tab3 << "<PointData" << point_data_optional_attr << ">" << endl;
    for (const auto &point_data : fields_)
    {
        file << tab4 << "<DataArray Name=\"" << point_data.name_
             << "\" type=\"" << string_Real_
             << "\" NumberOfComponents=\""<< point_data.num_components_
             << "\" format=\"ascii\">" << endl;

        file << tab5;
        for (const auto &v : *point_data.values_)
            file << v << " ";
        file << endl;

        file << tab4 << "</DataArray>" << endl;
    }
    file << tab3 << "</PointData>" << endl;
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // writing the <CellData> section
    string cell_data_optional_attr;
    if (!names_cell_data_scalar_.empty())
    {
        cell_data_optional_attr += " Scalars=\"";
        for (const string &name : names_cell_data_scalar_)
            cell_data_optional_attr += name + " ";
        cell_data_optional_attr+= "\"";
    }
    if (!names_cell_data_vector_.empty())
    {
        cell_data_optional_attr += " Vectors=\"";
        for (const string &name : names_cell_data_vector_)
            cell_data_optional_attr += name + " ";
        cell_data_optional_attr+= "\"";
    }
    if (!names_cell_data_tensor_.empty())
    {
        cell_data_optional_attr += " Tensors=\"";
        for (const string &name : names_cell_data_tensor_)
            cell_data_optional_attr += name + " ";
        cell_data_optional_attr+= "\"";
    }

    file << tab3 << "<CellData" << cell_data_optional_attr << ">" << endl;
    for (const auto &cell_data : cell_data_double_)
    {
        file << tab4 << "<DataArray Name=\"" << cell_data.name_
             << "\" type=\"" << string_Real_
             << "\" NumberOfComponents=\""<< cell_data.num_components_
             << "\" format=\"ascii\">" << endl;
        file << tab5;
        for (const double &v : *cell_data.values_)
            file << v << " ";
        file << endl;
        file << tab4 << "</DataArray>" << endl;
    }
    for (const auto &cell_data : cell_data_int_)
    {
        file << tab4 << "<DataArray Name=\"" << cell_data.name_
             << "\" type=\"" << string_int_
             << "\" NumberOfComponents=\""<< cell_data.num_components_
             << "\" format=\"ascii\">" << endl;
        file << tab5;
        for (const int &v : *cell_data.values_)
            file << v << " ";
        file << endl;
        file << tab4 << "</DataArray>" << endl;
    }
    file << tab3 << "</CellData>" << endl;
    //--------------------------------------------------------------------------

    file << tab2 << "</Piece>" << endl;
    file << tab1 << "</UnstructuredGrid>" << endl;
    file << "</VTKFile>";
}



template<int dim, int codim, class T>
void Writer<dim, codim, T>::save_appended(const string &filename,
                                          const std::vector< std::vector< std::array<T,3> > > &points_in_iga_elements,
                                          const std::vector< std::vector< std::array< int,n_vertices_per_vtk_element_> > >
                                          &vtk_elements_connectivity) const
{
    ofstream file(filename);
    file.setf(ios::scientific);
    file.precision(precision_);

    const string tab1("\t");
    const string tab2 = tab1 + tab1;
    const string tab3 = tab2 + tab1;
    const string tab4 = tab3 + tab1;
    const string tab5 = tab4 + tab1;

    int offset = 0;


    file << "<?xml version=\"1.0\"?>" << endl;
    file << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"" << byte_order_ << "\">" << endl;

    file << tab1 << "<UnstructuredGrid>" << endl;

    file << tab2 << "<Piece NumberOfPoints=\"" << to_string(n_vtk_points_) << "\" NumberOfCells=\""<< to_string(n_vtk_elements_) << "\">" << endl;

    file << tab3 << "<Points>" << endl;
    file << tab4 << "<DataArray type=\"" << string_Real_ << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;
    const int n_bytes_points = n_vtk_points_ * 3 * sizeof_Real_;
    offset += sizeof_int_ + n_bytes_points;
    file << tab3 << "</Points>" << endl;


    file << tab3 << "<Cells>" << endl;
    file << tab4 << "<DataArray Name=\"connectivity\" type=\"" << string_int_ << "\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;
    const int n_bytes_connectivity = n_vtk_elements_ * n_vertices_per_vtk_element_ * sizeof_int_;
    offset += sizeof_int_ + n_bytes_connectivity;


    file << tab4 << "<DataArray Name=\"offsets\" type=\"" << string_int_ << "\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;
    const int n_bytes_offsets = n_vtk_elements_ * sizeof_int_;
    offset += sizeof_int_ + n_bytes_offsets;


    file << tab4 << "<DataArray Name=\"types\" type=\"" << string_uchar_ << "\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;
    const int n_bytes_types = n_vtk_elements_ * sizeof_uchar_;
    offset += sizeof_int_ + n_bytes_types;
    file << tab3 << "</Cells>" << endl;


    //--------------------------------------------------------------------------
    // writing the <PointData> section
    string point_data_optional_attr;
    if (!names_point_data_scalar_.empty())
    {
        point_data_optional_attr += " Scalars=\"";
        for (const string &name : names_point_data_scalar_)
            point_data_optional_attr += name + " ";
        point_data_optional_attr+= "\"";
    }
    if (!names_point_data_vector_.empty())
    {
        point_data_optional_attr += " Vectors=\"";
        for (const string &name : names_point_data_vector_)
            point_data_optional_attr += name + " ";
        point_data_optional_attr+= "\"";
    }
    if (!names_point_data_tensor_.empty())
    {
        point_data_optional_attr += " Tensors=\"";
        for (const string &name : names_point_data_tensor_)
            point_data_optional_attr += name + " ";
        point_data_optional_attr+= "\"";
    }

    vector<int> n_bytes_point_data;
    file << tab3 << "<PointData" << point_data_optional_attr << ">" << endl;
    for (const auto &point_data : fields_)
    {

        file << tab4 << "<DataArray Name=\"" << point_data.name_
             << "\" type=\"" << string_Real_
             << "\" NumberOfComponents=\""<< point_data.num_components_
             << "\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;

        n_bytes_point_data.emplace_back(point_data.values_->size() * sizeof_Real_);
        offset += sizeof_int_ + n_bytes_point_data.back();
    }
    file << tab3 << "</PointData>" << endl;
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // writing the <CellData> section
    string cell_data_optional_attr;
    if (!names_cell_data_scalar_.empty())
    {
        cell_data_optional_attr += " Scalars=\"";
        for (const string &name : names_cell_data_scalar_)
            cell_data_optional_attr += name + " ";
        cell_data_optional_attr+= "\"";
    }
    if (!names_cell_data_vector_.empty())
    {
        cell_data_optional_attr += " Vectors=\"";
        for (const string &name : names_cell_data_vector_)
            cell_data_optional_attr += name + " ";
        cell_data_optional_attr+= "\"";
    }
    if (!names_cell_data_tensor_.empty())
    {
        cell_data_optional_attr += " Tensors=\"";
        for (const string &name : names_cell_data_tensor_)
            cell_data_optional_attr += name + " ";
        cell_data_optional_attr+= "\"";
    }

    file << tab3 << "<CellData" << cell_data_optional_attr << ">" << endl;

    vector<int> n_bytes_cell_data_double;
    for (const auto &cell_data : cell_data_double_)
    {
        file << tab4 << "<DataArray Name=\"" << cell_data.name_
             << "\" type=\"" << string_Real_
             << "\" NumberOfComponents=\""<< cell_data.num_components_
             << "\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;

        n_bytes_cell_data_double.emplace_back(cell_data.values_->size() * sizeof_Real_);
        offset += sizeof_int_ + n_bytes_cell_data_double.back();
    }

    vector<int> n_bytes_cell_data_int;
    for (const auto &cell_data : cell_data_int_)
    {
        file << tab4 << "<DataArray Name=\"" << cell_data.name_
             << "\" type=\"" << string_int_
             << "\" NumberOfComponents=\""<< cell_data.num_components_
             << "\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;

        n_bytes_cell_data_int.emplace_back(cell_data.values_->size() * sizeof_int_);
        offset += sizeof_int_ + n_bytes_cell_data_int.back();
    }
    file << tab3 << "</CellData>" << endl;
    //--------------------------------------------------------------------------



    file << tab2 << "</Piece>" << endl;

    file << tab1 << "</UnstructuredGrid>" << endl;


    file << tab1 << "<AppendedData encoding=\"raw\">" << endl;
    file << tab2 << "_";

    //--------------------------------------------------------------------------
    // writing the points coordinate
    file.write((char *) &n_bytes_points, sizeof_int_);
    for (const auto &point_in_iga_element : points_in_iga_elements)
        for (const auto &point : point_in_iga_element)
            file.write((char *) &point[0], 3 * sizeof_Real_);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // writing the element connectivity
    file.write((char *) &n_bytes_connectivity, sizeof_int_);
    for (const auto &iga_elem_connectivity : vtk_elements_connectivity)
        for (const auto &vtk_elem_connectivity : iga_elem_connectivity)
            file.write((char *) vtk_elem_connectivity.data(), n_vertices_per_vtk_element_ * sizeof_int_);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // writing the element offsets
    file.write((char *) &n_bytes_offsets, sizeof_int_);
    for (int i = 1; i <= n_vtk_elements_; ++i)
    {
        const int tmp = i * n_vertices_per_vtk_element_;
        file.write((char *) &tmp, sizeof_int_);
    }
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // writing the element types
    file.write((char *) &n_bytes_types, sizeof_int_);
    for (int vtk_elem_id = 1; vtk_elem_id <= n_vtk_elements_; ++vtk_elem_id)
        file.write((char *) &vtk_element_type_, sizeof_uchar_);
    //--------------------------------------------------------------------------



    //--------------------------------------------------------------------------
    // writing the point data
    const int n_point_data = fields_.size();
    for (int i = 0; i < n_point_data; ++i)
    {
        file.write((char *) &n_bytes_point_data[i], sizeof_int_);
        file.write((char *) fields_[i].values_->data(), n_bytes_point_data[i]);
    }
    //--------------------------------------------------------------------------



    //--------------------------------------------------------------------------
    // writing the cell data (double)
    const int n_cell_data_double = cell_data_double_.size();
    for (int i = 0; i < n_cell_data_double; ++i)
    {
        file.write((char *) &n_bytes_cell_data_double[i], sizeof_int_);

        const int n_values = cell_data_double_[i].values_->size();
        //here we convert the type double in CellData.values_ to type T
        vector<T> buffer(n_values);
        for (int j = 0; j < n_values; ++j)
            buffer[j] = (*cell_data_double_[i].values_)[j];

        file.write((char *) buffer.data(), n_bytes_cell_data_double[i]);
    }
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // writing the cell data (int)
    const int n_cell_data_int = cell_data_int_.size();
    for (int i = 0; i < n_cell_data_int; ++i)
    {
        file.write((char *) &n_bytes_cell_data_int[i], sizeof_int_);
        file.write((char *) cell_data_int_[i].values_->data(), n_bytes_cell_data_int[i]);
    }
    //--------------------------------------------------------------------------


    file << endl;
    file << tab1 << "</AppendedData>" << endl;

    file << "</VTKFile>";
}



template<int dim, int codim, class T>
void Writer<dim, codim, T>::
save(const string &filename, const string &format) const
{
    //--------------------------------------------------------------------------
    Assert(format == "ascii" || format == "appended",
           ExcMessage("Unsupported format."));
    //--------------------------------------------------------------------------

    std::vector< std::vector< std::array<T,3> > >
    points_in_iga_elements(n_iga_elements_, vector< array<T,3> >(n_points_per_iga_element_));

    std::vector< std::vector< std::array< int,n_vertices_per_vtk_element_> > >
    vtk_elements_connectivity(n_iga_elements_);
    for (auto &iga_elem_connectivity : vtk_elements_connectivity)
        iga_elem_connectivity.resize(n_vtk_elements_per_iga_element_);

    this->fill_points_and_connectivity(points_in_iga_elements, vtk_elements_connectivity);

    const string vtu_filename = filename + ".vtu";


    if (format == "ascii")
    {
        ofstream file(vtu_filename);
        file.setf(ios::scientific);
        file.precision(precision_);
        this->save_ascii(file, points_in_iga_elements, vtk_elements_connectivity);
    }
    else if (format == "appended")
    {
        this->save_appended(vtu_filename, points_in_iga_elements, vtk_elements_connectivity);
    }
}


template<int dim, int codim, class T>
void Writer<dim, codim, T>::print_info(LogStream &out) const
{

    std::vector< std::vector< std::array<T,3> > >
    points_in_iga_elements(n_iga_elements_, vector< array<T,3> >(n_points_per_iga_element_));

    std::vector< std::vector< std::array< int,n_vertices_per_vtk_element_> > >
    vtk_elements_connectivity(n_iga_elements_);
    for (auto &iga_elem_connectivity : vtk_elements_connectivity)
        iga_elem_connectivity.resize(n_vtk_elements_per_iga_element_);

    this->fill_points_and_connectivity(points_in_iga_elements, vtk_elements_connectivity);

    this->save_ascii(out, points_in_iga_elements, vtk_elements_connectivity);
}


IGA_NAMESPACE_CLOSE

#include <igatools/io/writer.inst>
