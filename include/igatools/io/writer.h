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

#ifndef __WRITER_H_
#define __WRITER_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/geometry/mapping.h>
#include <igatools/base/quadrature.h>

#include <boost/variant.hpp>

IGA_NAMESPACE_OPEN

template < class RefSpace, class PushForward >
class PhysicalSpace;

template< int dim_domain, int range, int rank >
class NURBSSpace;


//todo:add the add_celldata function

/**
 * This class is use to generate and save output in graphical
 * format.
 *
 * @todo This class needs to be documented, and explained with an example.
 */
template<int dim, int codim = 0, class T = double>
class Writer
{
private:
    using self_t = Writer<dim, codim, T>;
    using Grid = CartesianGrid<dim>;
    using Map  = Mapping<dim, codim>;

public:
    static const int space_dim = dim + codim;
    // TODO (pauletti, Jun 18, 2014): document this constructors
    //TODO(pauletti, Jun 21, 2014): should be a const Grid
    Writer(const std::shared_ptr<Grid> grid);

    Writer(const std::shared_ptr<Grid> grid,
           const Index num_points_direction);

    Writer(const std::shared_ptr<const Map> mapping,
           const Index num_points_direction);

    /**
     * This constructor builds a Writer object using a distribution for
     * the evaluation points given by the @p quadrature scheme.
     * \note Any field that will be added to the writer must refer to the same
     * CartesianGrid used here, otherwise an exception will be raised.
     * \note The number of points in each coordinate direction must be greater or equal than 2,
     * otherwise an exception will be raised.
     * \see add_field
     */
    Writer(const std::shared_ptr<const Map> mapping,
           const std::shared_ptr<const Quadrature<dim>> quadrature);


    /**
     * Default constructor. Not allowed to be used.
     */
    Writer() = delete;

    /**
     * Copy constructor. Not allowed to be used.
     */
    Writer(const self_t &writer) = delete;

    /**
     * Assignment operator. Not allowed to be used.
     */
    self_t &operator=(const self_t &writer) = delete;

    /**
     * \brief Add a field to the output file.
     * \param[in] space Space (i.e. basis function) that will be used for the field evaluation.
     * \param[in] coefs Coefficients of the field as linear combination of the set of basis functions specified by space.
     * \param[in] transformation_flag This flag specifies which kind of transformation will be used to push-forward the field.
     * It can be only one choice between {h_grad, h_div, h_curl}. For their meaning see the documentation in Mapping.
     * \param[in] name Name of the field.
     * \param[in] format - If 'ascii' the field will be written using ASCII characters,
     * if 'appended' the field will be written in binary format using.
     * See the VTK file format documentation for details. Default is 'appended'.
     * \note The number of coefficients must be equal to the number of basis functions specified
     * by the space, otherwise an exception will be raised.
     */
    template<class Space, LAPack la_pack = LAPack::trilinos>
    void add_field(std::shared_ptr<Space> space,
                   const Vector<la_pack> &coefs,
                   const std::string &name);


    void add_element_data(const vector<double> &element_data,
                          const std::string &name);

    void add_element_data(const vector<int> &element_data,
                          const std::string &name);


    /**
     * Save the data on a .vtu file.
     * \param[in] filename - Output file name.
     * \param[in] format - Output format. It can be "ascii" or "appended".
     * \note The .vtu extension should NOT part of the file name.
     */
    void save(const std::string &filename,
              const std::string &format = "ascii") const;

    /**
     * Writes the vtu into a LogStream, filtering it for uniform
     * output across different systems.
     * @note this function is only for testing purposes
     */
    void print_info(LogStream &out) const;

    /**
     * Returns the number of IGA elements handled by the Writer.
     */
    int get_num_iga_elements() const;

    /**
     * Returns the number of VTK elements handled by the Writer.
     */
    int get_num_vtk_elements() const;

    /**
     * Returns the number of VTK elements used for each IGA element.
     */
    int get_num_vtk_elements_per_iga_element() const;

    /**
     * Returns the number of evaluation points used for each IGA element.
     */
    int get_num_points_per_iga_element() const;


    /**
     * \brief Add data for every evaluation point to the output file.
     * \param[in] n_values_per_point
     * \param[in] type Type of add to be added. It can be scalar, vector
     * or tensor.
     * \param[in] name Name of the field.
     * \param[in] data_iga_elements Data to be added. The different levels of
     * the container are: the first vector level corresponds to the IGA
     * elements; the second one to the evaluation points inside an IGA element;
     * and the third one to the components of the data to be plotted.
     * \note The number of entries of @p data_iga_elements must be equal to the
     * number of IGA elements, otherwise an exception will be raised.
     * \note The number of entries of eachy entry of @p data_iga_elements must
     * be equal to the number of IGA elements, otherwise an exception will be
     * raised.
     * \note The number of values associated to every plot points that are
     * specified in @p data_iga_elements must be equal to @ n_values_per_point,
     * otherwise an exception will be raised.
     */
    void add_point_data(const int n_values_per_point,
                        const std::string &type,
                        const vector<vector<vector<T>>> &data_iga_elements,
                        const std::string &name);



private:
    std::string byte_order_;

    static const int n_vertices_per_vtk_element_ = UnitElement< dim >::vertices_per_element;

    const std::string filename_;

    std::shared_ptr<const Grid> grid_;

    std::shared_ptr<const Map> map_;

    /**
     * Unit element quadrature rule used for the plot.
     */
    Quadrature< dim > quad_plot_;


    TensorSize<dim> num_points_direction_;


    TensorSize<dim> num_subelements_direction_;

    /**
     * Number of VTK elements contained in each IGA element.
     */
    Size n_vtk_elements_per_iga_element_;

    /**
     * Number of VTK elements handled by the Writer.
     */
    Size n_vtk_elements_;


    /**
     * Number of IGA elements handled by the Writer.
     */
    Size n_iga_elements_;

    /**
     * Number of evaluation points in each IGA element.
     */
    Size n_points_per_iga_element_;

    /**
     * Number of VTK elements handled by the Writer.
     */
    Size n_vtk_points_;

    unsigned char vtk_element_type_;

//TODO(pauletti, Jul 8, 2014): this documentation is incorrect
    /**
     * This function take as input an @p iga_element_id and a set of points in
     * the [0,1]^dim domain, and maps
     * those points to the reference domain and to the physical domain defined
     * by the mapping used in the costructor.
     * Moreover it returns the connectivity of the points on the element.
     * @param[in] iga_element_id Element ID.
     * @param[in] elem_quad Evaluation points in the [0,1]^dim domain.
     * @param[out] element_connectivity Connectivity of the points defined on
     * the element.
     * @param[out] points_phys_iga_element Coordinate of the points in the
     * physical domain.
     * \note Due to the fact that VTK needs always points in 3D, when we have
     *  the dimension of the physical space less than 3,
     * we set the coordinate of the "missing dimension" to 0. In other words,
     * if the dimension of the physical domain is 2, then the
     * points are located on the plane z=0.
     * If the dimension of the physical domain is 1, then the points are
     * located on the line with y=0 and z=0.
     */
    void get_subelements(
        const typename Mapping< dim, codim>::ElementIterator elem,
        vector< std::array< int, n_vertices_per_vtk_element_ > > &vtk_elements_connectivity,
        vector< std::array<T,3> > &points_phys_iga_element) const;


private:
    struct PointData
    {

        PointData(
            const std::string &name,
            const std::string &type,
            const Size num_elements,
            const Size num_points_per_element,
            const Size num_components,
            std::shared_ptr< vector<T> > values)
            :
            name_(name),
            type_(type),
            num_elements_(num_elements),
            num_points_per_element_(num_points_per_element),
            num_components_(num_components),
            values_(values)
        {};

        const std::string name_;

        const std::string type_;

        Size num_elements_;

        Size num_points_per_element_;

        Size num_components_;


        std::shared_ptr< vector<T> > values_;
    };

    vector< PointData > fields_;

    vector<std::string> names_point_data_scalar_;
    vector<std::string> names_point_data_vector_;
    vector<std::string> names_point_data_tensor_;

    template<class data_type>
    struct CellData
    {
        CellData(
            const vector<data_type> &values,
            const std::string &name)
            :
            values_(new vector<data_type>(values)),
            name_(name),
            type_("scalar"),
            num_components_(1)
        {};

        std::shared_ptr< vector<data_type> > values_;

        const std::string name_;

        const std::string type_;

        iga::Size num_components_;
    };

    vector< CellData<double> > cell_data_double_;

    vector< CellData<int> > cell_data_int_;

    vector<std::string> names_cell_data_scalar_;
    vector<std::string> names_cell_data_vector_;
    vector<std::string> names_cell_data_tensor_;

    const int sizeof_Real_ = 0;
    std::string string_Real_;

    const int sizeof_int_  = 0;
    std::string string_int_;

    const int sizeof_uchar_  = 0;
    std::string string_uchar_;

    std::stringstream appended_data_;

    int offset_;

    int precision_;

    template<class Out>
    void save_ascii(Out &file,
                    const vector< vector< std::array<T,3> > > &points_in_iga_elements,
                    const vector< vector< std::array< int,n_vertices_per_vtk_element_> > >
                    &vtk_elements_connectivity) const;

    void save_appended(const std::string &filename,
                       const vector< vector< std::array<T,3> > > &points_in_iga_elements,
                       const vector< vector< std::array< int,n_vertices_per_vtk_element_> > >
                       &vtk_elements_connectivity) const;

    void fill_points_and_connectivity(
        vector< vector< std::array<T,3> > > &points_in_iga_elements,
        vector< vector< std::array< int,n_vertices_per_vtk_element_> > >
        &vtk_elements_connectivity) const;
};


IGA_NAMESPACE_CLOSE

#endif /* __WRITER_H_ */
