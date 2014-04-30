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
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/geometry/mapping.h>
#include <igatools/base/quadrature.h>

#include <boost/variant.hpp>

IGA_NAMESPACE_OPEN

template < class RefSpace, class PushForward >
class PhysicalSpace ;

template< int dim_domain, int range, int rank >
class NURBSSpace ;




//todo:add the add_celldata function

/**
 * @todo This class needs to be documented, and explained with an example.
 */
template< int dim_ref, int dim_phys=dim_ref, class T = double >
class Writer
{
public:

    static const int codim = dim_phys-dim_ref;

    Writer(const std::shared_ptr< CartesianGrid< dim_ref > > grid);

    Writer(const std::shared_ptr< CartesianGrid< dim_ref > > grid,
           const Index num_points_direction);

    Writer(const std::shared_ptr< const Mapping< dim_ref, codim > > mapping,
           const Index num_points_direction);


    //TODO: can we remove this constructor?
    /**
     * This constructor builds a Writer object with using a uniform distribution for the evaluation points.
     * \param[in] grid Grid.
     * \param[in] mapping Mapping used to transform the field(s) from the reference domain to the physical domain.
     * \param[in] num_points_direction Number of plotting points along the coordinate directions.
     * \param[in] eps This value is the amount of the displacement of the first and last point.
     * \param[in] map_is_double_precision If true the map values will be written using double precision data.
     * Default value is TRUE.
     * (in each coordinate direction) of the quad rule on the unit element.
     * This means that the element for the quad rule it will not [0,1]^dim but [eps,1.0-eps]^dim.
     * Valid values for eps are the ones in the interval [0.0,0.5).
     * \note Any field that will be added to the writer must refer to the same
     * CartesianGrid used here, otherwise an exception will be raised.
     * \note The number of points in each coordinate direction must be greater or equal than 2,
     * otherwise an exception will be raised.
     * \see add_field
     */
//    Writer(const std::shared_ptr< const CartesianGrid< dim_ref > > grid,
//           const std::shared_ptr< const Mapping< dim_ref, dim_phys > > mapping,
//           const std::array< Index, dim_ref > &num_points_direction) {}


    /**
     * This constructor builds a Writer object with using a distribution for the evaluation points
     * given by the @p quadrature scheme.
     * \note Any field that will be added to the writer must refer to the same
     * CartesianGrid used here, otherwise an exception will be raised.
     * \note The number of points in each coordinate direction must be greater or equal than 2,
     * otherwise an exception will be raised.
     * \see add_field
     */
    Writer(const std::shared_ptr< const Mapping< dim_ref, codim > > mapping,
           const std::shared_ptr< const Quadrature<dim_ref> > quadrature) ;


    /**
     * Default constructor. Not allowed to be used.
     */
    Writer() = delete ;


    /**
     * Copy constructor. Not allowed to be used.
     */
    Writer(const Writer< dim_ref, codim > &writer) = delete ;


    /**
     * Assignment operator. Not allowed to be used.
     */
    Writer< dim_ref, dim_phys > &
    operator=(const Writer< dim_ref, dim_phys > &writer) = delete ;


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
    template<class Space>
    void add_field(
        std::shared_ptr<Space> space,
        const Vector<LinearAlgebraPackage::trilinos> &coefs,
        const std::string &name) ;


    void add_element_data(
        const std::vector<double> &element_data,
        const std::string &name) ;

    void add_element_data(
        const std::vector<int> &element_data,
        const std::string &name) ;


    /**
     * Save the data on a .vtu file.
     * \param[in] filename - Output file name.
     * \param[in] format - Output format. It can be "ascii" or "appended".
     * \note The .vtu extension should NOT part of the file name.
     */
    void save(
        const std::string &filename,
        const std::string &format = "ascii") ;


    /**
     * Get the points (in the 3D physical domain) associated to each iga element.
     * The first index refers to the iga element,
     * the second index refers to the point within the iga element specified by the first index.
     * \code
     * auto points = writer.get_points_in_elements() ;
     * //points[i][j] is the j-th point in the i-th iga element
     * \endcode
     */
    const std::vector< std::vector< std::array<T,3> > > &get_points_in_iga_elements() const ;

    /**
     * Returns the number of IGA elements handled by the Writer.
     */
    int get_num_iga_elements() const ;

    /**
     * Returns the number of VTK elements handled by the Writer.
     */
    int get_num_vtk_elements() const ;

    /**
     * Returns the number of VTK elements used for each IGA element.
     */
    int get_num_vtk_elements_per_iga_element() const ;

    /**
     * Returns the number of evaluation points used for each IGA element.
     */
    int get_num_points_per_iga_element() const ;


    void add_point_data(const int n_iga_elements,
                        const int n_points_per_iga_element,
                        const int n_values_per_point,
                        const std::string &type,
                        const std::vector<std::vector<std::vector<T>>> &data_iga_elements,
                        const std::string &name) ;



private:
    std::string byte_order_ ;

    static const int n_vertices_per_vtk_element_ = UnitElement< dim_ref >::vertices_per_element ;

    const std::string filename_ ;

    std::shared_ptr< const CartesianGrid< dim_ref > > grid_ ;


    std::shared_ptr< const Mapping< dim_ref, codim > > map_ ;


    /**
     * Unit element quadrature rule used for the plot.
     */
    Quadrature< dim_ref > quad_plot_ ;


    TensorSize<dim_ref> num_points_direction_ ;


    TensorSize<dim_ref> num_subelements_direction_ ;

    /**
     * Number of VTK elements contained in each IGA element.
     */
    Size n_vtk_elements_per_iga_element_ ;

    /**
     * Number of VTK elements handled by the Writer.
     */
    Size n_vtk_elements_ ;


    /**
     * Number of IGA elements handled by the Writer.
     */
    Size n_iga_elements_ ;

    /**
     * Number of evaluation points in each IGA element.
     */
    Size n_points_per_iga_element_ ;

    /**
     * Number of VTK elements handled by the Writer.
     */
    Size n_vtk_points_ ;

    /**
     * Coordinates of the evaluation points associated to each iga element.
     * \note The point coordinates are referred to the physical domain extended to 3D space.
     */
    std::vector< std::vector< std::array<T,3> > > points_in_iga_elements_ ;

    unsigned char vtk_element_type_ ;

    /**
     * Connectivity of the vtk elements inside each iga element.
     */
    std::vector< std::vector< std::array< int,n_vertices_per_vtk_element_> > > vtk_elements_connectivity_ ;


    /**
     * This function take as input an @p iga_element_id and a set of points in the [0,1]^dim domain, and maps
     * those points to the reference domain and to the physical domain defined by the mapping used in the costructor.
     * Moreover it returns the connectivity of the points on the element.
     * @param[in] iga_element_id Element ID.
     * @param[in] elem_quad Evaluation points in the [0,1]^dim domain.
     * @param[out] element_connectivity Connectivity of the points defined on the element.
     * @param[out] points_phys_iga_element Coordinate of the points in the physical domain.
     * \note Due to the fact that VTK needs always points in 3D, when we have the dimension of the physical space less than 3,
     * we set the coordinate of the "missing dimension" to 0. In other words, if the dimension of the physical domain is 2, then the
     * points are located on the plane z=0.
     * If the dimension of the physical domain is 1, then the points are located on the line with y=0 and z=0.
     */
    void get_subelements(
        const typename Mapping< dim_ref, codim>::ElementIterator elem,
        std::vector< std::array< int, n_vertices_per_vtk_element_ > > &vtk_elements_connectivity,
        std::vector< std::array<T,3> > &points_phys_iga_element) const ;


private:
    struct PointData
    {

        PointData(
            const std::string &name,
            const std::string &type,
            const Size num_elements,
            const Size num_points_per_element,
            const Size num_components,
            std::shared_ptr< std::vector<T> > values)
            :
            name_(name),
            type_(type),
            num_elements_(num_elements),
            num_points_per_element_(num_points_per_element),
            num_components_(num_components),
            values_(values)
        {} ;

        const std::string name_ ;

        const std::string type_ ;

        Size num_elements_ ;

        Size num_points_per_element_ ;

        Size num_components_ ;


        std::shared_ptr< std::vector<T> > values_ ;
    } ;

    std::vector< PointData > fields_ ;

    std::vector<std::string> names_point_data_scalar_ ;
    std::vector<std::string> names_point_data_vector_ ;
    std::vector<std::string> names_point_data_tensor_ ;

    template<class data_type>
    struct CellData
    {
        CellData(
            const std::vector<data_type> &values,
            const std::string &name)
            :
            values_(new std::vector<data_type>(values)),
            name_(name),
            type_("scalar"),
            num_components_(1)
        {} ;

        std::shared_ptr< std::vector<data_type> > values_ ;

        const std::string name_ ;

        const std::string type_ ;

        iga::Size num_components_ ;
    } ;

    std::vector< CellData<double> > cell_data_double_ ;

    std::vector< CellData<int> > cell_data_int_ ;

    std::vector<std::string> names_cell_data_scalar_ ;
    std::vector<std::string> names_cell_data_vector_ ;
    std::vector<std::string> names_cell_data_tensor_ ;

    const int sizeof_Real_ = 0;
    std::string string_Real_ ;

    const int sizeof_int_  = 0;
    std::string string_int_ ;

    const int sizeof_uchar_  = 0;
    std::string string_uchar_ ;





    std::stringstream appended_data_ ;


    int offset_ ;


    int precision_ ;


    void save_ascii(const std::string &filename) const ;

    void save_appended(const std::string &filename) const ;

    void fill_points_and_connectivity() ;
} ;


IGA_NAMESPACE_CLOSE

#endif /* __WRITER_H_ */
