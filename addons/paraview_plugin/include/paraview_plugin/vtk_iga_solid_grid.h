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

#ifndef __VTK_IGA_SOLID_GRID_H_
#define __VTK_IGA_SOLID_GRID_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/safe_stl_vector.h>
#include <igatools/utils/tensor_size.h>

#include <igatools/base/objects_container.h>
#include <igatools/geometry/domain.h>
#include <igatools/functions/function.h>

#include <boost/mpl/lambda.hpp>
#include <boost/mpl/remove_if.hpp>

class vtkPointSet;
class vtkPoints;
class vtkPointData;
class vtkCellArray;
template<class T> class vtkSmartPointer;


IGA_NAMESPACE_OPEN

class ObjectsContainer;
template <int dim> class Quadrature;
template <int dim, int codim> class Domain;
struct VtkGridInformation;


template <class Domain>
class VtkIgaSolidGrid
{
private:

  /**
   * Dimension.
   */
  static const int dim = Domain::dim;

  struct PointsTopology
  {
  private:
      /// Alias for mesh grid information shared pointer.
      typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

      /// Quadrature container shared pointer type.
      typedef std::shared_ptr<Quadrature<dim>>  QuadPtr_;

      /// Connectivity container.
      typedef SafeSTLVector<SafeSTLVector<Index>>  Connectivity_;

      /// Points map container.
      typedef SafeSTLVector <SafeSTLVector <Index>>  Map_;

      /// Self type;
      typedef PointsTopology Self_;

      PointsTopology() = delete;
      PointsTopology(const PointsTopology &) = delete;
      PointsTopology(const PointsTopology &&) = delete;
      void operator=(const PointsTopology &) = delete;
      void operator=(const PointsTopology &&) = delete;

  public:
      PointsTopology(const std::shared_ptr<const Grid<dim>> cartesian_grid,
                     const GridInfoPtr_ grid_info);

  private:
      /// Connectivity for the vtk cells of a single Bezier element.
      const Connectivity_ connectivity_;

      /// Visualization quadrature.
      const QuadPtr_ quad_;

      /// Number of vtk cells per direction in each Bezier element.
      const TensorSize <dim> n_vis_elements_;

      /// Points map.
      Map_ map_;

      /// Points mask.
      SafeSTLVector <Index> mask_;

      /// Total number of points in the visualization.
      Size n_total_points_;

  public:
      /**
       * TODO: DOCUMENT.
       */
      Size get_num_pts_per_bezier_elem() const;

      /**
       * TODO: DOCUMENT.
       */
      Size get_num_bezier_elems() const;

      /**
       * TODO: DOCUMENT.
       */
      Size get_flat_num_cells_per_bezier_elem() const;

      /**
       * TODO: DOCUMENT.
       */
      Size get_num_total_pts() const;

      /**
       * TODO: DOCUMENT.
       */
      Size get_num_pts_per_single_vtk_cell() const;

      /**
       * TODO: DOCUMENT.
       */
      const SafeSTLVector<Index> &get_mask() const;

      /**
       * TODO: DOCUMENT.
       */
      const Connectivity_ &get_connectivity() const;

      /**
       * TODO: DOCUMENT.
       */
      Size get_num_vtk_cells_per_bezier_elem (const Index &dir) const;

      /**
       * TODO: DOCUMENT.
       */
      Map_::const_iterator map_cbegin() const;

      /**
       * TODO: DOCUMENT.
       */
      QuadPtr_ get_quadrature() const;

  private:
      /**
       * Todo: to document
       */
      void fill_points_map_mask(const std::shared_ptr<const Grid<dim>> cartesian_grid,
                                const GridInfoPtr_ grid_info);

      /**
       * Creates the connectivity of the VTK cells for a single
       * Bezier element.
       */
      static Connectivity_ create_element_connectivity(const GridInfoPtr_ grid_info);

      /**
       * Creates the connectivity of the VTK linear cells for a single
       * Bezier element.
       */
      static Connectivity_ create_linear_element_connectivity
      (const GridInfoPtr_ grid_info);

      /**
       * Creates the connectivity of the VTK quadratic cells for a single
       * Bezier element. For the 1D case.
       */
      template <int aux_dim>
      static Connectivity_ create_quadratic_element_connectivity
      (const GridInfoPtr_ grid_info,
       typename std::enable_if_t<aux_dim == 1> * = 0);

      /**
       * Creates the connectivity of the VTK quadratic cells for a single
       * Bezier element. For the 2D case.
       */
      template <int aux_dim>
      static Connectivity_ create_quadratic_element_connectivity
      (const GridInfoPtr_ grid_info,
       typename std::enable_if_t<aux_dim == 2> * = 0);

      /**
       * Creates the connectivity of the VTK quadratic cells for a single
       * Bezier element. For the 3D case.
       */
      template <int aux_dim>
      static Connectivity_ create_quadratic_element_connectivity
      (const GridInfoPtr_ grid_info,
       typename std::enable_if_t<aux_dim == 3> * = 0);

      /**
       * Creates the quadrature needed for the visualization.
       */
      static QuadPtr_ create_visualization_quadrature(const GridInfoPtr_ grid_info);
  };

  /**
   * Space dimension.
   */
  static const int space_dim = Domain::space_dim;

  /**
   * Self type.
   */
  typedef VtkIgaSolidGrid Self_;

  /**
   * Self shared poitner type.
   */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /**
   * Alias for a shared pointer of a domain type.
   */
  typedef std::shared_ptr<const Domain> DomainPtr_;

  /**
   * Alias for mesh grid information shared pointer.
   */
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /**
   * Alias for vtk grid object for visualization.
   */
  typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

  /**
   * Quadrature container shared pointer type.
   */
  typedef std::shared_ptr<Quadrature<dim>>  QuadPtr_t_;

  /**
   * Functions container shared pointer type.
   */
  typedef typename std::shared_ptr<ObjectsContainer> ObjContPtr_t_;

  template <class T>
  struct IsInValidFunction :
          boost::mpl::or_<
          boost::mpl::bool_<(T::dim != Domain::dim)>,
          boost::mpl::bool_<(T::space_dim != Domain::space_dim)>>
          {};

  template <class T>
  struct IsInValidGridFunction :
          boost::mpl::or_<
          boost::mpl::bool_<(T::dim != Domain::dim)>,
          boost::mpl::bool_<(Domain::space_dim != Domain::dim)>>
          {};

  /**
   * Valid functions.
   */
  using ValidFuncs_ = typename boost::fusion::result_of::as_vector<
          typename boost::mpl::transform<
            typename boost::mpl::remove_if<
              InstantiatedTypes::Functions,
              typename boost::mpl::lambda< IsInValidFunction< boost::mpl::_1 > >::type>::type,
            std::shared_ptr<boost::mpl::_1>>::type>::type;

  /**
   * Valid grid functions.
   */
  using ValidGridFuncs_ = typename boost::fusion::result_of::as_vector<
          typename boost::mpl::transform<
            typename boost::mpl::remove_if<
              InstantiatedTypes::GridFunctions,
              typename boost::mpl::lambda< IsInValidGridFunction< boost::mpl::_1 > >::type>::type,
            std::shared_ptr<boost::mpl::_1>>::type>::type;

  /**
   * Constructor, copy and assignment operators not allowed to be used.
   */
  VtkIgaSolidGrid() = delete;
  VtkIgaSolidGrid(const VtkIgaSolidGrid &) = delete;
  VtkIgaSolidGrid(const VtkIgaSolidGrid &&) = delete;
  void operator=(const VtkIgaSolidGrid &) = delete;
  void operator=(const VtkIgaSolidGrid &&) = delete;


public:

  /**
   * Creates and returns the vtk grid for the visualization.
   */
  static VtkGridPtr_ create(const DomainPtr_ domain,
                            const ObjContPtr_t_ obj_container,
                            const GridInfoPtr_ grid_info,
                            const bool is_physical);

private:

  /**
   * Creates and returns the structured vtk grid for the visualization.
   */
  static VtkGridPtr_ create_grid_vts(const DomainPtr_ domain,
                                     const ObjContPtr_t_ objs_container,
                                     const GridInfoPtr_ grid_info,
                                     const bool is_physical);

  /**
   * Creates and returns the unstructured vtk grid for the visualization
   */
  static VtkGridPtr_ create_grid_vtu(const DomainPtr_ domain,
                                     const ObjContPtr_t_ objs_container,
                                     const GridInfoPtr_ grid_info,
                                     const bool is_physical);

  /**
   * Creates the points for the visualization grids.
   * (Both, structured and unstructured).
   */
  static vtkSmartPointer<vtkPoints> create_points(const DomainPtr_ domain,
                                                  const PointsTopology &points_top);

  /**
   * Creates the point data associated to the mapping.
   * TODO: to document.
   */
  static void create_point_data_physical(const DomainPtr_ domain,
                                         const ObjContPtr_t_ objs_container,
                                         const PointsTopology &points_top,
                                         const VtkGridPtr_ vtk_grid);

  /**
   * Creates the point data associated to the mapping.
   * TODO: to document.
   */
  static void create_point_data_parametric(const DomainPtr_ domain,
                                           const ObjContPtr_t_ objs_container,
                                           const PointsTopology &points_top,
                                           const VtkGridPtr_ vtk_grid);

};

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_SOLID_GRID_H_
