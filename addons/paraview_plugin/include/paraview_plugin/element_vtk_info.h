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

#ifndef ELEMENT_VTK_INFO_H_
#define ELEMENT_VTK_INFO_H_

#if PARAVIEW_PLUGIN

#include <igatools/base/config.h>
// #include <igatools/utils/tensor_size.h>
// #include <igatools/base/tensor.h>


// template<class T> class vtkSmartPointer;
// class vtkMultiBlockDataSet;
// class vtkPointData;
// class vtkPoints;
// class vtkStructuredGrid;
// class vtkUnstructuredGrid;
// class vtkCellArray;
// class vtkPointSet;


IGA_NAMESPACE_OPEN

template <int dim> class Quadrature;
template <int dim> class TensorSize;

template <int dim>
class ElementVTKInfo
{
private:
    /*
     * Self type.
     */
    typedef ElementVTKInfo<dim> Self_;

    /*
     * Self shared pointer type.
     */
    typedef std::shared_ptr<Self_t_> SelfPtr_;

    /*
     * Alias for a shared pointer of a Quadrature type.
     */
    typedef std::shared_ptr<const Quadrature<dim>> QuadConstPtr_;

    /*
     * Alias for the points map.
     */
    typedef SafeSTLVector<SafeSTLVector<Index>> PointsMap_;

    /*
     * Alias for the points mask.
     */
    typedef SafeSTLVector<Index> PointsMask_;

    /** @name Constructors */
    ///@{
    /**
     * Constructs an instance of the class for the given
     * @p num_vis_element by side and the @p grid_type.
     */
    ElementVTKInfo(const TensorSize<dim> &num_vis_elements,
                   const vtkGridType &grid_type);

    /** Default constructor. Not allowed to be used. */
    Self_() = delete;

    /** Copy constructor. Not allowed to be used. */
    Self_(const Self_ &) = delete;

    /** Move constructor. Not allowed to be used. */
    Self_(const Self_ &&) = delete;

    ///@}

    /** @name Assignment operators */
    ///@{
    /** Copy assignment. Not allowed to be used. */
    void operator=(const Self_ &) = delete;

    /** Move assignment. Not allowed to be used. */
    void operator=(const Self_ &&) = delete;
    ///@}

public:

    /**
     * @name Creators.
     */
    ///@{
    /**
     * Builds and returns an instance of the class, wrapped by a shared pointer,
     * for the given @p num_vis_element by side and the @p grid_type.
     */
    static SelfPtr_ create(const TensorSize &num_vis_elements,
                           const vtkGridType &grid_type);
    ///@}

    /** Retrieves the number of visualization points by element. */
    Size get_num_points() const;

    /** Retrieves the quadrature used for visualization. */
    QuadConstPtr_ get_quadrature() const;

    /** Retrieves the points map. */
    const PointsMap_ &get_points_map() const;

    /** Retrieves the points mask. */
    const PointsMask_ &get_points_mask() const;

private:

  /** Quadrature used for visualization. */
  const QuadConstPtr_ quadrature_;

  /** Points map
   * @todo: to define here what is this map.
   */
  const PointsMap_ points_map_;

  /** Points mask
   * @todo: to define here what is this mask.
   */
  const PointsMask_ points_mask_;


  /**
  * Creates and returns the visualization quadrature wrapped into a shared
  * pointer.
  */
  template <int dim>
  static QuadConstPtr_
  create_quadrature(const TensorSize<dim> &num_vis_elements,
                    const vtkGridType &grid_type);

  template <int dim>
  static PointsMap_
  create_points_map(const TensorSize<dim> &num_vis_elements,
                    const vtkGridType &grid_type);
  static PointsMask_
  create_points_mask(const TensorSize<dim> &num_vis_elements,
                     const vtkGridType &grid_type);

};

IGA_NAMESPACE_CLOSE

#endif // PARAVIEW_PLUGIN

#endif // ELEMENT_VTK_INFO_H_