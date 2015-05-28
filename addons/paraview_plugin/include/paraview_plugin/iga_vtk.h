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

#ifndef IGA_VTK_H_
#define IGA_VTK_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_size.h>

#include <igatools/base/tensor.h>
#include <igatools/contrib/variant.h>
#include <memory>

namespace iga
{
  template<class T> class SafeSTLVector;
  template<class T, int N> class SafeSTLArray;
  template<int dim> class Quadrature;
  template<int dim> class CartesianGrid;
  class FunctionsContainer;
  template<int dim, int codim, int range, int rank>  class Function;
};

template<class T> class vtkSmartPointer;
class vtkStructuredGrid;
class vtkMultiBlockDataSet;
class vtkIdTypeArray;
class vtkPointData;


class IGAVTK
{
private:
  /*
   * Self type.
   */
  typedef IGAVTK Self_t_;

  /*
   * Self shared pointer type.
   */
  typedef std::shared_ptr<Self_t_> SelfPtr_t_;

  /*
   * Container type for the connectivity of a block of cells contained
   * in a single Bezier element.
   */
  template <int dim>
  using Connectivity_t_ = 
    iga::SafeSTLVector<iga::SafeSTLArray<iga::Index, iga::constexpr_pow(2, dim)>>;

  /*
   * Alias for a shared pointer of a Quadrature type.
   */
  template <int dim>
  using QuadPtr_ = std::shared_ptr<iga::Quadrature<dim>>;

  /*
   * Constructor, copy and assignement opertors not allowed to be used.
   */
  IGAVTK (const IGAVTK &) = delete;
  IGAVTK (const IGAVTK &&) = delete;
  void operator=(const IGAVTK &) = delete;
  void operator=(const IGAVTK &&) = delete;

public:
  /*
   * Default constructor.
   */
  IGAVTK ();

  /*
   * Set the number of visualization points.
   */
  void set_number_visualization_points (const int* const num_visualization_points);

  /*
   * Set the file name and path.
   */
  void set_file (const std::string& file_name, const std::string& file_path);

  /*
   * Clears the class, destroying the read information.
   */
  void clear ();

  /*
   * Parses the input file.
   */
  void parse_file ();

  /*
   * Create temporary geometries.
   */
  template <int dim>
  void create_geometries ();

  /*
   * Generates the VTK grids.
   */
  void generate_vtk_grids(const int& grid_type,
                          const bool& create_control_mesh,
                          const bool& create_parametric_mesh,
                          const bool& create_physical_mesh,
                          vtkMultiBlockDataSet* const mb) const;

private:
  /*
   * File name.
   */
  std::string file_name_;

  /*
   * File path.
   */
  std::string file_path_;

  /*
   * Number of visualization points per direction.
   */
  iga::TensorSize<3> num_visualization_points_;

  /*
   * Container for the mapping and field functions.
   */
  std::shared_ptr<iga::FunctionsContainer> funcs_container_;

  /*
   * Generates the control mesh vtk grids.
   */
  template <int dim, int codim>
  void generate_control_mesh_grids (vtkMultiBlockDataSet* const mb,
                                    unsigned int& id) const;

  /*
   * Generates the parametric vtk grids.
   */
  template <int dim, int codim>
  void generate_parametric_mesh_grids (vtkMultiBlockDataSet* const mb,
                                       unsigned int& id,
                                       const bool unstructured) const;

  /*
   * Generates the physical vtk grids.
   */
  template <int dim, int codim>
  void generate_solid_mesh_grids (vtkMultiBlockDataSet* const mb,
                                  unsigned int& id,
                                  const bool unstructured,
                                  const bool is_parametric) const;

  /*
   * Returns the namesof identity and mapped functions from the function
   * container.
   */
  iga::SafeSTLArray<iga::SafeSTLVector<std::string>, 3> get_map_names () const;

  /*
   * Returns true if the passed mapping is and identity mapping. False elsewhere.
   */
  template <int dim, int codim>
  bool is_identity_mapping (std::shared_ptr<iga::Function<dim, 0, dim+codim, 1>> map) const;

  /*
   * Create the cell ids container needed for defining vtk cells.
   */
  template <int dim>
  static vtkSmartPointer<vtkIdTypeArray>
  create_vtu_cell_ids (const iga::TensorSize<dim>& n_points_per_direction,
                       const iga::Size& n_bezier_elements);

  /*
   * Creates a map between the number of Bezier element and the number o point
   * inside the element, and the global number of the point in the vtk grid.
   */
  template <int dim>
  iga::SafeSTLVector<iga::SafeSTLVector<iga::Index>>
  create_points_numbering_map (const std::shared_ptr<const iga::CartesianGrid<dim>> grid,
                               const iga::TensorSize<dim>& n_points,
                               const bool is_unstructured) const;

  template <int dim, int codim>
  void
  create_point_data (const std::shared_ptr<iga::Function<dim, 0, dim + codim, 1>> mapping,
                     const iga::Quadrature<dim> &quad,
                     const iga::SafeSTLVector<iga::SafeSTLVector<iga::Index>>& points_map,
                     vtkPointData* const data) const;

  template <int dim, int codim, int range, int rank>
  void
  create_point_data (const std::shared_ptr<iga::Function<dim, 0, dim + codim, 1>> mapping,
                     const iga::Quadrature<dim> &quad,
                     const iga::SafeSTLVector<iga::SafeSTLVector<iga::Index>>& point_num_map,
                     vtkPointData* const data) const;


  void
  inline
  tensor_to_tuple(const iga::Tdouble t, iga::Real* const tuple, int& pos) const
  {
    *(tuple + pos++) = t;
  }


  template <class Tensor>
  inline
  void
  tensor_to_tuple (const Tensor& t, iga::Real* const tuple, int& pos) const
  {
    for (int i = 0; i < Tensor::size; ++i)
      tensor_to_tuple (t[i], tuple, pos);
  };

  template <class Tensor>
  inline
  void
  tensor_to_tuple (const Tensor& t, iga::Real* const tuple) const
  {
    int pos = 0;
    tensor_to_tuple (t, tuple, pos);
  };
};

#endif // IGA_VTK_H_