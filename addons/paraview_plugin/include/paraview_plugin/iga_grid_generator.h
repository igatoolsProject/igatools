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

#ifndef IGA_GRID_GENERATOR_H_
#define IGA_GRID_GENERATOR_H_

#include <igatools/base/config.h>

#include <igatools/contrib/variant.h>
#include <boost/property_tree/ptree.hpp>
#include <memory>

class vtkInformation;
template <class T> class vtkSmartPointer;
class vtkPoints;
class vtkCellArray;
class vtkIdTypeArray;

namespace iga
{
template <int dim, int codim, int range, int rank> class Function;
template <int dim> class QuadratureTensorProduct;
template <class T> class vector;
template <class T, int dim> class special_array;
template <int dim> class TensorSize;
}

using FunctionPtrVariant = Variant<
std::shared_ptr<iga::Function<1,0,1,1>>,
std::shared_ptr<iga::Function<1,0,2,1>>,
std::shared_ptr<iga::Function<1,0,3,1>>,
std::shared_ptr<iga::Function<2,0,1,1>>,
std::shared_ptr<iga::Function<2,0,2,1>>,
std::shared_ptr<iga::Function<2,0,3,1>>,
std::shared_ptr<iga::Function<3,0,1,1>>,
std::shared_ptr<iga::Function<3,0,3,1>>,
std::shared_ptr<iga::Function<0,0,1,1>>,
std::shared_ptr<iga::Function<0,0,2,1>>,
std::shared_ptr<iga::Function<0,0,3,1>>>;

using QuadraturePtrVariant = Variant<
std::shared_ptr<iga::QuadratureTensorProduct<1>>,
std::shared_ptr<iga::QuadratureTensorProduct<2>>,
std::shared_ptr<iga::QuadratureTensorProduct<3>>>;




class IGAVTKGridGenerator
{
private:
  /*
   * Self type.
   */
  typedef IGAVTKGridGenerator Self_t_;

  /*
   * Self shared pointer type.
   */
  typedef std::shared_ptr<Self_t_> SelfPtr_t_;

  /*
   * Container type for the connectivity of a block of cells contained
   * in a single Bezier element.
   */
  template <int dim>
  using Connectivity_t_ = iga::vector<iga::special_array<iga::Index, iga::constexpr_pow(2, dim)>>;

  /*
   * Constructor, copy and assignement opertors not allowed to be used.
   */
  IGAVTKGridGenerator () = delete;
  IGAVTKGridGenerator (const IGAVTKGridGenerator &) = delete;
  IGAVTKGridGenerator (const IGAVTKGridGenerator &&) = delete;
  void operator=(const IGAVTKGridGenerator &) = delete;
  void operator=(const IGAVTKGridGenerator &&) = delete;

  /*
   * Constructor based on the number of points per direction.
   */
  IGAVTKGridGenerator (const std::string& file_name,
                       const std::string& file_path,
                       const int* num_visualization_points);

public:
  /*
   * Creates a new object and returns it wrapped in a shared pointer.
   */
  static SelfPtr_t_ create (const std::string& file_name,
                            const std::string& file_path,
                            const int* num_visualization_points);

  /*
   * Fill the solid unstructured vtk grid output.
   */
  int fill_solid_output (vtkInformation *outInfo) const;

  /*
   * Fill the indentity map unstructured vtk grid output.
   */
  int fill_identity_output (vtkInformation *outInfo) const;

  /*
   * Fill the control mesh unstructured vtk grid output.
   */
  int fill_control_mesh_output (vtkInformation *outInfo) const;


private:

  /*
   * File name.
   */
  const std::string file_name_;

  /*
   * File path.
   */
  const std::string file_path_;

  /*
   * Dimension of the IgMapping;
   */
  int dim_;

  /*
   * Codimension of the IgMapping;
   */
  int codim_;

  /*
   * Variant for containing the function read from the xml file.
   */
  FunctionPtrVariant function_variant_;

  /*
   * Variant for containing the quadratur for the visualization.
   */
  QuadraturePtrVariant quadrature_variant_;

  /*
   * Retrieves the xml tree.
   */
  boost::property_tree::ptree parse_xml_file () const;

  /*
   * Retrieves the dimensions of the function contained in the xml tree.
   */
  std::pair<int, int>
    get_dimensions_from_xml (const boost::property_tree::ptree& tree) const;

  /*
   * Builds the function from the xml file and stores it in function_variant_
   */
  FunctionPtrVariant
    create_function_from_xml (const boost::property_tree::ptree& xml_tree);

  /*
   * Creates the quadrature given the points in every direction.
   */
  QuadraturePtrVariant create_quadrature (const iga::vector<iga::Size>& n_points) const;

  /**
   * Visitor used for filling points and cells of the solid grid.
   */
  struct FillSolidGridVisitor
  {
  public:
    typedef bool result_type;

    /*
     * This function fill the points and arrays of the solid grid.
     */
    template <class FunctionPtr_t_>
    result_type operator () (const FunctionPtr_t_ func_ptr,
                             const vtkSmartPointer<vtkPoints> points,
                             const vtkSmartPointer<vtkCellArray> cellsArray,
                             const QuadraturePtrVariant* const quad_var_ptr);

  };

  /**
   * Visitor used for filling points and cells of the identity mapping grid.
   */
  struct FillIdentityGridVisitor
  {
  public:
    typedef bool result_type;

    /*
     * This function fill the points and arrays of the solid grid.
     */
    template <class FunctionPtr_t_>
    result_type operator () (const FunctionPtr_t_ func_ptr,
                             const vtkSmartPointer<vtkPoints> points,
                             const vtkSmartPointer<vtkCellArray> cellsArray,
                             const QuadraturePtrVariant* const quad_var_ptr);

  };

  /**
   * Visitor used for filling points and cells of the control mesh grid.
   */
  struct FillControlMeshGridVisitor
  {
  public:
    typedef bool result_type;

    /*
     * This function fill the points and arrays of the control mesh grid.
     */
    template <class FunctionPtr_t_>
    result_type operator () (const FunctionPtr_t_ func_ptr,
                             const vtkSmartPointer<vtkPoints> points,
                             const vtkSmartPointer<vtkCellArray> cellsArray);

  };

  /*
   * Create the cell ids container needed for defining vtk cells.
   */
  template <int dim>
  static vtkSmartPointer<vtkIdTypeArray>
  create_cell_ids (const iga::TensorSize<dim>& n_points_per_direction,
                   const iga::Size& n_bezier_elements);

  /*
   * Creates a container for creating the connectivity of the patch.
   */
  template <int dim>
  static Connectivity_t_<dim>
  create_connectivity_base (const iga::TensorSize<dim>& n_points_per_direction);

};

#endif // IGA_GRID_GENERATOR_H_