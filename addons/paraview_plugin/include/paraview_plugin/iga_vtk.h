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


template<class T> class vtkSmartPointer;
class vtkMultiBlockDataSet;
class vtkPointData;
class vtkPoints;
class vtkStructuredGrid;
class vtkUnstructuredGrid;
class vtkCellArray;
class vtkPointSet;


IGA_NAMESPACE_OPEN

template <int dim> class Quadrature;
template <int dim> class CartesianGrid;
template <int dim> class TensorSize;
template <int dim, int codim, int range, int rank> class Function;
class FunctionsContainer;

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
     * Alias for a shared pointer of a Quadrature type.
     */
    template <int dim>
    using QuadPtr_ = std::shared_ptr<Quadrature<dim>>;

    /*
     * Alias for a shared pointer of a map function type.
     */
    template <int dim, int codim>
    using MapFunPtr_ = std::shared_ptr<Function<dim, 0, dim+codim, 1>>;

    /*
     * Constructor, copy and assignement opertors not allowed to be used.
     */
    IGAVTK(const IGAVTK &) = delete;
    IGAVTK(const IGAVTK &&) = delete;
    void operator=(const IGAVTK &) = delete;
    void operator=(const IGAVTK &&) = delete;

public:
    /*
     * Default constructor.
     */
    IGAVTK();

    /*
     * Set the number of visualization elements and the flag for quadratic
     * cells.
     */
    void set_visualization_element_properties(const int *const num_visualization_elements_physical,
                                              const int &grid_type_physical,
                                              const int *const num_visualization_elements_parametric,
                                              const int &grid_type_parametric);

    /*
     * Set the file name and path.
     */
    void set_file(const std::string &file_name, const std::string &file_path);

    /*
     * Clears the class, destroying the read information.
     */
    void clear();

    /*
     * Parses the input file.
     */
    void parse_file();

    /*
     * Create temporary geometries.
     */
    template <int dim>
    void create_geometries();

    /*
     * Generates the VTK grids.
     */
    void generate_vtk_grids(const bool &create_physical_mesh,
                            const bool &create_solid_physical_mesh,
                            const bool &create_control_physical_mesh,
                            const bool &create_knot_physical_mesh,
                            const bool &create_parametric_mesh,
                            const bool &create_solid_parametric_mesh,
                            const bool &create_knot_parametric_mesh,
                            vtkMultiBlockDataSet *const mb) const;

    /*
     * Fill VTK grids.
     */
    void fill_vtk_grids(vtkMultiBlockDataSet *const mb,
                        const Size &n_solid_mesh,
                        const Size &n_control_mesh,
                        const bool is_identity,
                        const bool create_solid_mesh,
                        const bool create_control_mesh,
                        const bool create_knot_mesh) const;

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
     * Number of visualization elements per direction.
     */
    TensorSize<3> num_visualization_elements_physical_;

    /*
     * Number of visualization elements per direction.
     */
    TensorSize<3> num_visualization_elements_parametric_;

    /*
     * Flag for the use of quadratic elements.
     */
    bool quadratic_cells_physical_;

    /*
     * Flag for the use of quadratic elements.
     */
    bool quadratic_cells_parametric_;

    /*
     * Flag for the use of VTK unstructured grids.
     */
    bool unstructured_grid_physical_;

    /*
     * Flag for the use of VTK unstructured grids.
     */
    bool unstructured_grid_parametric_;

    /*
     * Container for the mapping and field functions.
     */
    std::shared_ptr<FunctionsContainer> funcs_container_;

    /*
     * Generates the physical vtk grids.
     */
    template <int dim, int codim>
    void
    generate_solid_mesh_grids(const MapFunPtr_<dim, codim> mapping,
                              const bool is_identity,
                              const Index &vtk_block_id,
                              vtkMultiBlockDataSet *const vtk_block) const;

    /*
     * Generates the control mesh vtk grids.
     */
    template <int dim, int codim>
    void
    generate_control_mesh_grids(const MapFunPtr_<dim, codim> mapping,
                                const Index &vtk_block_id,
                                vtkMultiBlockDataSet *const vtk_block) const;

    /*
     * Generates the knot mesh vtk grids.
     */
    template <int dim, int codim>
    void
    generate_knot_mesh_grids(const MapFunPtr_<dim, codim> mapping,
                             const bool is_identity,
                             const Index &vtk_block_id,
                             vtkMultiBlockDataSet *const vtk_block) const;

    /*
     * Returns the number of identity functions (first entry of the array),
     * not identity functions (second entry in the array and ig function functions
     * (third entry).
     */
    SafeSTLArray<Size, 3> get_number_functions() const;

    /*
     * Create the cell ids container needed for defining vtk cells.
     */
    template <int dim>
    vtkSmartPointer<vtkCellArray>
    create_cells_solid_vtu_grid(const TensorSize<dim> &n_visualization_elements,
                                const Size &n_bezier_elements,
                                const bool quadratic_cells) const;

    /*
     * Creates a VTK unstructured grid for the solid block.
     */
    template <int dim, int codim>
    vtkSmartPointer<vtkUnstructuredGrid>
    create_solid_vtu_grid(const MapFunPtr_<dim, codim> mapping,
                          const TensorSize<dim> &n_vis_elements,
                          const bool quadratic_cells) const;

    /*
     * Creates a VTK structured grid for the solid block.
     * For the 1D case, a VTK unstructured grid is returned.
     */
    template <int dim, int codim>
    vtkSmartPointer<vtkPointSet>
    create_solid_vts_grid(const MapFunPtr_<dim, codim> mapping,
                          const TensorSize<dim> &n_vis_elements) const;

    /*
     * Creates a the points for the solid VTK grid.
     */
    template <int dim, int codim>
    vtkSmartPointer<vtkPoints>
    create_points_solid_vtk_grid(const MapFunPtr_<dim, codim> mapping,
                                 const TensorSize<dim> &n_vis_elements,
                                 const bool is_structured,
                                 const bool is_quadratic,
                                 vtkPointData *const point_data) const;

    template <int dim, int codim>
    void
    create_point_data_dim_codim(const std::shared_ptr<Function<dim, 0, dim + codim, 1>> mapping,
                                const Quadrature<dim> &quad,
                                const SafeSTLVector<SafeSTLVector<Index>> &points_map,
                                const SafeSTLVector<Index> &points_mask,
                                vtkPointData *const data,
                                typename std::enable_if_t<(dim == 1 && codim == 0)>* = 0) const;

    template <int dim, int codim>
    void
    create_point_data_dim_codim(const std::shared_ptr<Function<dim, 0, dim + codim, 1>> mapping,
                                const Quadrature<dim> &quad,
                                const SafeSTLVector<SafeSTLVector<Index>> &points_map,
                                const SafeSTLVector<Index> &points_mask,
                                vtkPointData *const data,
                                typename std::enable_if_t<!(dim == 1 && codim == 0)>* = 0) const;

    template <int dim, int codim, int range, int rank>
    void
    create_point_data(const std::shared_ptr<Function<dim, 0, dim + codim, 1>> mapping,
                      const Quadrature<dim> &quad,
                      const SafeSTLVector<SafeSTLVector<Index>> &point_num_map,
                      const SafeSTLVector<Index> &points_mask,
                      vtkPointData *const data) const;

    /*
     * Create the quadratures for VTK quadratic cells.
     */
    template <int dim>
    static std::shared_ptr<Quadrature<dim>>
                                         create_visualization_quadrature(const TensorSize<dim> &n_elements_per_direction,
                                                 const bool is_quadratic);

    /*
     * Creates the mapping between the number of the points local to the element
     * and the global points in the VTK grid.
     * It also returns the mask for choosing the points in the element.
     */
    template <int dim>
    static void
    create_points_numbering_map(const std::shared_ptr<const CartesianGrid<dim>> grid,
                                const std::shared_ptr<Quadrature<dim>> quad,
                                const bool is_structured,
                                const bool is_quadratic,
                                SafeSTLVector<SafeSTLVector<Index>> &points_map,
                                SafeSTLVector<Index> &points_mask,
                                Size &n_total_points);

    /*
     * Creates the connectivity of the VTK quadratic cells for a single
     * Bezier element.
     */
    template <int dim>
    static void create_VTK_quadratic_element_connectivity
    (const TensorSize<dim> &n_vis_elements,
     SafeSTLVector<SafeSTLVector<Index>> &connectivity,
     Size &n_points_per_bezier_element);


    void
    inline
    tensor_to_tuple(const Tdouble t, Real *const tuple, int &pos) const
    {
        *(tuple + pos++) = t;
    }


    template <class Tensor>
    inline
    void
    tensor_to_tuple(const Tensor &t, Real *const tuple, int &pos) const
    {
        for (int i = 0; i < Tensor::size; ++i)
            tensor_to_tuple(t[i], tuple, pos);
    };

    template <class Tensor>
    inline
    void
    tensor_to_tuple(const Tensor &t, Real *const tuple) const
    {
        int pos = 0;
        tensor_to_tuple(t, tuple, pos);
    };
};

IGA_NAMESPACE_CLOSE

#endif // IGA_VTK_H_
