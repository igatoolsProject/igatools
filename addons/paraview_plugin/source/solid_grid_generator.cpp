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

#include <paraview_plugin/solid_grid_generator.h>

#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <boost/range/irange.hpp>

#include <igatools/base/quadrature_lib.h>
#include <igatools/functions/function.h>
#include <igatools/functions/function_element.h>
#include <igatools/functions/functions_container.h>
#include <igatools/utils/multi_array_utils.h>
#include <paraview_plugin/grid_information.h>
using std::shared_ptr;

IGA_NAMESPACE_OPEN

template <int dim, int codim>
VtkIgaSolidGridGenerator<dim, codim>::
VtkIgaSolidGridGenerator (const MapFunPtr_ mapping,
                       const GridInfoPtr_ grid_info,
                       const FunContPtr_t_ func_container)
    :
                         map_fun_ (mapping),
                         grid_info_ (grid_info),
                         funcs_container_ (func_container)
{
    Assert (mapping != nullptr, ExcNullPtr());
    Assert (grid_info != nullptr, ExcNullPtr());

    this->init_points_info();
}



template <int dim, int codim>
auto
VtkIgaSolidGridGenerator<dim, codim>::
get_grid (const MapFunPtr_ mapping,
          const GridInfoPtr_ grid_info,
          const FunContPtr_t_ func_container) -> VtkGridPtr_
{
    VtkIgaSolidGridGenerator generator (mapping, grid_info, func_container);
    return generator.create_grid();
}



template <int dim, int codim>
auto
VtkIgaSolidGridGenerator<dim, codim>::
create_grid () const -> VtkGridPtr_
{
    VtkGridPtr_ grid;

    // 1D geometries must be unstructured.
    if (!grid_info_->is_structured () || dim == 1) // VTK unstructured grid.
        grid =  this->create_grid_vtu ();
    else // VTK structured grid.
        grid = this->create_grid_vts ();

    this->create_point_data_dim_codim<dim, codim> (grid->GetPointData ());

    return grid;
}



template <int dim, int codim>
auto
VtkIgaSolidGridGenerator<dim, codim>::
create_grid_vts () const ->
VtkGridPtr_
{
    AssertThrow (dim != 1, ExcMessage("For 1D is not possible to create vtk"
                                      " structured meshes."))
    auto vts_grid = vtkSmartPointer <vtkStructuredGrid>::New ();

    const auto points = this->create_points();

    const auto n_intervals = map_fun_->get_grid ()->get_num_intervals ();
    int grid_dim[3] = { 1, 1, 1 };
    for (const auto &dir : boost::irange (0, dim))
        grid_dim[dir] = (n_vis_elements_[dir] + 1) * n_intervals[dir];
    vts_grid->SetDimensions (grid_dim[0], grid_dim[1], grid_dim[2]);

    vts_grid->SetPoints (points);

    return vts_grid;
}



template <int dim, int codim>
auto
VtkIgaSolidGridGenerator<dim, codim>::
create_grid_vtu () const ->VtkGridPtr_
{
    const bool quad_cells = grid_info_->is_quadratic();

    auto vtu_grid = vtkSmartPointer <vtkUnstructuredGrid>::New ();

    const auto points = this->create_points();

    vtu_grid->SetPoints (points);

    const Size n_bezier_elements = points_map_.size ();

    // Number of cells per Bezier element.
    const Size n_cells_per_bezier = n_vis_elements_.flat_size ();

    // Total number of cells.
    const int n_total_cells = n_bezier_elements * n_cells_per_bezier;

    auto cell_ids = vtkSmartPointer <vtkIdTypeArray>::New ();

    const Size n_points_per_single_cell = connectivity_[0].size ();
    const Size tuple_size = n_points_per_single_cell + 1;
    cell_ids->SetNumberOfComponents (tuple_size);
    cell_ids->SetNumberOfTuples (n_total_cells);

    vtkIdType *tuple = new vtkIdType[tuple_size];
    tuple[0] = tuple_size - 1;

    Index cell_id = 0;
    for (int be = 0; be < n_bezier_elements; ++be)
    {
        const int vtk_vertex_id_offset = points_mask_.size() * be;

        // Iterating along the cells of one Bezier element.
        for (const auto &cc : connectivity_)
        {
            int point_id = 1;
            for (const auto &c : cc)
                tuple[point_id++] = c + vtk_vertex_id_offset;

            cell_ids->SetTupleValue (cell_id++, tuple);
        }
    }

    const auto cells = vtkSmartPointer <vtkCellArray>::New ();
    cells->SetCells (cell_ids->GetNumberOfTuples (), cell_ids);

    vtu_grid->Allocate (cells->GetNumberOfCells (), 0);
    const int vtk_enum_type =
            quad_cells ?
                    dim == 3 ? VTK_QUADRATIC_HEXAHEDRON :
                    dim == 2 ? VTK_QUADRATIC_QUAD :
                               VTK_QUADRATIC_EDGE
                    :
                    dim == 3 ? VTK_HEXAHEDRON :
                    dim == 2 ? VTK_QUAD :
                               VTK_LINE;

    vtu_grid->SetCells (vtk_enum_type, cells);

    return vtu_grid;
}



template <int dim, int codim>
vtkSmartPointer<vtkPoints>
VtkIgaSolidGridGenerator<dim, codim>::
create_points () const
{
    const auto points = vtkSmartPointer <vtkPoints>::New ();
    points->SetNumberOfPoints (n_total_points_);

    double point_tmp[3] = { 0.0, 0.0, 0.0 };

    const auto &map_fun = map_fun_;

    const auto flag = ValueFlags::point | ValueFlags::value;
    map_fun->reset (flag, *quad_);

    auto elem = map_fun->begin ();
    const auto end = map_fun->end ();

    const auto topology = Topology <dim> ();
    map_fun->init_cache (elem, topology);

    for (auto pm_el = points_map_.cbegin (); elem != end; ++elem, ++pm_el)
    {
        map_fun->fill_cache (elem, topology, 0);

        auto element_vertices_tmp =
                elem->template get_values <_Value, dim> (0);
        auto pm = pm_el->cbegin ();
        for (const auto &mask : points_mask_)
        {
            const auto &point = element_vertices_tmp[mask];
            for (const auto &dir : boost::irange (0, space_dim))
                point_tmp[dir] = point[dir];
            points->SetPoint (*pm++, point_tmp);
        }
    }

    return points;
}



template <int dim, int codim>
void
VtkIgaSolidGridGenerator<dim, codim>::
init_points_info ()
{
    const auto &num_visualization_elements =
            grid_info_->get_num_cells_per_element ();

    for (int dir = 0; dir < dim; ++dir)
    {
        n_vis_elements_[dir] = num_visualization_elements[dir];
        Assert (n_vis_elements_[dir] > 0,
                ExcMessage("The number of visualization elements must be > 0."))
    }

    this->create_visualization_quadrature();

    points_map_.clear ();
    points_mask_.clear ();

    const auto cartesian_grid = map_fun_->get_grid ();
    const Size n_bezier_elements = cartesian_grid->get_num_all_elems ();

    AssertThrow(n_bezier_elements > 0,
                ExcMessage ("0 Bezier elements found."));

    points_map_.resize (n_bezier_elements);

    if (!grid_info_->is_structured() || dim == 1) // VTK unstructured grid
    {
        if (grid_info_->is_quadratic()) // VTK quadratic elements
        {
            // Number of visualization points per Bezier element in each direction.
            // const auto n_pts_dir_per_bezier_elem = quad->get_num_coords_direction();

            // Iteration along all the quadrature point.
            // Only the points in an edge are added to the mask.
            const auto &quad_points_1d = quad_->get_points_1d ();
            const Size n_quad_points = quad_->get_num_points ();
            for (int i_pt = 0; i_pt < n_quad_points; ++i_pt)
            {
                const auto tensor_id = quad_points_1d.flat_to_tensor (i_pt);

                // If the number of odd values is > 1. This means that the point is
                // not in an edge.
                // For the case 1D, this conditions is never fulfilled, no
                // points will be removed (this is what is desired).
                Size n_odd_values = 0;
                for (const auto &dir : boost::irange (0, dim))
                    n_odd_values += tensor_id[dir] % 2;

                if (n_odd_values < 2) // It's on an edge.
                    points_mask_.push_back (i_pt);
            }

            // Total number of visualization points per Bezier element.
            const Size n_pts_per_bezier_elem = points_mask_.size ();

            Index point_id = 0;
            for (auto &pm_el : points_map_)
            {
                pm_el.resize (n_pts_per_bezier_elem);
                for (auto &pm : pm_el)
                    pm = point_id++;
            } // points_map_
            n_total_points_ = point_id;

            this->create_quadratic_element_connectivity<dim>();

        } // VTK quadratic elements

        else // VTK linear elements
        {
            // Total number of visualization points per Bezier element.
            const Size n_pts_per_bezier_elem = quad_->get_num_points ();

            points_mask_.resize (n_pts_per_bezier_elem);
            Index id = 0;
            for (auto &pm : points_mask_)
                pm = id++;

            Index point_id = 0;
            for (auto &pm_el : points_map_)
            {
                pm_el.resize (n_pts_per_bezier_elem);
                for (auto &pm : pm_el)
                    pm = point_id++;
            } // points_map_
            n_total_points_ = point_id;

            this->create_linear_element_connectivity();

        } // VTK linear elements
    }
    else // VTK structured grid
    {
        // Total number of visualization points per Bezier element.
        const Size n_pts_per_bezier_elem = quad_->get_num_points ();
        // Number of visualization points per Bezier element in each direction.
        const auto n_pts_dir_per_bezier_elem =
                quad_->get_num_coords_direction ();

        points_mask_.resize (n_pts_per_bezier_elem);
        Index id = 0;
        for (auto &pm : points_mask_)
            pm = id++;

        n_total_points_ = n_pts_per_bezier_elem * n_bezier_elements;

        TensorSize <dim> n_pts_per_mesh; // Number of points per direction of
        // VTK structured grid.
        const auto n_intervals = cartesian_grid->get_num_intervals ();
        for (const auto &dir : boost::irange (0, dim))
            n_pts_per_mesh[dir] = n_intervals[dir]
                                              * n_pts_dir_per_bezier_elem[dir];

        TensorIndex <dim> elem_t_id; // Tensorial index of the Bezier element.
        TensorIndex <dim> pt_mesh_t_offset; // Tensorial index of the first point in
        // Bezier element.
        TensorIndex <dim> pt_mesh_t_id;     // Tensorial index of the point.
        TensorIndex <dim> pt_elem_t_id; // Tensorial index of the point referred
        // to the number of points in a
        // single element.

        const auto w_elem_pts = MultiArrayUtils <dim>::compute_weight (
                n_pts_dir_per_bezier_elem);
        const auto w_mesh_pts = MultiArrayUtils <dim>::compute_weight (
                n_pts_per_mesh);

        for (const auto &i_el : boost::irange (0, n_bezier_elements))
        {
            auto &pmi = points_map_[i_el];
            pmi.resize (n_pts_per_bezier_elem);

            // Computing the tensor index of the first point of the element.
            elem_t_id = cartesian_grid->flat_to_tensor (i_el);
            for (const auto &dir : boost::irange (0, dim))
                pt_mesh_t_offset[dir] = elem_t_id[dir]
                                                  * n_pts_dir_per_bezier_elem[dir];

            Index point_id = 0;
            for (auto &pm : pmi)
            {
                // Computing the tensor index of the point referred to the number
                // of points into a single Bezier element.
                pt_elem_t_id = MultiArrayUtils <dim>::flat_to_tensor_index (
                        point_id++, w_elem_pts);

                // Computing the tensor index of the point referred to the number
                // of points into the whole mesh.
                for (const auto &dir : boost::irange (0, dim))
                    pt_mesh_t_id[dir] = pt_mesh_t_offset[dir]
                                                         + pt_elem_t_id[dir];

                pm = MultiArrayUtils <dim>::tensor_to_flat_index (
                        pt_mesh_t_id, w_mesh_pts);
            }
        }
    }
}



template <int dim, int codim>
void
VtkIgaSolidGridGenerator<dim, codim>::
create_visualization_quadrature()
{
    const Size n_pts_per_cell_dir =  grid_info_->is_quadratic() ? 3 : 2;
    TensorSize <dim> n_points;
    for (const auto &dir : boost::irange (0, dim))
        n_points[dir] = n_vis_elements_[dir] * (n_pts_per_cell_dir - 1) + 1;

    quad_ = QUniform <dim>::create (n_points);

#if 0
    // TODO: Right now is not possible to use non tensor product quadrature
    // for cache elements in igatools. It will be in the future.
    // The code below defines a non tensor product quadrature for VTK quadratic
    // cells.

    const auto &uniform_array_points = uniform_quad.get_points_1d();
    const Size n_total_points = uniform_quad.get_num_points();

    // The points that are not in an edge are removed from the quadrature points.
    // The condition to be in an edge is that just one or zero directions are
    // active. (This also works for the 1D case, where no point is removed.)
    SafeSTLVector<Index> point_indices;
    for (int i_pt = 0; i_pt < n_total_points; ++i_pt)
    {
        const auto t_id = uniform_array_points.flat_to_tensor(i_pt);

        int active_directions = 0;
        for (int dir = 0; dir < dim; ++dir)
            if ((t_id[dir] > 0) && (t_id[dir] < (n_points[dir] - 1)))
                ++active_directions;

        if (active_directions < 2)
            point_indices.push_back(i_pt);
    }

    ValueVector<Points<dim>> quad_points(point_indices.size());
    auto qp = quad_points.begin();
    for (const auto &id : point_indices)
        *qp++ = uniform_quad.get_point(id);

    BBox<dim> box;

    return shared_ptr<Quadrature<dim>> (new Quadrature<dim> (quad_points,
                                                             uniform_quad.get_weights_1d(), box));
#endif
}



template <int dim, int codim>
void
VtkIgaSolidGridGenerator<dim, codim>::
create_linear_element_connectivity ()
{
    // Number of vertices in a dim-dimensional square.
    static constexpr int n_vertices =
            UnitElement <dim>::template num_elem <0> ();

    TensorSize <dim> n_elem_bound_per_dir;
    for (const auto &dir : boost::irange (0, dim))
        n_elem_bound_per_dir[dir] = n_vis_elements_[dir] + 1;

    // This grid is going to help in building the connectivity.
    // Every element of the grid refers to a cell.
    const auto cells_grid = CartesianGrid <dim>::create (
            n_elem_bound_per_dir);

    const Size n_cells_per_bezier = n_vis_elements_.flat_size ();
    connectivity_.resize (n_cells_per_bezier);

    const Size n_points_per_single_cell = n_vertices;

    // Creating the connectivity ---------------------------------------------//

    // Building the offsets container. According to the vtk elements connectivity.
    using T_ = SafeSTLArray < SafeSTLArray<int, dim>, n_vertices>;
    const T_ delta_idx =
            dim == 1 ? T_ ({{ 0 },
                            { 1 } }) :   // dim = 1

            dim == 2 ? T_ ({{ 0, 0 },
                            { 1, 0 },
                            { 1, 1 },
                            { 0, 1 } }) :   // dim = 2

                       T_ ({
                            { 0, 0, 0 },
                            { 1, 0, 0 },
                            { 1, 1, 0 },
                            { 0, 1, 0 },
                            { 0, 0, 1 },
                            { 1, 0, 1 },
                            { 1, 1, 1 },
                            { 0, 1, 1 } }); // dim = 3

    const TensorIndex <dim> weight_points =
        MultiArrayUtils <dim>::compute_weight (n_elem_bound_per_dir);

    TensorIndex <dim> vtk_vertex_tensor_idx;

    auto conn_el = connectivity_.begin ();
    auto cell = cells_grid->begin ();
    const auto end = cells_grid->end ();
    for (; cell != end; ++cell, ++conn_el)
    {
        conn_el->resize (n_points_per_single_cell);

        SafeSTLArray <Index, dim> vtk_elem_tensor_idx =
                cell->get_tensor_index ();

        auto conn = conn_el->begin ();
        for (int iVertex = 0; iVertex < n_points_per_single_cell;
                ++iVertex, ++conn)
        {
            for (int i = 0; i < dim; ++i)
                vtk_vertex_tensor_idx[i] = vtk_elem_tensor_idx[i]
                                                               + delta_idx[iVertex][i];

            *conn = MultiArrayUtils <dim>::tensor_to_flat_index
                    (vtk_vertex_tensor_idx, weight_points);
        }
    }
    //--------------------------------------------------------------------------//

}



template <int dim, int codim>
template <int aux_dim>
void
VtkIgaSolidGridGenerator<dim, codim>::
create_quadratic_element_connectivity (
        typename std::enable_if_t<aux_dim == 1>*)
{
    //  This is the connectivity pattern of the 1D VTK quadratic element.
    //
    //         0 -- 2 -- 1

    // Number of cells per Bezier element.
    const Size n_cells_per_bezier = n_vis_elements_.flat_size ();

    connectivity_.clear ();
    connectivity_.resize (n_cells_per_bezier);

    Index point_id = 0;
    for (auto &conn_el : connectivity_)
    {
        conn_el =
                {   point_id, point_id+2, point_id+1};
        point_id += 2;
    }
}


template <int dim, int codim>
template <int aux_dim>
void
VtkIgaSolidGridGenerator<dim, codim>::
create_quadratic_element_connectivity (
        typename std::enable_if_t<aux_dim == 2>*)
{
    //  This is the connectivity pattern of the 2D VTK quadratic element.
    //
    //         3 -- 6 -- 2
    //         |         |
    //         7         5
    //         |         |
    //         0 -- 4 -- 1

    // Number of cells per Bezier element.
    const Size n_cells_per_bezier = n_vis_elements_.flat_size ();

    connectivity_.clear ();
    connectivity_.resize (n_cells_per_bezier);

    TensorSize <aux_dim> n_full_points_per_dir;
    for (int dir = 0; dir < aux_dim; ++dir)
        n_full_points_per_dir[dir] = n_vis_elements_[dir] * 2 + 1;

    static const int n_points_per_single_cell = 8;

    TensorSize <aux_dim> n_elem_bound_per_dir;
    for (int dir = 0; dir < aux_dim; ++dir)
        n_elem_bound_per_dir[dir] = n_vis_elements_[dir] + 1;

    // This grid is going to help in building the connectivity.
    // Every element of the grid refers to a cell.
    const auto cells_grid = CartesianGrid <aux_dim>::create (
            n_elem_bound_per_dir);

    // This array constaints the offsets of the points along the
    // first (u) direction.
    //   Offset for the line:  3 -- 6 -- 2 -> offset = 2
    //   Offset for the line:  7         5 -> offset = 1
    //   Offset for the line:  0 -- 4 -- 1 -> offset = 2

    const SafeSTLArray <Index, 8> offsets_u =
            {
                    { 2, 2, 2, 2, 2, 1, 2, 1 } };

    // This container if for the offsets of the points along the
    // second (v) direction.
    const Index offsets_v = n_full_points_per_dir[0] * 2
            - n_vis_elements_[0];

    SafeSTLVector <Index> vtk_vertex_id_0 (n_points_per_single_cell);
    vtk_vertex_id_0[0] = 0;
    vtk_vertex_id_0[4] = 1;
    vtk_vertex_id_0[1] = 2;

    vtk_vertex_id_0[7] = n_full_points_per_dir[0];
    vtk_vertex_id_0[5] = vtk_vertex_id_0[7] + 1;

    vtk_vertex_id_0[3] = vtk_vertex_id_0[7] + n_full_points_per_dir[0]
                                                                    - n_vis_elements_[0];
    vtk_vertex_id_0[6] = vtk_vertex_id_0[3] + 1;
    vtk_vertex_id_0[2] = vtk_vertex_id_0[3] + 2;

    auto conn_el = connectivity_.begin ();
    auto cell = cells_grid->begin ();
    const auto end = cells_grid->end ();
    for (; cell != end; ++cell, ++conn_el)
    {
        conn_el->resize (n_points_per_single_cell);

        SafeSTLArray <Index, aux_dim> vtk_elem_tensor_idx =
                cell->get_tensor_index ();
        auto conn = conn_el->begin ();
        for (int i_pt = 0; i_pt < n_points_per_single_cell; ++i_pt, ++conn)
            *conn = vtk_vertex_id_0[i_pt]
                                    + offsets_u[i_pt] * vtk_elem_tensor_idx[0]
                                                                            + offsets_v * vtk_elem_tensor_idx[1];
    }
}



template <int dim, int codim>
template <int aux_dim>
void
VtkIgaSolidGridGenerator<dim, codim>::
create_quadratic_element_connectivity (
        typename std::enable_if_t<aux_dim == 3>*)
{
    //  This is the connectivity pattern of the 3D VTK quadratic element.
    //
    //         w = 0             w = 1             w = 2
    //
    //     3 -- 10 --  2    19 -------- 18     7 -- 14 --  6
    //     |           |     |           |     |           |
    //    11           9     |           |    15          13
    //     |           |     |           |     |           |
    //     0 --  8 --  1    16 -------- 17     4 -- 12 --  5
    //

    // Number of cells per Bezier element.
    const Size n_cells_per_bezier = n_vis_elements_.flat_size ();

    connectivity_.clear ();
    connectivity_.resize (n_cells_per_bezier);

    TensorSize <aux_dim> n_full_points_per_dir;
    for (int dir = 0; dir < aux_dim; ++dir)
        n_full_points_per_dir[dir] = n_vis_elements_[dir] * 2 + 1;

    static const int n_points_per_single_cell = 20;

    TensorSize <aux_dim> n_elem_bound_per_dir;
    for (int dir = 0; dir < aux_dim; ++dir)
        n_elem_bound_per_dir[dir] = n_vis_elements_[dir] + 1;

    // This grid is going to help in building the connectivity.
    // Every element of the grid refers to a cell.
    const auto cells_grid = CartesianGrid <dim>::create (
            n_elem_bound_per_dir);

    // This array constaints the offsets of the points along the
    // first (u) direction.
    //   Offset for the line:  3 -- 10 --  2   -> offset = 2
    //   Offset for the line: 11 --------  9   -> offset = 1
    //   Offset for the line:  0 --  8 --  1   -> offset = 2

    //   Offset for the line: 19 -------- 18   -> offset = 1
    //   Offset for the line: 16 -------- 17   -> offset = 1

    //   Offset for the line:  7 -- 14 --  6   -> offset = 2
    //   Offset for the line: 15 -------- 13   -> offset = 1
    //   Offset for the line:  4 -- 12 --  5   -> offset = 2

    const SafeSTLArray <Index, 20> offsets_u =
            {
                    { 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 2,
                            1, 1, 1,
                            1, 1 } };

    TensorSize <aux_dim> n_full_points_per_dir_red;
    for (int dir = 0; dir < aux_dim; ++dir)
        n_full_points_per_dir_red[dir] = n_full_points_per_dir[dir]
                                                               - n_vis_elements_[dir];

    // This array constaints the offsets of the points along the
    // second (v) direction.
    //   Offset for the line:  3 -- 10 --  2   -> offset = ov0
    //   Offset for the line: 11 --------  9   -> offset = ov0
    //   Offset for the line:  0 --  8 --  1   -> offset = ov0

    //   Offset for the line: 19 -------- 18   -> offset = ov1
    //   Offset for the line: 16 -------- 17   -> offset = ov1

    //   Offset for the line:  7 -- 14 --  6   -> offset = ov0
    //   Offset for the line: 15 -------- 13   -> offset = ov0
    //   Offset for the line:  4 -- 12 --  5   -> offset = ov0

    const Index ov0 = n_full_points_per_dir[0]
                                            + n_full_points_per_dir_red[0];
    const Index ov1 = n_full_points_per_dir_red[0];
    const SafeSTLArray <Index, 20> offsets_v =
            {
                    {
                            ov0, ov0, ov0, ov0, ov0, ov0, ov0, ov0,
                            ov0, ov0, ov0, ov0, ov0, ov0,
                            ov0, ov0, ov1, ov1, ov1, ov1
                    }
            };

    // Offset along the third (w) direction.
    const Size ow0 = n_full_points_per_dir[0] * n_full_points_per_dir_red[1]
                                                                          + n_full_points_per_dir_red[0] * n_vis_elements_[1];
    const Size ow1 = n_full_points_per_dir_red[0]
                                               * n_full_points_per_dir_red[1];

    const Index offsets_w = ow0 + ow1;

    SafeSTLVector <Index> vtk_vertex_id_0 (n_points_per_single_cell);
    vtk_vertex_id_0[0] = 0;
    vtk_vertex_id_0[8] = 1;
    vtk_vertex_id_0[1] = 2;
    vtk_vertex_id_0[11] = vtk_vertex_id_0[0] + n_full_points_per_dir[0];
    vtk_vertex_id_0[9] = vtk_vertex_id_0[11] + 1;
    vtk_vertex_id_0[3] = vtk_vertex_id_0[11] + n_full_points_per_dir_red[0];
    vtk_vertex_id_0[10] = vtk_vertex_id_0[3] + 1;
    vtk_vertex_id_0[2] = vtk_vertex_id_0[3] + 2;

    vtk_vertex_id_0[16] = vtk_vertex_id_0[0] + ow0;
    vtk_vertex_id_0[17] = vtk_vertex_id_0[16] + 1;
    vtk_vertex_id_0[19] = vtk_vertex_id_0[16]
                                          + n_full_points_per_dir_red[0];
    vtk_vertex_id_0[18] = vtk_vertex_id_0[19] + 1;

    vtk_vertex_id_0[4] = vtk_vertex_id_0[16] + ow1;
    vtk_vertex_id_0[12] = vtk_vertex_id_0[4] + 1;
    vtk_vertex_id_0[5] = vtk_vertex_id_0[4] + 2;
    vtk_vertex_id_0[15] = vtk_vertex_id_0[4] + n_full_points_per_dir[0];
    vtk_vertex_id_0[13] = vtk_vertex_id_0[15] + 1;
    vtk_vertex_id_0[7] = vtk_vertex_id_0[15] + n_full_points_per_dir_red[0];
    vtk_vertex_id_0[14] = vtk_vertex_id_0[7] + 1;
    vtk_vertex_id_0[6] = vtk_vertex_id_0[7] + 2;

    auto conn_el = connectivity_.begin ();
    auto cell = cells_grid->begin ();
    const auto end = cells_grid->end ();
    for (; cell != end; ++cell, ++conn_el)
    {
        conn_el->resize (n_points_per_single_cell);

        SafeSTLArray <Index, aux_dim> vtk_elem_tensor_idx =
                cell->get_tensor_index ();
        auto conn = conn_el->begin ();
        for (int i_pt = 0; i_pt < n_points_per_single_cell; ++i_pt, ++conn)
            *conn = vtk_vertex_id_0[i_pt]
                                    + offsets_u[i_pt] * vtk_elem_tensor_idx[0]
                                                                            + offsets_v[i_pt] * vtk_elem_tensor_idx[1]
                                                                                                                    + offsets_w * vtk_elem_tensor_idx[2];
    }
}


template <int dim, int codim>
template <int aux_dim, int aux_codim>
void
VtkIgaSolidGridGenerator<dim, codim>::
create_point_data_dim_codim (vtkPointData * const point_data,
        typename std::enable_if_t <(aux_dim == 1 && aux_codim == 0)>*) const
{
    this->template create_point_data <1, 1> (point_data);
}



template <int dim, int codim>
template <int aux_dim, int aux_codim>
void
VtkIgaSolidGridGenerator<dim, codim>::
create_point_data_dim_codim (vtkPointData * const point_data,
        typename std::enable_if_t <!(aux_dim == 1 && aux_codim == 0)>*) const
{
    this->template create_point_data <1, 1> (point_data);
    this->template create_point_data <dim + codim, 1> (point_data);
}


template <int dim, int codim>
template <int range, int rank>
void
VtkIgaSolidGridGenerator<dim, codim>::
create_point_data (vtkPointData * const point_data) const
{

    // Getting the functions associated with the map.
    const auto &funcs_map = funcs_container_->template
            get_functions_associated_to_mapping <dim, codim, range, rank> (
                    map_fun_);

    using Value = typename Function<dim, codim, range, rank>::Value;

    static constexpr Size n_comp = Value::size;
    Real tuple[n_comp];

    const auto flag = ValueFlags::value | ValueFlags::point;
    for (const auto &it : funcs_map)
    {
        const auto &fun = it.second;
        const auto &name = fun->get_name ();

        const int n_bezier_elements = points_map_.size ();
        const vtkIdType n_tuples = n_bezier_elements * points_mask_.size();

        auto arr = vtkSmartPointer <vtkDoubleArray>::New ();

        arr->SetName (name.c_str ());
        arr->SetNumberOfComponents (n_comp);
        arr->SetNumberOfTuples (n_tuples);

        fun->reset (flag, *quad_);

        auto elem = fun->begin ();
        const auto end = fun->end ();

        const auto topology = Topology <dim> ();
        fun->init_cache (elem, topology);

        auto pnm_it = points_map_.cbegin ();
        for (; elem != end; ++elem, ++pnm_it)
        {
            fun->fill_cache (elem, topology, 0);

            auto pnm = pnm_it->cbegin ();
            auto values = elem->template get_values <_Value, dim> (0);
            for (const auto &pm : points_mask_)
            {
                this->template tensor_to_tuple <Value> (values[pm], tuple);
                arr->SetTuple (*pnm++, tuple);
            }
        }

        point_data->AddArray (arr.Get ());
    }

}



template <int dim, int codim>
void
VtkIgaSolidGridGenerator<dim, codim>::
tensor_to_tuple(const Tdouble t, Real *const tuple, int &pos) const
{
    *(tuple + pos++) = t;
}



template <int dim, int codim>
template <class Tensor>
void
VtkIgaSolidGridGenerator<dim, codim>::
tensor_to_tuple(const Tensor &t, Real *const tuple, int &pos) const
{
    for (int i = 0; i < Tensor::size; ++i)
        tensor_to_tuple(t[i], tuple, pos);
}



template <int dim, int codim>
template <class Tensor>
void
VtkIgaSolidGridGenerator<dim, codim>::
tensor_to_tuple(const Tensor &t, Real *const tuple) const
{
    int pos = 0;
    tensor_to_tuple(t, tuple, pos);
}


template class VtkIgaSolidGridGenerator<1, 0>;
template class VtkIgaSolidGridGenerator<1, 1>;
template class VtkIgaSolidGridGenerator<1, 2>;
template class VtkIgaSolidGridGenerator<2, 0>;
template class VtkIgaSolidGridGenerator<2, 1>;
template class VtkIgaSolidGridGenerator<3, 0>;


IGA_NAMESPACE_CLOSE
