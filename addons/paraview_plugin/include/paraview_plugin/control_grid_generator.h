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

#ifndef VTK_CONTROL_GRID_GENERATOR_H_
#define VTK_CONTROL_GRID_GENERATOR_H_

#include <igatools/base/config.h>

class vtkPoints;
class vtkPointSet;
template<class T> class vtkSmartPointer;


IGA_NAMESPACE_OPEN

template <int dim, int codim, int range, int rank> class IgFunction;
template <int dim, int codim, int range, int rank> class Function;
struct VtkControlGridInformation;


template <int dim, int codim>
class VtkIgaControlGridGenerator
{
private:
    /**
     * Space dimension.
     */
    static const int space_dim = dim + codim;

    /**
     * Self type.
     */
    typedef VtkIgaControlGridGenerator Self_;

    /**
     * Self shared poitner type.
     */
    typedef std::shared_ptr<Self_> SelfPtr_;

    /**
     * Alias for a shared pointer of a map function type.
     */
    typedef std::shared_ptr<Function<dim, 0, space_dim, 1>> MapFunPtr_;

    /**
     * Alias for mesh grid information shared pointer.
     */
    typedef std::shared_ptr<VtkControlGridInformation> ControlGridInfoPtr_;

    /**
     * Alias for vtk grid object for visualization.
     */
    typedef vtkSmartPointer<vtkPointSet> VtkGridPtr_;

    /**
     * Alias for a ig function type.
     */
    typedef IgFunction<dim, 0, space_dim, 1> IgFun_;

    /**
     * Alias for a shared pointer of a ig function type.
     */
    typedef std::shared_ptr<IgFun_> IgFunPtr_;

    /**
     * Constructor.
     */
    VtkIgaControlGridGenerator(const MapFunPtr_ mapping,
                            const ControlGridInfoPtr_ grid_info);

    /**
     * Constructor, copy and assignment operators not allowed to be used.
     */
    VtkIgaControlGridGenerator() = delete;
    VtkIgaControlGridGenerator(const VtkIgaControlGridGenerator &) = delete;
    VtkIgaControlGridGenerator(const VtkIgaControlGridGenerator &&) = delete;
    void operator=(const VtkIgaControlGridGenerator &) = delete;
    void operator=(const VtkIgaControlGridGenerator &&) = delete;


public:

    /**
     * Creates and returns the vtk grid for the visualization.
     */
    static VtkGridPtr_ get_grid (const MapFunPtr_ mapping,
                                 const ControlGridInfoPtr_ grid_info);

private:

    /**
     * Creates and returns the vtk grid for the visualization.
     */
    VtkGridPtr_ create_grid () const;

    /**
     * Shared pointer of the mapping function representing the geometry.
     */
    const IgFunPtr_ ig_fun_;

    /**
     * Shared pointer of the control grid information for representing the
     * geometry.
     */
    const ControlGridInfoPtr_ grid_info_;

    /**
     * Creates and returns the vtk structured grid for the visualization.
     */
    VtkGridPtr_
    create_grid_vts(const vtkSmartPointer<vtkPoints> points) const;

    /**
     * Creates and returns the vtk unstructured grid for the visualization.
     */
    VtkGridPtr_
    create_grid_vtu(const vtkSmartPointer<vtkPoints> points) const;

};

IGA_NAMESPACE_CLOSE

#endif // VTK_CONTROL_GRID_GENERATOR_H_
