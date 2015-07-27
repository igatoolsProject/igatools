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

#ifndef VTK_GRID_INFORMATION_H_
#define VTK_GRID_INFORMATION_H_

#include <igatools/base/config.h>
#include <igatools/utils/tensor_size.h>

#include <paraview_plugin/types.h>

IGA_NAMESPACE_OPEN


struct VtkGridInformation
{
private:
    /**
     * Container type for the number of vtk cells per Bezier element.
     */
    typedef TensorSize<3> NumCellsContainer_;

    /**
     * Type of the class.
     */
    typedef VtkGridInformation Self_;

    /**
     * Shared pointer type of the class.
     */
    typedef std::shared_ptr<Self_> SelfPtr_;

    /**
     * Constructor.
     */
    VtkGridInformation (const NumCellsContainer_ &num_cells,
                        const vtkGridType &grid_type);

    /**
     * Default and move constructor and move assignment operators not allowed to be used.
     */
    VtkGridInformation() = delete;
    VtkGridInformation(const VtkGridInformation & grid_info) = delete;
    VtkGridInformation(const VtkGridInformation &&) = delete;
    void operator=(const VtkGridInformation & grid_info) = delete;
    void operator=(const VtkGridInformation &&) = delete;

public:

    /**
     * Creates and returns a shared pointer with a new instance of the class.
     */
    static SelfPtr_ create(const NumCellsContainer_ &num_cells,
                             const vtkGridType &grid_type);

    /**
     * Updates the informations.
     * Returns a boolean, the value is true if the information
     * was updated, false elsewhere.
     */
    bool update(SelfPtr_ grid_info);

    /**
     * Retrieves the grid type.
     */
    const vtkGridType& get_grid_type () const;

    /**
     * Returns true is the grid type is structured, false elsewhere.
     */
    bool is_structured () const;

    /**
     * Returns true is the grid type is quadratic, false elsewhere.
     */
    bool is_quadratic () const;

    /**
     * Returns the number of cells per element.
     */
    const NumCellsContainer_& get_num_cells_per_element () const;

private:

    /**
     * Vtk grid type.
     */
    vtkGridType grid_type_;

    /**
     * Number of cells by Bezier element.
     */
    NumCellsContainer_ cells_per_element_;
};



struct VtkControlGridInformation
{
private:

    /**
     * Type of the class.
     */
    typedef VtkControlGridInformation Self_;

    /**
     * Shared pointer type of the class.
     */
    typedef std::shared_ptr<Self_> SelfPtr_;


    /**
     * Constructor.
     */
    VtkControlGridInformation (const bool structured);

public:
    /**
     * Default constructor.
     */
    VtkControlGridInformation () : VtkControlGridInformation(true) {};
private:


    /**
     * Copy and move constructors and assignment operators not allowed to be used.
     */
    VtkControlGridInformation(const VtkControlGridInformation &) = delete;
    VtkControlGridInformation(const VtkControlGridInformation &&) = delete;
    void operator=(const VtkControlGridInformation &) = delete;
    void operator=(const VtkControlGridInformation &&) = delete;


public:
    /**
     * Creates and returns a shared pointer with a new instance of the class.
     */
    static SelfPtr_ create(const bool structured);

    /**
     * Updates the informations.
     * Returns a boolean, the value is true if the information
     * was updated, false elsewhere.
     */
    bool update(SelfPtr_ grid_info);

    /**
     * Returns true is the grid type is structured, false elsewhere.
     *
     * @Note: this function recalls the function is_structured() from the
     * base class.
     */
    bool is_structured () const;

    /**
     * Retrieves the grid type.
     */
    const vtkGridType& get_grid_type () const;

private:
    /**
     * Vtk grid type.
     */
    vtkGridType grid_type_;

};



IGA_NAMESPACE_CLOSE

#endif // VTK_GRID_INFORMATION_H_
