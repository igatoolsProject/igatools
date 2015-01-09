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


#ifndef __GRID_WRAPPER_H_
#define __GRID_WRAPPER_H_

#include <igatools/base/config.h>


#include <array>
#include <memory>
#include <boost/signals2.hpp>

IGA_NAMESPACE_OPEN



//TODO(pauletti, Mar 4, 2014): the purpose of this class is not well explained
/**
 * @brief This class wraps a grid of type @p GridType with a std::shared_ptr.
 *
 * It is used as base class for classes that are based on a grid, e.g. BSplineSpace,
 * PhysicalSpace, Mapping, etc.
 *
 * It also provide an handler for the h-refinement mechanism based
 * on the signal/slot technology implemented by the
 * <a href="http://www.boost.org/doc/libs/1_55_0/doc/html/signals2.html">Boost.Signals2</a> library.
 *
 * @ingroup h_refinement
 */
template<class Grid_>
class GridWrapper
{
public:

    using GridType = Grid_;

    static constexpr std::array<Size, Grid_::dim> dims = Grid_::dims;

    /** @name Constructor and destructor. */
    ///@{
    /** Default constructor. Not allowed to be used. */
    GridWrapper() = delete;

    /** Construct a GridWrapper copying the pointer @p grid. */
    GridWrapper(std::shared_ptr<GridType> grid);

    /** Copy constructor. */
    GridWrapper(const GridWrapper<GridType> &grid) = default;

    /** Move constructor. Not allowed to be used. */
    GridWrapper(GridWrapper<GridType> &&grid) = delete;

    /** Destructor. */
    ~GridWrapper();
    ///@}

    /** @name Assignment operator. */
    ///@{
    /** Copy assignment operator. Not allowed to be used. */
    GridWrapper<GridType> &operator=(const GridWrapper<GridType> &grid) = delete ;

    /** Move assignment operator. Not allowed to be used. */
    GridWrapper<GridType> &operator=(GridWrapper<GridType> &&grid) = delete ;
    ///@}


    /** @name Getting the grid */
    ///@{

    /** Returns the grid (non-const version). */
    std::shared_ptr<GridType> get_grid();

    /** Returns the grid (const version). */
    std::shared_ptr<const GridType> get_grid() const;
    ///@}


public:
    /**
     * @name Functions for performing h-refinement
     */
    ///@{

    /**
     * Perform the h-refinement of the grid
     * along the directions specified by the true values in the entries of the
     * array of bools @p refinement_directions,
     * and with a number of subdivisions for each interval (along each direction)
     * specified by @p n_subdivisions.
     *
     * @note If the i-th direction is not active for the refinement
     * (i.e. <tt>refinement_directions[i] == false</tt>),
     *  then the corresponding value <tt>n_subdivisions[i]</tt> will be ignored.
     *
     * @ingroup h_refinement
     */
    void refine_h_directions(
        const std::array<bool,GridType::dim> &refinement_directions,
        const std::array<Size,GridType::dim> &n_subdiv_directions);

    /**
     * Perform a uniform h-refinement of the grid along the @p direction_id direction,
     * dividing each interval into @p n_subdivisions intervals.
     * @param[in] direction_id Direction along which the refinement is performed.
     * @param[in] n_subdivisions Number of subdivision in which each interval in the grid
     * along the specified direction is divided. This value must be >= 2.
     *
     * @ingroup h_refinement
     */
    void refine_h_direction(const int direction_id, const Size n_subdivisions);


    /**
     * Perform the h-refinement of the grid in all the directions.
     *
     * Each interval in the unrefined grid is uniformly divided in @p n_subdivisions
     * sub-intervals.
     *
     * @ingroup h_refinement
     */
    void refine_h(const Size n_subdivisions = 2);
    ///@}

protected:
    /**
     * Connect the function @p subscriber to the h-refinement signal in the grid object and create
     * the relative connection.
     *
     * @ingroup h_refinement
     */
    void
    connect_refinement_h_function(const typename GridType::SignalRefineSlot &subscriber);

private:
    /** Grid object. */
    std::shared_ptr<GridType> grid_ ;


    /**
     * Connection to the signal for the h-refinement.
     */
    boost::signals2::connection refine_h_connection_ ;
};





IGA_NAMESPACE_CLOSE

#endif /* __GRID_WRAPPER_H_ */
