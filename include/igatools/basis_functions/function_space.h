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

#ifndef __FUNCTION_SPACE_H_
#define __FUNCTION_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/geometry/grid_wrapper.h>

IGA_NAMESPACE_OPEN


/**
 * This class represent the "concept" of isogeometric function space
 * that is built on a certain type of grid.
 *
 * The main feature of this class is that contains a space identifier that is unique to the object,
 * i.e. two different objects will have two different space id.
 *
 * @author M.Martinelli, 2013, 2014.
 */
template<class Grid_>
class FunctionSpaceOnGrid :
    public GridWrapper<Grid_>
{

private:
    using self_t = FunctionSpaceOnGrid<Grid_>;

public:
    using typename GridWrapper<Grid_>::GridType;

    using Topology = typename Grid_::Topology;


protected:
    /** @name Constructor and destructor. */
    ///@{
    /** Default constructor. Not allowed to be used. */
    FunctionSpaceOnGrid() = delete;

    /** Construct the object from the @p grid on which the function space will be built upon. */
    FunctionSpaceOnGrid(std::shared_ptr<GridType> grid);

    /** Copy constructor. */
    FunctionSpaceOnGrid(const self_t &grid) = default;

    /** Move constructor. */
    FunctionSpaceOnGrid(self_t &&grid) = default;

    /** Destructor. */
    ~FunctionSpaceOnGrid() = default;
    ///@}

    /** @name Assignment operator. */
    ///@{
    /** Copy assignment operator. Not allowed to be used. */
    self_t &operator=(const self_t &) = delete;

    /** Move assignment operator. Not allowed to be used. */
    self_t &operator=(self_t &&) = delete;
    ///@}

public:
    /**
     * Returns the space id.
     */
    Index get_space_id() const;

protected:
    Index space_id_ = 0;
};


IGA_NAMESPACE_CLOSE

#endif // __FUNCTION_SPACE_H_
