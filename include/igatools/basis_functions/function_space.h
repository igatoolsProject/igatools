//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
 * @brief This class represent the "concept" of function space.
 * @note Its purpose is to use a function trait in order to check if a class is a function
 * space.
 *
 * @sa Is_function_space()
 * @author M.Martinelli, 2013
 */
class FunctionSpace
{};


/**
 * This class represent the "concept" of isogeometric function space
 * that is built on a grid.
 *
 * A function space holds functions
 * \f$\phi:A\to B\f$, where A is called domain and B range.
 *
 * In our particular setting A is an dim dimensional set in
 * R<sup>space_dim</sup> with codim = space_dim - dim non negative.
 *
 * B is a space of tensor \f$\mathcal{T}^{rank}(R^{range})\f$, i.e.
 * the tensor of rank @p rank over the R<sup>range</sup>.
 *
 * @note This class gives a common set of members expected to be used with
 * meta-programming techniques.
 *
 * In particular any function space should provide the following interface
 *
 * @code
 *  // Type of the push-forward.
 *  using PushForwardType = ;
 *
 *  // Type of the reference space.
 *  using RefSpace = ;
 *
 *  // Type of the grid on which the space is built upon.
 *  using GridType = typename PushForwardType::GridType;
 *
 *  // Dimension of the functions domain.
 *  static const int dim = PushForwardType::dim;
 *
 *  // Co-dimension of the functions domain.
 *  static const int codim = PushForwardType::codim;
 *
 *  // Dimension of the functions domain embedding space.
 *  static const int space_dim = PushForwardType::space_dim;
 *
 *  // Dimension of the range.
 *  static const int range = PushForwardType::template PhysRange<RefSpace::range>::value;
 *
 *  // Rank of the range.
 *  static const int rank = RefSpace::rank;
 *
 *  // Number of scalar components of the range.
 *  static constexpr int n_components = constexpr_pow(dim_range, rank);
 * @endcode
 *
 *
 * @author M.Martinelli, 2013
 */
template<class Grid_>
class FunctionSpaceOnGrid :
    public FunctionSpace,
    public GridWrapper<Grid_>
{

private:
    using self_t = FunctionSpaceOnGrid<Grid_>;

public:
    using typename GridWrapper<Grid_>::GridType;
    static constexpr std::array<Size, Grid_::dim> dims = Grid_::dims;

public:
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

    Index get_id() const
    {
        return id_;
    }

    void set_id(const Index id)
    {
        Assert(id >= 0,ExcLowerRange(id,0));
        id_ = id;
    }

protected:
    Index id_ = 0;
};


/**
 * Return true if the class T is derived from FunctionSpace.
 * @relates FunctionSpace
 *
 */
template<class T>
constexpr bool is_function_space()
{
    return std::is_base_of<FunctionSpace,T>::value;
}

IGA_NAMESPACE_CLOSE

#endif // __FUNCTION_SPACE_H_
