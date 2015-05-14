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

#ifndef __SPACE_H_
#define __SPACE_H_

#include <igatools/base/config.h>
#include <igatools/geometry/grid_wrapper.h>
#include <igatools/geometry/cartesian_grid.h>

IGA_NAMESPACE_OPEN


/**
 * @brief This class represent the "concept" of isogeometric function space.
 * It is used as base class of ReferenceSpace and PhysicalSpace.
 *
 * The main feature of this class is that contains a space identifier that is unique to the object,
 * i.e. two different objects will have two different space id.
 *
 * @author M.Martinelli, 2013, 2014, 2015.
 *
 * @ingroup serializable
 */
template<int dim_>
class Space :
    public GridWrapper<CartesianGrid<dim_>>
{

private:
	using base_t = GridWrapper<CartesianGrid<dim_>>;
    using self_t = Space<dim_>;

public:
    using GridType = CartesianGrid<dim_>;

    using GridElement = typename GridType::ElementAccessor;

    virtual void get_element_dofs(
        const GridElement &element,
        SafeSTLVector<Index> &dofs_global,
        SafeSTLVector<Index> &dofs_local_to_patch,
        SafeSTLVector<Index> &dofs_local_to_elem,
        const std::string &dofs_property = DofProperties::active) const = 0;


    /** @name Constructor and destructor. */
    ///@{
protected:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    Space() = default;

    /** Construct the object from the @p grid on which the function space will be built upon. */
    Space(std::shared_ptr<GridType> grid);

    /** Copy constructor. */
    Space(const self_t &grid) = default;

    /** Move constructor. */
    Space(self_t &&grid) = default;

public:
    /** Destructor. */
    virtual ~Space() = default;
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


private:

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif // SERIALIZATION
};


IGA_NAMESPACE_CLOSE

#endif // __SPACE_H_
