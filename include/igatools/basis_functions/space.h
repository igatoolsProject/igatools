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
//#include <igatools/basis_functions/space_element.h>


IGA_NAMESPACE_OPEN

template <int> class CartesianGridElementAccessor;
template <int,int,int,int> class SpaceElement;
template <int,int,int,int> class SpaceElementHandler;



/**
 * @brief This is an auxiliary class used represent the "concept" of isogeometric function space, defined
 * over <tt>dim</tt>-dimensional parametric domain.
 * It is used as base class of Space and its purpose it is isolate the methods that depends only
 * on the topological informations (in order to do not perform unnecessary instantiations).
 *
 * The main feature of this class is that contains a space identifier that is unique to the object,
 * i.e. two different objects will have two different space id.
 *
 * @author M.Martinelli, 2013, 2014, 2015.
 *
 * @ingroup serializable
 */
template<int dim_>
class SpaceBase :
    public GridWrapper<CartesianGrid<dim_>>
{
private:
    using base_t = GridWrapper<CartesianGrid<dim_>>;
    using self_t = SpaceBase<dim_>;


public:
    virtual void get_element_dofs(
        const Index element_id,
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
    SpaceBase() = default;

    /** Construct the object from the @p grid on which the function space will be built upon. */
    SpaceBase(std::shared_ptr<CartesianGrid<dim_>> grid);

    /** Copy constructor. */
    SpaceBase(const self_t &) = default;

    /** Move constructor. */
    SpaceBase(self_t &&) = default;

public:
    /** Destructor. */
    virtual ~SpaceBase() = default;
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




template <int,int,int> class DofDistribution;


/**
 * @brief This class represent the "concept" of isogeometric function space.
 * It is used as base class of ReferenceSpace and PhysicalSpace.
 *
 * @author M.Martinelli, 2015.
 *
 * @ingroup serializable
 */
template <int dim_,int codim_,int range_,int rank_>
class Space
    : public SpaceBase<dim_>
{
private:
    using base_t = SpaceBase<dim_>;
    using self_t = Space<dim_,codim_,range_,rank_>;

protected:
    /**
     * Inheriting the constructors from the base class.
     */
    using base_t::SpaceBase;

public:

    static const int dim = dim_;

    virtual std::shared_ptr<const DofDistribution<dim_,range_,rank_> >
    get_dof_distribution() const = 0;

    /**
     * Create and element (defined on this space) with a given flat_index
     */
    virtual std::shared_ptr<SpaceElement<dim_,codim_,range_,rank_>>
            create_element(const Index flat_index) const = 0;

    virtual std::shared_ptr< SpaceElementHandler<dim_,codim_,range_,rank_> >
    get_elem_handler() const = 0;

    virtual void print_info(LogStream &out) const = 0;

    /**
     * Default destructor.
     */
    virtual ~Space() = default;


    using ElementAccessor = SpaceElement<dim_,codim_,range_,rank_>;
    using ElementIterator = CartesianGridIterator<ElementAccessor>;

    /** @name Functions involving the element iterator */
    ///@{
    /**
     * Returns a element iterator to the first element of the patch
     * with the property @p element_property.
     */
    ElementIterator begin(const std::string &element_property = ElementProperties::none) const;

    /**
     * Returns a element iterator to the last element of the patch
     * with the property @p element_property.
     */
    ElementIterator last(const std::string &element_property = ElementProperties::none) const;

    /**
     * Returns a element iterator to one-pass the end of patch
     * with the property @p element_property.
     */
    ElementIterator end(const std::string &element_property = ElementProperties::none) const;
    ///@}

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
