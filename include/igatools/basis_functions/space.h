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
#include <igatools/utils/shared_ptr_constness_handler.h>
#include <igatools/geometry/cartesian_grid.h>
//#include <igatools/base/function.h>

//#include <igatools/basis_functions/space_element.h>



IGA_NAMESPACE_OPEN

template <int,int,int,int> class Function;
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
class SpaceBase
{
private:
//    using base_t = GridWrapper<dim_>;
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

    /** Construct the object from the (const) @p grid on which the function space will be built upon. */
    SpaceBase(const std::shared_ptr<const CartesianGrid<dim_>> &grid);

    /** Construct the object from the (non-const) @p grid on which the function space will be built upon. */
    SpaceBase(const std::shared_ptr<CartesianGrid<dim_>> &grid);

    /** Copy constructor. */
    SpaceBase(const self_t &) = delete;

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
     * Returns the unique identifier associated to each object instance.
     */
    Index get_object_id() const;


    std::shared_ptr<CartesianGrid<dim_>> get_ptr_grid();

    std::shared_ptr<const CartesianGrid<dim_>> get_ptr_const_grid() const;

    /**
     * Get the name associated to the object instance.
     */
    const std::string &get_name() const;

    /**
     * Set the name associated to the object instance.
     */
    void set_name(const std::string &name);



#ifdef MESH_REFINEMENT

    /**
     * Perform the h-refinement of the space in all the directions.
     *
     * Each interval in the unrefined grid is uniformly divided in @p n_subdivisions
     * sub-intervals.
     *
     * @ingroup h_refinement
     */
    void refine_h(const Size n_subdivisions = 2);

#endif // MESH_REFINEMENT

protected:

    /**
     * Unique identifier associated to each object instance.
     */
    Index object_id_ = 0;

    /**
     * Name associated to the object instance.
     */
    std::string name_;

private:

//    std::shared_ptr<CartesianGrid<dim_> > grid_;

    SharedPtrConstnessHandler<CartesianGrid<dim_> > grid_;

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
    :
    public std::enable_shared_from_this<Space<dim_,codim_,range_,rank_> >,
    public SpaceBase<dim_>
{
private:
    using base_t = SpaceBase<dim_>;
    using self_t = Space<dim_,codim_,range_,rank_>;


protected:
    using MapFunc = Function<dim_,0,dim_+codim_,1>;


    /** @name Constructor and destructor. */
    ///@{
protected:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    Space() = default;

    /**
     * Construct the object from the (non-const) @p grid on which the function space will be built upon
     * and the function representing the mapping.
     *
     * @pre The shared_pointer <tt>map_func</tt> must be unique.
     *
     * @warning After the object construction the state of <tt>map_func</tt> will be no longer valid.
     */
    Space(const std::shared_ptr<CartesianGrid<dim_>> &grid,const std::shared_ptr<MapFunc> &map_func);

    /**
     * Construct the object from the (const) @p grid on which the function space will be built upon
     * and the function representing the mapping.
     *
     * @pre The shared_pointer <tt>map_func</tt> must be unique.
     *
     * @warning After the object construction the state of <tt>map_func</tt> will be no longer valid.
     */
    Space(const std::shared_ptr<const CartesianGrid<dim_>> &grid,const std::shared_ptr<MapFunc> &map_func);

    /** Copy constructor. */
    Space(const self_t &) = delete;

    /** Move constructor. */
    Space(self_t &&) = default;

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


    static const int dim = dim_;
    static const int codim = codim_;
    static const int range = range_;
    static const int rank = rank_;



    std::shared_ptr<MapFunc> get_ptr_map_func()
    {
        return map_func_.get_ptr_data();
    }

    std::shared_ptr<const MapFunc> get_ptr_const_map_func() const
    {
        return map_func_.get_ptr_const_data();
    }

    virtual std::shared_ptr<const DofDistribution<dim_,range_,rank_> >
    get_ptr_const_dof_distribution() const = 0;


    virtual std::shared_ptr<DofDistribution<dim_,range_,rank_> >
    get_ptr_dof_distribution() = 0;


    /**
     * Returns the dofs that are on the interior of the <tt>dim</tt>-dimensional hypercube
     * (i.e. the dofs that are not on the boundary).
     */
    std::set<Index> get_interior_dofs() const;

    using topology_variant = TopologyVariants<dim_>;

    std::set<Index> get_boundary_dofs(const int s_id, const topology_variant &topology) const;

    template<int k>
    std::set<Index> get_boundary_dofs(const int s_id) const
    {
        return this->get_boundary_dofs(s_id,Topology<k>());
    }

    /** @name Functions for retrieving information about the number of basis function. */
    ///@{
    Size get_num_basis() const;

    Size get_num_basis(const int comp) const;

    Size get_num_basis(const int comp, const int dir) const;

    Size get_elem_num_basis() const;

    /**
     * This function returns the global dof id corresponding to the basis function
     * with tensor index <p>tensor_index</p> on the @p comp component of the space.
     */
    Index
    get_global_dof_id(const TensorIndex<dim> &tensor_index,
                      const Index comp) const;
    ///@}


    /**
     * Return the maximum value of the polynomial degree, for each component, for each direction;
     */
    virtual int get_max_degree() const = 0;


    /**
     * Create and element (defined on this space) with a given flat_index
     */
    virtual std::shared_ptr<SpaceElement<dim_,codim_,range_,rank_> >
    create_element(const Index flat_index) const = 0;


    virtual std::shared_ptr< SpaceElementHandler<dim_,codim_,range_,rank_> >
    get_elem_handler() const = 0;


    virtual void print_info(LogStream &out) const = 0;



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


#ifdef MESH_REFINEMENT

    virtual std::shared_ptr<const self_t> get_space_previous_refinement() const = 0;

#endif


private:

    SharedPtrConstnessHandler<MapFunc>  map_func_;

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
