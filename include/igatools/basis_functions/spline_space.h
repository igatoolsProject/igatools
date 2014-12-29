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

#ifndef SPLINE_SPACE_H_
#define SPLINE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/geometry/grid_wrapper.h>
#include <igatools/base/function.h>
#include <igatools/base/function_element.h>


IGA_NAMESPACE_OPEN

template <int,int,int> class DofDistribution;

template <class T, int dim>
inline
vector<T>
unique_container(std::array <T, dim> a)
{
    auto it = std::unique(a.begin(), a.end());
    return vector<T>(a.begin(), it);
}

enum class EndBehaviour
{
    /**
     * Interpolatory basis functions at knots bounday (i.e. open knot vector).
     */
    interpolatory,

    periodic,

    end_knots
};

/**
 * @brief Tensor product spline space
 *
 * A one dimensional polynomial spline space is determined by:
 * - a domain, the interval [a,b]
 * - the polynomial order
 * - a partition of [a,b], the knots
 * - the interior knots smoothness (given by the interior multiplicity)
 *
 * This class provides the realization of a scalar, vector or k-tensor
 * value tensor product spline space.
 *
 * @note This is independent of the basis functions one may wish to use
 * for the given space.
 *
 * @author pauletti, 2014
 *
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class SplineSpace : public GridWrapper<CartesianGrid<dim_> >
{
public:
	using self_t = SplineSpace<dim_,range_,rank_>;


    /**
     *  Class to manage the component quantities with the knowledge of
     * uniform range spaces
     */
    template<class T>
    class ComponentContainer : public StaticMultiArray<T,range_,rank_>
    {
        using base_t = StaticMultiArray<T,range_,rank_>;
    public:
        /** Type of the iterator. */
        using iterator =  MultiArrayIterator<ComponentContainer<T>>;

        /** Type of the const iterator. */
        using const_iterator =  MultiArrayConstIterator<ComponentContainer<T>>;
    public:
        using base_t::n_entries;

        using ComponentMap = std::array <Index, n_entries>;

        ComponentContainer(const ComponentMap &comp_map =
                               sequence<n_entries>())
            :
            base_t(),
            comp_map_(comp_map),
            active_components_id_(unique_container<Index, n_entries>(comp_map)),
            inactive_components_id_(n_entries)
        {
            auto all = sequence<n_entries>();
            auto it=std::set_difference(all.begin(), all.end(),
                                        active_components_id_.begin(),active_components_id_.end(),
                                        inactive_components_id_.begin());

            inactive_components_id_.resize(it-inactive_components_id_.begin());
        }


        ComponentContainer(const ComponentMap &comp_map, const T &val)
            :
            base_t(),
            comp_map_(comp_map),
            active_components_id_(unique_container<Index, n_entries>(comp_map)),
            inactive_components_id_(n_entries)
        {
            auto all = sequence<n_entries>();
            auto it=std::set_difference(all.begin(), all.end(),
                                        active_components_id_.begin(),active_components_id_.end(),
                                        inactive_components_id_.begin());

            inactive_components_id_.resize(it-inactive_components_id_.begin());

            for (auto i : active_components_id_)
                base_t::operator[](i) = val;
        }


        /**
         * Construct a homogenous range table with val value
         */
        ComponentContainer(const T &val)
            :
            comp_map_(filled_array<Index, n_entries>(0)),
            active_components_id_(1,0),
            inactive_components_id_(n_entries-1)
        {
            for (int i=1; i<n_entries; ++i)
                inactive_components_id_[i-1] = i;

            base_t::operator[](0) = val;
        }


        ComponentContainer(std::initializer_list<T> list)
            :
            base_t(list),
            comp_map_(sequence<n_entries>()),
            active_components_id_(unique_container<Index, n_entries>(comp_map_))
        {};


        const_iterator
        cbegin() const
        {
            return const_iterator(*this,0);
        }

        const_iterator
        cend() const
        {
            return const_iterator(*this,IteratorState::pass_the_end);
        }


        const_iterator
        begin() const
        {
            return cbegin();
        }

        const_iterator
        end() const
        {
            return cend();
        }

        iterator
        begin()
        {
            return iterator(*this,0);
        }

        iterator
        end()
        {
            return iterator(*this,IteratorState::pass_the_end);
        }

        /**
         *  Flat index access operator (non-const version).
         */
        T &operator[](const Index i)
        {
            return base_t::operator[](comp_map_[i]);
        }

        /**
         *  Flat index access operator (const version).
         */
        const T &operator[](const Index i) const
        {
            return base_t::operator[](comp_map_[i]);
        }

        const Index active(const Index i) const
        {
            return comp_map_[i];
        }

        const vector<Index> &get_active_components_id() const
        {
            return active_components_id_;
        }

        const vector<Index> &get_inactive_components_id() const
        {
            return inactive_components_id_;
        }

        void
        print_info(LogStream &out) const
        {
            out.begin_item("Raw componets: ");
            base_t::print_info(out);
            out.end_item();
            out.begin_item("Active componets ids: ");
            active_components_id_.print_info(out);
            out.end_item();
            out.begin_item("Inactive componets ids: ");
            inactive_components_id_.print_info(out);
            out.end_item();
        }

        const std::array <Index, n_entries> &get_comp_map() const
        {
            return comp_map_;
        }

    private:
        /** For each component return the index of the active component */
        ComponentMap comp_map_;

        /** list of the active components */
        vector<Index> active_components_id_;

        /** list of the inactive components */
        vector<Index> inactive_components_id_;
    };


    /**
     * Component holding the number of basis functions
     */
    class SpaceDimensionTable : public ComponentContainer<TensorSize<dim_> >
    {
        using base_t = ComponentContainer<TensorSize<dim_>>;
    public:
        //using base_t::ComponentContainer;

        SpaceDimensionTable() = default;

        SpaceDimensionTable(const base_t &n_basis)
            :
            base_t(n_basis),
            comp_dimension(n_basis.get_comp_map()),
            total_dimension(0)
        {
            for (auto comp : this->get_active_components_id())
            {
                auto size = (*this)[comp].flat_size();
                comp_dimension[comp] = size;
            }
            for (auto size : comp_dimension)
                total_dimension += size;
        }

        //TODO(pauletti, Sep 8, 2014): make this private and write some getters
        ComponentContainer<Size> comp_dimension;
        Size total_dimension;
    };

private:
    using GridType = CartesianGrid<dim_>;
//    using typename GridType::GridType;


public:
//    using GridSpace::dims;

    using Func = Function<dim_, 0, range_, rank_>;

public:
    template <int order>
    using Derivative = typename Func::template Derivative<order>;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Div   = typename Func::Div;

public:
//    using RefSpace::n_components;
//    static constexpr int n_components = RefSpace::template ComponentContainer<Size>::n_entries;
//    static const std::array<Size, n_components> components;


public:
    using Degrees  = TensorIndex<dim_>;
    using DegreeTable = ComponentContainer<Degrees>;

    using KnotCoordinates = typename GridType::KnotCoordinates;
    using BoundaryKnots = std::array<CartesianProductArray<Real,2>, dim_>;
    using Multiplicity = CartesianProductArray<Size, dim_>;

    using MultiplicityTable = ComponentContainer<Multiplicity>;
    using BoundaryKnotsTable = ComponentContainer<BoundaryKnots>;
    using KnotsTable = ComponentContainer<KnotCoordinates>;

    using IndexSpaceTable = ComponentContainer<DynamicMultiArray<Index,dim_>>;
    using IndexSpaceMarkTable = Multiplicity;

    static constexpr int n_components = ComponentContainer<Size>::n_entries;

    using ComponentMap = typename ComponentContainer<Size>::ComponentMap;

    /**
     * Type alias for the boundary conditions on each face of each scalar component of the space.
     */
    using BCTable = ComponentContainer<std::array<BoundaryConditionType,UnitElement<dim_>::n_faces>>;

public:
    static const std::array<Size, n_components> components;


public:

    using EndBehaviourTable = ComponentContainer<std::array<EndBehaviour, dim_> >;

    // For the interior multiplicities
    // maximum regularity
    // minimul regularity discontinous
    enum class InteriorReg
    {
        maximum, minimun
    };

protected:

    SplineSpace() = delete;

    /**
     * Most general constructor
     */
    explicit SplineSpace(const DegreeTable &deg,
                         std::shared_ptr<GridType> knots,
                         std::shared_ptr<const MultiplicityTable> interior_mult,
                         const EndBehaviourTable &end_behaviour = EndBehaviourTable(filled_array<EndBehaviour,dim_>(EndBehaviour::interpolatory)));


    explicit SplineSpace(const DegreeTable &deg,
                         std::shared_ptr<GridType> knots,
                         const InteriorReg &interior_mult,
                         const EndBehaviourTable &ebt = EndBehaviourTable(filled_array<EndBehaviour,dim_>(EndBehaviour::interpolatory)))
        :SplineSpace(deg, knots, fill_max_regularity(deg, knots), ebt)
    {}

public:
    static std::shared_ptr<self_t> create(const DegreeTable &deg,
                         std::shared_ptr<GridType> knots,
                         const InteriorReg &interior_mult,
                         const EndBehaviourTable &ebt = EndBehaviourTable(filled_array<EndBehaviour,dim_>(EndBehaviour::interpolatory)));

    static std::shared_ptr<self_t> create(
    		const DegreeTable &deg,
   	        std::shared_ptr<GridType> knots,
    	    std::shared_ptr<const MultiplicityTable> interior_mult,
    	    const EndBehaviourTable &end_behaviour = EndBehaviourTable(filled_array<EndBehaviour,dim_>(EndBehaviour::interpolatory)));

    virtual ~SplineSpace() = default;

    const DegreeTable &get_degree() const
    {
        return deg_;
    }


    const ComponentMap &get_components_map() const
    {
        return interior_mult_->get_comp_map();
    }

    /** @name Getting information about the space */
    ///@{
    /**
     * Total number of basis functions. This is the dimensionality
     * of the space.
     */
    Size get_num_basis() const
    {
        return space_dim_.total_dimension;
    }

    /**
     * Total number of basis functions
     * for the comp space component.
     */
    Size get_num_basis(const int comp) const
    {
        return space_dim_.comp_dimension[comp];
    }

    /**
     *  Total number of basis functions for the comp space component
     *  and the dir direction.
     */
    Size get_num_basis(const int comp, const int dir) const
    {
        return  space_dim_[comp][dir];
    }

    /**
     * Component-direction indexed table with the number of basis functions
     * in each direction and component
     */
    const SpaceDimensionTable &get_num_basis_table() const
    {
        return space_dim_;
    }

    //TODO (MM, Dec 17,2014): this function is misleading because refers to an element and not to the space. REMOVE IT!
    /**
     * Component table with the offset of basis functions
     * in each component of an element.
     */
    SpaceDimensionTable get_num_all_element_basis() const
    {
        ComponentContainer<TensorSize<dim_>> n_basis(deg_.get_comp_map());
        for (auto comp : deg_.get_active_components_id())
            n_basis[comp] = TensorSize<dim_>(deg_[comp]+1);

        return SpaceDimensionTable(n_basis);
    }

    /**
     * Component table with the offset of basis functions
     * in each component of the space.
     */
    ComponentContainer<Size> get_basis_offset() const
    {
        ComponentContainer<Size> offset;
        offset[0] = 0;
        for (int comp = 1; comp < n_components; ++comp)
            offset[comp] = offset[comp-1] + space_dim_.comp_dimension[comp];

        return offset;
    }

    ///@}

    template<int k>
    using SubSpace = SplineSpace<k, range_, rank_>;

    template<int k>
    std::shared_ptr<typename SubSpace<k>::MultiplicityTable>
    get_sub_space_mult(const Index s_id) const;

    template<int k>
    typename SubSpace<k>::DegreeTable
    get_sub_space_degree(const Index s_id) const;

    template<int k>
    typename SubSpace<k>::EndBehaviourTable
    get_sub_space_end_b(const Index s_id) const;


public:

    KnotsTable compute_knots_with_repetition(const BoundaryKnotsTable &boundary_knots) const;

    KnotsTable compute_knots_with_repetition(const EndBehaviourTable &ends) const;

    /**
     * For each element and for each component there is an initial
     * tensor index in the Index space from where all non-zero basis
     * function can be determined.
     */
    MultiplicityTable accumulated_interior_multiplicities() const;

    void print_info(LogStream &out) const;


private:
    /**
     * Fill the multiplicy for the maximum possible regularity
     *  of the given number of knots
     */
    std::shared_ptr<MultiplicityTable> fill_max_regularity(const DegreeTable &deg, std::shared_ptr<const GridType> grid);

#if 0
    BoundaryKnotsTable interpolatory_end_knots() const;
#endif

    CartesianProductArray<Real,2> interpolatory_end_knots(const int comp_id,const int dir) const;


private:

    std::shared_ptr<const MultiplicityTable> interior_mult_;

    DegreeTable deg_;

    /** Table with the dimensionality of the space in each component and direction */
    SpaceDimensionTable space_dim_;

    EndBehaviourTable end_behaviour_;

    /**
     * Boundary conditions on each face of each scalar component of the space.
     */
    BCTable boundary_conditions_table_;

public:

    /** Returns the multiplicity of the internal knots that defines the space. */
    std::shared_ptr<const MultiplicityTable> get_interior_mult() const
    {
        return interior_mult_;
    }

    const EndBehaviourTable &get_end_behaviour() const
    {
        return end_behaviour_;
    }


    // TODO (pauletti, Dec 12, 2014): boundary condition is not a general property
    // of the space, rather specific flag for some application, this should be
    // done in some other layer

    /**
     * Returns a const-reference to the table containing
     * the boundary conditions on each face of each scalar component of the space.
     *
     * For example, with the code
     * @code{.cpp}
       const auto &bc_table = space.get_boundary_conditions_table();

       BoundaryConditionType bc_id = bc_table[1][3]; // boundary condition on face 3 of space's component 1
       @endcode
     * we copy to the variable <tt>bc_id</tt> the value of the boundary condition
     * on the face 3 of the space component 1.
     *
     * @sa BoundaryConditionType
     */
    const BCTable &get_boundary_conditions_table() const
    {
        return boundary_conditions_table_;
    }

    /**
     * Returns a reference to the table containing
     * the boundary conditions on each face of each scalar component of the space.
     *
     * For example, with the code
     * @code{.cpp}
       auto &bc_table = space.get_boundary_conditions_table();

       bc_table[1][3] = BoundaryConditionType::DirichletHomogeneous; // setting Dirichlet homogeneous boundary condition on face 3 of space's component 1
       @endcode
     * we assign the value <tt>BoundaryConditionType::DirichletHomogeneous</tt> to the
     * boundary condition on the face 3 of the space component 1.
     *
     * @sa BoundaryConditionType
     */
    BCTable &get_boundary_conditions_table()
    {
        return boundary_conditions_table_;
    }



    /**
     * Refines the function space after a grid uniform refinement.
     *
     * @param[in] refinement_directions Directions along which the refinement is performed.
     * @param[in] grid_old Grid before the refinement.
     *
     * @pre Before invoking this function, must be invoked the function grid_->refine().
     * @note This function is connected to the CartesianGrid's signal for the refinement, and
     * it is necessary in order to avoid infinite loops in the refine() function calls.
     *
     * @ingroup h_refinement
     */
    void refine_h_after_grid_refinement(
        const std::array<bool,dim_> &refinement_directions,
        const GridType &grid_old) ;

    void create_connection_for_h_refinement(std::shared_ptr<self_t> space);

    std::shared_ptr<const SplineSpace<dim_,range_,rank_> > spline_space_previous_refinement_;

public:
    std::shared_ptr<const SplineSpace<dim_,range_,rank_> >
    get_spline_space_previous_refinement() const
    {
        return spline_space_previous_refinement_;
    }

protected:

    /** This function initialize the member variables from the constructor
     * arguments or after an h-refinement. */
    void init();

public:

    vector<Index> get_loc_to_global(const CartesianGridElement<dim_> &element) const
    {
        Assert(false,ExcMessage("This class should not have this function."))
        return vector<Index>();
    }

    vector<Index> get_loc_to_patch(const CartesianGridElement<dim_> &element) const
    {
        Assert(false,ExcMessage("This class should not have this function."))
        return vector<Index>();
    }


    /** Returns the container with the global dof distribution (const version). */
    const DofDistribution<dim_, range_, rank_> &
    get_dof_distribution_global() const
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return *reinterpret_cast<const DofDistribution<dim_,range_,rank_> *>(this);
    }

    /** Returns the container with the global dof distribution (non const version). */
    DofDistribution<dim_, range_, rank_> &
    get_dof_distribution_global()
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return *reinterpret_cast<DofDistribution<dim_,range_,rank_> *>(this);
    }

    /** Returns the container with the patch dof distribution (const version). */
    const DofDistribution<dim_, range_, rank_> &
    get_dof_distribution_patch() const
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return *reinterpret_cast<const DofDistribution<dim_,range_,rank_> *>(this);
    }


    /** Returns the container with the patch dof distribution (non const version). */
    DofDistribution<dim_, range_, rank_> &
    get_dof_distribution_patch()
    {
        Assert(false,ExcMessage("This class should not have this function."));
        return *reinterpret_cast<DofDistribution<dim_,range_,rank_> *>(this);
    }

};









IGA_NAMESPACE_CLOSE

#endif
