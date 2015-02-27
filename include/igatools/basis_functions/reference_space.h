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

#ifndef REFERENCE_SPACE_H_
#define REFERENCE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/array_utils.h>
#include <igatools/base/function_element.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/basis_functions/function_space.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/basis_functions/space_element.h>
#include <igatools/basis_functions/spline_space.h>

IGA_NAMESPACE_OPEN

class SpaceManager;


template<Transformation,int, int> class PushForward;

template <int, int, int ,int,Transformation> class PhysicalSpace;

template <int, int, int> class ReferenceElement;
template <int,int,int> class ReferenceElementHandler;

template <int, int, int> class BSplineSpace;
template <int, int, int> class NURBSSpace;


template <int,int,int> class DofDistribution;

template<int dim_, int range_ = 1, int rank_ = 1>
class ReferenceSpace : public FunctionSpaceOnGrid<CartesianGrid<dim_>>
{
public:
    static const int dim       = dim_;
    static const int codim     = 0;
    static const int space_dim = dim_;
    static const int range     = range_;
    static const int rank      = rank_;
    static const bool is_physical_space = false;

    /**
     * See documentation in \ref FunctionSpaceOnGrid
     *
     * @see FunctionSpaceOnGrid
     */
    using PushForwardType = PushForward<Transformation::h_grad, dim, codim>;

    using RefSpace = ReferenceSpace<dim_,range_,rank_>;

    using Func = Function<dim, 0, range, rank>;

    template <int order>
    using Derivative = typename Func::template Derivative<order>;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Div   = typename Func::Div;
    using RefPoint = Point;


    using GridSpace = FunctionSpaceOnGrid<CartesianGrid<dim>>;
    using typename GridSpace::GridType;



    /** Type for the element accessor. */
    using ElementAccessor = ReferenceElement<dim,range,rank>;

    /** Type for iterator over the elements.  */
    using ElementIterator = CartesianGridIterator<ElementAccessor>;

    using ElementHandler = ReferenceElementHandler<dim_, range_, rank_>;

    using SpaceData = SplineSpace<dim_,range_,rank_>;

    using Degrees = typename SpaceData::Degrees;
    using Multiplicity = typename SpaceData::Multiplicity;
    using EndBehaviour = typename SpaceData::EndBehaviour;
    using Periodicity = typename SpaceData::Periodicity;

    using KnotsTable = typename SpaceData::KnotsTable;
    using DegreeTable = typename SpaceData::DegreeTable;
    using MultiplicityTable = typename SpaceData::MultiplicityTable;
    using SpaceDimensionTable = typename SpaceData::SpaceDimensionTable;
    using PeriodicTable = typename SpaceData::PeriodicTable;
    using EndBehaviourTable = typename SpaceData::EndBehaviourTable;

    using BCTable = typename SpaceData::BCTable;

    template <class T>
    using ComponentContainer = typename SpaceData::template ComponentContainer<T>;

    using ComponentMap = typename SpaceData::template ComponentContainer<int>::ComponentMap;

    static const auto n_components = SpaceData::n_components;

protected:

    ReferenceSpace() = delete;

    explicit ReferenceSpace(const std::shared_ptr<SpaceData> space_data);

public:
    virtual ~ReferenceSpace() = default;

    /**
     * Create and element (defined on this space) with a given flat_index
     */
    virtual std::shared_ptr<ElementAccessor>
    create_element(const Index flat_index) const = 0;

    template <int k>
    using InterGridMap = typename GridType::template InterGridMap<k>;

    template <int k>
    using InterSpaceMap = vector<Index>;

    template <int k>
    using SubRefSpace = ReferenceSpace<k, range, rank>;

    template <int k>
    using SubSpace = PhysicalSpace<k,range,rank, dim-k, Transformation::h_grad>;

    virtual bool is_bspline() const = 0;

    virtual vector<Index> get_loc_to_global(const CartesianGridElement<dim> &element) const = 0;

    virtual vector<Index> get_loc_to_patch(const CartesianGridElement<dim> &element) const = 0;

    /**
     * Returns the degree of the BSpline space for each component and for each coordinate direction.
     * \return The degree of the BSpline space for each component and for each coordinate direction.
     * The first index of the returned object is the component id, the second index is the direction id.
     */
    const DegreeTable &get_degree() const
    {
        return space_data_->get_degree();
    }

    /** @name Functions for retrieving information about the number of basis function. */
    ///@{
    SpaceDimensionTable get_num_all_element_basis() const
    {
        return space_data_->get_num_all_element_basis();
    }

    const SpaceDimensionTable &get_num_basis_table() const
    {
        return space_data_->get_num_basis_table();
    }

    Size get_num_basis() const
    {
        return space_data_->get_num_basis();
    }


    Size get_num_basis(const int comp) const
    {
        return space_data_->get_num_basis(comp);
    }

    Size get_num_basis(const int comp, const int dir) const
    {
        return space_data_->get_num_basis(comp,dir);
    }

    ComponentContainer<Size> get_basis_offset() const
    {
        return space_data_->get_basis_offset();
    }
    ///@}


    const auto &get_active_components_id() const
    {
        return space_data_->get_active_components_id();
    }

    const ComponentMap &get_components_map() const
    {
        return space_data_->get_components_map();
    }

    static std::array<Size,SpaceData::n_components> get_components()
    {
        return SpaceData::components;
    }

    const BCTable &get_boundary_conditions_table() const
    {
        return space_data_->get_boundary_conditions_table();
    }

    BCTable &get_boundary_conditions_table()
    {
        return space_data_->get_boundary_conditions_table();
    }

    /**
     * Returns a reference to the end behaviour table of the BSpline space.
     */
    virtual EndBehaviourTable &get_end_behaviour_table() = 0;

    /**
     * Returns a const reference to the end behaviour table of the BSpline space.
     */
    virtual const EndBehaviourTable &get_end_behaviour_table() const = 0;


    /** Returns the container with the global dof distribution (const version). */
    virtual const DofDistribution<dim, range, rank> &
    get_dof_distribution_global() const = 0;

    /** Returns the container with the global dof distribution (non const version). */
    virtual DofDistribution<dim, range, rank> &
    get_dof_distribution_global() = 0;

    /** Returns the container with the patch dof distribution (const version). */
    virtual const DofDistribution<dim, range, rank> &
    get_dof_distribution_patch() const = 0;


    /** Returns the container with the patch dof distribution (non const version). */
    virtual DofDistribution<dim, range, rank> &
    get_dof_distribution_patch() = 0;


    //TODO (MM, Dec 22, 2014): implement ReferenceSpace::get_space_manager()
    // instead of the two implementation in BSplineSpace and NURBSSpace
    virtual std::shared_ptr<SpaceManager> get_space_manager() = 0;
    virtual std::shared_ptr<const SpaceManager> get_space_manager() const = 0;


    /** @name Functions involving the element iterator */
    ///@{
    /**
     * Returns a element iterator to the first element of the patch
     */
    virtual ElementIterator begin() const = 0;

    /**
     * Returns a element iterator to the last element of the patch
     */
    virtual ElementIterator last() const = 0;


    /**
     * Returns a element iterator to one-pass the end of patch.
     */
    virtual ElementIterator end() const = 0;
    ///@}

    template<int k>
    std::shared_ptr< SubRefSpace<k> >
    get_ref_sub_space(const int s_id,
                      InterSpaceMap<k> &dof_map,
                      std::shared_ptr<CartesianGrid<k>> sub_grid = nullptr) const;

    template<int k>
    std::shared_ptr<SubSpace<k> >
    get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
                  std::shared_ptr<CartesianGrid<k>> sub_grid,
                  std::shared_ptr<InterGridMap<k>> elem_map) const;

    virtual void print_info(LogStream &out) const = 0;

    virtual std::shared_ptr<ElementHandler> create_elem_handler() const = 0;

protected:
    std::shared_ptr<SpaceData > space_data_;


public:
    //TODO (pauletti, Feb 26, 2015): the use of this function may be a design problem
    std::shared_ptr<SpaceData> get_space_data() const;
};

IGA_NAMESPACE_CLOSE

#endif // #ifndef REFERENCE_SPACE_H_
