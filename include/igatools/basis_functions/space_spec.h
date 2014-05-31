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

#ifndef SPACE_SPEC_H_
#define SPACE_SPEC_H_

#include <igatools/base/config.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/basis_functions/function_space.h>
#include <igatools/geometry/cartesian_grid.h>
#include <algorithm>

IGA_NAMESPACE_OPEN

template <class T, int dim>
inline
std::array<T,dim>
iota_array(const int init=0)
{
    std::array<T,dim> res;
    std::iota(res.begin(), res.end(), init);
    return res;
}



/**
 * @brief Tensor product spline space specification class
 *
 * A polynomial spline space is determined by:
 * - the interval
 * - the order
 * - a partition (interior knots without repetition)
 * - the interior knots multiplicity (or smoothness a the knots)
 *
 * This is independent of the basis functions one may wish to use
 * for the given space.
 *
 * @author pauletti, 2013-2014
 *
 */
template<int dim, int range = 1, int rank = 1>
class SplineSpace :
        public FunctionSpaceOnGrid<CartesianGrid<dim> >
{

private:
    using GridType  = CartesianGrid<dim>;
    using GridSpace = FunctionSpaceOnGrid<GridType>;

public:
    template<class> class ComponentContainer;
    static const int n_components = ComponentContainer<int>::n_entries;

public:
    using Knots = CartesianProductArray<Real, dim>;
    using BoundaryKnots = std::array<CartesianProductArray<Real,2>, dim>;
    using Degrees  = TensorIndex<dim>;
    using Multiplicity = CartesianProductArray<Size, dim>;

    using DegreeTable = ComponentContainer<Degrees>;
    using MultiplicityTable = ComponentContainer<Multiplicity>;
    using BoundaryKnotsTable = ComponentContainer<BoundaryKnots>;
    using KnotsTable = ComponentContainer<Knots>;
    using PeriodicTable = ComponentContainer<std::array<bool, dim> >;

    using IndexSpaceTable = ComponentContainer<DynamicMultiArray<Index,dim>>;
    using IndexSpaceMarkTable = Multiplicity;

    class SpaceDimensionTable : public ComponentContainer<TensorSize<dim> >
    {
    public:
        ComponentContainer<Size> comp_dimension;
        Size total_dimension;
    };

    // For the boundary kntos types
    // interpolatory (open knot)
    enum class EndBehaviour
    {
        interpolatory
    };

    // For the interior multiplicities
    // maximum regularity
    // minimul regularity discontinous
    enum class InteriorReg
    {
        maximum, minimun
    };

public:
    /**
     * Most general constructor
     */
    explicit SplineSpace(const DegreeTable &deg,
                         std::shared_ptr<GridType> knots,
                         std::shared_ptr<const MultiplicityTable> interior_mult,
                         const PeriodicTable periodic = PeriodicTable(filled_array<bool,dim>(false)));

    explicit SplineSpace(const DegreeTable &deg,
                         std::shared_ptr<GridType> knots,
                         const InteriorReg interior_mult,
                         const PeriodicTable periodic = PeriodicTable(filled_array<bool,dim>(false)))
        :SplineSpace(deg, knots, fill_max_regularity(knots), periodic)
    {}

    const DegreeTable &get_degree() const
    {
        return deg_;
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
        return space_dim_.comp_dimension(comp);
    }

    /**
     *  Total number of basis functions for the comp space component
     *  and the dir direction.
     */
    Size get_num_basis(const int comp, const int dir) const
    {
        return  space_dim_(comp)[dir];
    }

    /**
     * Component-direction indexed table with the number of basis functions
     * in each direction and component
     */
    const SpaceDimensionTable &get_num_basis_table() const
    {
        return space_dim_;
    }

    const SpaceDimensionTable &get_num_basis_per_element_table() const
    {
        return elem_n_basis_;
    }
    /**
     * Returns the number of dofs per element.
     */
    Size get_num_basis_per_element() const
    {
        return elem_n_basis_.total_dimension;
    }

    /**
     *  Return the number of dofs per element for the i-th space component.
     */
    Size get_num_basis_per_element(int i) const
    {
        return elem_n_basis_.comp_dimension(i);
    }

    ///@}


    KnotsTable compute_knots_with_repetition(const BoundaryKnotsTable &boundary_knots);

    KnotsTable compute_knots_with_repetition(const EndBehaviour type)
    {
        return compute_knots_with_repetition(interpolatory_end_knots());
    }

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
    std::shared_ptr<MultiplicityTable> fill_max_regularity(std::shared_ptr<const GridType> grid);

    BoundaryKnotsTable interpolatory_end_knots();



private:
    std::shared_ptr<const MultiplicityTable> interior_mult_;

    DegreeTable deg_;

    /** Table with the dimensionality of the space in each component and direction */
    SpaceDimensionTable space_dim_;

    /** Table with the number of element non zero basis in each component and direction */
    SpaceDimensionTable elem_n_basis_;

    PeriodicTable periodic_;

public:
    /**
     *  Class to manage the component quantities with the knowledge of
     * uniform range spaces
     */
    template<class T>
    class ComponentContainer : public StaticMultiArray<T,range,rank>
    {
        using base_t = StaticMultiArray<T,range,rank>;

    public:
        using base_t::n_entries;

        using  ComponentMap = std::array <Index, n_entries>;

        ComponentContainer(const ComponentMap &comp_map =
                iota_array<Index, n_entries>());

        ComponentContainer(const T &val);

        ComponentContainer(std::initializer_list<T> list);

        /**
         *  Flat index access operator (non-const version).
         */
        T &operator()(const Index i);

        /**
         *  Flat index access operator (const version).
         */
        const T &operator()(const Index i) const;

        const Index active(const Index i) const
        {
        	return comp_map_[i];
        }

        const std::vector<Index> &get_active_components() const
        {return active_components_;}

        const std::vector<Index> &get_inactive_components() const
        {return inactive_components_;}
    private:
        /** For each component return de index of the active component */
        std::array <Index, n_entries> comp_map_;

        /** list of the active components */
        std::vector<Index> active_components_;

        /** list of the inactive components */
        std::vector<Index> inactive_components_;
    };

};


template <class T, int dim>
inline
std::vector<T>
unique_container(std::array <T, dim> a)
{
    std::unique(a.begin(), a.end());
    return std::vector<T>(a.begin(), a.end());
}



template<int dim, int range, int rank>
template<class T>
SplineSpace<dim, range, rank>::ComponentContainer<T>::
ComponentContainer(const ComponentMap &comp_map)
    :
    base_t(),
    comp_map_(comp_map),
    active_components_(unique_container<Index, n_entries>(comp_map)),
    inactive_components_(n_entries)
{
	auto all = iota_array<Index, n_entries>();
	auto it=std::set_difference (all.begin(), all.end(), active_components_.begin(),
			active_components_.end(), inactive_components_.begin());

	inactive_components_.resize(it-inactive_components_.begin());
}



template<int dim, int range, int rank>
template<class T>
SplineSpace<dim, range, rank>::ComponentContainer<T>::
ComponentContainer(std::initializer_list<T> list)
    :
    base_t(list),
    comp_map_(iota_array<Index, n_entries>()),
    active_components_(unique_container<Index, n_entries>(comp_map_))
{}



template<int dim, int range, int rank>
template<class T>
SplineSpace<dim, range, rank>::ComponentContainer<T>::
ComponentContainer(const T &val)
    :
    comp_map_(filled_array<Index, n_entries>(0)),
    active_components_(unique_container<Index, n_entries>(comp_map_))
{
    base_t::operator()(0) = val;
}



template<int dim, int range, int rank>
template<class T>
T &
SplineSpace<dim, range, rank>::ComponentContainer<T>::
operator()(const Index i)
{
    return base_t::operator()(comp_map_[i]);
}


template<int dim, int range, int rank>
template<class T>
const T &
SplineSpace<dim, range, rank>::ComponentContainer<T>::
operator()(const Index i) const
{
    return base_t::operator()(comp_map_[i]);
}

IGA_NAMESPACE_CLOSE

#endif
