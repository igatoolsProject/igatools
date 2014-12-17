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
#include <igatools/base/function.h>
#include <igatools/base/function_element.h>
#include <igatools/utils/cartesian_product_array.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/basis_functions/function_space.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/basis_functions/space_element.h>

IGA_NAMESPACE_OPEN



template <class T, int dim>
inline
vector<T>
unique_container(std::array <T, dim> a)
{
    auto it = std::unique(a.begin(), a.end());
    return vector<T>(a.begin(), it);
}


template <class,int,Transformation> class PhysicalSpace;

template <int, int, int> class ReferenceElement;

template <int,int,int> class ReferenceElementHandler;

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

    using Func = Function<dim, 0, range, rank>;

    template <int order>
    using Derivative = typename Func::template Derivative<order>;
    using Point = typename Func::Point;
    using Value = typename Func::Value;
    using Div   = typename Func::Div;
    using RefPoint = Point;


    using GridSpace = FunctionSpaceOnGrid<CartesianGrid<dim>>;
    using typename GridSpace::GridType;

    using GridSpace::dims;

    using GridSpace::GridSpace;

    /** Type for the element accessor. */
    using ElementAccessor = ReferenceElement<dim,range,rank>;

    /** Type for iterator over the elements.  */
    using ElementIterator = CartesianGridIterator<ElementAccessor>;

    using ElementHandler = ReferenceElementHandler<dim_, range_, rank_>;


    virtual ~ReferenceSpace() = default;

    /**
     *  Class to manage the component quantities with the knowledge of
     * uniform range spaces
     */
    template<class T>
    class ComponentContainer : public StaticMultiArray<T,range,rank>
    {
        using base_t = StaticMultiArray<T,range,rank>;
    public:
        /** Type of the iterator. */
        using iterator =  MultiArrayIterator<ComponentContainer<T>>;

        /** Type of the const iterator. */
        using const_iterator =  MultiArrayConstIterator<ComponentContainer<T>>;
    public:
        using base_t::n_entries;

        using  ComponentMap = std::array <Index, n_entries>;

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
        std::array <Index, n_entries> comp_map_;

        /** list of the active components */
        vector<Index> active_components_id_;

        /** list of the inactive components */
        vector<Index> inactive_components_id_;


    };


    /**
     * Component holding the number of basis functions
     */
    class SpaceDimensionTable : public ComponentContainer<TensorSize<dim> >
    {
        using base_t = ComponentContainer<TensorSize<dim>>;
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

    static constexpr int n_components = ComponentContainer<Size>::n_entries;
    static const std::array<Size, n_components> components;


    using Degrees  = TensorIndex<dim>;
    using DegreeTable = ComponentContainer<Degrees>;

    /**
     * Type alias for the boundary conditions on each face of each scalar component of the space.
     */
    using BCTable = ComponentContainer<std::array<BoundaryConditionType,UnitElement<dim>::n_faces>>;



    template <int k>
    using InterGridMap = typename GridType::template InterGridMap<k>;

    template <int k>
    using InterSpaceMap = vector<Index>;

    template <int k>
    using SubRefSpace = ReferenceSpace<k, range, rank>;

    template <int k>
    using SubSpace = PhysicalSpace<ReferenceSpace<k>, dim-k, Transformation::h_grad>;



    virtual bool is_bspline() const = 0;

    virtual vector<Index> get_loc_to_global(const CartesianGridElement<dim> &element) const = 0;
    virtual vector<Index> get_loc_to_patch(const CartesianGridElement<dim> &element) const = 0;

    virtual SpaceDimensionTable get_num_all_element_basis() const = 0 ;

    virtual const SpaceDimensionTable &get_num_basis_table() const = 0;

    virtual Size get_num_basis() const = 0;

    virtual Size get_num_basis(const int comp) const = 0;
    virtual Size get_num_basis(const int comp, const int dir) const = 0;

    virtual const std::array<Index,n_components> &get_components_map() const = 0;

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


    virtual void print_info(LogStream &out) const = 0;


    template<int k>
    std::shared_ptr<SubRefSpace<k> >
    get_ref_sub_space(const int sub_elem_id,
                      InterSpaceMap<k> &dof_map,
                      std::shared_ptr<CartesianGrid<k>> sub_grid = nullptr) const
    {
        Assert(false,ExcNotImplemented());
        return nullptr;
    }

    template<int k>
    std::shared_ptr<SubSpace<k> >
    get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
                  std::shared_ptr<CartesianGrid<k>> sub_grid,
                  std::shared_ptr<typename GridType::template InterGridMap<k>> elem_map) const
    {
        Assert(false,ExcNotImplemented());
        return nullptr;
    }


    virtual std::shared_ptr<ElementHandler> create_elem_handler() const = 0;
};







IGA_NAMESPACE_CLOSE

#endif // #ifndef REFERENCE_SPACE_H_
