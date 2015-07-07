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

#ifndef __FUNCTION_H_
#define __FUNCTION_H_

#include <igatools/base/config.h>
#include <igatools/base/tensor.h>
#include <igatools/utils/value_vector.h>
#include <igatools/geometry/cartesian_grid_iterator.h>
#include <igatools/geometry/grid_element_handler.h>
#include <igatools/base/quadrature.h>

IGA_NAMESPACE_OPEN



template <int, int, int, int> class FunctionElement;

/**
 * Function Class
 *
 * @ingroup serializable
 */
template<int dim_, int codim_ = 0, int range_ = 1, int rank_ = 1>
class Function
    : public GridElementHandler<dim_>,
      public std::enable_shared_from_this<Function<dim_,codim_,range_,rank_> >
{
private:
    using base_t = Function<dim_, codim_, range_, rank_>;
    using self_t = Function<dim_, codim_, range_, rank_>;
    using parent_t = GridElementHandler<dim_>;


public:
    using GridType = const CartesianGrid<dim_>;

    using topology_variant = TopologyVariants<dim_>;
    using eval_pts_variant = SubElemVariants<Quadrature,dim_>;

    using ElementAccessor = FunctionElement<dim_, codim_, range_, rank_>;
    using ElementIterator = CartesianGridIterator<ElementAccessor>;

    static const int space_dim = dim_ + codim_;
    static const int dim       = dim_;
    static const int codim     = codim_;
    static const int range     = range_;
    static const int rank      = rank_;

    /** Types for the input/output evaluation arguments */
    ///@{
    using RefPoint = Points<dim>;

    /**
     * Type for the input argument of the function.
     */
    using Point = Points<space_dim>;

    /**
     * Type for the return of the function.
     */
    using Value = Values<space_dim, range_, rank_>;

    /**
     * Type for the derivative of the function.
     */
    template <int order>
    using Derivative = Derivatives<space_dim, range_, rank_, order>;

    /**
     * Type for the gradient of the function.
     */
    using Gradient = Derivative<1>;

    /**
     * Type for the hessian of the function.
     */
    using Hessian = Derivative<2>;

    /**
     * Type for the divergence of function.
     */
    using Div = Values<space_dim, range_, rank_-1>;
    ///@}

    /** @name Constructors and destructor. */
    ///@{
protected:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    Function() = default;

    /** Constructor */
    Function(std::shared_ptr<GridType> grid);

    /**
     * Copy constructor.
     */
    Function(const self_t &func);

public:
    /** Destructor */
    virtual ~Function() = default;
    ///@}


    virtual std::shared_ptr<base_t> clone() const = 0;

    virtual void reset(const ValueFlags &flag, const eval_pts_variant &quad);

    virtual void init_cache(ElementAccessor &elem, const topology_variant &k);

    void init_cache(ElementIterator &elem, const topology_variant &k);

    void init_element_cache(ElementAccessor &elem);

    void init_element_cache(ElementIterator &elem);

    virtual void fill_cache(ElementAccessor &elem, const topology_variant &k,const int j);

    void fill_cache(ElementIterator &elem, const topology_variant &k, const int j);

    void fill_element_cache(ElementAccessor &elem);


    void fill_element_cache(ElementIterator &elem);

    std::shared_ptr<ElementAccessor> create_element(const Index flat_index) const;

    ElementIterator begin() const;

    ElementIterator end() const;



    virtual void print_info(LogStream &out) const;

    /**
     * Returns the unique identifier associated to each object instance.
     */
    Index get_object_id() const;

    /**
     * Get the name associated to the object instance.
     */
    const std::string &get_name() const;

    /**
     * Set the name associated to the object instance.
     */
    void set_name(const std::string &name);

private:

    /**
     * Unique identifier associated to each object instance.
     */
    Index object_id_;

    /**
     * Name associated to the object instance.
     */
    std::string name_;


    struct ResetDispatcher : boost::static_visitor<void>
    {
        ResetDispatcher(const ValueFlags flag_in,
                        GridElementHandler<dim_> &grid_elem_handler,
                        SafeSTLArray<ValueFlags, dim_ + 1> &flags)
            :
            flag_in_(flag_in),
            grid_elem_handler_(grid_elem_handler),
            flags_(flags)
        {}

        template<int sub_elem_dim>
        void operator()(const Quadrature<sub_elem_dim> &quad)
        {
            flags_[sub_elem_dim] = flag_in_;

            grid_elem_handler_.template reset<sub_elem_dim>(flag_in_, quad);
        }

        const ValueFlags flag_in_;
        GridElementHandler<dim_> &grid_elem_handler_;
        SafeSTLArray<ValueFlags, dim_ + 1> &flags_;
    };

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        FillCacheDispatcher(const int sub_elem_id,
                            const GridElementHandler<dim_> &grid_elem_handler,
                            ElementAccessor &func_elem)
            :
            sub_elem_id_(sub_elem_id),
            grid_elem_handler_(grid_elem_handler),
            func_elem_(func_elem)
        {}

        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem)
        {
            grid_elem_handler_.template fill_cache<sub_elem_dim>(func_elem_, sub_elem_id_);
        }

        int sub_elem_id_;
        const GridElementHandler<dim_> &grid_elem_handler_;
        ElementAccessor &func_elem_;
    };

    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        InitCacheDispatcher(const GridElementHandler<dim_> &grid_elem_handler,
                            const SafeSTLArray<ValueFlags, dim_ + 1> &flags,
                            ElementAccessor &func_elem)
            :
            grid_elem_handler_(grid_elem_handler),
            flags_(flags),
            func_elem_(func_elem)
        {}


        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem)
        {
            grid_elem_handler_.template init_cache<sub_elem_dim>(func_elem_);

            auto &cache = func_elem_.all_sub_elems_cache_;
            if (cache == nullptr)
            {
                using Cache = typename ElementAccessor::CacheType;
                cache = std::make_shared<Cache>();
            }

            const auto n_pts = grid_elem_handler_.template get_num_points<sub_elem_dim>();
            for (const auto s_id: UnitElement<dim_>::template elems_ids<sub_elem_dim>())
            {
                auto &s_cache = cache->template get_sub_elem_cache<sub_elem_dim>(s_id);
                s_cache.resize(flags_[sub_elem_dim],n_pts);
            }
        }

        const GridElementHandler<dim_> &grid_elem_handler_;
        const SafeSTLArray<ValueFlags, dim_ + 1> &flags_;
        ElementAccessor &func_elem_;
    };



protected:
    std::shared_ptr<typename ElementAccessor::CacheType>
    &get_cache(ElementAccessor &elem);



    /**
     * One flag for each possile subdim
     */
    SafeSTLArray<ValueFlags, dim_ + 1> flags_;


    std::shared_ptr<CartesianGrid<dim_> > grid_;

    std::shared_ptr<self_t> function_previous_refinement_;

#ifdef MESH_REFINEMENT

public:



    const std::shared_ptr<self_t> &get_function_previous_refinement() const
    {
//        std::cout << "Function::get_function_previous_refinement()    Counter = " << function_previous_refinement_.use_count() << std::endl;

        return function_previous_refinement_;
    }


#endif // MESH_REFINEMENT


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
    /*
    {
        ar &boost::serialization::make_nvp("grid_elem_handler_",
                                           boost::serialization::base_object<GridElementHandler<dim_>>(*this));

        ar &boost::serialization::make_nvp("flags_",flags_);
    }
    //*/
    ///@}
#endif // SERIALIZATION
};


template<int dim, int space_dim>
using MapFunction = Function<dim, 0, space_dim, 1>;

IGA_NAMESPACE_CLOSE




#endif
