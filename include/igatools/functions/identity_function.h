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

#ifndef IDENTITY_FUNCTIONS_H
#define IDENTITY_FUNCTIONS_H

#include <igatools/functions/function.h>
#include <igatools/base/value_types.h>

IGA_NAMESPACE_OPEN

template<int dim,int space_dim>
auto
create_id_tensor()
{
    typename Function<dim, 0, space_dim, 1>::Gradient res;
    for (int i=0; i<dim; ++i)
        res[i][i] = 1.;
    return res;
}


/**
 * The identity function from R^dim to R^spacedim,
 * if dim < space_dim, the last space_dim-dim  coordinates are zero.
 *
 * @note this function is not inherited from formula function because
 * we want to optimize its computation
 *
 *
 * @ingroup serializable
 *
 * @author martinelli 2015
 * @author pauletti 2015
 */
template<int dim, int space_dim = dim>
class IdentityFunction :
    public Function<dim, 0, space_dim, 1>
{
private:
    using parent_t = Function<dim, 0, space_dim, 1>;
    using self_t = IdentityFunction<dim,space_dim>;

protected:
    using typename parent_t::GridType;

public:
    using typename parent_t::topology_variant;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::Hessian;

    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;



    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    IdentityFunction(std::shared_ptr<GridType> grid);

    virtual ~IdentityFunction() = default;

    static std::shared_ptr<parent_t>
    create(std::shared_ptr<GridType> grid);

    std::shared_ptr<parent_t> clone() const override final;


    void fill_cache(ElementAccessor &elem, const topology_variant &k,
                    const int j) override final;

    virtual void print_info(LogStream &out) const override final;

private:
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    IdentityFunction() = default;

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        FillCacheDispatcher(const int sub_elem_id,self_t &function,ElementAccessor &elem)
            :
            sub_elem_id_(sub_elem_id),
            function_(function),
            elem_(elem)
        {}

        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem)
        {
            auto &local_cache = function_.get_cache(elem_);
            auto &cache = local_cache->template get_sub_elem_cache<sub_elem_dim>(sub_elem_id_);

            if (!cache.fill_none())
            {
                if (cache.template status_fill<_Point>() || cache.template status_fill<_Value>())
                {
                    const auto points =
                        elem_.CartesianGridElement<dim>::template get_points<sub_elem_dim>(sub_elem_id_);

                    if (cache.template status_fill<_Point>())
                    {
                        auto &cache_pts = cache.template get_data<_Point>();
                        cache_pts = points;

                        cache.template set_status_filled<_Point>(true);
                    }
                    if (cache.template status_fill<_Value>())
                    {
                        const auto n_pts = points.get_num_points();

                        auto &values = cache.template get_data<_Value>();
                        for (int pt = 0 ; pt < n_pts ; ++pt)
                            for (int i = 0 ; i < dim ; ++i)
                                values[pt][i] = points[pt][i];

                        cache.template set_status_filled<_Value>(true);
                    }
                }
                if (cache.template status_fill<_Gradient>())
                {
                    // TODO (pauletti, Apr 17, 2015): this can be static const
                    const auto identity = create_id_tensor<dim,space_dim>();
                    cache.template get_data<_Gradient>().fill(identity);

                    cache.template set_status_filled<_Gradient>(true);
                }
                if (cache.template status_fill<_Hessian>())
                {
                    // TODO (pauletti, Apr 17, 2015): this can be static const
                    Hessian zero;
                    cache.template get_data<_Hessian>().fill(zero);

                    cache.template set_status_filled<_Hessian>(true);
                }
            }
            cache.set_filled(true);
        }

        const int sub_elem_id_;
        self_t &function_;
        ElementAccessor &elem_;
    };

    friend struct FillCacheDispatcher;


#ifdef MESH_REFINEMENT

    void create_connection_for_insert_knots(std::shared_ptr<self_t> &ig_function);

    void rebuild_after_insert_knots(
        const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
        const CartesianGrid<dim> &old_grid);

#endif // MESH_REFINEMENT

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
        ar &boost::serialization::make_nvp("IdentityFunction_base_t",
                                           boost::serialization::base_object<parent_t>(*this));
    }
    ///@}
#endif // SERIALIZATION

};



IGA_NAMESPACE_CLOSE

#endif
