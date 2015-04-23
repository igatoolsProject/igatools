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

#include <igatools/base/function.h>
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



template<int dim,int space_dim = dim>
class IdentityFunction :
    public Function<dim, 0, space_dim, 1>,
    public std::enable_shared_from_this<IdentityFunction<dim,space_dim> >
{
private:
    using parent_t = Function<dim, 0, space_dim, 1>;
    using self_t = IdentityFunction<dim,space_dim>;

    std::shared_ptr<const parent_t> shared_from_derived() const override final
    {
        return this->shared_from_this();
    }

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

    static std::shared_ptr<parent_t>
    create(std::shared_ptr<GridType> grid)
    {
        return std::shared_ptr<parent_t>(new self_t(grid));
    }

    std::shared_ptr<parent_t> clone() const override final
    {

        return std::make_shared<self_t>(self_t(*this));
    }


    void fill_cache(ElementAccessor &elem, const topology_variant &k, const int j) override;

    virtual void print_info(LogStream &out) const override final
    {
        using std::to_string;
        out.begin_item("IdentityFunction<" + to_string(dim) + "," + to_string(space_dim) +">");
        parent_t::print_info(out);
        out.end_item();
    }


private:

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<int sub_elem_dim>
        void operator()(const Topology<sub_elem_dim> &sub_elem)
        {
            auto &local_cache = function->get_cache(*elem);
            auto &cache = local_cache->template get_value_cache<sub_elem_dim>(j);
            auto &flags = cache.flags_handler_;

            if (!flags.fill_none())
            {

                if (flags.template fill<_Point>() || flags.template fill<_Value>())
                {
                    const auto points =
                        elem->CartesianGridElement<dim>::template get_points<sub_elem_dim>(j);

                    if (flags.template fill<_Point>())
                    {
                        auto &cache_pts = cache.template get_der<_Point>();
                        cache_pts = points;

                        flags.template set_filled<_Point>(true);
                    }
                    if (flags.template fill<_Value>())
                    {
                        const auto n_pts = points.get_num_points();

                        auto &values = cache.template get_der<_Value>();
                        for (int pt = 0 ; pt < n_pts ; ++pt)
                            for (int i = 0 ; i < dim ; ++i)
                                values[pt][i] = points[pt][i];

                        flags.template set_filled<_Value>(true);
                    }
                }
                if (flags.template fill<_Gradient>())
                {
                    // TODO (pauletti, Apr 17, 2015): this can be static const
                    auto identity = create_id_tensor<dim,space_dim>();
                    cache.template get_der<_Gradient>().fill(identity);
//                    std::get<1>(cache.values_).fill(identity);

                    flags.template set_filled<_Gradient>(true);
                }
                if (flags.template fill<_Hessian>())
                {
                    // TODO (pauletti, Apr 17, 2015): this can be static const
                    Hessian zero;
                    cache.template get_der<_Hessian>().fill(zero);
                }
            }
            cache.set_filled(true);
        }

        int j;
        self_t *function;
        ElementAccessor *elem;
        std::array<FunctionFlags, dim + 1> *flags_;
    };

    FillCacheDispatcher fill_cache_impl;
    friend class FillCacheDispatcher;
};



IGA_NAMESPACE_CLOSE

#endif
