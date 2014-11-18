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

#ifndef IDENTITY_FUNCTIONS_H
#define IDENTITY_FUNCTIONS_H

#include <igatools/base/new_function.h>

IGA_NAMESPACE_OPEN

template<int dim>
auto
create_id_tensor()
{
    typename NewFunction<dim, 0, dim, 1>::Gradient res;
    for (int i=0; i<dim; ++i)
        res[i][i] = 1.;
    return res;
}



template<int dim>
class IdentityFunction : public NewFunction<dim, 0, dim, 1>
{
private:
    using parent_t = NewFunction<dim, 0, dim, 1>;
    using self_t = IdentityFunction<dim>;
protected:
    using typename parent_t::GridType;
public:
    using typename parent_t::variant_1;
    using typename parent_t::variant_2;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    using parent_t::space_dim;

    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    IdentityFunction(std::shared_ptr<GridType> grid);

    static std::shared_ptr<parent_t>
    create(std::shared_ptr<GridType> grid)
    {
        return std::shared_ptr<parent_t>(new self_t(grid));
    }

    std::shared_ptr<parent_t> clone() const override
    {

        return std::make_shared<self_t>(self_t(*this));
    }


    void fill_cache(ElementAccessor &elem, const int j, const variant_2 &k) override;


private:

    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T &quad)
        {
            auto &local_cache = function->get_cache(*elem);
            auto &cache = local_cache->template get_value_cache<T::k>(j);
            auto &flags = cache.flags_handler_;

            if (!flags.fill_none())
            {
                const auto points =
                    elem->CartesianGridElement<dim>::template get_points<T::k>(j);

                if (flags.fill_points())
                    cache.points_ = points;
                if (flags.fill_values())
                    std::get<0>(cache.values_) = points;
                if (flags.fill_gradients())
                {
                    auto identity = create_id_tensor<dim>();
                    std::get<1>(cache.values_).fill(identity);
                }
                if (flags.fill_hessians())
                {
                    Assert(false, ExcNotImplemented());
                    //std::get<2>(cache.values_) = 0.;
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
