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
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
//    using parent_t::space_dim;

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
        template<class T>
        void operator()(const T &quad)
        {
            auto &local_cache = function->get_cache(*elem);
            auto &cache = local_cache->template get_value_cache<T::k>(j);
            auto &flags = cache.flags_handler_;

            if (!flags.fill_none())
            {

                if (flags.fill_points() || flags.fill_values())
                {
                    const auto points =
                        elem->CartesianGridElement<dim>::template get_points<T::k>(j);

                    if (flags.fill_points())
                    {
                        cache.points_ = points;
                    }
                    if (flags.fill_values())
                    {
                        const auto n_pts = points.get_num_points();

                        auto &values = std::get<0>(cache.values_);
                        for (int pt = 0 ; pt < n_pts ; ++pt)
                            for (int i = 0 ; i < dim ; ++i)
                                values[pt][i] = points[pt][i];
                    }
                }
                if (flags.fill_gradients())
                {
                    auto identity = create_id_tensor<dim,space_dim>();
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
