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

#ifndef SUB_MAPPING_H_
#define SUB_MAPPING_H_

#include <igatools/base/new_function.h>
#include <igatools/base/function_element.h>

IGA_NAMESPACE_OPEN

/**
 *
 * @author pauletti 2014
 */
template<int k, int dim, int spacedim>
class SubFunction : public NewFunction<dim - k, 0, spacedim, 1>
{
public:
    using base_t = NewFunction<dim - k, 0, spacedim, 1>;
    using SupFunc = NewFunction<dim, 0, spacedim, 1>;

    using typename base_t::GridType;

    using SuperGrid = typename SupFunc::GridType;
    template <int j>
    using InterGridMap = typename SuperGrid::template InterGridMap<j>;

public:

    SubFunction(std::shared_ptr<GridType> grid,
                std::shared_ptr<SupFunc> func,
                const int s_id,
                InterGridMap<k> &elem_map)
:
    base_t(grid),
    func_(std::make_shared<SupFunc>(*func)),
    s_id_(s_id),
    elem_map_(elem_map)
{}


 //   void reset(const NewValueFlags &flag, const variant_1& quad) override;

//    void init_cache(ElementAccessor &elem, const variant_2& k) override;
//
//    void fill_cache(ElementAccessor &elem, const int j, const variant_2& k) override;


private:
    std::shared_ptr<SupFunc> func_;
    const int s_id_;
    InterGridMap<k> &elem_map_;

    typename SupFunc::ElementIterator elem_;


private:
    struct ResetDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T& quad)
        {
            (*flags_)[T::dim] = flag;
            auto sup_quad = extend_sub_elem_quad<T::dim, dim>(quad, s_id);
            func_->template reset<T::dim+k>(flag, sup_quad);
        }

        NewValueFlags flag;
        SupFunc *func_;
        std::array<FunctionFlags, dim - k + 1> *flags_;
        int s_id;
    };

#if 0

    struct InitCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T& quad)
        {
            func_->template init_cache<T::k>(*sup_elem);
        }

        SupFunc *func_;
        typename SupFunc::ElementIterator *sup_elem;


    };


    struct FillCacheDispatcher : boost::static_visitor<void>
    {
        template<class T>
        void operator()(const T& quad)
        {
            space_handler_->template fill_cache<T::k>(*space_elem, j);

            auto &local_cache = function->get_cache(*func_elem);
            auto &cache = local_cache->template get_value_cache<T::k>(j);
            auto &flags = cache.flags_handler_;

            if (flags.fill_values())
                std::get<0>(cache.values_) =
                        space_elem->template linear_combination<0, T::k>(*loc_coeff, j);
            if (flags.fill_gradients())
                std::get<1>(cache.values_) =
                        space_elem->template linear_combination<1, T::k>(*loc_coeff, j);
            if (flags.fill_hessians())
                std::get<2>(cache.values_) =
                        space_elem->template linear_combination<2, T::k>(*loc_coeff, j);

            cache.set_filled(true);
        }

        int j;
        self_t *function;
        typename Space::ElementHandler *space_handler_;
        ElementAccessor *func_elem;
        typename Space::ElementAccessor *space_elem;
        vector<Real> *loc_coeff;
    };
#endif


    ResetDispatcher reset_impl;
//    InitCacheDispatcher init_cache_impl;
//    FillCacheDispatcher fill_cache_impl;

};



IGA_NAMESPACE_CLOSE

#endif /* MAPPING_SLICE_H_ */

