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

/*
 *  Test for Function class, as a prototype for an spline function
 *  author: pauletti
 *  date: Oct 11, 2014
 */

#include "../tests.h"

#include <igatools/base/ig_function.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/base/function_element.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
//#include <igatools/linear_algebra/distributed_vector.h>

#if 0
template<class Space>
class IgFunction : public NewFunction<Space::dim, Space::range, Space::rank>
{
public:
    static const int dim = Space::dim;
    static const int range = Space::range;
    static const int rank = Space::rank;

    using parent_t = NewFunction<dim, range>;
    using typename NewFunction<dim, range>::Point;
    using typename NewFunction<dim, range>::Value;
    using typename NewFunction<dim, range>::Gradient;
    using typename NewFunction<dim, range>::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;

    using CoeffType = Vector<LAPack::trilinos>;

    IgFunction(const ValueFlags &flag, const Quadrature<dim> &quad,
               std::shared_ptr<const Space> space,
               const CoeffType &coeff)
        :
        parent_t::NewFunction(space->get_grid(), flag, quad),
        flag_(flag),
        quad_(quad),
        space_(space),
        coeff_(coeff),
        elem_(space_->begin()),
        space_filler_(space_, flag, quad)
    {}


    void init_element(ElementIterator &elem)
    {
        auto &el    = elem.get_accessor();
        GridUniformQuadCache<dim>::init_element_cache(el);
        auto &cache = this->get_cache(elem);
        if (cache == nullptr)
        {
            using Cache = typename ElementAccessor::CacheType;
            cache = shared_ptr<Cache>(new Cache);
        }
        cache->resize(flag_, quad_.get_num_points());

        space_filler_.init_element_cache(elem_);
    }

    void fill_element(ElementIterator &elem)
    {
        auto &el    = elem.get_accessor();
        GridUniformQuadCache<dim>::fill_element_cache(el);
        space_filler_.fill_element_cache(elem_);
        auto &cache = this->get_cache(elem);

        elem_.move_to(elem->get_flat_index());
        const auto loc_coeff = coeff_.get_local_coefs(elem_->get_local_to_global());
        if (flag_.fill_values())
            cache->values_ = elem_->evaluate_field(loc_coeff);
        if (flag_.fill_gradients())
            std::get<1>(cache->derivatives_) = elem_->evaluate_field_gradients(loc_coeff);
        if (flag_.fill_hessians())
            std::get<2>(cache->derivatives_) = elem_->evaluate_field_hessians(loc_coeff);
    }

private:
    ValueFlagsHandler flag_;
    Quadrature<dim> quad_;
    std::shared_ptr<const Space> space_;
    const CoeffType coeff_;
    typename Space::ElementIterator elem_;
    typename Space::UniformQuadCache space_filler_;
};


#endif

template<int dim, int range>
void test()
{
    using Space = BSplineSpace<dim>;
    using Function = IgFunction<Space>;


    auto flag = ValueFlags::value| ValueFlags::gradient |  ValueFlags::hessian;
//    auto flag = ValueFlags::point | ValueFlags::value | ValueFlags::gradient |
//                ValueFlags::hessian;
    auto quad = QGauss<dim>(2);
    auto grid = CartesianGrid<dim>::create(3);
    const int deg = 1;
    auto space = Space::create(deg, grid);
    typename Function::CoeffType coeff(space->get_num_basis());
    coeff(0) = 1.;
    Function F(flag, quad, space, coeff);


    GridForwardIterator<FunctionElement<dim,0,range,1>> elem(grid, 0);
    GridForwardIterator<FunctionElement<dim,0,range,1>> end(grid, IteratorState::pass_the_end);

    F.init_element(elem);
    for (; elem != end; ++elem)
    {
        F.fill_element(elem);
        //elem->get_points().print_info(out);
        //out << endl;
        elem->get_values().print_info(out);
        out << endl;
        elem->get_gradients().print_info(out);
        out << endl;
        elem->get_hessians().print_info(out);
        out << endl;
    }
}


int main()
{
    test<2,1>();
//    test<3,3>();

    return 0;
}

