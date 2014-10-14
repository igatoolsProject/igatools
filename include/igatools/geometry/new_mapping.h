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

#ifndef NEW_MAPPING_H_
#define NEW_MAPPING_H_

#include <igatools/base/config.h>
#include <igatools/base/new_function.h>

IGA_NAMESPACE_OPEN

//Forward declaration to avoid including header file.
template <int, int> class MappingElement;

/**
 * @brief The mapping is a deformation \f$ F : \hat\Omega \to \Omega\f$
 * which maps the reference domain \f$\hat\Omega \in \mathbb{R}^{dim}\f$ to the
 * physical domain \f$\Omega \in \mathbb{R}^{dim+codim}\f$.
 *
 * In igatools the mapping is a grid-like container, so that on top
 * of being a deformation is aware of an underlying grid structure.
 *
 * Most of the methods in this class are virtual,
 * because we want to use Mapping as general interface for user-defined Mapping
 * specializations (unknown at compile time).
 *
 * @ingroup containers
 *
 * @author pauletti 2014
 */
template<int dim, int codim = 0>
class NewMapping
{
private:
    using self_t = NewMapping<dim, codim>;
    using FuncType = NewFunction<dim, codim, dim+codim>;
public:
    using ElementAccessor = MappingElement<dim, codim>;
    using ElementIterator = GridForwardIterator<ElementAccessor>;

    static const int space_dim = FuncType::space_dim;

private:
    /** Type for the given order derivatives of the
     *  the mapping. */
    template<int order>
    using Derivative = typename FuncType::template Derivative<order>;

    /** Type for the diferent order derivatives of the inverse of
     * the mapping
     */
    template<int order>
    using InvDerivative = Derivatives<space_dim, dim, 1, order>;

public:
    /** Type of the mapping evaluation point. */
    using Point = typename FuncType::Point;

    /** Type of the mapping return value. */
    using Value = typename FuncType::Value;

    /** Type of the mapping gradient. */
    using Gradient = typename FuncType::Gradient;

    /** Typedef for the mapping hessian. */
    using Hessian = typename FuncType::Hessian;

public:

    NewMapping(std::shared_ptr<FuncType> F,
               const ValueFlags flag,
               const Quadrature<dim> &quad);

    NewMapping() = delete;

    ~NewMapping();

    void init_element(ElementIterator &elem);

    void fill_element(ElementIterator &elem);

private:
    std::shared_ptr<FuncType> F_;

    friend ElementAccessor;
};


template<int dim, int codim>
NewMapping<dim, codim>::
NewMapping(std::shared_ptr<FuncType> F,
           const ValueFlags flag,
           const Quadrature<dim> &quad)
    :
    F_(F)
{}



template<int dim, int codim>
NewMapping<dim, codim>::
~NewMapping()
{}



template<int dim, int codim>
auto
NewMapping<dim, codim>::
init_element(ElementIterator &elem) -> void
{
    auto &el = elem.get_accessor();
    F_->init_elem(el);
//    auto &cache = this->get_cache(elem);
//    if (cache == nullptr)
//    {
//        using Cache = typename ElementAccessor::CacheType;
//        cache = shared_ptr<Cache>(new Cache);
//    }
//    cache->resize(flag_, quad_.get_num_points());
}



template<int dim, int codim>
auto
NewMapping<dim, codim>::
fill_element(ElementIterator &elem) -> void
{
    auto &el    = elem.get_accessor();
    F_->fill_elem(el);
    //const auto points = el.CartesianGridElement<dim>::get_points();
//    auto &cache = this->get_cache(elem);
//    if (flag_.fill_points())
//        this->parametrization(points, cache->points_);
//    if (flag_.fill_values())
//        this->evaluate_0(cache->points_, cache->values_);
//    if (flag_.fill_gradients())
//        this->evaluate_1(cache->points_, std::get<1>(cache->derivatives_));
//    if (flag_.fill_hessians())
//        this->evaluate_2(cache->points_, std::get<2>(cache->derivatives_));
}

IGA_NAMESPACE_CLOSE

#endif
