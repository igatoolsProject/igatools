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

#ifndef NEW_PUSHFORWARD_H_
#define NEW_PUSHFORWARD_H_

#include <igatools/base/config.h>
#include <igatools/geometry/new_mapping.h>
#include <igatools/geometry/new_mapping_element_accessor.h>

IGA_NAMESPACE_OPEN

//Forward declaration to avoid including header file.
template <Transformation, int, int> class PushForwardElement;

template<Transformation type, int dim, int codim = 0>
class NewPushForward : public NewMapping<dim, codim>
{
private:
    using self_t = NewPushForward<type, dim, codim>;
    using MapType = NewMapping<dim, codim>;
    using typename MapType::FuncType;
public:
    using ElementAccessor = PushForwardElement<type, dim, codim>;
    using ElementIterator = GridForwardIterator<ElementAccessor>;
#if 0
    static const int space_dim = dim + codim;

private:
    /** Type for the given order derivatives of the
     *  the mapping. */
    template<int order>
    using Derivative = typename MapType::template Derivative<order>;

public:
    /** Type for the diferent order derivatives of the inverse of
     * the mapping
     */
    template<int order>
    using InvDerivative = Derivatives<space_dim, dim, 1, order>;


    /** Type of the mapping evaluation point. */
    using Point = typename MapType::Point;

    /** Type of the mapping return value. */
    using Value = typename MapType::Value;

    /** Type of the mapping gradient. */
    using Gradient = typename MapType::Gradient;

    /** Typedef for the mapping hessian. */
    using Hessian = typename MapType::Hessian;
#endif
public:

    NewPushForward(std::shared_ptr<FuncType> F,
                   const ValueFlags flag,
                   const Quadrature<dim> &quad);

   // NewMapping() = delete;

    ~NewPushForward() = default;

    void init_element(ElementIterator &elem);

    void fill_element(ElementIterator &elem);

//protected:
//    std::shared_ptr<typename ElementAccessor::CacheType>
//    &get_cache(ElementAccessor &elem);

private:
//    std::shared_ptr<MapType> F_;
//    MappingElemValueFlagsHandler flag_;
//    Quadrature<dim> quad_;
    //friend ElementAccessor;
};


template<Transformation type, int dim, int codim>
NewPushForward<type, dim, codim>::
NewPushForward(std::shared_ptr<FuncType> F,
               const ValueFlags flag,
               const Quadrature<dim> &quad)
               :
               MapType::NewMapping(F, flag, quad)
               {}


template<Transformation type, int dim, int codim>
auto
NewPushForward<type, dim, codim>::
init_element(ElementIterator &elem) ->void
{
    auto &el = elem.get_accessor();
    MapType::init_element(el);
}



template<Transformation type, int dim, int codim>
auto
NewPushForward<type, dim, codim>::
fill_element(ElementIterator &elem) ->void
{
    auto &el = elem.get_accessor();
    MapType::fill_element(el);
}

IGA_NAMESPACE_CLOSE

#endif
