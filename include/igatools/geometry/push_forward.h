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

#ifndef NEW_PUSHFORWARD_H_
#define NEW_PUSHFORWARD_H_

#include <igatools/base/config.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_element.h>

IGA_NAMESPACE_OPEN
constexpr
int physical_range(const int ref_range, const int space_dim, const Transformation type)
{
    return type == Transformation::h_grad ? ref_range : space_dim;
}

//Forward declaration to avoid including header file.
template <Transformation, int, int> class PushForwardElement;

/**
 *
 * @ingroup containers
 */
template<Transformation type_, int dim_, int codim_ = 0>
class PushForward : public Mapping<dim_, codim_>
{
private:
    using self_t = PushForward<type_, dim_, codim_>;
    using MapType = Mapping<dim_, codim_>;
    using typename MapType::FuncType;
public:
    using Topology = UnitElement<dim_>;

    using ElementAccessor = PushForwardElement<type_, dim_, codim_>;
    using ElementIterator = CartesianGridIterator<ElementAccessor>;
    using MapType::space_dim;
    static const int dim   = dim_;
    static const int codim = codim_;
    static const Transformation type = type_;

    template<int ref_range>
    struct PhysRange
    {
        static const int value = physical_range(ref_range, space_dim, type_);
    };

    template <int range, int rank>
    using RefValue = Values<dim, range, rank>;

    template <int range, int rank, int order>
    using RefDerivative = Derivatives<dim, range, rank, order>;

    template <int range, int rank>
    using PhysValue = Values<space_dim, PhysRange<range>::value,rank>;

    template <int range, int rank, int order>
    using PhysDerivative = Derivatives<space_dim, PhysRange<range>::value, rank, order>;


public:

    PushForward(std::shared_ptr<FuncType> F);

    ~PushForward() = default;

    template<int k>
    void reset(const ValueFlags flag, const Quadrature<k> &eval_pts);

    std::shared_ptr<ElementAccessor> create_element(const Index flat_index) const
    {
        auto elem = std::shared_ptr<ElementAccessor>(
                        new ElementAccessor(this->get_function(),flat_index));
        Assert(elem != nullptr,ExcNullPtr());

        return elem;
    }

    auto begin()  const -> ElementIterator
    {
        return ElementIterator(this->create_element(0),ElementProperties::none);
    }

    auto end() const -> ElementIterator
    {
        return ElementIterator(this->create_element(IteratorState::pass_the_end),ElementProperties::none);
    }

};


IGA_NAMESPACE_CLOSE

#endif
