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

#ifndef NEW_PUSH_FORWARD_ELEMENT_ACCESSOR_H_
#define NEW_PUSH_FORWARD_ELEMENT_ACCESSOR_H_

#include <igatools/geometry/new_push_forward.h>
#include <igatools/geometry/mapping_element.h>

IGA_NAMESPACE_OPEN


template<Transformation type_, int dim_, int codim_ = 0>
class PushForwardElement
    : public MappingElement<dim_, codim_>
{
private:
    using self_t  = PushForwardElement<type_, dim_, codim_>;
    using parent_t = MappingElement<dim_, codim_>;
    using PF      = NewPushForward<type_, dim_, codim_>;
public:

    using parent_t::dim;
    using parent_t::codim;
    using parent_t::space_dim;

    static const Transformation type = type_;
    // using ContainerType = const PF;
    using parent_t::MappingElement;

    template <int range, int rank>
    using RefValue = typename PF::template RefValue<range, rank>;

    template <int range, int rank>
    using PhysValue = typename PF::template PhysValue<range, rank>;

    template <int range, int rank, int order>
    using RefDerivative = typename PF::template RefDerivative<range, rank, order>;

    template <int range, int rank, int order>
    using PhysDerivative = typename PF::template PhysDerivative<range, rank, order>;


public:

    template <int range, int rank, Transformation ttype=type_>
    void
    transform_0(const ValueContainer<RefValue<range, rank>> &v_hat,
                ValueContainer< PhysValue<range, rank> > &v,
                EnableIf<ttype == Transformation::h_grad> *= 0) const
    {
        v = v_hat;
    }


    template <int range, int rank, Transformation ttype=type_>
    void
    transform_1(const std::tuple<
                ValueContainer<RefValue<range, rank>>,
                ValueContainer<RefDerivative<range, rank, 1>>> &ref_values,
                const ValueContainer<PhysValue<range, rank>>   &phys_values,
                ValueContainer<PhysDerivative<range, rank, 1>> &Dv,
                EnableIf<ttype == Transformation::h_grad> *= 0) const
    {
        const auto &Dv_hat = std::get<1>(ref_values);

        const int n_func   = Dv_hat.get_num_functions();
        const int n_points = Dv_hat.get_num_points();
        auto Dv_it     = Dv.begin();
        auto Dv_hat_it = Dv_hat.cbegin();

        const auto &DF_inv = this->get_inverse_gradients();
        for (int i_fn = 0; i_fn < n_func; ++i_fn)
            for (Index j_pt = 0; j_pt < n_points; ++j_pt)
            {
                (*Dv_it) = compose((*Dv_hat_it), DF_inv[j_pt]);
                ++Dv_hat_it;
                ++Dv_it;
            }
    }


    template <int range, int rank, Transformation ttype=type_>
    void
    transform_2(const std::tuple<
                ValueContainer<RefValue<range, rank>>,
                ValueContainer<RefDerivative<range, rank, 1>>,
                ValueContainer<RefDerivative<range, rank, 2>> > &ref_values,
                const std::tuple<
                ValueContainer<PhysValue<range, rank>>,
                ValueContainer<PhysDerivative<range, rank, 1>>> &phys_values,
                ValueContainer<PhysDerivative<range, rank, 2>> &D2v,
                EnableIf<ttype == Transformation::h_grad> *= 0) const
    {
        const auto &D2v_hat = std::get<2>(ref_values);
        const auto &D1v     = std::get<1>(phys_values);

        const int n_func   = D2v_hat.get_num_functions();
        const int n_points = D2v_hat.get_num_points();
        auto D2v_it     = D2v.begin();
        auto D1v_it     = D1v.cbegin();
        auto D2v_hat_it = D2v_hat.cbegin();
        const auto D2F     = this->get_hessians();
        const auto &DF_inv = this->get_inverse_gradients();

        for (int i_fn = 0; i_fn < n_func; ++i_fn)
            for (Index j_pt = 0; j_pt < n_points; ++j_pt)
            {
                for (int u=0; u<dim; ++u)
                {
                    (*D2v_it)[u] = compose(
                                       (*D2v_hat_it)[u] - compose((*D1v_it),D2F[j_pt][u]),
                                       DF_inv[j_pt]);
                }
                ++D2v_hat_it;
                ++D1v_it;
                ++D2v_it;
            }

    }


#if 0
    template <int order, int range, int rank, Transformation ttype=type >
    void
    transform_0(const ValueContainer<RefValue<range, rank>> &D0v_hat,
                ValueContainer< PhysValue<range, rank> > &D0v,
                EnableIf<ttype == Transformation::h_grad> *= 0) const;



    template < int range, int rank, Transformation ttype=type >
    void
    transform_values(
        const ValueContainer< RefValue<range, rank> > &D0v_hat,
        ValueContainer< PhysValue<range, rank> > &D0v,
        EnableIf<ttype == Transformation::h_div> *= 0) const;
#endif

private:
    template <typename Accessor> friend class GridForwardIterator;
    friend class NewPushForward<type_, dim_, codim_>;

};

IGA_NAMESPACE_CLOSE

#endif
