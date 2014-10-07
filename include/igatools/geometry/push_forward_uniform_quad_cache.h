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

#ifndef __PUSH_FORWARD_UNIFORM_QUAD_CACHE_H_
#define __PUSH_FORWARD_UNIFORM_QUAD_CACHE_H_

#include <igatools/base/config.h>
#include <igatools/geometry/push_forward.h>
#include<igatools/geometry/mapping_uniform_quad_cache.h>

//#include <igatools/utils/value_table.h>
//#include <igatools/utils/value_vector.h>
//#include <igatools/base/quadrature.h>
//#include <igatools/geometry/mapping_element_accessor.h>

IGA_NAMESPACE_OPEN

template <class PushForward_>
class PushFowardUniformQuadCache :
    public MappingUniformQuadCache<PushForward_::dim, PushForward_::codim>
{
    using base_t = MappingUniformQuadCache<PushForward_::dim, PushForward_::codim>;
    using PF = PushForward_;
    static const Transformation transformation_type = PushForward_::transformation_type;
    using ElementIterator = typename PushForward_::ElementIterator;
protected:
    using ElementAccessor = typename PushForward_::ElementAccessor;
    void init_element_cache(ElementAccessor &elem)
    {
        base_t::init_element_cache(elem);
    }
    void fill_element_cache(ElementAccessor &elem)
    {
        base_t::fill_element_cache(elem);
        elem.fill_cache();
    }
    //  void fill_face_cache(ElementAccessor &elem, const int face);

public:
    static const int dim = PF::dim;

    //Allocates and fill the (global) cache
    PushFowardUniformQuadCache(std::shared_ptr<const PF> pf,
                               const ValueFlags flag,
                               const Quadrature<dim> &quad)
        :
        base_t(pf->get_mapping(), value_to_mapping_flag(flag), quad),
        flags_(value_to_mapping_flag(flag)),
        quad_(quad)
    {}

    //Allocates the ElementIterator element_cache
    void init_element_cache(ElementIterator &elem)
    {
        init_element_cache(elem.get_accessor());
    }
    //Fill the ElementIterator element_cache
    void fill_element_cache(ElementIterator &elem)
    {
        fill_element_cache(elem.get_accessor());
    }

    /**
     * Fills the ElementIterator face_cache
     * element dependent part
     */
    //void fill_face_cache(ElementIterator &elem, const int face);

    void print_info(LogStream &out) const
    {
        base_t::print_info(out);
    }
private:
    ValueFlags flags_;
    Quadrature<dim> quad_;

    auto value_to_mapping_flag(
        const ValueFlags v_flag) const -> ValueFlags
    {
        const ValueFlags common_flag =
            ValueFlags::point|ValueFlags::map_gradient|ValueFlags::map_hessian|
            ValueFlags::w_measure|ValueFlags::face_point|ValueFlags::map_face_gradient|
            ValueFlags::map_face_hessian|ValueFlags::face_w_measure|ValueFlags::face_normal;

        /*
         * For each MappingValueFlags there is an if that checks for all
         * ValueFlags that activate the given value flag.
         */
        ValueFlags fill_flag = common_flag & v_flag;

        if (contains(v_flag, ValueFlags::point))
            fill_flag |= ValueFlags::point;

        if (contains(v_flag, ValueFlags::w_measure))
            fill_flag |= (ValueFlags::measure |
                          ValueFlags::map_gradient);

        if (contains(v_flag, ValueFlags::face_point))
            fill_flag |= ValueFlags::face_point;

        if (contains(v_flag, ValueFlags::face_w_measure))
            fill_flag |= (ValueFlags::face_measure |
                          ValueFlags::map_face_gradient);

        if (contains(v_flag, ValueFlags::face_normal))
            fill_flag |= (ValueFlags::map_face_inv_gradient |
                          ValueFlags::map_face_gradient);




        if (transformation_type == Transformation::h_grad)
        {
            auto flag = v_flag;
            if (contains(v_flag,ValueFlags::tran_hessian))
            {
                flag |= ValueFlags::tran_gradient;
                fill_flag |= (ValueFlags::map_hessian|
                              ValueFlags::map_inv_gradient|
                              ValueFlags::map_face_hessian|
                              ValueFlags::map_face_inv_hessian |
                              ValueFlags::map_face_inv_gradient);
            }
            if (contains(flag,ValueFlags::tran_value))
                fill_flag |= (ValueFlags::point | ValueFlags::face_point);

            if (contains(flag,ValueFlags::tran_gradient))
                fill_flag |= (ValueFlags::map_gradient |
                              ValueFlags::map_inv_gradient|
                              ValueFlags::map_face_gradient |
                              ValueFlags::map_face_inv_gradient);


        }
        else if (transformation_type == Transformation::h_div)
        {
            if (contains(v_flag,ValueFlags::tran_value))
                fill_flag |= (ValueFlags::map_gradient |
                              ValueFlags::map_face_gradient);
            if (contains(v_flag,ValueFlags::tran_gradient))
                fill_flag |= (ValueFlags::map_gradient |
                              ValueFlags::map_hessian |
                              ValueFlags::map_face_gradient |
                              ValueFlags::map_face_hessian);
            if (contains(v_flag,ValueFlags::tran_hessian))
                AssertThrow(false,ExcNotImplemented());
        }
        else if (transformation_type == Transformation::h_curl)
        {
            AssertThrow(false,ExcNotImplemented());
            if (contains(v_flag,ValueFlags::tran_value))
                fill_flag |= (ValueFlags::map_gradient |
                              ValueFlags::map_face_gradient);
            if (contains(v_flag,ValueFlags::tran_gradient))
                fill_flag |= (ValueFlags::map_gradient |
                              ValueFlags::map_hessian |
                              ValueFlags::map_face_gradient |
                              ValueFlags::map_face_hessian);
            if (contains(v_flag,ValueFlags::tran_hessian))
                AssertThrow(false,ExcNotImplemented());
        }
        else if (transformation_type == Transformation::l_2)
        {
            AssertThrow(false,ExcNotImplemented());
            if (contains(v_flag,ValueFlags::tran_value))
                AssertThrow(false,ExcNotImplemented());
            if (contains(v_flag,ValueFlags::tran_gradient))
                AssertThrow(false,ExcNotImplemented());
            if (contains(v_flag,ValueFlags::tran_hessian))
                AssertThrow(false,ExcNotImplemented());
        }



        // We fill extra stuff as the computation is performed anyways
        if (contains(fill_flag , ValueFlags::measure))
            fill_flag |= (ValueFlags::map_gradient |
                          ValueFlags::map_face_gradient);

        if (contains(fill_flag , ValueFlags::map_inv_gradient))
            fill_flag |= (ValueFlags::map_gradient |
                          ValueFlags::measure |
                          ValueFlags::map_face_gradient |
                          ValueFlags::face_measure);

        if (contains(fill_flag , ValueFlags::map_inv_hessian))
            fill_flag |= (ValueFlags::map_hessian |
                          ValueFlags::map_face_hessian);

        return fill_flag;
    }

};

IGA_NAMESPACE_CLOSE

#endif
