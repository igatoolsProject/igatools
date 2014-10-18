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


#include <igatools/geometry/new_push_forward.h>


using std::array;
using std::shared_ptr;

IGA_NAMESPACE_OPEN

namespace
{
auto
pushforward_to_mapping_flag(const Transformation type, const ValueFlags v_flag)
-> ValueFlags
{
//    const ValueFlags common_flag =
//        ValueFlags::value|
//        ValueFlags::gradient|ValueFlags::map_hessian|
//        ValueFlags::map_inv_gradient|ValueFlags::map_inv_hessian|
//        ValueFlags::w_measure|ValueFlags::face_point|
//        ValueFlags::map_face_gradient|ValueFlags::map_face_hessian|
//        ValueFlags::map_face_inv_gradient|ValueFlags::map_face_inv_hessian|
//        ValueFlags::face_w_measure|ValueFlags::face_normal;

//    /*
//     * For each MappingValueFlags there is an if that checks for all
//     * ValueFlags that activate the given value flag.
//     */
//    ValueFlags fill_flag = common_flag & v_flag;

    ValueFlags fill_flag = ValueFlags::none;

//    if (contains(v_flag, ValueFlags::point))
//        fill_flag |= ValueFlags::point;
//
//    if (contains(v_flag, ValueFlags::w_measure))
//        fill_flag |= (ValueFlags::measure |
//                      ValueFlags::map_gradient);
//
//    if (contains(v_flag, ValueFlags::face_point))
//        fill_flag |= ValueFlags::face_point;
//
//    if (contains(v_flag, ValueFlags::face_w_measure))
//        fill_flag |= (ValueFlags::face_measure |
//                      ValueFlags::map_face_gradient);
//
//    if (contains(v_flag, ValueFlags::face_normal))
//        fill_flag |= (ValueFlags::map_face_inv_gradient |
//                      ValueFlags::map_face_gradient);

    if (type == Transformation::h_grad)
    {
        if (contains(v_flag, ValueFlags::tran_hessian))
        {
            fill_flag |= (ValueFlags::hessian|
            ValueFlags::map_inv_gradient);
        }

        if (contains(v_flag, ValueFlags::tran_value))
        {}

        if (contains(v_flag, ValueFlags::tran_gradient))
            fill_flag |= (ValueFlags::map_inv_gradient);


    }
//    else if (type == Transformation::h_div)
//    {
//        if (contains(v_flag,ValueFlags::tran_value))
//            fill_flag |= (ValueFlags::map_gradient |
//                          ValueFlags::map_face_gradient);
//        if (contains(v_flag,ValueFlags::tran_gradient))
//            fill_flag |= (ValueFlags::map_gradient |
//                          ValueFlags::map_hessian |
//                          ValueFlags::map_face_gradient |
//                          ValueFlags::map_face_hessian);
//        if (contains(v_flag,ValueFlags::tran_hessian))
//            AssertThrow(false,ExcNotImplemented());
//    }
//    else if (type == Transformation::h_curl)
//    {
//        AssertThrow(false,ExcNotImplemented());
//        if (contains(v_flag,ValueFlags::tran_value))
//            fill_flag |= (ValueFlags::map_gradient |
//                          ValueFlags::map_face_gradient);
//        if (contains(v_flag,ValueFlags::tran_gradient))
//            fill_flag |= (ValueFlags::map_gradient |
//                          ValueFlags::map_hessian |
//                          ValueFlags::map_face_gradient |
//                          ValueFlags::map_face_hessian);
//        if (contains(v_flag,ValueFlags::tran_hessian))
//            AssertThrow(false,ExcNotImplemented());
//    }
//    else if (type == Transformation::l_2)
//    {
//        AssertThrow(false,ExcNotImplemented());
//        if (contains(v_flag,ValueFlags::tran_value))
//            AssertThrow(false,ExcNotImplemented());
//        if (contains(v_flag,ValueFlags::tran_gradient))
//            AssertThrow(false,ExcNotImplemented());
//        if (contains(v_flag,ValueFlags::tran_hessian))
//            AssertThrow(false,ExcNotImplemented());
//    }
//
//
//
//    // We fill extra stuff as the computation is performed anyways
//    if (contains(fill_flag , ValueFlags::measure))
//        fill_flag |= (ValueFlags::map_gradient |
//                      ValueFlags::map_face_gradient);
//
//    if (contains(fill_flag , ValueFlags::map_inv_gradient))
//        fill_flag |= (ValueFlags::map_gradient |
//                      ValueFlags::measure |
//                      ValueFlags::map_face_gradient |
//                      ValueFlags::face_measure);
//
//    if (contains(fill_flag , ValueFlags::map_inv_hessian))
//        fill_flag |= (ValueFlags::map_hessian |
//                      ValueFlags::map_face_hessian);

    return fill_flag;
}
}



template<Transformation type, int dim, int codim>
NewPushForward<type, dim, codim>::
NewPushForward(std::shared_ptr<FuncType> F,
               const ValueFlags flag,
               const Quadrature<dim> &quad)
    :
    MapType::NewMapping(F, pushforward_to_mapping_flag(type, flag), quad)
{}



template<Transformation type, int dim, int codim>
auto
NewPushForward<type, dim, codim>::
init_element(ElementAccessor &elem) ->void
{
    MapType::init_element(elem);
}



template<Transformation type, int dim, int codim>
auto
NewPushForward<type, dim, codim>::
fill_element(ElementAccessor &elem) ->void
{
    MapType::fill_element(elem);
}

template<Transformation type, int dim, int codim>
auto
NewPushForward<type, dim, codim>::
init_element(ElementIterator &elem) ->void
{
    init_element(elem.get_accessor());
}



template<Transformation type, int dim, int codim>
auto
NewPushForward<type, dim, codim>::
fill_element(ElementIterator &elem) ->void
{
    fill_element(elem.get_accessor());
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/new_push_forward.inst>


