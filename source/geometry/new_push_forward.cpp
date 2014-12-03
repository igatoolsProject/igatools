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
pushforward_to_mapping_flag(const Transformation type, const NewValueFlags flags)
-> NewValueFlags
{
    NewValueFlags transfer_flag =
    NewValueFlags::measure |
    NewValueFlags::w_measure |
    NewValueFlags::outer_normal|
    NewValueFlags::boundary_normal|
    NewValueFlags::point|
    NewValueFlags::value|
    NewValueFlags::gradient|
    NewValueFlags::hessian;

    NewValueFlags map_flag = flags & transfer_flag;

    if (type == Transformation::h_grad)
    {
        if (contains(flags, NewValueFlags::tran_value))
        {}

        if (contains(flags, NewValueFlags::tran_gradient))
            map_flag|= (NewValueFlags::inv_gradient);

        if (contains(flags, NewValueFlags::tran_hessian))
        {
            map_flag |= (NewValueFlags::hessian | NewValueFlags::inv_gradient);
        }
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

    return map_flag;
}
}


template<Transformation type, int dim, int codim>
NewPushForward<type, dim, codim>::
NewPushForward(std::shared_ptr<FuncType> F)
    :
    MapType::NewMapping(F)
{}



template<Transformation type, int dim, int codim>
template<int k>
auto
NewPushForward<type, dim, codim>::
reset(const NewValueFlags flag, const Quadrature<k> &quad) -> void
{
    MapType::template reset<k>(pushforward_to_mapping_flag(type, flag), quad);
}



IGA_NAMESPACE_CLOSE

#include <igatools/geometry/new_push_forward.inst>


