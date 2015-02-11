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

#include <igatools/geometry/mapping.h>
#include <igatools/base/exceptions.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>


using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

template<int dim_, int codim_>
void
Mapping<dim_, codim_>::
evaluate(ValueVector<Value> &values) const
{
    Assert(false, ExcNotImplemented());
    AssertThrow(false, ExcNotImplemented());
}



template<int dim_, int codim_>
void
Mapping<dim_,codim_>::
evaluate_gradients(ValueVector<Gradient> &gradients) const
{
    Assert(false, ExcNotImplemented());
    AssertThrow(false, ExcNotImplemented());
}



template<int dim_, int codim_>
void
Mapping<dim_,codim_>::
evaluate_hessians(ValueVector<Hessian> &hessians) const
{
    Assert(false, ExcNotImplemented());
    AssertThrow(false, ExcNotImplemented());
}



template<int dim_, int codim_>
void
Mapping<dim_, codim_>::
evaluate_face(const Index face_id, ValueVector<Value> &values) const
{
    Assert(false,ExcNotImplemented()) ;
    AssertThrow(false, ExcNotImplemented());
}



template<int dim_, int codim_>
void
Mapping<dim_,codim_>::
evaluate_face_gradients(const Index face_id, ValueVector<Gradient> &gradients) const
{
    Assert(false,ExcNotImplemented()) ;
    AssertThrow(false, ExcNotImplemented());
}



template<int dim_, int codim_>
void
Mapping<dim_,codim_>::
evaluate_face_hessians(const Index face_id, ValueVector<Hessian> &hessians) const
{
    Assert(false,ExcNotImplemented()) ;
    AssertThrow(false, ExcNotImplemented());
}


template<int dim_, int codim_>
void
Mapping<dim_,codim_>::
evaluate_at_points(const ValueVector<Point> &points, ValueVector<Value> &values) const
{
    Assert(false,ExcMessage("This function must be implemented by the concrete class derived from Mapping")) ;
    AssertThrow(false, ExcMessage("This function must be implemented by the concrete class derived from Mapping"));
}

template<int dim_, int codim_>
void
Mapping<dim_,codim_>::
evaluate_gradients_at_points(const ValueVector<Point> &points, ValueVector<Gradient> &gradients) const
{
    Assert(false,ExcMessage("This function must be implemented by the concrete class derived from Mapping")) ;
    AssertThrow(false, ExcMessage("This function must be implemented by the concrete class derived from Mapping"));
}

template<int dim_, int codim_>
void
Mapping<dim_,codim_>::
evaluate_hessians_at_points(const ValueVector<Point> &points, ValueVector<Hessian> &hessians) const
{
    Assert(false,ExcMessage("This function must be implemented by the concrete class derived from Mapping")) ;
    AssertThrow(false, ExcMessage("This function must be implemented by the concrete class derived from Mapping"));
}

template<int dim_, int codim_>
Mapping<dim_,codim_>::
Mapping(const shared_ptr<GridType> grid)
    :
    GridWrapper<GridType>(grid)
{}



template<int dim_, int codim_>
Mapping<dim_,codim_>::
~Mapping()
{}



template<int dim_, int codim_>
Mapping<dim_,codim_>::
Mapping(const Mapping<dim_,codim_> &map)
    :
    GridWrapper<GridType>(make_shared<GridType>(GridType(*(map.get_grid()))))
{}


template<int dim_, int codim_>
ValueFlags
Mapping<dim_,codim_>::
required_flags() const
{
    return ValueFlags::none;
}



template<int dim_, int codim_>
auto
Mapping<dim_,codim_>::
values() const -> ValueVector<Value>
{
    Assert(false, ExcNotImplemented());
    AssertThrow(false, ExcNotImplemented());
    return (ValueVector<Value>());
}



template<int dim_, int codim_>
auto
Mapping<dim_,codim_>::
gradients() const -> ValueVector<Gradient>
{
    Assert(false, ExcNotImplemented());
    AssertThrow(false, ExcNotImplemented());
    return (ValueVector<Gradient>());
}



template<int dim_, int codim_>
auto
Mapping<dim_,codim_>::
hessians() const -> ValueVector<Hessian>
{
    Assert(false, ExcNotImplemented());
    AssertThrow(false, ExcNotImplemented());
    return (ValueVector<Hessian>());
}



template<int dim_, int codim_>
auto
Mapping<dim_,codim_>::
begin() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(), 0);
}



template<int dim_, int codim_>
auto
Mapping<dim_,codim_>::
last() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           this->get_grid()->get_num_active_elems() - 1);
}



template<int dim_, int codim_>
auto
Mapping<dim_,codim_>::
end() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           IteratorState::pass_the_end);
}



template<int dim_, int codim_>
void
Mapping<dim_,codim_>::
print_info(LogStream &out) const
{
    AssertThrow(false,
                ExcMessage("This function is called from an abstract class. Try to call a specialization of it from a derived (concrete) class."));
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/mapping.inst>
