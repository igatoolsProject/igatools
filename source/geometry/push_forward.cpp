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


#include <igatools/geometry/push_forward.h>


using std::array;
using std::shared_ptr;

IGA_NAMESPACE_OPEN


template< Transformation transformation_type_, int dim_, int codim_>
PushForward<transformation_type_, dim_, codim_>::
PushForward(const std::shared_ptr< Map > map)
    :
    map_(map)
{};


template< Transformation transformation_type_, int dim_, int codim_>
auto
PushForward<transformation_type_, dim_, codim_>::
create(const std::shared_ptr< Map > map) -> std::shared_ptr<Self>
{
    return std::shared_ptr<Self>(new Self(map));
}


template< Transformation transformation_type_, int dim_, int codim_>
PushForward<transformation_type_, dim_, codim_>::
PushForward(const Self &push_forward)
    :
    map_(push_forward.map_)
{}



template< Transformation transformation_type_, int dim_, int codim_>
auto
PushForward<transformation_type_, dim_, codim_>::
begin() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(), 0);
}



template< Transformation transformation_type_, int dim_, int codim_>
auto
PushForward<transformation_type_, dim_, codim_>::
end() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           IteratorState::pass_the_end);
}



template< Transformation transformation_type_, int dim_, int codim_>
void
PushForward<transformation_type_, dim_, codim_>::
reset_map(const std::shared_ptr<  Map > map)
{
    Assert(map != nullptr, ExcNullPtr());
    map_= map;
};

template< Transformation transformation_type_, int dim_, int codim_>
auto
PushForward<transformation_type_, dim_, codim_>::
clone() const -> std::shared_ptr<Self>
{
    return std::shared_ptr<Self>(new Self(*this));
}


template< Transformation transformation_type_, int dim_, int codim_>
auto
PushForward<transformation_type_, dim_, codim_>::
get_mapping() const -> shared_ptr< const Map >
{
    return std::shared_ptr< const Map >(map_);
}


template< Transformation transformation_type_, int dim_, int codim_>
void
PushForward<transformation_type_, dim_, codim_>::
print_info(LogStream &out) const
{
    out << "Transformation type: " << int(transformation_type) << std::endl;

    out.begin_item("Mapping:");
    map_->print_info(out);
    out.end_item();
}

template< Transformation transformation_type_, int dim_, int codim_>
void
PushForward<transformation_type_, dim_, codim_>::
print_memory_info(LogStream &out) const
{
    using std::endl;
    out << "PushForward memory info" << endl;
    out << "this address = " << this << endl;

    out.push("\t");
    out << "map_ memory address = " << map_ << endl;
    out.pop();
}

IGA_NAMESPACE_CLOSE


#include <igatools/geometry/push_forward.inst>


