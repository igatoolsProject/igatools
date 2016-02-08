//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

#include <igatools/utils/safe_stl_container.h>
#include <igatools/base/exceptions.h>

#include <set>
#include <map>

IGA_NAMESPACE_OPEN
#if 0
template <class STLContainer>
template <class A>
EnableIf<(!std::is_pointer<STLContainer>::value &&!std::is_arithmetic<A>::value &&!std::is_pointer<A>::value), void>
SafeSTLContainer<STLContainer>::
t_print_info(LogStream &out) const
{
  int entry_id = 0;
  for (auto &entry : *this)
  {
    out.begin_item("Entry id: " + std::to_string(entry_id++));
    entry.print_info(out);
    out.end_item();
  }
}

template <class STLContainer>
template <class A>
EnableIf<(!std::is_pointer<STLContainer>::value &&std::is_arithmetic<A>::value &&!std::is_pointer<A>::value), void>
SafeSTLContainer<STLContainer>::
t_print_info(LogStream &out) const
{
  out << "[ ";
  for (auto &entry : *this)
    out << entry << " ";
  out << "]";
}


template <class STLContainer>
template <class A>
EnableIf<(!std::is_pointer<STLContainer>::value &&!std::is_arithmetic<A>::value &&std::is_pointer<A>::value), void>
SafeSTLContainer<STLContainer>::
t_print_info(LogStream &out) const
{
  int entry_id = 0;
  for (auto &entry : *this)
  {
    out.begin_item("Entry id: " + std::to_string(entry_id++));
    entry->print_info(out);
    out.end_item();
  }
}

template <class STLContainer>
template <class A>
EnableIf<(!std::is_pointer<STLContainer>::value &&std::is_arithmetic<A>::value &&std::is_pointer<A>::value), void>
SafeSTLContainer<STLContainer>::
t_print_info(LogStream &out) const
{
  out << "[ ";
  for (auto &entry : *this)
    out << *entry << " ";
  out << "]";
}





template <class STLContainer>
template <class A>
EnableIf<(std::is_pointer<STLContainer>::value &&!std::is_arithmetic<A>::value &&!std::is_pointer<A>::value), void>
SafeSTLContainer<STLContainer>::
t_print_info(LogStream &out) const
{
  int entry_id = 0;
  for (auto &entry : *(*this))
  {
    out.begin_item("Entry id: " + std::to_string(entry_id++));
    entry.print_info(out);
    out.end_item();
  }
}

template <class STLContainer>
template <class A>
EnableIf<(std::is_pointer<STLContainer>::value &&std::is_arithmetic<A>::value &&!std::is_pointer<A>::value), void>
SafeSTLContainer<STLContainer>::
t_print_info(LogStream &out) const
{
  out << "[ ";
  for (auto &entry : *(*this))
    out << entry << " ";
  out << "]";
}


template <class STLContainer>
template <class A>
EnableIf<(std::is_pointer<STLContainer>::value &&!std::is_arithmetic<A>::value &&std::is_pointer<A>::value), void>
SafeSTLContainer<STLContainer>::
t_print_info(LogStream &out) const
{
  int entry_id = 0;
  for (auto &entry : *(*this))
  {
    out.begin_item("Entry id: " + std::to_string(entry_id++));
    entry->print_info(out);
    out.end_item();
  }
}

template <class STLContainer>
template <class A>
EnableIf<(std::is_pointer<STLContainer>::value &&std::is_arithmetic<A>::value &&std::is_pointer<A>::value), void>
SafeSTLContainer<STLContainer>::
t_print_info(LogStream &out) const
{
  out << "[ ";
  for (auto &entry : *(*this))
    out << *entry << " ";
  out << "]";
}

template <class STLContainer>
void
SafeSTLContainer<STLContainer>::
print_info(LogStream &out) const
{
//  t_print_info<typename STLContainer::value_type>(out);
  Assert(false,ExcNotImplemented());
}
#endif

IGA_NAMESPACE_CLOSE

#include <igatools/utils/safe_stl_container.inst>
