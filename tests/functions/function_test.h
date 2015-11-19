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


template<int dim, int codim, int range, int rank>
void
function_values(shared_ptr<const Function<dim, codim, range, rank>> func)
{
  const int sdim = dim;
  using Flags = function_element::Flags;
  auto flag = Flags::D0 | Flags::D1;
  auto handler = func->create_cache_handler();

  handler->template set_flags<sdim>(flag);
  auto quad   = QGauss<sdim>::create(2);

  auto elem = func->cbegin();
  auto end  = func->cend();
  handler->init_cache(elem, quad);

  for (; elem != end; ++elem)
  {
    handler->template fill_cache<dim>(elem, 0);
    elem->template get_values_from_cache<function_element::template _D<0>, dim>(0).print_info(out);
    out << endl;
    elem->template get_values_from_cache<function_element::template _D<1>, dim>(0).print_info(out);
    out << endl;
//    elem->template get_values<function_element::_D2, dim>(0).print_info(out);
//    out << endl;
  }
}
