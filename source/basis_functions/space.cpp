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

#include <igatools/basis_functions/space.h>
#include <igatools/utils/unique_id_generator.h>

IGA_NAMESPACE_OPEN

template <int dim_>
Space<dim_>::
Space(std::shared_ptr<GridType> grid)
    :
    base_t(grid),
    space_id_(UniqueIdGenerator::get_unique_id())
{};


template <int dim_>
Index
Space<dim_>::
get_space_id() const
{
    return space_id_;
}

#ifdef SERIALIZATION
template <int dim_>
template<class Archive>
void
Space<dim_>::
    serialize(Archive &ar, const unsigned int version)
    {
        ar &boost::serialization::make_nvp("Space_base_t",
                                           boost::serialization::base_object<base_t>(*this));

        ar &boost::serialization::make_nvp("space_id_",space_id_);
    }
    ///@}
#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/space.inst>
