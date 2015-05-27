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

#include <igatools/functions/functions_container.h>

IGA_NAMESPACE_OPEN


void
FunctionsContainer::
print_info(LogStream &out) const
{
    boost::fusion::for_each(data_varying_dim_,
                            [&](const auto & type_and_data_same_dim)
    {
        using Type_Value = typename std::remove_reference<decltype(type_and_data_same_dim)>::type;
        using Type = typename Type_Value::first_type;

        out.begin_item("Dim : " + std::to_string(Type::value));
        type_and_data_same_dim.second.print_info(out);
        out.end_item();

    } // end lambda function
                           );
}

#ifdef SERIALIZATION
template<class Archive>
void
FunctionsContainer::
serialize(Archive &ar, const unsigned int version)
{
    boost::fusion::for_each(data_varying_dim_,
                            [&](auto & type_and_data_same_dim)
    {
        using Type_Value = typename std::remove_reference<decltype(type_and_data_same_dim)>::type;
        using Type = typename Type_Value::first_type;

        const std::string tag_name = "data_dim_" + std::to_string(Type::value);
        ar &boost::serialization::make_nvp(tag_name.c_str(),type_and_data_same_dim.second);
    } // end lambda function
                           );
}
#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE

#ifdef SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(iga::FunctionsContainer)

template void iga::FunctionsContainer::serialize(OArchive &, const unsigned int);
template void iga::FunctionsContainer::serialize(IArchive &, const unsigned int);
#endif // SERIALIZATION

