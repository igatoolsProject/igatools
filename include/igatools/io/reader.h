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



#ifndef READER_H_
#define READER_H_

#include <igatools/base/config.h>
#include <igatools/geometry/mapping.h>

#include <memory>



IGA_NAMESPACE_OPEN

//TODO: This function has to be documented
//TODO: Document the xml format

/**
 * Reads an IgMapping from an xml file.
 *
 * @note The reference space for the IgMapping can be either BSplineSpace or NURBSSpace,
 * because this function returns a shared_ptr to the base class Mapping.
 *
 * @todo document the XML file formats (version 1.0 and 2.0) for IgMapping
 *
 * @author M. Martinelli
 * @date 04 Mar 2014
 */
template <int dim, int codim = 0>
std::shared_ptr< Mapping<dim,codim> >
ig_mapping_reader(const std::string &filename);


/**
 * Reads an IgMapping from an xml file in which the IgMapping
 * is described with the format version 2.0.
 *
 * @note The reference space for the IgMapping can be either BSplineSpace or NURBSSpace,
 * because this function returns a shared_ptr to the base class Mapping.
 *
 * @todo document the XML file format version 1.0 for IgMapping
 *
 * @author M. Martinelli
 * @date 04 Mar 2014
 */
template <int dim, int codim = 0>
std::shared_ptr< Mapping<dim,codim> >
ig_mapping_reader_version_1_0(const std::string &filename);


/**
 * Reads an IgMapping from an xml file in which the IgMapping
 * is described with the format version 2.0.
 *
 * @note The reference space for the IgMapping can be either BSplineSpace or NURBSSpace,
 * because this function returns a shared_ptr to the base class Mapping.
 *
 * @todo document the XML file format version 2.0 for IgMapping
 *
 * @author M. Martinelli
 * @date 04 Mar 2014
 */
template <int dim, int codim = 0>
std::shared_ptr< Mapping<dim,codim> >
ig_mapping_reader_version_2_0(const std::string &filename);


IGA_NAMESPACE_CLOSE


#endif

