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



#ifndef READER_H_
#define READER_H_

#include <igatools/base/config.h>
#include <igatools/functions/function.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/nurbs_space.h>

#include <boost/property_tree/ptree.hpp>

#include <memory>

IGA_NAMESPACE_OPEN

/**
 * Collection of useful functions for reading xml trees.
 */
namespace reader_utils
{
/**
 * Returns the XML tree contained in the file @p filename.
 */
boost::property_tree::ptree
get_xml_tree(const std::string &filename);

/**
 * Extracts from the XML @p element tree, the subtrees corresponding to the tag attributes.
 * @note If any element is present in the tree, an assertion will be raised
 * (in Debug and Release mode).
 */
boost::property_tree::ptree
get_xml_element_attributes(const boost::property_tree::ptree &element);


/**
 * Returns the number of nodes in the XML @p tree that have the tag @p tag_name.
 */
Size
count_xml_elements_same_tag(const boost::property_tree::ptree &tree,
                            const std::string &tag_name);

/**
 * Returns true if the XML @p tree has at least one node with the tag @p tag_name.
 */
bool
xml_element_is_present(const boost::property_tree::ptree &tree,
                       const std::string &tag_name);

/**
 * Returns true if the XML @p tree has exactly one node with the tag @p tag_name.
 * @note If no node with the tag @p tag_name are present,
 * an assertion (in Debug and Release mode) will be raised.
 */
bool
xml_element_is_unique(const boost::property_tree::ptree &tree,
                      const std::string &tag_name);

/**
 * Extracts from the XML @tree, the subtrees corresponding to the tag @p tag_name.
 * @note If any element is present in the tree, an assertion will be raised
 * (in Debug and Release mode).
 */
SafeSTLVector< boost::property_tree::ptree >
get_xml_element_vector(const boost::property_tree::ptree &tree,
                       const std::string &tag_name);

/**
 * Extracts from the XML @tree, the unique subtree corresponding to the tag @p tag_name.
 * @note The must be only one element in the @p tree with the given @p tag_name,
 * otherwise (more than one element on no element present) an assertion will be raised
 * (in Debug and Release mode).
 */
boost::property_tree::ptree
get_xml_element(const boost::property_tree::ptree &tree,
                const std::string &tag_name);


/**
 * Extracts a vector of scalars from the XML @p tree.
 * The type of scalars is determined by the template parameter @p ScalarType.
 */
template <class ScalarType>
SafeSTLVector<ScalarType>
get_vector_data_from_xml(const boost::property_tree::ptree &tree);

/**
 * Extracts the multiplicity of the internal knots from the XML @p tree.
 * The number of multiplicity vectors is determined by the template parameter
 * @p dim.
 */
template <int dim>
CartesianProductArray<Size,dim>
get_interior_multiplicity_from_xml(const boost::property_tree::ptree &tree);

} // end of namespace reader_utils.

/**
 * Returns a string containing the format of the igatools XML input file.
 *
 * @ingroup input
 * @author M. Martinelli
 * @date 04 Mar 2014
 */
std::string
get_xml_input_file_format(const std::string &filename);

/**
 * Reads an MapFunction from an xml file.
 *
 * @note Due to the virtual nature of the class MapFunction, this function can allocate
 * different kind of functions:
 *    - FormulaFunction and its specializations (not implemented)
 *    - IgFunction (implemented)
 *
 * @todo document the XML file formats (version 1.0 and 2.0) for IgMapping
 *
 * @ingroup input
 * @author M. Martinelli
 * @date 04 Mar 2014
 */
template <int dim, int codim = 0>
std::shared_ptr< MapFunction<dim,dim+codim> >
get_mapping_from_file(const std::string &filename);


/**
 * Returns a MapFunction object (wrapped by a std::shared_ptr) from a Boost XML tree
 * containing exactly one node with the a tag describing a supported mapping.
 * @warning The kind of instantiated MapFunction object depends
 * on the XML tag describing the mapping class.
 * Currently, only the IgFunction is supported (for which the specifying tag is "IgMapping").
 * @note An assertion will be raised (in Debug and Release mode)
 * if no node or more than one node with a tag with a supported mapping are present in XML tree.
 * @ingroup input_v2
 * @author M. Martinelli
 * @date 04 Mar 2014
 */
template <int dim, int codim = 0>
std::shared_ptr< MapFunction<dim,dim+codim> >
get_mapping_from_xml(const boost::property_tree::ptree &tree);



/**
 * Returns a CartesianGrid object (wrapped by a std::shared_ptr) from a Boost XML tree
 * containing exactly one node with the tag "CartesianGrid".
 * @note An assertion will be raised (in Debug and Release mode)
 * if no node or more than one node with the tag "CartesianGrid" are present in XML tree.
 * @ingroup input_v2
 * @author M. Martinelli
 * @date 04 Mar 2014
 */
template <int dim>
std::shared_ptr< CartesianGrid<dim> >
get_cartesian_grid_from_xml(const boost::property_tree::ptree &tree);





/**
 * Returns a BSplineSpace object (wrapped by a std::shared_ptr) from a Boost XML tree
 * containing exactly one node with the tag "BSplineSpace".
 * @note An assertion will be raised (in Debug and Release mode)
 * if no node or more than one node with the tag "BSplineSpace" are present in XML tree.
 * @ingroup input_v2
 * @author M. Martinelli
 * @date 04 Mar 2014
 */
template <int dim, int range, int rank>
std::shared_ptr< BSplineSpace<dim,range,rank> >
get_bspline_space_from_xml(const boost::property_tree::ptree &tree);

#ifdef NURBS
/**
 * Returns a NURBSSpace object (wrapped by a std::shared_ptr) from a Boost XML tree
 * containing exactly one node with the tag "NURBSSpace".
 * @note An assertion will be raised (in Debug and Release mode)
 * if no node or more than one node with the tag "NURBSSpace" are present in XML tree.
 * @ingroup input_v2
 * @author M. Martinelli
 * @date 04 Mar 2014
 */
template <int dim, int range, int rank>
std::shared_ptr< NURBSSpace<dim,range,rank> >
get_nurbs_space_from_xml(const boost::property_tree::ptree &tree);
#endif

IGA_NAMESPACE_CLOSE

#endif // #ifdef READER_H_


