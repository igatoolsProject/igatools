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

#include <igatools/io/xml_element.h>

#ifdef XML_IO

#include <igatools/utils/safe_stl_vector.h>

#include <xercesc/util/XMLString.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMLSSerializer.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>

using std::string;
using std::shared_ptr;
using std::unique_ptr;
using xercesc::XMLString;
using xercesc::DOMElement;
using xercesc::DOMNode;
using xercesc::DOMNodeList;
using xercesc::DOMText;

IGA_NAMESPACE_OPEN


XMLElement::
XMLElement(const DOMElemPtr_ dom_elem)
  :
  root_elem_(dom_elem)
{
  Assert(root_elem_ != nullptr, ExcNullPtr());
}



auto
XMLElement::
create(const DOMElemPtr_ dom_elem) ->
SelfPtr_
{
  return SelfPtr_(new XMLElement(dom_elem));
}



auto
XMLElement::
get_children_elements() const ->
SafeSTLVector<SelfPtr_>
{
  DOMNodeList *elems = root_elem_->getChildNodes();

  const Size n_children = elems->getLength();

  SafeSTLVector<SelfPtr_> children;
  for (int i = 0; i < n_children; ++i)
  {
    DOMNode *n = elems->item(i);

    if (n->getNodeType() && // true is not NULL
    n->getNodeType() == DOMNode::ELEMENT_NODE)  // is element
    {
      const auto elem_ptr = dynamic_cast<DOMElemPtr_>(n);
      Assert(elem_ptr != nullptr, ExcNullPtr());
      children.push_back(Self_::create(elem_ptr));
    }
  }

  return children;
}



auto
XMLElement::
get_children_elements(const string &name) const ->
SafeSTLVector<SelfPtr_>
{
  SafeSTLVector<SelfPtr_> children;

  DOMNodeList *elems = root_elem_->getChildNodes();
  const Size n_children = elems->getLength();

  for (int i = 0; i < n_children; ++i)
  {
    DOMNode *n = elems->item(i);
    if (n->getNodeType() && // true is not NULL
    XMLString::transcode(n->getNodeName()) == name &&
    n->getNodeType() == DOMNode::ELEMENT_NODE)  // is element
    {
      const auto elem_ptr = dynamic_cast<DOMElemPtr_>(n);
      Assert(elem_ptr != nullptr, ExcNullPtr());
      children.push_back(Self_::create(elem_ptr));
    }
  }

  return children;
}



string
XMLElement::
get_name() const
{
  return XMLString::transcode(root_elem_->getNodeName());
}



bool
XMLElement::
has_attribute(const string &name) const
{
  return root_elem_->hasAttribute(XMLString::transcode(name.c_str()));
}



bool
XMLElement::
has_element(const string &name) const
{
  for (const auto &el : this->get_children_elements())
  {
    const auto el_name = el->get_name();
    if (el_name == name)
      return true;
  }
  return false;
}



DOMText *
XMLElement::
get_single_text_element() const
{
  DOMNodeList *elems = root_elem_->getChildNodes();
  const Size n_children = elems->getLength();

  DOMText *element;

  for (int i = 0; i < n_children; ++i)
  {
    DOMNode *n = elems->item(i);

    if (n->getNodeType() && // true is not NULL
        n->getNodeType() == DOMNode::TEXT_NODE)  // is element
    {
      element = dynamic_cast<DOMText *>(n);
      break;
    }
  }

  // if there is more than, or less than, one element, an error is thrown.
#ifndef NDEBUG
  Size n_elems = 0;
  for (int i = 0; i < n_children; ++i)
  {
    DOMNode *n = elems->item(i);

    if (n->getNodeType() && // true is not NULL
        n->getNodeType() == DOMNode::TEXT_NODE) // is element
      ++n_elems;
  }
  Assert(n_elems == 1, ExcDimensionMismatch(n_elems, 1));
#endif

  return element;
}



#if 0
template <class T>
T
XMLElement::
get_value() const
{
  try
  {
    const auto text_elem = this->get_single_text_element();
    return XMLString::transcode(text_elem->getWholeText());
  }
  catch (...)
  {
    AssertThrow(false, ExcMessage("Error parsing text element."));
    return "Error"; // For avoiding the warning.
  }
}
#endif



template <>
string
XMLElement::
get_value<string> () const
{
  try
  {
    const auto text_elem = this->get_single_text_element();
    string s = XMLString::transcode(text_elem->getWholeText());
    // Trimming (removing spaces, tabs, jumps ...) both ends of the string.
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
  }
  catch (...)
  {
    AssertThrow(false, ExcMessage("Error parsing text element."));
    return "Error"; // For avoiding the warning.
  }
}



template <>
Index
XMLElement::
get_attribute<Index> (const string &name) const
{
  Assert(this->has_attribute(name), ExcMessage("Attribute not present."));
  const string str = XMLString::transcode(
                       root_elem_->getAttribute(XMLString::transcode(name.c_str())));

  try
  {
    return std::stoul(str);
  }
  catch (...)
  {
    AssertThrow(false, ExcMessage("Value has not type Index."));
    return 0;
  }

}



template <>
Real
XMLElement::
get_attribute<Real> (const string &name) const
{
  Assert(this->has_attribute(name), ExcMessage("Attribute not present."));
  const string str = XMLString::transcode(
                       root_elem_->getAttribute(XMLString::transcode(name.c_str())));
  try
  {
    return std::stod(str);
  }
  catch (...)
  {
    AssertThrow(false, ExcMessage("Value has not type Real."));
    return 0.0;
  }
}



template <>
string
XMLElement::
get_attribute<string> (const string &name) const
{
  Assert(this->has_attribute(name), ExcMessage("Attribute not present."));
  return XMLString::transcode(
           root_elem_->getAttribute(XMLString::transcode(name.c_str())));
}



template <>
bool
XMLElement::
get_attribute<bool> (const string &name) const
{
  Assert(this->has_attribute(name), ExcMessage("Attribute not present."));
  string str = XMLString::transcode(
                 root_elem_->getAttribute(XMLString::transcode(name.c_str())));

  try
  {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
  }
  catch (...)
  {
    AssertThrow(false, ExcMessage("Value has not type boolean."));
    return false;
  }
}



template <>
SafeSTLVector<Real>
XMLElement::
get_values_vector() const
{
  SafeSTLVector<Real> data;
  const auto text_elem = this->get_single_text_element();

  const string str = XMLString::transcode(text_elem->getWholeText());

  try
  {
    Real v;
    std::stringstream line_stream(str);
    while (line_stream >> v)
      data.push_back(v);

  }
  catch (...)
  {
    AssertThrow(false, ExcMessage("Impossible to parse vector."));
  }

  return data;
}



template <>
SafeSTLVector<Index>
XMLElement::
get_values_vector() const
{
  SafeSTLVector<Index> data;
  const auto text_elem = this->get_single_text_element();

  const string str = XMLString::transcode(text_elem->getWholeText());

  try
  {
    Index v;
    std::stringstream line_stream(str);
    while (line_stream >> v)
      data.push_back(v);

  }
  catch (...)
  {
    AssertThrow(false, ExcMessage("Impossible to parse vector."));
  }

  return data;
}



template <>
SafeSTLVector<bool>
XMLElement::
get_values_vector() const
{
  SafeSTLVector<bool> data;
  const auto text_elem = this->get_single_text_element();

  const string str = XMLString::transcode(text_elem->getWholeText());

  try
  {
    bool v;
    std::stringstream line_stream(str);
    while (line_stream  >> std::boolalpha >> v)
      data.push_back(v);
  }
  catch (...)
  {
    AssertThrow(false, ExcMessage("Impossible to parse vector."));
  }

  return data;
}



template <>
SafeSTLVector<string>
XMLElement::
get_values_vector() const
{
  SafeSTLVector<string> data;
  const auto text_elem = this->get_single_text_element();

  const string str = XMLString::transcode(text_elem->getWholeText());

  try
  {
    string v;
    std::stringstream line_stream(str);
    while (line_stream  >> v)
      data.push_back(v);
  }
  catch (...)
  {
    AssertThrow(false, ExcMessage("Impossible to parse vector."));
  }

  return data;
}



auto
XMLElement::
get_single_element() -> SelfPtr_
{
  DOMNodeList *children = root_elem_->getChildNodes();

  const Size n_children = children->getLength();

  SelfPtr_ element;

  for (int i = 0; i < n_children; ++i)
  {
    DOMNode *n = children->item(i);

    if (n->getNodeType() && // true is not NULL
    n->getNodeType() == DOMNode::ELEMENT_NODE)  // is element
    {
      element = Self_::create(dynamic_cast<DOMElemPtr_>(n));
      break;
    }
  }

  // if there is more than one element, an error is thrown.
#ifndef NDEBUG
  Size n_children_elems = 0;
  for (int i = 0; i < n_children; ++i)
  {
    auto *n = children->item(i);

    if (n->getNodeType() && // true is not NULL
        n->getNodeType() == DOMNode::ELEMENT_NODE)  // is element
      ++n_children_elems;
  }
  Assert(n_children_elems == 1, ExcDimensionMismatch(n_children_elems, 1));
  Assert(element != nullptr, ExcNullPtr());
#endif

  return element;
}



auto
XMLElement::
get_single_element(const string &name) -> SelfPtr_
{
  // Getting all the elements with the given name.
  auto *children = root_elem_->getChildNodes();
  const Size n_children = children->getLength();

  Size n_matching_childs = 0;
  SelfPtr_ element;
  for (int c = 0; c < n_children; ++c)
  {
    auto *elem = dynamic_cast<DOMElemPtr_>(children->item(c));
    if (elem != nullptr && XMLString::transcode(elem->getNodeName()) == name)
    {
      element = Self_::create(elem);
      ++n_matching_childs;
    }
  }

  Assert(n_matching_childs == 1, ExcDimensionMismatch(n_matching_childs, 1));
  Assert(element != nullptr, ExcNullPtr());

  return element;
}



template <class T>
void
XMLElement::
add_attribute(const string &name,
              const T &value)
{
    XMLCh* name_ch = XMLString::transcode(name.c_str());
    const auto value_str = std::to_string(value);
    XMLCh* value_ch = XMLString::transcode(value_str.c_str());
    root_elem_->setAttribute(name_ch, value_ch);
    XMLString::release(&name_ch);
    XMLString::release(&value_ch);
}



void
XMLElement::
add_attribute(const string &name,
              const string &value)
{
    XMLCh* name_ch = XMLString::transcode(name.c_str());
    XMLCh* value_ch = XMLString::transcode(value.c_str());
    root_elem_->setAttribute(name_ch, value_ch);
    XMLString::release(&name_ch);
    XMLString::release(&value_ch);
}



void
XMLElement::
add_attribute(const char *name,
              const char *value)
{
    XMLCh* name_ch = XMLString::transcode(name);
    XMLCh* value_ch = XMLString::transcode(value);
    root_elem_->setAttribute(name_ch, value_ch);
    XMLString::release(&name_ch);
    XMLString::release(&value_ch);
}



void
XMLElement::
append_child_element (const SelfPtr_ xml_elem)
{
    root_elem_->appendChild(xml_elem->root_elem_);
}



void
XMLElement::
print_info(LogStream &out) const
{

  // Creating XML writer and writing the DOM element.
  try
  {
      xercesc::DOMImplementation *impl = xercesc::
              DOMImplementationRegistry::getDOMImplementation(XMLString::transcode("LS"));
      xercesc::DOMLSSerializer *writer = ((xercesc::DOMImplementationLS *)impl)
                  ->createLSSerializer();

      const auto *xmlch_output = writer->writeToString(root_elem_);
      const auto output_string = XMLString::transcode(xmlch_output);

      out.begin_item("XMLElement:");
      out << output_string;
      out.end_item();

      delete xmlch_output;
      delete writer;
  }
  catch(const xercesc::XMLException &ex)
  {
      char *msg = XMLString::transcode(ex.getMessage());
      AssertThrow(false, ExcXMLError("An Exception occurred when "
              + string("writing element: ") + msg, 0, 0));
      XMLString::release(&msg);
  }
  catch(const xercesc::DOMException &ex)
  {
      char *msg = XMLString::transcode(ex.getMessage());
      AssertThrow(false, ExcXMLError("An Exception occurred when "
              + string("writing element: ") + msg, 0, 0));
      XMLString::release(&msg);
  }
  catch (const xercesc::OutOfMemoryException& ex)
  {
      char *msg = XMLString::transcode(ex.getMessage());
      AssertThrow(false, ExcXMLError("An Exception occurred when "
              + string("writing element: ") + msg, 0, 0));
      XMLString::release(&msg);
  }
  catch (...)
  {
      AssertThrow(false, ExcXMLError("Unknown Exception occurred when "
              "writing element.", 0, 0));
  }
}

template void XMLElement::add_attribute<Index> (const string &, const Index &);
template void XMLElement::add_attribute<Real> (const string &, const Real &);

IGA_NAMESPACE_CLOSE

#endif // XML_IO
