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

#include <igatools/io/xml_document.h>

#ifdef XML_IO

#include <igatools/io/xml_element.h>
#include <igatools/io/xml_parser_error_handler.h>

#include <igatools/base/logstream.h>

//#include <igatools/utils/safe_stl_vector.h>


#undef Assert
#include <xercesc/util/XMLString.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMLSSerializer.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/dom/DOMLSOutput.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>


//#include <xercesc/util/PlatformUtils.hpp>
//#include <xercesc/dom/DOM.hpp>

#include <sys/stat.h>

using std::shared_ptr;
using std::string;
using xercesc::XMLString;
using xercesc::DOMNode;
using xercesc::DOMDocument;

IGA_NAMESPACE_OPEN


XMLDocument::
XMLDocument(const char *name)
{
  this->initialize_xml();

  // Creating a new DOM document with a single XML element with the given
  // name.
  try
  {
    // Creating main XML DOM document.
    XMLCh* name_ch = XMLString::transcode(name);
    xml_doc_ = dom_impl_->createDocument(0, name_ch, 0);
    XMLString::release(&name_ch);
#ifndef NDEBUG
    AssertThrow (xml_doc_ != nullptr, ExcNullPtr());
#endif
  }
  catch (xercesc::XMLException &exception)
  {
    char *error_msg = xercesc::XMLString::transcode(exception.getMessage());
    AssertThrow(false, ExcXMLError(error_msg, 0, 0));
    xercesc::XMLString::release(&error_msg);
  }
}



XMLDocument::
XMLDocument(const string &file_path)
{
  this->initialize_xml();

  this->check_file(file_path);

  // Creating and configuring DOM parser
  const auto parser = new xercesc::XercesDOMParser;

  const auto error_handler = XMLParserErrorHandler::create();

  parser->setErrorHandler(error_handler.get());

  parser->setDoNamespaces(true);
  parser->setDoSchema(true);
  parser->setLoadExternalDTD(false);
  parser->setValidationSchemaFullChecking(true);
  parser->setValidationConstraintFatal(true);
  parser->useCachedGrammarInParse(true);
  // parser_->setHandleMultipleImports(true);

  parser->parse(file_path.c_str());
  xml_doc_ = parser->getDocument();

  delete parser;
}



XMLDocument::
XMLDocument(const string &file_path,
            const string &grammar_definition)
{
  this->initialize_xml();

  this->check_file(file_path);

  // Creating and configuring DOM parser
  const auto parser = new xercesc::XercesDOMParser;

  const auto error_handler = XMLParserErrorHandler::create();

  parser->setErrorHandler(error_handler.get());

  parser->setDoNamespaces(true);
  parser->setDoSchema(true);
  parser->setLoadExternalDTD(false);
  parser->setValidationSchemaFullChecking(true);
  parser->setValidationConstraintFatal(true);
  parser->useCachedGrammarInParse(true);
  // parser_->setHandleMultipleImports(true);

  // Activation of the validation.
  parser->setValidationScheme(xercesc::XercesDOMParser::Val_Always);

  // Loading grammar
  xercesc::MemBufInputSource grammar(
    reinterpret_cast<const XMLByte *>(grammar_definition.c_str()),
    grammar_definition.size(), "/igatools_xml_schema.xsd");
  parser->loadGrammar(grammar, xercesc::Grammar::SchemaGrammarType, true);

  // Parsing the file.
  parser->parse(file_path.c_str());

  xml_doc_ = parser->getDocument();

  delete parser;
}



void
XMLDocument::
initialize_xml()
{
  try
  {
    xercesc::XMLPlatformUtils::Initialize();
    XMLCh* core = XMLString::transcode("core");
    dom_impl_ = xercesc::DOMImplementationRegistry::getDOMImplementation(core);
    XMLString::release(&core);
#ifndef NDEBUG
    AssertThrow (dom_impl_ != nullptr, ExcNullPtr());
#endif
  }
  catch (xercesc::XMLException &exception)
  {
    char *error_msg = xercesc::XMLString::transcode(exception.getMessage());
    AssertThrow(false, ExcXMLError(error_msg, 0, 0));
    xercesc::XMLString::release(&error_msg);
  }
}



XMLDocument::
~XMLDocument()
{
  xml_doc_->release();
  delete xml_doc_;
  delete dom_impl_;

  try
  {
    xercesc::XMLPlatformUtils::Terminate();
  }
  catch (xercesc::XMLException &exception)
  {
    char *error_msg = xercesc::XMLString::transcode(exception.getMessage());
    AssertThrow(false, ExcXMLError(error_msg, 0, 0));
    xercesc::XMLString::release(&error_msg);
  }
}



auto
XMLDocument::
create_void_document(const string &name) ->
SelfPtr_
{
  return SelfPtr_(new XMLDocument(name.c_str()));
}



auto
XMLDocument::
parse_from_file(const string &file_path) ->
SelfPtr_
{
  return SelfPtr_(new XMLDocument(file_path));
}



auto
XMLDocument::
parse_from_file(const string &file_path,
                const string &grammar_definition) ->
SelfPtr_
{
  return SelfPtr_(new XMLDocument(file_path, grammar_definition));
}



shared_ptr<XMLElement>
XMLDocument::
get_document_element () const
{
  return XMLElement::create(xml_doc_->getDocumentElement());
}



shared_ptr<XMLElement>
XMLDocument::
create_new_element(const string &name) const
{
  XMLCh* name_ch = XMLString::transcode(name.c_str());
  const auto xml_elem = XMLElement::create(xml_doc_->createElement(name_ch));
  XMLString::release(&name_ch);
  return xml_elem;
}



void
XMLDocument::
write_to_file(const string &file_path) const
{
    // Creating XML writer and writting the DOM document.
    try
    {
        xercesc::DOMLSOutput *xml_stream = dom_impl_->createLSOutput();
        XMLCh* output_file = XMLString::transcode(file_path.c_str());
        xml_stream->setSystemId(output_file);
        XMLString::release(&output_file);

        xercesc::DOMLSSerializer *writer = dom_impl_->createLSSerializer();

        xercesc::DOMConfiguration* dc = writer->getDomConfig();
        dc->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);

        writer->write(xml_doc_, xml_stream);

        writer->release();
        xml_stream->release();
    }
    catch(const xercesc::XMLException &ex)
    {
        char *msg = XMLString::transcode(ex.getMessage());
        AssertThrow(false, ExcXMLError("An Exception occurred when "
                    "parsing file " + file_path + ": " + msg, 0, 0));
        XMLString::release(&msg);
    }
    catch(const xercesc::DOMException &ex)
    {
        char *msg = XMLString::transcode(ex.getMessage());
        AssertThrow(false, ExcXMLError("An Exception occurred when "
                    "parsing file " + file_path + ": " + msg, 0, 0));
        XMLString::release(&msg);
    }
    catch (const xercesc::OutOfMemoryException& ex)
    {
        char *msg = XMLString::transcode(ex.getMessage());
        AssertThrow(false, ExcXMLError("An Exception occurred when "
                    "parsing file " + file_path + ": " + msg, 0, 0));
        XMLString::release(&msg);
    }
    catch (...)
    {
        AssertThrow(false, ExcXMLError("Unknown Exception occurred when "
                    "parsing file " + file_path + ".", 0, 0));
    }
}



void
XMLDocument::
print_info(LogStream &out) const
{
  const auto root_elem = this->get_document_element();
  out.begin_item("XMLDocument:");
  root_elem->print_info(out);
  out.end_item();
}



void
XMLDocument::
check_file(const string &file_path)
{
  // Checking if the file exists.
  errno = 0;
  struct stat buffer;
  if (stat(file_path.c_str(), &buffer) == -1) // == 0 ok; == -1 error
  {
    string error_msg = "Parsing file path " + file_path + " : ";
    if (errno == ENOENT)   // errno declared by include file errno.h
      error_msg += "Path file does not exist, or path is an empty string.";
    else if (errno == ENOTDIR)
      error_msg += "A component of the path is not a directory.";
    else if (errno == ELOOP)
      error_msg += "Too many symbolic links encountered while traversing the path.";
    else if (errno == EACCES)
      error_msg += "Permission denied.";
    else if (errno == ENAMETOOLONG)
      error_msg += "File can not be read";
    else
      error_msg += "An unknown problem was encountered.";
    AssertThrow(false, ExcXMLError(error_msg, 0, 0));
  }
}

IGA_NAMESPACE_CLOSE

#endif // XML_IO
