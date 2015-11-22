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

#include <igatools/io/xml_file_parser.h>

#ifdef XML_IO

#include <igatools/io/xml_parser_error_handler.h>

#undef Assert // Notice this!!
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/HandlerBase.hpp>

#include <sys/stat.h>

using namespace xercesc;
using std::string;

IGA_NAMESPACE_OPEN


XMLFileParser::
XMLFileParser(const string &file_path)
  :
  file_path_(file_path)
{
  this->check_file(file_path);

  try
  {
    XMLPlatformUtils::Initialize();
  }
  catch (XMLException &exception)
  {
    char *error_msg = XMLString::transcode(exception.getMessage());
    AssertThrow (false, ExcMessage(error_msg));
    XMLString::release(&error_msg);
  }

  parser_ = new XercesDOMParser;
}



auto
XMLFileParser::
create(const string &file_path) -> SelfPtr_
{
  return SelfPtr_(new Self_(file_path));
}



XMLFileParser::
~XMLFileParser()
{
  delete parser_;

  try
  {
    XMLPlatformUtils::Terminate();
  }
  catch (XMLException &exception)
  {
    char *error_msg = XMLString::transcode(exception.getMessage());
    AssertThrow (false, ExcMessage(error_msg));
    XMLString::release(&error_msg);
  }
}



void
XMLFileParser::
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
    AssertThrow (false, ExcXMLError(error_msg, 0, 0));
  }
}



DOMDocument *
XMLFileParser::
parse(void (*load_grammar)(XercesDOMParser *const))
{
  load_grammar(parser_);

  // Configuring DOM parser
  const auto error_handler = XMLParserErrorHandler::create();
  parser_->setErrorHandler(error_handler.get());
  parser_->setValidationScheme(XercesDOMParser::Val_Always);
  parser_->setDoNamespaces(true);
  parser_->setDoSchema(true);
  parser_->setLoadExternalDTD(false);
  parser_->setValidationSchemaFullChecking(true);
  parser_->setValidationConstraintFatal(true);
  parser_->useCachedGrammarInParse(true);
  //    parser_->setHandleMultipleImports(true);

  parser_->parse(file_path_.c_str());

  // no need to free this pointer - owned by the parent parser object
  return parser_->getDocument();
}

IGA_NAMESPACE_CLOSE

#endif // XML_IO
