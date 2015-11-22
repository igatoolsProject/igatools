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

#include <igatools/io/xml_parser_error_handler.h>

#ifdef XML_IO

#include <xercesc/sax/HandlerBase.hpp>

using namespace xercesc;

IGA_NAMESPACE_OPEN

auto
XMLParserErrorHandler::
create () -> SelfPtr_
{
    return SelfPtr_ (new Self_());
}


void
XMLParserErrorHandler::
warning(const SAXParseException &ex)
{
  char *msg = XMLString::transcode(ex.getMessage());
  AssertThrow (false, ExcXMLWarning (msg, ex.getLineNumber(), ex.getColumnNumber()));
  XMLString::release(&msg);
}



void
XMLParserErrorHandler::
error(const SAXParseException &ex)
{
  char *msg = XMLString::transcode(ex.getMessage());
  AssertThrow (false, ExcXMLError (msg, ex.getLineNumber(), ex.getColumnNumber()));
  XMLString::release(&msg);
}



void
XMLParserErrorHandler::
fatalError(const SAXParseException &ex)
{
  char *msg = XMLString::transcode(ex.getMessage());
  AssertThrow (false, ExcXMLError (msg, ex.getLineNumber(), ex.getColumnNumber()));
  XMLString::release(&msg);
}



void
XMLParserErrorHandler::
resetErrors()
{}


IGA_NAMESPACE_CLOSE

#endif // XML_IO
