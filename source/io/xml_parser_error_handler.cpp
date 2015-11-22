//-+--------------------------------------------------------------------
// Isolde (Isogeometric Solver Environment) is a software for
// analyzing continuum mechanics problems by means of Isogeometric Methods.
// Copyright (C) 2012-2014 by the isolde authors (see authors.txt).
//
// This file is part of the isolde software.
//
// Isolde is property of the University of Pavia and IMATI / CNR,
// Italy. It can not be neither redistributed nor modify without
// an express authorization of the authors.
//
// This software is based on igatools library, that is distributed
// under GNU General Public License version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//-+--------------------------------------------------------------------

#include <igatools/io/xml_parser_error_handler.h>

#if XML_IO

#include <xercesc/sax/HandlerBase.hpp>

using namespace xercesc;

IGA_NAMESPACE_OPEN


auto
XMLParserErrorHandler::
create() -> SelfPtr_
{
  return SelfPtr_(new Self_());
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
