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

#include <xercesc/dom/DOMDocument.hpp>

//
//#undef Assert
//#include <xercesc/parsers/XercesDOMParser.hpp>
//
//#include <igatools/io/xml_parser_error_handler.h>
//#include <xercesc/sax/HandlerBase.hpp>
//#include <xercesc/framework/MemBufInputSource.hpp>
//#include <fstream>
//#include <streambuf>
//#include <sys/stat.h>
//
//using std::string;
using std::shared_ptr;

IGA_NAMESPACE_OPEN


XMLDocument::
XMLDocument(const shared_ptr<xercesc::DOMDocument> dom_doc)
{
}



auto
XMLDocument::
create(const shared_ptr<xercesc::DOMDocument> dom_doc) ->
SelfPtr_
{
    return SelfPtr_ (new XMLDocument (dom_doc));
}

IGA_NAMESPACE_CLOSE

#endif // XML_IO
