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

#include <igatools/io/objects_container_parser.h>

#ifdef XML_IO

#include <igatools/io/xml_file_parser.h>
#include <igatools/base/objects_container.h>

#include <xercesc/dom/DOMDocument.hpp>

//#undef Assert // Notice this!!
//#include <xercesc/parsers/XercesDOMParser.hpp>
//#include <xercesc/sax/HandlerBase.hpp>
//#include <xercesc/framework/MemBufInputSource.hpp>
//#include <fstream>
//#include <streambuf>
//
//#include <sys/stat.h>
//
//using namespace xercesc;
using std::string;
using std::shared_ptr;

IGA_NAMESPACE_OPEN


ObjectsContainerParser::
ObjectsContainerParser(const string &file_path)
  :
  file_parser_(XMLFileParser::create(file_path))
{}



shared_ptr<ObjectsContainer>
ObjectsContainerParser::
parse(const string &schema_file) const
{
    const auto xml_document = file_parser_->parse(schema_file);
    const auto container = ObjectsContainer::create();
    return container;
}

IGA_NAMESPACE_CLOSE

#endif // XML_IO
