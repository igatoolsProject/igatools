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

#include <igatools/io/xml_file_parser.h>

#if XML_IO

#include <igatools/io/xml_parser_error_handler.h>

//#include <isolde/base/exceptions.h>
#undef Assert
//
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/HandlerBase.hpp>
//
//#include <sys/stat.h>



using namespace xercesc;
using std::string;
using std::to_string;
using std::shared_ptr;

IGA_NAMESPACE_OPEN


XMLFileParser::
XMLFileParser(const std::shared_ptr<XMLParserErrorHandler> eh,
              const string &file_path)
  :
  file_path_(file_path),
  eh_(eh)
{
//  IsoldeAssert(eh_ != nullptr, ExcNullPtr());

  this->check_file();

  try
  {
    XMLPlatformUtils::Initialize();
  }
  catch (XMLException &exception)
  {
    char *error_msg = XMLString::transcode(exception.getMessage());
//    eh_->error(false, error_msg);
    XMLString::release(&error_msg);
  }

  parser_ = new XercesDOMParser;
}



auto
XMLFileParser::
create(const std::shared_ptr<XMLParserErrorHandler> eh,
       const string &file_path) -> SelfPtr_
{
  return SelfPtr_(new Self_(eh, file_path));
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
//    eh_->error(false, error_msg);
    XMLString::release(&error_msg);
  }
}



void
XMLFileParser::
check_file() const
{
#if 0
  // Checking if the file exists.
  errno = 0;
  struct stat buffer;
  if (stat(file_path_.c_str(), &buffer) == -1) // ==0 ok; ==-1 error
  {
    string error_msg = "Parsing file path " + file_path_ + " : ";
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
    eh_->error(error_msg);
  }
#endif
}



DOMDocument *
XMLFileParser::
parse(void (*load_grammar)(XercesDOMParser *const))
{
  load_grammar(parser_);

  // Configuring DOM parser

//  parser_->setErrorHandler(eh_.get());
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
