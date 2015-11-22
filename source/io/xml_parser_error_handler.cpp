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

//#include <isolde/base/exceptions.h>
//#include <isolde/io/input_tree_log.h>
//#include <isolde/io/log_manager.h>

#include <xercesc/sax/HandlerBase.hpp>


using namespace iga;
using namespace xercesc;
using std::string;
using std::endl;
using std::to_string;
using std::shared_ptr;

IGA_NAMESPACE_OPEN


//XMLFileParserHandler::
//XMLFileParserHandler(const shared_ptr<InputTreeLog> log,
//                   const std::shared_ptr<LogManager> log_manager,
//                   const string &file_path)
//  :
//  logs_(1, log),
//  log_manager_(log_manager),
//  file_paths_(1, file_path)
//{
////  IsoldeAssert(log != nullptr, ExcNullPtr());
////  IsoldeAssert(log_manager_ != nullptr, ExcNullPtr());
//}

XMLFileParserHandler::
XMLFileParserHandler(const std::shared_ptr<LogManager> log_manager,
                    const string &file_path)
  :
//  logs_(1, log),
  log_manager_(log_manager),
  file_paths_(1, file_path)
{
//  IsoldeAssert(log != nullptr, ExcNullPtr());
//  IsoldeAssert(log_manager_ != nullptr, ExcNullPtr());
}
//
//
//
//auto
//XMLFileParserHandler::
//create(const shared_ptr<InputTreeLog> log,
//       const std::shared_ptr<LogManager> log_manager,
//       const string &file_path) -> SelfPtr_
//{
//  return SelfPtr_(new Self_(log, log_manager, file_path));
//}



auto
XMLFileParserHandler::
create(const std::shared_ptr<LogManager> log_manager,
       const string &file_path) -> SelfPtr_
{
  return SelfPtr_(new Self_(log_manager, file_path));
}


void
XMLFileParserHandler::
warning(const SAXParseException &ex)
{
  char *msg = XMLString::transcode(ex.getMessage());
  this->error(msg, ex.getLineNumber(), ex.getColumnNumber());
  XMLString::release(&msg);
}



void
XMLFileParserHandler::
error(const SAXParseException &ex)
{
  char *msg = XMLString::transcode(ex.getMessage());
  this->error(msg, ex.getLineNumber(), ex.getColumnNumber());
  XMLString::release(&msg);
}



void
XMLFileParserHandler::
fatalError(const SAXParseException &ex)
{
  char *msg = XMLString::transcode(ex.getMessage());
  this->error(msg, ex.getLineNumber(), ex.getColumnNumber());
  XMLString::release(&msg);
}



void
XMLFileParserHandler::
resetErrors()
{}



void
XMLFileParserHandler::
error(const bool cond,
      const string &msg,
      const int &line_number,
      const int &column_number) const
{
//  IsoldeAssert(logs_.size() == file_paths_.size(),
//               ExcDimensionMismatch(logs_.size(), file_paths_.size()));

  const string line_str = line_number == -1 ?
                          "unknown" : to_string(line_number);
  const string col_str = column_number == -1 ?
                         "unknown" : to_string(column_number);

//  auto logs_it = logs_.crbegin();
//  auto files_it = file_paths_.crbegin();
//  auto files_end = file_paths_.crend();
//
//  string error_msg =
//    "An error occurred parsing file: \"" + *files_it +
//    "\" at line " + line_str  + " column " + col_str + ":\n\t" +
//    (*logs_it)->get_tree_as_string("In XML tree tag ", ":\n\n");
//
//  ++files_it;
//  ++logs_it;
//
//  for (; files_it != files_end; ++files_it, ++logs_it)
//    error_msg +=
//      "Called from file: \"" + *files_it + ":\n\t" +
//      (*logs_it)->get_tree_as_string("In XML tree tag ", ":\n\n");
//
//  if (!cond)
//    log_manager_->throw_error(error_msg + "\t" + msg);

}




void
XMLFileParserHandler::
error(const string &msg) const
{
  this->error(false, msg);
}



void
XMLFileParserHandler::
error(const string &msg,
      const int &line_number,
      const int &column_number) const
{
  this->error(false, msg, line_number, column_number);
}



void
XMLFileParserHandler::
add_id(const Index &id)
{
//  const auto &log = logs_.back();
//  log->add_id(id);
}



void
XMLFileParserHandler::
add_tag_name(const string &tag_name)
{
//  const auto &log = logs_.back();
//  log->add_tag_name(tag_name);
}




void
XMLFileParserHandler::
pop_tag_name()
{
//  const auto &log = logs_.back();
//  log->pop_tag_name();
}



void
XMLFileParserHandler::
add_file(const string &file_path, const string &tag_name)
{
//  IsoldeAssert(logs_.size() == file_paths_.size(),
//               ExcDimensionMismatch(logs_.size(), file_paths_.size()));
  file_paths_.push_back(file_path);
//  logs_.push_back(InputTreeLog::create(tag_name));
}



void
XMLFileParserHandler::
pop_file()
{
//  IsoldeAssert(logs_.size() == file_paths_.size(),
//               ExcDimensionMismatch(logs_.size(), file_paths_.size()));
  file_paths_.pop_back();
//  logs_.pop_back();
}



void
XMLFileParserHandler::
print_info(LogStream &out) const
{
//  IsoldeAssert(logs_.size() == file_paths_.size(),
//               ExcDimensionMismatch(logs_.size(), file_paths_.size()));

//  Index level_id = 0;
//  auto logs_it = logs_.cbegin();
//  auto files_it = file_paths_.cbegin();
//  auto files_end = file_paths_.cend();
//  for (; files_it != files_end; ++files_it, ++logs_it, ++level_id)
//  {
//    out.begin_item("Level " + to_string(level_id));
//    out << "Input file path:  " << *files_it << endl;
//    (*logs_it)->print_info(out);
//    out.end_item();
//  }
}


IGA_NAMESPACE_CLOSE

#endif // XML_IO
