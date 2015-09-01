//-+--------------------------------------------------------------------
// This file is part of igatools, a general purpose Isogeometric analysis
// library. It was copied from the deal.II project where it is licensed
// under the LGPL (see http://www.dealii.org/).
// It has been modified by the igatools authors to fit the igatools framework.
//-+--------------------------------------------------------------------



#include <igatools/contrib/table_handler.h>
#include <igatools/contrib/table.h>

#include <sstream>
#include <iostream>
#include <iomanip>


IGA_NAMESPACE_OPEN


/*---------------------------------------------------------------------*/

// inline and template functions
namespace internal
{
TableEntry::TableEntry()
{}


double TableEntry::get_numeric_value() const
{
  // we don't quite know the data type in 'value', but
  // it must be one of the ones in the type list of the
  // boost::variant. Go through this list and return
  // the value if this happens to be a number
  //
  // first try with int
  try
  {
    return boost::get<int>(value);
  }
  catch (...)
  {}


  // ... then with int...
  try
  {
    return boost::get<int>(value);
  }
  catch (...)
  {}

  // ...and finally with double precision:
  try
  {
    return boost::get<double>(value);
  }
  catch (...)
  {
    Assert(false, ExcMessage("The number stored by this element of the "
                             "table is not a number."))
  }

  return 0;
}


namespace Local
{
// see which type we can cast to, then use this type to create
// a default constructed object
struct GetDefaultValue : public boost::static_visitor<>
{
  template <typename T>
  void operator()(T &operand) const
  {
    operand = T();
  }
};
}

TableEntry TableEntry::get_default_constructed_copy() const
{
  TableEntry new_entry = *this;
  boost::apply_visitor(Local::GetDefaultValue(), new_entry.value);

  return new_entry;
}


}

/* ------------------------------------------------ */

TableHandler::Column::Column(const std::string &tex_caption)
  :
  tex_caption(tex_caption),
  tex_format("c"),
  precision(4),
  scientific(0),
  flag(0)
{}



TableHandler::Column::Column()
  :
  tex_caption(),
  tex_format("c"),
  precision(4),
  scientific(0),
  flag(0)
{}



void
TableHandler::Column::pad_column_below(const uint size)
{
  // we should never have a column that is completely
  // empty and that needs to be padded
  Assert(entries.size() > 0, ExcInternalError());

  // add as many elements as necessary
  while (entries.size() < size)
    entries.push_back(entries.back().get_default_constructed_copy());
}


/*---------------------------------------------------------------------*/


TableHandler::TableHandler()
  :
  auto_fill_mode(false)
{}



void
TableHandler::set_auto_fill_mode(const bool state)
{
  auto_fill_mode = state;
}


void TableHandler::add_column_to_supercolumn(const std::string &key,
                                             const std::string &superkey)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));

  if (!supercolumns.count(superkey))
  {
    std::pair<std::string, std::vector<std::string> >
    new_column(superkey, std::vector<std::string>());
    supercolumns.insert(new_column);
    // replace key in column_order
    // by superkey
    for (auto column_order_entry : column_order)
      if (column_order_entry==key)
      {
        column_order_entry=superkey;
        break;
      }
  }
  else
  {
    // remove key from column_order
    // for erase we need an iterator
    for (std::vector<std::string>::iterator order_iter=column_order.begin();
         order_iter!=column_order.end(); ++order_iter)
      if (*order_iter==key)
      {
        column_order.erase(order_iter);
        break;
      }
  }

  if (supercolumns.count(superkey))
  {
    supercolumns[superkey].push_back(key);
    // By default set the
    // tex_supercaption to superkey
    std::pair<std::string, std::string> new_tex_supercaption(superkey, superkey);
    tex_supercaptions.insert(new_tex_supercaption);
  }
  else
    Assert(false, ExcInternalError());
}



void TableHandler::set_column_order(const std::vector<std::string> &new_order)
{
  for (auto new_order_entry : new_order)
    Assert(supercolumns.count(new_order_entry) || columns.count(new_order_entry),
           ExcColumnOrSuperColumnNotExistent(new_order_entry));

  column_order=new_order;
}


void TableHandler::set_tex_caption(const std::string &key,
                                   const std::string &tex_caption)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  columns[key].tex_caption=tex_caption;
}



void TableHandler::set_tex_table_caption(const std::string &table_caption)
{
  tex_table_caption=table_caption;
}



void TableHandler::set_tex_table_label(const std::string &table_label)
{
  tex_table_label=table_label;
}



void TableHandler::set_tex_supercaption(const std::string &superkey,
                                        const std::string &tex_supercaption)
{
  Assert(supercolumns.count(superkey), ExcSuperColumnNotExistent(superkey));
  Assert(tex_supercaptions.count(superkey), ExcInternalError());
  tex_supercaptions[superkey]=tex_supercaption;
}



void TableHandler::set_tex_format(const std::string &key,
                                  const std::string &tex_format)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  Assert(tex_format=="l" || tex_format=="c" || tex_format=="r",
         ExcUndefinedTexFormat(tex_format));
  columns[key].tex_format=tex_format;
}



void TableHandler::set_precision(const std::string &key,
                                 const int precision)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  columns[key].precision=precision;
}


void TableHandler::set_scientific(const std::string &key,
                                  const bool scientific)
{
  Assert(columns.count(key), ExcColumnNotExistent(key));
  columns[key].scientific=scientific;
}


void TableHandler::write_text(std::ostream &out,
                              const TextOutputFormat format) const
{
  AssertThrow(out, ExcIO());

  // first pad the table from below if necessary
  if (auto_fill_mode == true)
  {
    uint max_rows = 0;
    for (std::map<std::string, Column>::const_iterator p = columns.begin();
         p != columns.end(); ++p)
      max_rows = std::max<uint>(max_rows, p->second.entries.size());

    for (std::map<std::string, Column>::iterator p = columns.begin();
         p != columns.end(); ++p)
      p->second.pad_column_below(max_rows);
  }

  std::vector<std::string> sel_columns;
  get_selected_columns(sel_columns);

  const int nrows  = n_rows();
  const int n_cols = sel_columns.size();

  // first compute the widths of each
  // entry of the table, in order to
  // have a nicer alignement
  Table<2,uint> entry_widths(nrows, n_cols);
  for (int i=0; i<nrows; ++i)
    for (int j=0; j<n_cols; ++j)
    {
      // get key and entry here
      std::string key = sel_columns[j];
      const std::map<std::string, Column>::const_iterator
      col_iter = columns.find(key);
      Assert(col_iter!=columns.end(), ExcInternalError());

      const Column &column = col_iter->second;

      // write it into a dummy
      // stream, just to get its
      // size upon output
      std::ostringstream dummy_out;

      dummy_out << std::setprecision(column.precision);

      if (col_iter->second.scientific)
        dummy_out.setf(std::ios::scientific, std::ios::floatfield);
      else
        dummy_out.setf(std::ios::fixed, std::ios::floatfield);

      dummy_out << column.entries[i].value;

      // get size. as a side note, if the
      // text printed is in fact the empty
      // string, then we get into a bit of
      // trouble since IgReaders would skip
      // over the resulting whitespaces. as
      // a consequence, we'll print ""
      // instead in that case.
      const uint size = dummy_out.str().length();
      if (size > 0)
        entry_widths[i][j] = size;
      else
        entry_widths[i][j] = 2;
    }

  // next compute the width each row
  // has to have to suit all entries
  std::vector<uint> column_widths(n_cols, 0);
  for (int i=0; i<nrows; ++i)
    for (int j=0; j<n_cols; ++j)
      column_widths[j] = std::max(entry_widths[i][j],
                                  column_widths[j]);

  // write the captions
  for (std::size_t j=0; j<column_order.size(); ++j)
  {
    const std::string &key = column_order[j];

    switch (format)
    {
      case table_with_headers:
      {
        // if the key of this column is
        // wider than the widest entry,
        // then adjust
        if (key.length() > column_widths[j])
          column_widths[j] = key.length();

        // now write key. try to center
        // it somehow
        const int front_padding = (column_widths[j]-key.length())/2,
                  rear_padding  = (column_widths[j]-key.length()) -
                                  front_padding;
        for (int i=0; i<front_padding; ++i)
          out << ' ';
        out << key;
        for (int i=0; i<rear_padding; ++i)
          out << ' ';

        // finally column break
        out << ' ';

        break;
      }

      case table_with_separate_column_description:
      {
        // print column key with column number. enumerate
        // columns starting with 1
        out << "# " << j+1 << ": " << key << std::endl;
        break;
      }

      default:
        Assert(false, ExcInternalError());
    }
  }
  if (format == table_with_headers)
    out << std::endl;

  for (int i=0; i<nrows; ++i)
  {
    for (int j=0; j<n_cols; ++j)
    {
      std::string key=sel_columns[j];
      const std::map<std::string, Column>::const_iterator
      col_iter=columns.find(key);
      Assert(col_iter!=columns.end(), ExcInternalError());

      const Column &column=col_iter->second;

      out << std::setprecision(column.precision);

      if (col_iter->second.scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);
      out << std::setw(column_widths[j]);

      // get the string to write into the
      // table. we could simply << it
      // into the stream but we have to
      // be a bit careful about the case
      // where the string may be empty,
      // in which case we'll print it as
      // "". note that ultimately we
      // still just << it into the stream
      // since we want to use the flags
      // of that stream, and the first
      // test is only used to determine
      // whether the size of the string
      // representation is zero
      {
        std::ostringstream text;
        text << column.entries[i].value;
        if (text.str().size() > 0)
          out << column.entries[i].value;
        else
          out << "\"\"";
      }

      // pad after this column
      out << " ";
    }
    out << std::endl;
  }
}


void TableHandler::write_tex(std::ostream &out, const bool with_header) const
{
  AssertThrow(out, ExcIO());
  if (with_header)
    out << "\\documentclass[10pt]{report}" << std::endl
        << "\\usepackage{float}" << std::endl << std::endl << std::endl
        << "\\begin{document}" << std::endl;

  out << "\\begin{table}[H]" << std::endl
      << "\\begin{center}" << std::endl
      << "\\begin{tabular}{|";

  // first pad the table from below if necessary
  if (auto_fill_mode == true)
  {
    uint max_rows = 0;
    for (std::map<std::string, Column>::const_iterator p = columns.begin();
         p != columns.end(); ++p)
      max_rows = std::max<uint>(max_rows, p->second.entries.size());

    for (std::map<std::string, Column>::iterator p = columns.begin();
         p != columns.end(); ++p)
      p->second.pad_column_below(max_rows);
  }

  std::vector<std::string> sel_columns;
  get_selected_columns(sel_columns);

  // write the column formats
  for (auto key : column_order)
  {
    // avoid `supercolumns[key]'
    const std::map<std::string, std::vector<std::string> >::const_iterator
    super_iter=supercolumns.find(key);

    if (super_iter!=supercolumns.end())
    {
      const int n_subcolumns=super_iter->second.size();
      for (int k=0; k<n_subcolumns; ++k)
      {
        // avoid `columns[supercolumns[key]]'
        const std::map<std::string, Column>::const_iterator
        col_iter=columns.find(super_iter->second[k]);
        Assert(col_iter!=columns.end(), ExcInternalError());

        out << col_iter->second.tex_format << "|";
      }
    }
    else
    {
      // avoid `columns[key]';
      const std::map<std::string, Column>::const_iterator
      col_iter=columns.find(key);
      Assert(col_iter!=columns.end(), ExcInternalError());
      out << col_iter->second.tex_format << "|";
    }
  }
  out << "} \\hline" << std::endl;

  // write the caption line of the table

  for (size_t j=0; j<column_order.size(); ++j)
  {
    std::string key=column_order[j];
    const std::map<std::string, std::vector<std::string> >::const_iterator
    super_iter=supercolumns.find(key);

    if (super_iter!=supercolumns.end())
    {
      const int n_subcolumns=super_iter->second.size();
      // avoid use of `tex_supercaptions[key]'
      std::map<std::string,std::string>::const_iterator
      tex_super_cap_iter=tex_supercaptions.find(key);
      out << std::endl << "\\multicolumn{" << n_subcolumns << "}{|c|}{"
          << tex_super_cap_iter->second << "}";
    }
    else
    {
      // col_iter->second=columns[col];
      const std::map<std::string, Column>::const_iterator
      col_iter=columns.find(key);
      Assert(col_iter!=columns.end(), ExcInternalError());
      out << col_iter->second.tex_caption;
    }
    if (j<column_order.size()-1)
      out << " & ";
  }
  out << "\\\\ \\hline" << std::endl;

  // write the n rows
  const int nrows=n_rows();
  for (int i=0; i<nrows; ++i)
  {
    const int n_cols=sel_columns.size();

    for (int j=0; j<n_cols; ++j)
    {
      std::string key=sel_columns[j];
      // avoid `column[key]'
      const std::map<std::string, Column>::const_iterator
      col_iter=columns.find(key);
      Assert(col_iter!=columns.end(), ExcInternalError());

      const Column &column=col_iter->second;

      out << std::setprecision(column.precision);

      if (col_iter->second.scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

      out << column.entries[i].value;

      if (j<n_cols-1)
        out << " & ";
    }
    out << "\\\\ \\hline" << std::endl;
  }

  out   << "\\end{tabular}" << std::endl
        << "\\end{center}" << std::endl;
  if (tex_table_caption!="")
    out << "\\caption{"  << tex_table_caption << "}" << std::endl;
  if (tex_table_label!="")
    out << "\\label{"   << tex_table_label << "}" << std::endl;
  out   << "\\end{table}" << std::endl;
  if (with_header)
    out << "\\end{document}" << std::endl;
}


int TableHandler::n_rows() const
{
  if (columns.size() == 0)
    return 0;

  std::map<std::string, Column>::const_iterator col_iter = columns.begin();
  auto n = col_iter->second.entries.size();
  std::string first_name=col_iter->first;

  for (++col_iter; col_iter!=columns.end(); ++col_iter)
    Assert(col_iter->second.entries.size()==n,
           ExcWrongNumberOfDataEntries(col_iter->first,
                                       col_iter->second.entries.size(),
                                       first_name, n));

  return n;
}


void TableHandler::get_selected_columns(std::vector<std::string> &sel_columns) const
{
  sel_columns.clear();

  for (auto key : column_order)
  {
    const std::map<std::string, std::vector<std::string> >::const_iterator
    super_iter=supercolumns.find(key);

    if (super_iter!=supercolumns.end())
    {
      // i.e. key is a supercolumn key
      for (const auto subkey : super_iter->second)
      {
        Assert(columns.count(subkey), ExcInternalError());
        sel_columns.push_back(subkey);
      }
    }
    else
    {
      Assert(columns.count(key), ExcInternalError());
      // i.e. key is a column key
      sel_columns.push_back(key);
    }
  }
}


IGA_NAMESPACE_CLOSE
