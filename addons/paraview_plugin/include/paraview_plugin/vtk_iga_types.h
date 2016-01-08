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

#ifndef __VTK_IGA_TYPES_H_
#define __VTK_IGA_TYPES_H_

#include <igatools/base/config.h>

IGA_NAMESPACE_OPEN

/**
* Bit field flags for specifying which VTK grid type must be used.
*/
enum class VtkGridType : std::int64_t
{
  /** VTK structured grid */
  Structured             = 1 << 0,

  /** VTK unstructured grid with linear cells */
  UnstructuredLinear     = 1 << 1,

  /** VTK unstructured grid with quadratic cells */
  UnstructuredQuadratic  = 1 << 2,

  /** None */
  None                   = 1 << 3,
};


/**
 * TODO: to ducment
 * When throwing this exception,
 * you can give a message as a
 * <tt>std::string</tt> as argument to the
 * exception that is then
 * displayed. The argument can, of
 * course, be constructed at run-time,
 * for example including the name of a
 * file that can't be opened, or any
 * other text you may want to assemble
 * from different pieces.
 *
 * This is exception is intended to be
 * used as vtk runtime warning to be
 * caught by ParaView
 */
class ExcVtkWarning: public std::exception
{
public:

    /**
     * Constructor
     * @param message The error message.
     */
    explicit ExcVtkWarning(const std::string& message):
      error_msg_(message)
    {}

    /** Destructor.
     * Virtual to allow for subclassing.
     */
    virtual ~ExcVtkWarning() throw (){};

    /**
     *  Returns a pointer to the (constant) error description.
     *  @return A pointer to a \c const \c char*. The underlying memory
     *          is in posession of the \c ExcVtkWarning object. Callers \a must
     *          not attempt to free the memory.
     */
    virtual const char* what() const throw ()
    {
       return error_msg_.c_str();
    }

protected:
    /** Error message.
     */
    std::string error_msg_;
};

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_TYPES_H_
