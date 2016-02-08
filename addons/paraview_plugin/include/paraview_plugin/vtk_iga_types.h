//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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

namespace paraview_plugin
{

/**
 * Bit field flags for specifying which VTK grid type must be used.
 *
 * The VTK grid can be structured or unstructured. And for unstructured
 * meshes, each of the cells of the mesh can be linear in each direction,
 * or quadratic.
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
enum class VtkGridType : std::int64_t
{
  /// VTK structured grid.
  Structured             = 1 << 0,

  /// VTK unstructured grid with linear cells.
  UnstructuredLinear     = 1 << 1,

  /// VTK unstructured grid with quadratic cells.
  UnstructuredQuadratic  = 1 << 2,

  /// None.
  None                   = 1 << 3,
};



/**
 * @brief Exception for throwing errors that are caught by ParaView.
 *
 * This class derives directly from the standard @p exception class.
 *
 * It just overrides the @ref what method, who retrieves a message
 * explaining the error.
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
class ExcVtkError: public std::exception
{
public:

  /**
   * @brief Constructor.
   * @param message The error message.
   */
  explicit ExcVtkError(const std::string &message):
    error_msg_(message)
  {}

  /**
   * @brief Destructor.
   *
   * This is virtual to allow for subclassing.
   */
  virtual ~ExcVtkError() throw () {};

  /**
   *  @brief Returns a pointer to the (constant) error description.
   *  @return A pointer to a \c const \c char*. The underlying memory
   *          is in possession of the \c ExcVtkError object.
   *          Callers \a must not attempt to free the memory.
   */
  virtual const char *what() const throw ()
  {
    return error_msg_.c_str();
  }

private:
  /// Error message.
  std::string error_msg_;
};



/**
 * @brief Macro to popping up error messages in the ParaView log window.
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
#define VtkIgaErrorMacro(x)                                       \
  {                                                                 \
    if (vtkObject::GetGlobalWarningDisplay())                     \
    {                                                             \
      vtkOStreamWrapper::EndlType endl;                         \
      vtkOStreamWrapper::UseEndl(endl);                         \
      vtkOStrStreamWrapper vtkmsg;                              \
      vtkmsg << "IGATOOLS PLUGIN ERROR!!:\n" x << "\n\n";       \
      vtkOutputWindowDisplayErrorText(vtkmsg.str());            \
      vtkmsg.rdbuf()->freeze(0);                                \
      vtkObject::BreakOnError();                                \
    }                                                             \
  }



/**
 * @brief Macro to popping up warnings in the ParaView log window.
 *
 * @author P. Antolin, 2016.
 *
 * @ingroup paraview_plugin
 */
#define VtkIgaWarningMacro(x)                                     \
  {                                                                 \
    if (vtkObject::GetGlobalWarningDisplay())                     \
    {                                                             \
      vtkOStreamWrapper::EndlType endl;                         \
      vtkOStreamWrapper::UseEndl(endl);                         \
      vtkOStrStreamWrapper vtkmsg;                              \
      vtkmsg << "IGATOOLS PLUGIN WARNING!!:\n" x << "\n\n";     \
      vtkOutputWindowDisplayErrorText(vtkmsg.str());            \
      vtkmsg.rdbuf()->freeze(0);                                \
      vtkObject::BreakOnError();                                \
    }                                                             \
  }

}; // namespace paraview_plugin

IGA_NAMESPACE_CLOSE

#endif // __VTK_IGA_TYPES_H_
