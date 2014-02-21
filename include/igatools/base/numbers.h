//-+--------------------------------------------------------------------
// Igatools a general purpose Isogeometric analysis library.
// Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

#ifndef __IGA_NUMBERS_H_
#define __IGA_NUMBERS_H_

#include <igatools/base/config.h>
#include <igatools/base/types.h>

IGA_NAMESPACE_OPEN

/**
  @brief Namespace for the declaration of universal constants.

  Since the
  availability in <tt>math.h</tt> is not always guaranteed, we put
  them here. Since this file is included by <tt>base/config.h</tt>,
  they are available to the whole library.

  The constants defined here are a subset of the <tt>M_XXX</tt> constants
  sometimes declared in the system include file <tt>math.h</tt>, but without
  the prefix <tt>M_</tt>.

 */
namespace numbers
{

//TODO(pauletti, Feb 19, 2014): remove commented lines if not used
//static const Real DegToRad = 0.0174532925199433;
//static const Real RadToDeg = 57.295779513082323;
/**
 * e
 */
static const Real  E       = 2.7182818284590452354;

/**
 * log_2 e
 */
static const Real  LOG2E   = 1.4426950408889634074;

/**
 * log_10 e
 */
static const Real  LOG10E  = 0.43429448190325182765;

/**
 * log_e 2
 */
static const Real  LN2     = 0.69314718055994530942;

/**
 * log_e 10
 */
static const Real  LN10    = 2.30258509299404568402;

/**
 * pi
 */
static const Real  PI      = 3.14159265358979323846;

}

IGA_NAMESPACE_CLOSE

#endif /* __IGA_NUMBERS_H_ */
