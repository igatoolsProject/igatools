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

#ifndef __IG_FUNCTION_H
#define __IG_FUNCTION_H

#include <igatools/base/value_types.h>
#include <igatools/functions/function.h>
#include <igatools/functions/ig_coefficients.h>

#include <igatools/basis_functions/spline_space.h>
#include <igatools/linear_algebra/epetra_vector.h>
//#include <igatools/basis_functions/bspline_space.h>
//#include <igatools/basis_functions/nurbs_space.h>


#include <boost/fusion/include/filter_if.hpp>
//#include <boost/fusion/include/iterator.hpp>
#include <boost/fusion/include/tag_of.hpp>
#include <boost/fusion/include/key_of.hpp>
#include <boost/mpl/not_equal_to.hpp>
#include <boost/fusion/include/begin.hpp>
IGA_NAMESPACE_OPEN



template <int,int,int,int>
class PhysicalSpace;


template <int,int,int,int>
class SpaceElementHandler;

template <int,int,int>
class BSplineElementHandler;

template <int,int,int>
class NURBSElementHandler;

template <int,int,int,int>
class PhysSpaceElementHandler;

template <int,int,int,int>
class SpaceElement;

template <int,int,int>
class BSplineElement;

template <int,int,int>
class NURBSElement;

template <int,int,int,int>
class PhysicalSpaceElement;


/**
 *
 * @ingroup serializable
 */
template<int dim,int codim,int range,int rank>
class IgFunction :
  public Function<dim,codim,range,rank>
{
public:

  using CoeffType = IgCoefficients;


private:
  using base_t = Function<dim,codim,range,rank>;
  using parent_t = Function<dim,codim,range,rank>;
  using self_t = IgFunction<dim,codim,range,rank>;
  using Sp = PhysicalSpace<dim,range,rank,codim>;

public:
  //TODO (pauletti, Mar 23, 2015): should we make this private?
  IgFunction(const SharedPtrConstnessHandler<Sp> &space,
             const EpetraTools::Vector &coeff,
             const std::string &property = DofProperties::active);

  IgFunction(const SharedPtrConstnessHandler<Sp> &space,
             const IgCoefficients &coeff,
             const std::string &property = DofProperties::active);



  virtual ~IgFunction() = default;

//  using typename parent_t::topology_variant;
//  using typename parent_t::eval_pts_variant;
  using typename parent_t::Point;
  using typename parent_t::Value;
  using typename parent_t::Gradient;
  using typename parent_t::ElementIterator;
  using typename parent_t::ElementAccessor;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

  using typename parent_t::DomainType;

public:

  std::unique_ptr<typename parent_t::ElementHandler>
  create_cache_handler() const override final;


  static std::shared_ptr<const parent_t>
  const_create(const std::shared_ptr<const Sp> &space,
               const EpetraTools::Vector &coeff,
               const std::string &property = DofProperties::active);

  static std::shared_ptr<const parent_t>
  const_create(const std::shared_ptr<const Sp> &space,
               const IgCoefficients &coeff,
               const std::string &property = DofProperties::active);

  static std::shared_ptr<parent_t>
  create(const std::shared_ptr<Sp> &space,
         const EpetraTools::Vector &coeff,
         const std::string &property = DofProperties::active);

  static std::shared_ptr<parent_t>
  create(const std::shared_ptr<Sp> &space,
         const IgCoefficients &coeff,
         const std::string &property = DofProperties::active);



  std::shared_ptr<const Sp> get_ig_space() const;

  const CoeffType &get_coefficients() const;

  const std::string &get_property() const;

  self_t &operator +=(const self_t &fun);

  virtual void print_info(LogStream &out) const override final;



protected:
  /**
   * Default constructor. It does nothing but it is needed for the
   * serialization mechanism.
   */
  IgFunction() = default;

private:

  SharedPtrConstnessHandler<Sp> space_;

  CoeffType coeff_;

  const std::string property_;

private:

#ifdef MESH_REFINEMENT

  void create_connection_for_insert_knots(std::shared_ptr<self_t> ig_function);

  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid);

#endif // MESH_REFINEMENT

#if 0
#ifdef SERIALIZATION
  /**
   * @name Functions needed for boost::serialization
   * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
   */
  ///@{
  friend class boost::serialization::access;

  template<class Archive>
  void
  serialize(Archive &ar, const unsigned int version);
  ///@}
#endif // SERIALIZATION
#endif
};



IGA_NAMESPACE_CLOSE

#ifdef SERIALIZATION

#include <igatools/functions/ig_function.serial>

#endif // SERIALIZATION

#endif


