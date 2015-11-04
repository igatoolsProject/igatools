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

#ifndef __IG_FUNCTION_HANDLER_H_
#define __IG_FUNCTION_HANDLER_H_

#include <igatools/functions/function_handler.h>
#include <igatools/functions/ig_function.h>
#include <igatools/basis_functions/space_element_handler.h>


IGA_NAMESPACE_OPEN



/**
 * @ingroup handlers
 */
template<int dim,int codim,int range,int rank>
class IgFunctionHandler :
  public FunctionHandler<dim,codim,range,rank>
{
private:
  using parent_t = FunctionHandler<dim,codim,range,rank>;
  using self_t = IgFunctionHandler<dim,codim,range,rank>;
protected:
//  using typename parent_t::GridType;
public:
  using FunctionType = const IgFunction<dim,codim,range,rank>;
  using typename parent_t::ElementAccessor;
  using typename parent_t::Flags;
  using typename parent_t::topology_variant;
  using typename parent_t::eval_pts_variant;

  IgFunctionHandler(const std::shared_ptr<FunctionType> &ig_function);


  virtual ~IgFunctionHandler() = default;

  void set_flags(const topology_variant &sdim,
                 const Flags &flag) override final;


  void fill_cache(const topology_variant &sdim,
                  ElementAccessor &elem,
                  const int s_id) const override final;

#if 0
  std::shared_ptr<FunctionType>
  get_ig_function() const;
#endif

private:
  struct FillCacheDispatcher : boost::static_visitor<void>
  {
//    template <int order>
//    using _D = typename ElementAccessor::template _D<order>;

    FillCacheDispatcher(const self_t &ig_function_handler,
                        ElementAccessor &ig_function_elem,
                        const int s_id)
      :
      ig_function_handler_(ig_function_handler),
      ig_function_elem_(ig_function_elem),
      s_id_(s_id)
    {}


    template<int sdim>
    void operator()(const Topology<sdim> &sub_elem);


    const self_t &ig_function_handler_;
    ElementAccessor &ig_function_elem_;
    const int s_id_;
  };

//  friend struct FillCacheDispatcher;

  using IgBasisHandler = SpaceElementHandler<dim,codim,range,rank>;
  std::unique_ptr<IgBasisHandler> ig_basis_handler_;

  std::shared_ptr<FunctionType> ig_function_;

};

IGA_NAMESPACE_CLOSE

#endif
