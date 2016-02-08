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

#include <igatools/functions/formula_function.h>


template <int dim,int codim>
Real
compute_measure(const Domain<dim,codim> &domain)
{
  auto handler = domain.create_cache_handler();

  handler->set_element_flags(domain_element::Flags::w_measure);

  auto elem = domain.begin();
  auto elem_end = domain.end();

  auto quad = QGauss<dim>::const_create(10);

  handler->init_element_cache(elem,quad);

  Real measure = 0.0;
  for (; elem != elem_end ; ++elem)
  {
    handler->fill_element_cache(elem);

    const auto w_meas = elem->get_element_w_measures();

    for (const auto &w :w_meas)
      measure += w;
  }

  return measure;
}



template<int dim, int codim>
class CustomFunction
  : public FormulaFunction<dim,codim,1,1>
{

public:
  using base_t = Function<dim, codim, 1, 1>;
  using parent_t = FormulaFunction<dim, codim, 1, 1>;
  using self_t = CustomFunction<dim, codim>;
  using typename base_t::DomainType;
  using typename parent_t::Point;
  using typename parent_t::Value;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  static std::shared_ptr<self_t>
  create(const std::shared_ptr<DomainType> &domain)
  {
    return std::shared_ptr<self_t>(new
                                   self_t(SharedPtrConstnessHandler<DomainType>(domain)));
  }

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const DomainType> &domain)
  {
    return std::shared_ptr<self_t>(new
                                   self_t(SharedPtrConstnessHandler<DomainType>(domain)));
  }

  CustomFunction(const self_t &) = default;

  virtual ~CustomFunction() = default;

  virtual void print_info(LogStream &out) const override final
  {
    Assert(false,ExcNotImplemented());
  }


protected:
  CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain)
    :
    parent_t(domain)
  {}

private:
  void evaluate_0(const ValueVector<Point> &points,
                  ValueVector<Value> &values) const override
  {
    const int sp_dim = dim+codim;

    int pt = 0;
    for (auto &val : values)
    {
      const auto &point = points[pt];

      val[0] = 0.0;
      for (int i = 0 ; i < sp_dim ; ++i)
        val[0] += point[i];

      ++pt;
    }
  }

  void evaluate_1(const ValueVector<Point> &points,
                  ValueVector<Derivative<1>> &values) const override
  {
    const int sp_dim = dim+codim;

    for (auto &val : values)
    {

      for (int i = 0 ; i < sp_dim ; ++i)
        val[i][0] = 1.0;
    }
  }

  void evaluate_2(const ValueVector<Point> &points,
                  ValueVector<Derivative<2>> &values) const override
  {
    //The hessian is zero.
  }
};



template<int dim, int codim, int range, int rank>
void
function_values(const Function<dim, codim, range, rank> &func)
{
  auto quad = QGauss<dim>::create(2);

  auto func_handler = func.create_cache_handler();

  using Flags = function_element::Flags;
  auto flag = Flags::D0 | Flags::D1 | Flags::D2;
  func_handler->set_element_flags(flag);

  auto elem = func.begin();
  auto end  = func.end();

  func_handler->init_cache(*elem,quad);

  for (; elem != end; ++elem)
  {
    func_handler->fill_element_cache(*elem);

    out.begin_item("FunctionElement:");
    elem->get_index().print_info(out);


    out.begin_item("Values:");
    elem->template get_values_from_cache<function_element::_D<0>, dim>(0).print_info(out);
    out.end_item();

    out.begin_item("Gradients:");
    elem->template get_values_from_cache<function_element::_D<1>, dim>(0).print_info(out);
    out.end_item();

    out.begin_item("Hessians:");
    elem->template get_values_from_cache<function_element::_D<2>, dim>(0).print_info(out);
    out.end_item();

    out.end_item();
  }
#if 0
  const int sdim = dim;
  using Flags = function_element::Flags;
  auto flag = Flags::D0 | Flags::D1;
  auto handler = func->create_cache_handler();

  handler->template set_flags<sdim>(flag);
  auto quad   = QGauss<sdim>::create(2);

  auto elem = func->cbegin();
  auto end  = func->cend();
  handler->init_cache(elem, quad);

  for (; elem != end; ++elem)
  {
    handler->template fill_cache<dim>(elem, 0);
    elem->template get_values_from_cache<function_element::template _D<0>, dim>(0).print_info(out);
    out << endl;
    elem->template get_values_from_cache<function_element::template _D<1>, dim>(0).print_info(out);
    out << endl;
//    elem->template get_values<function_element::_D2, dim>(0).print_info(out);
//    out << endl;
  }
#endif
}
