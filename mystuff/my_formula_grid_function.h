
#include <igatools/geometry/formula_grid_function.h>

IGA_NAMESPACE_OPEN

namespace grid_functions
{

//------------------------------------------------------------------------------
/**
 * F(x) = A * x + b
 */
template<int dim, int space_dim>
class MyFormulaGridFunction :
  public FormulaGridFunction<dim,space_dim>
{
  using base_t = GridFunction<dim,space_dim>;
  using parent_t = FormulaGridFunction<dim,space_dim>;
  using self_t = MyFormulaGridFunction<dim,space_dim>;
  using typename base_t::GridType;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  static std::shared_ptr<base_t>
  create(const std::shared_ptr<GridType> &domain) {
  
    auto func = std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<GridType>(domain)));
    return func;
  };

  static std::shared_ptr<const base_t>
  const_create(const std::shared_ptr<const GridType> &domain) { 
    return std::shared_ptr<const self_t>(new self_t(SharedPtrConstnessHandler<GridType>(domain)));
  }

  MyFormulaGridFunction(const self_t &) = default;

  virtual ~MyFormulaGridFunction() = default;

  virtual void print_info(LogStream &out) const override final;

protected:
  MyFormulaGridFunction(
    const SharedPtrConstnessHandler<GridType> &domain)
  :
  parent_t(domain)

{};

private:
  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const {
    auto point = points.begin();
    for (auto &val : values ) {
      val *= 1.0;
      ++point;
    }
 
  };

  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const override;

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const override;

  const Derivative<1> A_;
  const Value    b_;

};


} // of namespace functions.

IGA_NAMESPACE_CLOSE
