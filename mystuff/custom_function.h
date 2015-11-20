
#include <igatools/functions/formula_function.h>

IGA_NAMESPACE_OPEN

namespace functions
{

template<int dim, int space_dim>
class CustomFunction : public FormulaFunction<dim,0,space_dim> {
  // some useful aliases
  using base_t   = Function<dim,0,space_dim>;
  using parent_t = FormulaFunction<dim,0,space_dim>;
  using self_t   = CustomFunction<dim,space_dim>;
  using typename base_t::DomainType;
public:
  using typename parent_t::Value;
  using typename parent_t::Point;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  // defaulted public constructor and destructor
  CustomFunction(const self_t &) = default;
  virtual ~CustomFunction() = default;

  // creators
  static std::shared_ptr<self_t>
  create(const std::shared_ptr<DomainType> &domain);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const DomainType> & domain,
               Value (*f_D0)(const Point));

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const DomainType> & domain,
               Value (*f_D0)(const Point),
               Derivative<1> (*f_D1)(const Point));

  // info printer method
  virtual void print_info(LogStream &out) const override final;

protected:
  CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain);
  CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain,
                     Value (*f_D0)(const Point));
  CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain,
                     Value (*f_D0)(const Point),
                     Derivative<1> (*f_D1)(const Point));

public:
  // function's functions (thanks anglosaxons for this this beautifully ambiguous expression)
  //std::array<Value (*)(const Point),dim> funct;
  Value (*funct_D0)(const Point);
  Derivative<1> (*funct_D1)(const Point);
  //std::array<array<void (*)(double),dim>,dim> grads;
  void test_custom_function(const Point x) {
    //std::cout << " function value is " << funct[0](x) << std::endl;
  }

private:
  void evaluate_0(const ValueVector<Point> &points, ValueVector<Value> &values) const;
  void evaluate_1(const ValueVector<Point> &points, ValueVector<Derivative<1>> &values) const;
  void evaluate_2(const ValueVector<Point> &points, ValueVector<Derivative<2>> &values) const;

#ifdef MESH_REFINEMENT
  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Domain<dim> &old_domain) {};
#endif // MESH_REFINEMENT
};

// ----------------------------------------------------------------------------
//  CONSTRUCTOR
// ----------------------------------------------------------------------------
template<int dim, int space_dim>
CustomFunction<dim,space_dim>::CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain) 
  :
  parent_t(domain)
{};

template<int dim, int space_dim>
CustomFunction<dim,space_dim>::CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain,
                                                      Value (*f_D0)(const Point))
  :
  parent_t(domain)
{funct_D0 = f_D0;};

template<int dim, int space_dim>
CustomFunction<dim,space_dim>::CustomFunction(const SharedPtrConstnessHandler<DomainType> &domain,
                                                      Value (*f_D0)(const Point),
                                                      Derivative<1> (*f_D1)(const Point))
  :
  parent_t(domain)
{funct_D0 = f_D0; funct_D1 = f_D1;};

// ----------------------------------------------------------------------------
//  CREATORS
// ----------------------------------------------------------------------------
template<int dim, int space_dim> // non const creator
auto CustomFunction<dim,space_dim>::create(const std::shared_ptr<DomainType> &domain) -> std::shared_ptr<self_t> {
  auto func = std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<DomainType>(domain)));
#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif // MESH_REFINEMENT
  return func;
};

template<int dim, int space_dim> // const creator with function
auto CustomFunction<dim,space_dim>::const_create(const std::shared_ptr<const DomainType> &domain,
                                                     Value (*f_D0)(const Point)) -> std::shared_ptr<const self_t> {
  return std::shared_ptr<const self_t>(new self_t(SharedPtrConstnessHandler<DomainType>(domain),f_D0));
};

template<int dim, int space_dim> // const creator with function and derivatives
auto CustomFunction<dim,space_dim>::const_create(const std::shared_ptr<const DomainType> &domain,
                                                     Value (*f_D0)(const Point),
                                                     Derivative<1> (*f_D1)(const Point)) -> std::shared_ptr<const self_t> {
  return std::shared_ptr<const self_t>(new self_t(SharedPtrConstnessHandler<DomainType>(domain),f_D0,f_D1));
};

// ----------------------------------------------------------------------------
//  EVALUATORS
// ----------------------------------------------------------------------------
template<int dim, int space_dim> // evaluate values
auto CustomFunction<dim,space_dim>::evaluate_0(const ValueVector<Point> &points, ValueVector<Value> &values) const -> void {
  auto point = points.begin();
  for (auto &val : values ) {
    val = funct_D0(*point);
    ++point;
  }
};

template<int dim, int space_dim> // evaluate first derivatives
auto CustomFunction<dim,space_dim>::evaluate_1(const ValueVector<Point> &points, ValueVector<Derivative<1>> &values) const -> void {
  auto point = points.begin();
  for (auto &val : values ) {
    val = funct_D1(*point);
    ++point;
  }
};

template<int dim, int space_dim> // evaluate second derivatives
auto CustomFunction<dim,space_dim>::evaluate_2(const ValueVector<Point> &points, ValueVector<Derivative<2>> &values) const ->void {
  std::cout << " not implemented yet!" << std::endl;
};

// ----------------------------------------------------------------------------
//   INFO PRINTER
// ----------------------------------------------------------------------------
template<int dim, int space_dim>
void CustomFunction<dim,space_dim>::print_info(LogStream &out) const {
  std::cout << " print_info not implemented yet!" << std::endl;
};

} // of namespace functions.

IGA_NAMESPACE_CLOSE

