
#include <igatools/geometry/formula_grid_function.h>

IGA_NAMESPACE_OPEN

namespace grid_functions
{

#define PI 3.14159265358979323846

template<int dim, int space_dim>
class CustomGridFunction : public FormulaGridFunction<dim,space_dim> {
  // some useful aliases
  using base_t   = GridFunction<dim,space_dim>;
  using parent_t = FormulaGridFunction<dim,space_dim>;
  using self_t   = CustomGridFunction<dim,space_dim>;
  using typename base_t::GridType;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  // defaulted public constructor and destructor
  CustomGridFunction(const self_t &) = default;
  virtual ~CustomGridFunction() = default;

  // creators
  static std::shared_ptr<self_t>
  create(const std::shared_ptr<GridType> &domain);

  static std::shared_ptr<const self_t>
  const_create(const std::shared_ptr<const GridType> &domain);

  // info printer method
  virtual void print_info(LogStream &out) const override final;

protected:
  CustomGridFunction(const SharedPtrConstnessHandler<GridType> &domain);

public:
  // function's functions (thanks anglosaxons for this this beautifully ambiguous expression)
  std::array<double (*)(const GridPoint),dim> funct;
  //std::array<array<void (*)(double),dim>,dim> grads;

private:
  void evaluate_0(const ValueVector<GridPoint> &points, ValueVector<Value> &values) const;
  void evaluate_1(const ValueVector<GridPoint> &points, ValueVector<Derivative<1>> &values) const;
  void evaluate_2(const ValueVector<GridPoint> &points, ValueVector<Derivative<2>> &values) const;

#ifdef MESH_REFINEMENT
  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) {};
#endif // MESH_REFINEMENT
};

// ----------------------------------------------------------------------------
//  CONSTRUCTOR
// ----------------------------------------------------------------------------
template<int dim, int space_dim>
CustomGridFunction<dim,space_dim>::CustomGridFunction(const SharedPtrConstnessHandler<GridType> &domain) 
  :
  parent_t(domain)
{};

// ----------------------------------------------------------------------------
//  CREATORS
// ----------------------------------------------------------------------------
template<int dim, int space_dim> // non const creator
auto CustomGridFunction<dim,space_dim>::create(const std::shared_ptr<GridType> &domain) -> std::shared_ptr<self_t> {
  auto func = std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<GridType>(domain)));
#ifdef MESH_REFINEMENT
  func->create_connection_for_insert_knots(func);
#endif // MESH_REFINEMENT
  return func;
};

template<int dim, int space_dim> // const creator
auto CustomGridFunction<dim,space_dim>::const_create(const std::shared_ptr<const GridType> &domain) -> std::shared_ptr<const self_t> {
  return std::shared_ptr<const self_t>(new self_t(SharedPtrConstnessHandler<GridType>(domain)));
};

// ----------------------------------------------------------------------------
//  EVALUATORS
// ----------------------------------------------------------------------------
template<int dim, int space_dim> // evaluate values
auto CustomGridFunction<dim,space_dim>::evaluate_0(const ValueVector<GridPoint> &points, ValueVector<Value> &values) const -> void {
  auto point = points.begin();
  for (auto &val : values ) {
    val = funct[0](*point);
    for (int idim=0; idim<dim; idim++)
      val *= sin( (*point)[idim] * PI );
    ++point;
  }
};

template<int dim, int space_dim> // evaluate first derivatives
auto CustomGridFunction<dim,space_dim>::evaluate_1(const ValueVector<GridPoint> &points, ValueVector<Derivative<1>> &values) const -> void {
  auto point = points.begin();
  for (auto &val : values ) {
    for (int idim=0; idim<dim; idim++) {
      val[idim] = PI; 
      for (int jdim=0; jdim<idim; jdim++)
        val[idim] *= sin((*point)[idim] * PI);
      val[idim] *= cos((*point)[idim] * PI);
      for (int jdim=idim+1; jdim<dim; jdim++)
        val[idim] *= sin((*point)[idim] * PI);
    }
  }
};

template<int dim, int space_dim> // evaluate second derivatives
auto CustomGridFunction<dim,space_dim>::evaluate_2(const ValueVector<GridPoint> &points, ValueVector<Derivative<2>> &values) const ->void {
  std::cout << " not implemented yet!" << std::endl;
};

template<int dim, int space_dim>
void CustomGridFunction<dim,space_dim>::print_info(LogStream &out) const {
  std::cout << " print_info method tested!" << std::endl;
};

} // of namespace functions.

IGA_NAMESPACE_CLOSE
