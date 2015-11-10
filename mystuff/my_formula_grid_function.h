
#include <igatools/geometry/formula_grid_function.h>

IGA_NAMESPACE_OPEN

namespace grid_functions
{

#define PI 3.14159265358979323846

template<int dim, int space_dim>
class MyFormulaGridFunction : public FormulaGridFunction<dim,space_dim> {
  // some useful aliases
  using base_t   = GridFunction<dim,space_dim>;
  using parent_t = FormulaGridFunction<dim,space_dim>;
  using self_t   = MyFormulaGridFunction<dim,space_dim>;
  using typename base_t::GridType;
public:
  using typename parent_t::Value;
  using typename parent_t::GridPoint;
  template <int order>
  using Derivative = typename parent_t::template Derivative<order>;

public:
  // creators
  static std::shared_ptr<base_t>
  create(const std::shared_ptr<GridType> &domain) {
    return std::shared_ptr<self_t>(new self_t(SharedPtrConstnessHandler<GridType>(domain)));
  }


  static std::shared_ptr<const base_t>
  const_create(const std::shared_ptr<const GridType> &domain) { 
    std::cout << " myformula grid function has been successfully created!" << std::endl;
    return std::shared_ptr<const self_t>(new self_t(SharedPtrConstnessHandler<GridType>(domain)));
  }

  MyFormulaGridFunction(const self_t &) = default;

  virtual ~MyFormulaGridFunction() = default;

  virtual void print_info(LogStream &out) const override final;
  
  void grid_loop();

  //  my experiment
  //const Value funct(Point pt);
  //array<dim,const Value funct_der(Point pt)> funct_grad;

protected:
  MyFormulaGridFunction(
    const SharedPtrConstnessHandler<GridType> &domain)
  :
  parent_t(domain)

{};

private:
  void evaluate_0(const ValueVector<GridPoint> &points,
                  ValueVector<Value> &values) const {
    //std::cout << "evaluation of function" << std::endl;
    auto point = points.begin();
    for (auto &val : values ) {
      val = 1.0;
      for (int idim=0; idim<dim; idim++)
        val *= sin( (*point)[idim] * PI );
      ++point;
    }
  };

  void evaluate_1(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<1>> &values) const {
    //std::cout << "evaluation of function derivatives" << std::endl;
    auto point = points.begin();
    for (auto &val : values ) {
      for (int idim=0; idim<dim; idim++) {
        val[idim] = PI; 
        for (int jdim=0; jdim<idim; jdim++)     val[idim] *= sin((*point)[idim] * PI);
	                                        val[idim] *= cos((*point)[idim] * PI);
	for (int jdim=idim+1; jdim<dim; jdim++) val[idim] *= sin((*point)[idim] * PI);
      }
    }
  };

  void evaluate_2(const ValueVector<GridPoint> &points,
                  ValueVector<Derivative<2>> &values) const {
    std::cout << "evaluation of second derivatives" << std::endl;
  };


#ifdef MESH_REFINEMENT
  void rebuild_after_insert_knots(
    const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
    const Grid<dim> &old_grid) {};
#endif // MESH_REFINEMENT

};

template<int dim, int space_dim>
void MyFormulaGridFunction<dim,space_dim>::print_info(LogStream &out) const {
  std::cout << " print_info method tested!" << std::endl;
};

} // of namespace functions.

IGA_NAMESPACE_CLOSE
