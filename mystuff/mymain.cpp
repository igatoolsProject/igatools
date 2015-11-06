#include <igatools/basis_functions/bspline.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/basis_functions/bspline_element_handler.h>
#include <igatools/base/quadrature_lib.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/base/logstream.h>
#include <igatools/functions/ig_grid_function.h>
#include <igatools/io/writer.h>

using namespace iga;
using namespace std;
LogStream out;

// ----------------------------------------------------------------------------
//    MY SUPERCOOL CLASS TURBO SPECIAL CONTAINER
// ----------------------------------------------------------------------------
template<int dim_>
class MySpace {

  private:
    using self_t = MySpace<dim_>; // creating the alias
    static int count;

    // constructors
  public:
    MySpace() = delete;
    MySpace(const Size nel, const Size deg);
    MySpace(const TensorSize<dim_> nel, const TensorIndex<dim_> deg);
    static const int dim = dim_;

  public:
      // data
    shared_ptr<const Grid<dim_>>        grid;
    shared_ptr<const SplineSpace<dim_>> space;
    shared_ptr<const BSpline<dim_>>     basis;

      // methods: creators
    static std::shared_ptr<self_t>
    create(const Size nel,const Size deg);

    static std::shared_ptr<const self_t> 
    const_create(const Size nel, const Size deg);

    /*static std::shared_ptr<const self_t> MySpace<dim_>::const_create(
                     const TensorSize<dim_> nel;
		     const TensorIndex<dim_> deg);*/
      // methods:
    void how_are_you_doin() const;
    void stiffness_assemble() const;
};
template<int dim_>
int MySpace<dim_>::count=0;

// ----------------------------------------------------------------------------
// constructors: the simple one
template<int dim_>
MySpace<dim_>::MySpace(const Size nel, const Size deg)
  : grid  {Grid<dim_>::const_create(nel+1)}
  , space {SplineSpace<dim_>::const_create(deg,grid)}
  , basis {BSpline<dim_>::const_create(space)}
{count++;}
// constructors: the cool one
template<int dim_>
MySpace<dim_>::MySpace(const TensorSize<dim_> nel, const TensorIndex<dim_> deg) {
  TensorSize<dim_> nknt;
  for (int idim=0; idim<dim_; idim++) nknt[idim]=nel[idim]+1;
  grid  = Grid<dim_>::const_create(nknt);
  space = SplineSpace<dim_>::const_create(deg,grid);
  basis = BSpline<dim_>::const_create(space);
  count++;
}

// ----------------------------------------------------------------------------
// creators: the simple ones in const/non const versions
template<int dim_>
auto MySpace<dim_>::create(const Size nel, const Size deg) -> shared_ptr<self_t> {
  return shared_ptr<self_t>(new self_t(nel,deg));
} 

template<int dim_>
auto MySpace<dim_>::const_create(const Size nel, const Size deg) -> shared_ptr<const self_t> {
  return create(nel,deg);
}

// ----------------------------------------------------------------------------
// methods
template<int dim_>
void MySpace<dim_>::how_are_you_doin() const {
  out << "regardless copyright infringements, this is the MySpace class instance number " << count << ":" << endl;
  out << "        elements: " << grid->get_num_elements() << " = " << grid->get_num_intervals() << endl;
  out << "         degrees: " << space->get_degree_table()[0] << endl;
  TensorSize<dim_> nbf;
  for (int idim=0; idim<dim_; idim++) nbf[idim]=space->get_num_basis(0,idim);
  out << " basis functions: " << space->get_num_basis() << " = " << nbf << endl;
}

template<int dim_>
void MySpace<dim_>::stiffness_assemble() const {

  // starting up the machinery for the cache
  auto       el  = basis->begin();
  const auto eld = basis->end();
  auto cache_handler = BSpline<dim_>::ElementHandler::create(basis);

/*  flag = ValueFlag::Flags::value |
              ValueFlag::Flags::gradient |
	      ValueFlag::Flags::w_measure;

  TensorSize<dim> nqn;
  typename SplineSpace<dim_>::DegreeTable deg = space->get_degree_table();
  for (int idim=0; idim<dim_; idim++) nqn[idim] = deg[0][idim]+1;
  QGauss<dim_> quad(nqn);

  cache_handler->set_element_flags(flag);
  cache_handler->init_element_cache(el,quad);

  // retrieving the last datum and then starting the loop
  auto Nbf = el->get_num_basis(DofProperties::active);
  for (; el!=eld; ++el) {
    cache_handler->fill_element_cache(el);
    // preparing some stuff
    DenseMatrix lA(Nbf,Nbf);
    lA=0.0;
    DenseVector lb(Nbf);
    lb=0.0;
    auto values = el->get_element_values(0,DofProperties::active);
    auto grad   = el->get_element_gradients(0,DofProperties::active);
    auto w_meas = el->get_w_measures(0);
  }*/

}

int main() {

  const int dim  = 2;
  const int nel  = 4;
  const int deg  = 1;

  // cool constructor for everything: grid, space, basis
  TensorSize<dim>  mnel; for (int idim=0; idim<dim; idim++) mnel[idim]=nel+idim;
  TensorIndex<dim> mdeg; for (int idim=0; idim<dim; idim++) mdeg[idim]=deg+idim;
  // testing the simple constructor
  MySpace<dim> space1(nel,deg);
  space1.how_are_you_doin();
  // testing the cool constructor
  MySpace<dim> space2(mnel,mdeg);
  space2.how_are_you_doin();
  // testing the simple creator
  auto space3 = MySpace<dim>::create(nel,deg+2);
  space3->how_are_you_doin();
  // testing the simple constant creator
  auto space4 = MySpace<dim>::const_create(nel,deg+3);
  space4->how_are_you_doin();
  space4->stiffness_assemble();

  return 0;
}
