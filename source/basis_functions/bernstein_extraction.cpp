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

#include <igatools/basis_functions/bernstein_extraction.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using std::endl;


using std::shared_ptr;
using std::make_shared;
using boost::numeric::ublas::matrix_row;

IGA_NAMESPACE_OPEN

auto
BernsteinOperator::
scale_action(const Real scale, const Values &b_values) const -> Values
{
  return scale * prec_prod(*this, b_values);
}

#ifdef SERIALIZATION
template<class Archive>
void
BernsteinOperator::
serialize(Archive &ar)
{
  ar &make_nvp("DenseMatrix",base_class<DenseMatrix>(this));
}
#endif // SERIALIZATION



template<int dim, int range, int rank>
auto
BernsteinExtraction<dim, range, rank>::
get_operator(const int dir, const int inter, const int comp) const -> const Operator &
{
  return ext_operators_[comp][dir][inter];
}


template<int dim, int range, int rank>
auto
BernsteinExtraction<dim, range, rank>::
get_element_operators(TensorIndex<dim> idx) const -> ElemOperTable
{
  ElemOperTable result(ext_operators_.get_comp_map());
  for (auto comp : result.get_active_components_id())
  {
    const auto &ext_operators_comp = ext_operators_[comp];
    auto &result_comp = result[comp];
    for (int dir = 0 ; dir < dim ; ++dir)
      result_comp[dir] = &(ext_operators_comp[dir][idx[dir]]);
  }
  return result;
}


template<int dim, int range, int rank>
auto
BernsteinExtraction<dim, range, rank>::
compute(const Operator &M_j_1,
        typename SafeSTLVector<Real>::const_iterator  y,
        const Real a,
        const Real b) -> Operator
{
  const int j = M_j_1.size1() + 1;
  Operator M_j(j,j);

  SafeSTLVector<Real> alpha(j);
  SafeSTLVector<Real> one_alpha(j,1);
  SafeSTLVector<Real> beta(j, b-a);

  for (int k = 0; k < j; ++k)
  {
    alpha[k] = (y[k+j] - a) / (y[k+j]-y[k]);
    one_alpha[k] -= alpha[k];

    beta[k] /= (y[k+j]-y[k]);
  }

  for (int l = 0; l < j-1; ++l)
  {
    //k = 0
    M_j(0, l) = alpha[0] * M_j_1(0, l);
    //k = 1,...,j-2
    for (int k = 1; k < j-1; ++k)
    {
      M_j(k, l) = alpha[k] * M_j_1(k, l) + one_alpha[k] * M_j_1(k-1, l);
    }
    //k = j-1
    M_j(j-1, l) = one_alpha[j-1] * M_j_1(j-2, l);
  }


  const int l = j-1;
  //k = 0
  M_j(0, l) = M_j(0, l-1) - beta[0] * M_j_1(0, l-1);
  //k = 1,...,j-2
  for (int k = 1; k < j-1; ++k)
  {
    M_j(k, l) = M_j(k, l-1) + beta[k] * (M_j_1(k-1, l-1) - M_j_1(k, l-1));
  }
  //k = j-1
  M_j(j-1, l) = M_j(j-1, l-1) + beta[j-1] * M_j_1(j-2, j-2);

  return M_j;
}



template<int dim, int range, int rank>
void
BernsteinExtraction<dim, range, rank>::
print_info(LogStream &out) const
{
  int c=0;
  for (const auto &comp : ext_operators_)
  {
    out << "Component[" << c++ << "]: " << endl;
    for (int j = 0; j < dim; ++j)
    {
      out << "Direction[" << j << "]:" << endl;
      for (const auto &M : comp[j])
        out << M << endl;
    }
  }
}



template<int dim, int range, int rank>
auto
BernsteinExtraction<dim, range, rank>::
fill_extraction(const int m,
                const SafeSTLVector<Real>  &knots,
                const SafeSTLVector<Real>  &rep_knots,
                const SafeSTLVector<Index> &acum_mult) -> SafeSTLVector<Operator>
{
  const int n_elem = knots.size()-1;

  SafeSTLVector<Operator>  operators(n_elem, Operator(m,m));
  // const auto &x = knots;
  const auto &y = rep_knots;

  auto x = knots;
  x[0] = y[m-1];
  x[n_elem] = *(y.end()-m);
  for (int n=0; n < n_elem; ++n)
  {
    const auto a = x[n];
    const auto b = x[n+1];

    Operator M(1,1);
    M(0,0) = 1/(b-a);
    for (int k = m-2; k>=0; --k)
    {
      const int s = acum_mult[n] + k;

      auto M1 = compute(M, y.begin()+s, a, b);
      M.assign_temporary(M1);
    }

    //Normalized
    auto M2(M);
    const int s = acum_mult[n];
    for (int k = 0; k < m; ++k)
    {
      matrix_row<Operator> mr(M2, k);
      mr *= (y[s+k+m]-y[s+k]);
    }
    operators[n] = M2;
  }
  return operators;
}



template<int dim, int range, int rank>
BernsteinExtraction<dim, range, rank>::
BernsteinExtraction(const Grid<dim> &grid,
                    const KnotsTable &rep_knots,
                    const MultiplicityTable &acum_mult,
                    const DegreeTable &deg)
  :
  ext_operators_(rep_knots.get_comp_map())
{
  for (const int comp : rep_knots.get_active_components_id())
  {
    const auto &deg_comp = deg[comp];
    const auto &rep_knots_comp = rep_knots[comp];
    const auto &acum_mult_comp = acum_mult[comp];
    auto &ext_operators_comp = ext_operators_[comp];
    for (int dir = 0; dir < dim; ++dir)
    {
      const int deg_comp_dir = deg_comp[dir] + 1;
      ext_operators_comp[dir] =
        fill_extraction(deg_comp_dir,
                        grid.get_knot_coordinates(dir),
                        rep_knots_comp[dir],
                        acum_mult_comp[dir]);
    }
  }
}

template<int dim, int range, int rank>
BernsteinExtraction<dim, range, rank>::
BernsteinExtraction(const Space &space_data,
                    const EndBehaviourTable &end_b)
  :
  BernsteinExtraction(*space_data.get_grid(),
                     space_data.compute_knots_with_repetition(end_b),
                     space_data.accumulated_interior_multiplicities(),
                     space_data.get_degree_table())
{}

#ifdef SERIALIZATION

template<int dim, int range, int rank>
template<class Archive>
void
BernsteinExtraction<dim, range, rank>::
serialize(Archive &ar)
{
  ar &make_nvp("ext_operators_",ext_operators_);
}

#endif // SERIALIZATION

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bernstein_extraction.inst>

