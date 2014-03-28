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

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/geometry/identity_mapping.h>

using std::endl;
using std::array;
using std::vector;
using std::shared_ptr;
using std::make_shared;

IGA_NAMESPACE_OPEN

namespace
{
/*
 * Fills a vector of size corresponding to the
 * number of intervals with the one dimensional Bezier extraction operator
 * given by the knot_values and multiplicities.
 *
 */
void evaluate_extraction_operators(
    const int degree,
    const vector< Real >    &knot_values,
    const vector< Index > &knot_multiplicities,
    vector< DenseMatrix >  &extraction_operators_data
)
{
    Assert(! knot_values.empty(), ExcEmptyObject()) ;
    Assert(! knot_multiplicities.empty(), ExcEmptyObject()) ;
    Assert(knot_values.size() == knot_multiplicities.size(),
           ExcDimensionMismatch(knot_values.size(), knot_multiplicities.size())) ;

    const int num_elements = knot_values.size() - 1 ;
    const int num_funcs = degree + 1 ;


    //----------------------------------------------------------------------------------------------
    // resizing the container for the extraction operators and initialization of the coefficient for
    // the Bernstein's basis
    extraction_operators_data.resize(num_elements) ;

    DenseMatrix identity_matrix = boost::numeric::ublas::identity_matrix<Real>(num_funcs);

    fill(extraction_operators_data.begin(), extraction_operators_data.end(), identity_matrix) ;
    //----------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------
    vector< Real > U ;

    const int num_unique_values = knot_values.size() ;
    for (int iValue = 0 ; iValue < num_unique_values ; iValue++)
    {
        const int multiplicity = knot_multiplicities[ iValue ] ;
        const Real value = knot_values[ iValue ] ;

        for (int iMult = 0 ; iMult < multiplicity ; iMult++)
        {
            U.push_back(value) ;
        }
    }
    //----------------------------------------------------------------------------------------------




    //----------------------------------------------------------------------------------------------
    // compute a vector containing the partial sum of the multiplicities
    vector< int > partial_sum_multiplicities(num_unique_values) ;
    partial_sum(knot_multiplicities.begin(), knot_multiplicities.end(), partial_sum_multiplicities.begin()) ;
    //----------------------------------------------------------------------------------------------




    //----------------------------------------------------------------------------------------------
    const int p = degree ;
    int k, j ;

    vector< Real > alphas(p) ;


    const int m = U.size() - 1 ;


    int iElement = 0 ;
    while (iElement < num_elements)
    {
        DenseMatrix &M = extraction_operators_data[ iElement ] ;


        // Count multiplicity of the knot at location b
        int mul_b = knot_multiplicities[ iElement + 1 ] ;
        int idx_b = partial_sum_multiplicities[ iElement + 1 ] - 1 ;


        if (mul_b < p)
        {
            const Real u_a = knot_values[ iElement ] ;
            const Real u_b = knot_values[ iElement + 1 ] ;

            Real numer = u_b - u_a ;


            for (j = p - mul_b ; j > 0 ; j--)
            {
                Real u_j = 0.0 ;

                for (int idx_j = iElement + 2 ; idx_j <= num_elements ; idx_j++)
                {
                    //TODO: check the correctness of this if condition!
                    if (j < partial_sum_multiplicities[ idx_j ] - idx_b)
                        //  if ( j <= partial_sum_multiplicities[ idx_j ] -  partial_sum_multiplicities[ iElement + 1 ] )
                    {
                        u_j = knot_values[ idx_j ] ;
                        Assert(u_j == U[ idx_b + j ], ExcMessage("Wrong knot value.")) ;
                        break ;
                    }
                }
                //*/

                alphas[ j - 1 ] = numer / (u_j - u_a) ;
            }

            int r = p - mul_b ; // Insert knot r times

            // Update the matrix coefficients for r new knots
            for (j = 1 ; j <= r ; j++)
            {

                int save = r - j ;
                int    s = mul_b + j ; // This many new points
                for (k = p ; k >= s ; k--)
                {
                    Real alpha = alphas[ k - s ] ;

                    for (int row = 0 ; row <= p ; row++)
                    {
                        M(row, k) = alpha * M(row, k) + (1.0 - alpha) * M(row, k - 1) ;
                    }
                }

                if (idx_b < m)
                {
                    DenseMatrix &M_next = extraction_operators_data[ iElement+1 ] ;

                    for (k = 0 ; k <= j ; k++)
                    {
                        M_next(save + k, save) = M(p - j + k, p) ;
                    }
                }
            }
        } // end if ( mul < p )

        iElement++ ; // Finished with the current operator
    }
}



/**
 * Determine the knot span index.
 *
 * @return The knot span index of the value @p u in the knot vector @p U.
 * @param[in] p Degree.
 * @param[in] u Knot values for which the span is requested.
 * @param[in] U Knot vector with repeated values.
 *
 * @note The implementation of this function is based on "The NURBS Book" Algorithm A2.1
 */
int find_span(
    const int p,
    const Real u,
    const vector<Real> &U)
{
    const int m = U.size()-1;
    const int n = m-p;

    // treat special case
    if (u == U[n+1])
        return n;

    // do binary search
    int low = p;
    int high = n+1;
    int mid = (low+high)/2;
    while (u < U[mid] || u >= U[mid+1])
    {
        if (u < U[mid])
            high = mid;
        else
            low = mid;

        mid = (low+high)/2;
    }
    return mid;
}

/**
 * Refine a knot vector
 *
 * @param[in] n Number of basis functions minus one.
 * @param[in] p Degree.
 * @param[in] U Knots vector with repetitions before the refinement.
 * @param[in] X Knots values that defines the refinement.
 * @param[in] Pw Control points before the refinement.
 * @param[out] Ubar Knots vector with repetitions after the refinement.
 * @param[out] alpha
 * @param[in] Qw Control points after the refinement.
 *
 * @note The implementation of this function is based on "The NURBS Book" Algorithm A5.4
 */
template <class T>
void refine_knot_vector(
    const int p,
    const vector<Real> &U,
    const vector<Real> &X,
    const vector<T> &Pw,
    vector<Real> &Ubar,
    vector<T> &Qw)
{
    const int m = U.size()-1;
    const int r = X.size()-1;
    const int a = find_span(p,X[0],U);
    const int b = find_span(p,X[r],U)+1;

//  LogStream out ;
//    out << a << "  " << b << endl;


    const int n = m-p;
    Assert(Pw.size() == n+1, ExcDimensionMismatch(Pw.size(),n+1));
    Assert(Qw.size() == n+r+2, ExcDimensionMismatch(Qw.size(),n+r+2));


    for (int j = 0 ; j <= a-p ; ++j)
        Qw[j] = Pw[j];

    for (int j = b-1 ; j <= n ; ++j)
        Qw[j+r+1] = Pw[j];
    //*/


    Ubar.resize(m+r+2);

    for (int j = 0 ; j <= a ; ++j)
        Ubar[j] = U[j];

    for (int j = b+p; j <= m ; ++j)
        Ubar[j+r+1] = U[j];

    int i = b + p - 1;
    int k = b + p + r;

//  alfa = vector<vector<Real>>(r+1,vector<Real>(p));
//    Assert(int(alfa.size()) == r+1, ExcDimensionMismatch(alfa.size(),r+1));
    for (int j = r ; j >= 0 ; --j)
    {
        while (X[j] <= U[i] && i > a)
        {
            Qw[k-p-1] = Pw[i-p-1];
            Ubar[k] = U[i];
            k = k-1;
            i = i-1;
        }
        Qw[k-p-1] = Qw[k-p];

//        Assert(int(alfa[j].size()) == p, ExcDimensionMismatch(alfa[j].size(),p));
        for (int l = 1 ; l <= p ; ++l)
        {
            int ind = k-p+1;

            Real alfa = Ubar[k+1] - X[j];
            if (abs(alfa) == 0.0)
            {
                Qw[ind-1] = Qw[ind];
            }
            else
            {
                alfa = alfa / (Ubar[k+1] - U[i-p+l]);
                Qw[ind-1] = alfa * Qw[ind-1] + (1.0-alfa) * Qw[ind];
            }
        }
        Ubar[k] = X[j];
        k = k-1;
    }
}

};



template<int dim_, int dim_range_, int rank_>
BSplineSpace<dim_, dim_range_, rank_>::
BSplineSpace(shared_ptr<GridType> cartesian_grid, const int degree)
    :
    BaseSpace(cartesian_grid),
    degree_(TensorIndex<dim>(degree)),
    homogeneous_range_(true)
{
    Assert((rank!=0) || ((rank==0) && (dim_range==1)),
           ExcScalarRange(dim_range));
    Multiplicity<dim> mult(this->get_grid()->get_num_knots_dim());
    mult.fill_max_regularity(degree);
    mult_.fill(mult) ;

    init() ;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
create(shared_ptr< GridType > knots, int degree) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(knots, degree));
}



template<int dim_, int dim_range_, int rank_>
BSplineSpace<dim_, dim_range_, rank_>::
BSplineSpace(shared_ptr<GridType> knots, const TensorIndex<dim> &degree)
    :
    BaseSpace(knots),
    degree_(degree),
    homogeneous_range_(true)
{
    Multiplicity<dim> mult(this->get_grid()->get_num_knots_dim());
    mult.fill_max_regularity(degree);
    for (int i = 0 ; i < n_components ; ++i)
        mult_(i) = mult;

    init() ;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
create(shared_ptr<GridType> knots,
       const TensorIndex<dim> &degree) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(knots, degree));
}



template<int dim_, int dim_range_, int rank_>
BSplineSpace<dim_, dim_range_, rank_>::
BSplineSpace(shared_ptr<GridType> knots,
             const StaticMultiArray<TensorIndex<dim>,dim_range,rank> &degree)
    :
    BaseSpace(knots),
    degree_(degree),
    homogeneous_range_(false)
{
    Multiplicity<dim> mult(this->get_grid()->get_num_knots_dim());
    for (int i = 0 ; i < n_components ; ++i)
    {
        mult.fill_max_regularity(degree(i));
        mult_(i) = mult;
    }

    init() ;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
create(shared_ptr<GridType> knots,
       const StaticMultiArray<TensorIndex<dim>,dim_range,rank> &degree) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(knots, degree));
}



template<int dim_, int dim_range_, int rank_>
BSplineSpace<dim_, dim_range_, rank_>::
BSplineSpace(shared_ptr<GridType> knots,
             const ComponentTable<Multiplicity<dim>> &mult_vectors,
             const StaticMultiArray<TensorIndex<dim>,dim_range,rank> &degree)
    :
    BaseSpace(knots),
    degree_(degree),
    mult_(mult_vectors),
    homogeneous_range_(false)
{
    init() ;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
create(shared_ptr<GridType> knots,
       const ComponentTable<Multiplicity<dim>> &mult_vectors,
       const StaticMultiArray<TensorIndex<dim>,dim_range,rank> &degree) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(knots,mult_vectors,degree));
}



template<int dim_, int dim_range_, int rank_>
void BSplineSpace<dim_, dim_range_, rank_>::
fill_num_dof_per_element()
{
    num_dofs_per_element_ = 0;
    for (int i = 0 ; i < n_components ; ++i)
    {
        int num_dofs_comp = 1 ;
        for (int j = 0 ; j < dim ; ++j)
            num_dofs_comp *= (degree_(i)[j] + 1) ;
        num_dofs_per_element_ += num_dofs_comp ;
    }
}



template<int dim_, int dim_range_, int rank_>
void BSplineSpace<dim_, dim_range_, rank_>::
fill_bezier_extraction_operator()
{

    if (homogeneous_range_)
    {
        map_component_to_active_data_.fill(0) ;
        num_active_components_ = 1;
    }
    else
    {
        for (int i = 0 ; i < n_components ; i++)
            map_component_to_active_data_(i) = i ;

        num_active_components_ = n_components;
    }

    // build the knots with repetitions
    for (int iComp = 0 ; iComp < n_components ; ++iComp)
    {
        for (int jDim = 0 ; jDim < dim ; ++jDim)
        {
            const vector<Real> &knot_values = this->get_grid()->get_knot_coordinates(jDim) ;
            const vector<Index> &knot_multiplicities = mult_(iComp).get_data_direction(jDim);

            Assert(knot_values.size() == knot_multiplicities.size(),
                   ExcDimensionMismatch(knot_values.size(), knot_multiplicities.size()));

            //--------------------------------------------------------------------------------------
            // filling the knots with repetitions

            const int num_unique_values = knot_values.size() ;

            vector<Real> knt_with_reps;

            for (int iValue = 0 ; iValue < num_unique_values ; iValue++)
            {
                const int multiplicity = knot_multiplicities[ iValue ] ;
                const Real value = knot_values[ iValue ] ;

                for (int iMult = 0 ; iMult < multiplicity ; ++iMult)
                    knt_with_reps.push_back(value) ;

            } // end loop iValue
            //--------------------------------------------------------------------------------------
            knots_with_repetitions_(iComp).copy_data_direction(jDim,knt_with_reps);
        } // end loop jDim
    } // end loop iComp
    //--------------------------------------------------------------------------



    //--------------------------------------------------------------------------
    // build the Bezier operator data
    for (int iComp = 0 ; iComp < num_active_components_ ; ++iComp)
    {
        for (int jDim = 0 ; jDim < dim ; ++jDim)
        {

            const int degree = degree_(iComp)[ jDim ] ;
            const vector<Real> &knot_values = this->get_grid()->get_knot_coordinates(jDim);
            const vector<Index> &knot_multiplicities = mult_(iComp).get_data_direction(jDim);

            Assert(knot_values.size() == knot_multiplicities.size(),
                   ExcDimensionMismatch(knot_values.size(), knot_multiplicities.size()));


            //TODO: this vector should be size here and not resize inside
            //evaluate_extraction_operator
            //TODO: we also have to think of a better way not to copy
            vector<DenseMatrix> bezier_op;

            evaluate_extraction_operators(
                degree, knot_values, knot_multiplicities, bezier_op) ;
            bezier_op_data_(iComp).copy_data_direction(jDim,bezier_op);

        }
    }

    //--------------------------------------------------------------------------
    // assign the Bezier operators data to the proper component/interval through
    // their memory address
    for (int iComp = 0 ; iComp < n_components ; iComp++)
    {
        const int component_data_id = map_component_to_active_data_(iComp) ;

        for (int jDim = 0 ; jDim < dim ; jDim++)
        {
            const auto &B_data_vec = bezier_op_data_(component_data_id).get_data_direction(jDim);
            const int num_operators = B_data_vec.size() ;

            vector< const DenseMatrix * > B_vec(num_operators);

            //TODO: avoid this temporary copy
            for (int iOp = 0 ; iOp < num_operators ; iOp++)
                B_vec[iOp] = &B_data_vec[iOp] ;

            bezier_op_(iComp).copy_data_direction(jDim,B_vec);
        }
    }
}



template<int dim_, int dim_range_, int rank_>
void BSplineSpace<dim_, dim_range_, rank_>::
fill_index_space_standard_policy()
{

    for (int iComp = 0 ; iComp < n_components ; iComp++)
    {
        index_space_offset_(iComp) =
            mult_(iComp).compute_index_space_offset(degree_(iComp));
        TensorSize<dim> index_space_size;
        for (int jDim = 0 ; jDim < dim ; ++jDim)
            index_space_size(jDim) = index_space_offset_(iComp).get_data_direction(jDim).back();
        index_space_(iComp).resize(index_space_size);

    }

    ComponentTable<Index> dof_offset_;
    dof_offset_(0) = 0;
    for (int i = 0 ; i < n_components-1 ; ++i)
        dof_offset_(i+1) = dof_offset_(i) + index_space_(i).flat_size();

    // Fill with the standard dof distribution policy
    for (int i = 0 ; i < n_components ; ++i)
        index_space_(i).fill_progression(dof_offset_(i));

}

template<int dim_, int dim_range_, int rank_>
void BSplineSpace<dim_, dim_range_, rank_>::
fill_element_dofs_from_index_space()
{
    const Index n_elements = this->get_grid()->get_num_elements();

    ComponentTable<TensorIndex<dim>> element_n_basis(degree_);
    for (auto &comp : element_n_basis)
        comp += 1;

//    for (int comp = 0; comp<n_components; ++comp)
//        for (int dir = 0; dir <dim; ++dir)
//            element_n_basis(comp)[dir] = degree_(comp)[dir] + 1;


    ComponentTable<int> element_n_basis_comp;
    element_n_basis_comp.fill(1);
    for (int comp = 0; comp<n_components; ++comp)
        for (int dir = 0; dir <dim; ++dir)
            element_n_basis_comp(comp) *= element_n_basis(comp)[dir];


    ComponentTable<Index> element_comp_offset;
    element_comp_offset(0) = 0;
    for (int comp = 1; comp<n_components; ++comp)
        element_comp_offset(comp) = element_n_basis_comp(comp-1);
    for (int i = 1; i < n_components-1; ++i)
        element_comp_offset(i+1) += element_comp_offset(i);

    //Allocate the space

    Index c_n_basis = 0;
    for (int comp = 0; comp<n_components; ++comp)
        c_n_basis += element_n_basis_comp(comp);
    const vector<Index> v0(c_n_basis);
    element_global_dofs_.resize(n_elements, v0);

    for (auto element : *(this->get_grid()))
    {
        const int elem_index  = element.get_flat_index();

        for (int comp = 0; comp < n_components; ++comp)
        {
            const auto comp_offset = element_comp_offset(comp);
            const auto origin      = index_space_offset_(comp).
                                     cartesian_product(element.get_tensor_index());
            const auto increment = element_n_basis(comp);
            const auto end = origin + increment;

            const auto comp_dofs = index_space_(comp).get_sub_array(origin,end).get_data();

            copy(comp_dofs.begin(),comp_dofs.end(),
                 element_global_dofs_[elem_index].begin()+comp_offset);

        }
    }
}



template<int dim_, int dim_range_, int rank_>
void
BSplineSpace<dim_, dim_range_, rank_>::
init_dofs()
{
    /*
     * This function is called from the constructors after the degree,
     * knots and multiplicity have been set.
     */
    fill_num_dof_per_element();
    fill_bezier_extraction_operator();
    fill_index_space_standard_policy();
    fill_element_dofs_from_index_space();

    // Quantities required for flat to tensor
    for (int iComp = 0 ; iComp < n_components ; iComp++)
    {
        for (int jDim = 0 ; jDim < dim ; ++jDim)
            num_dofs_(iComp)[jDim] =
                index_space_offset_(iComp).get_data_direction(jDim).back();
    }
}


template<int dim_, int dim_range_, int rank_>
void
BSplineSpace<dim_, dim_range_, rank_>::
init()
{
    init_dofs();

    //----------------------------------
    // create a signal and a connection for the grid refinement
    this->connect_refinement_h_function(
        std::bind(&self_t::refine_h_after_grid_refinement,
                  this,
                  std::placeholders::_1,std::placeholders::_2));
    //----------------------------------
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
get_reference_space() const -> shared_ptr<const BSplineSpace<dim,dim_range,rank>>
{
    return this->shared_from_this();
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
get_component_num_basis() const -> array<Size,n_components>
{
    array<Size,n_components> n_basis_components;
    for (uint comp_id = 0; comp_id < n_components; ++comp_id)
        n_basis_components[comp_id] = this->get_component_num_basis(comp_id);

    return n_basis_components;
}



template<int dim_, int dim_range_, int rank_>
Index
BSplineSpace<dim_, dim_range_, rank_>::
get_component_num_basis(int iComp) const
{
    //TODO: implement something similar in CartesianProductArray?
    Assert(iComp >= 0 && iComp < n_components, ExcIndexRange(iComp, 0, n_components)) ;

    return index_space_(iComp).flat_size();
}



template<int dim_, int dim_range_, int rank_>
Index
BSplineSpace<dim_, dim_range_, rank_>::
get_component_dir_num_basis(int comp, int dir) const
{
    return num_dofs_(comp)[dir];
}



template<int dim_, int dim_range_, int rank_>
Index
BSplineSpace<dim_, dim_range_, rank_>::
get_component_num_basis_per_element(int iComp) const
{
    //TODO: implement something similar in CartesianProductArray?
    Assert(iComp >= 0 && iComp < n_components, ExcIndexRange(iComp, 0, n_components)) ;


    const auto &degree_component = degree_(iComp) ;

    Index num_dofs_per_element_component = 1 ;
    for (const auto & p : degree_component)
        num_dofs_per_element_component *= (p + 1) ;

    return (num_dofs_per_element_component) ;
}



template<int dim_, int dim_range_, int rank_>
Index
BSplineSpace<dim_, dim_range_, rank_>::
get_num_basis() const
{
    Index result = 0;
    for (int iComp = 0 ; iComp < n_components ; iComp++)
        result += get_component_num_basis(iComp) ;

    return  result;
}



template<int dim_, int dim_range_, int rank_>
Size
BSplineSpace<dim_, dim_range_, rank_>::
get_num_basis_per_element() const
{
    return num_dofs_per_element_;
}



template<int dim_, int dim_range_, int rank_>
bool
BSplineSpace<dim_, dim_range_, rank_>::
is_range_homogeneous() const
{
    return homogeneous_range_;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::begin() const -> ElementIterator
{
    return ElementIterator(
               const_cast< BSplineSpace< dim, dim_range,rank > & >(*this),
               0) ;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::last() const -> ElementIterator
{
    return ElementIterator(
               const_cast< BSplineSpace< dim, dim_range,rank > & >(*this),
               this->get_grid()->get_num_elements() - 1) ;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::end() const -> ElementIterator
{
    return ElementIterator(const_cast<self_t &>(*this),
                           IteratorState::pass_the_end);
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
flat_to_tensor(const Index index, const Index comp) const -> TensorIndex<dim>
{
    return index_space_(comp).flat_to_tensor(index);
}



template<int dim_, int dim_range_, int rank_>
Index
BSplineSpace<dim_, dim_range_, rank_>::
tensor_to_flat(const TensorIndex<dim> &tensor_index,
               const Index comp) const
{
    return index_space_(comp).tensor_to_flat(tensor_index);
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
get_multiplicities() const -> const ComponentTable<Multiplicity<dim>> &
{
    return mult_ ;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
get_push_forward() -> shared_ptr<PushForwardType>
{

    return
    PushForwardType::create(IdentityMapping<dim>::create(this->get_grid()));
}

template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
get_push_forward() const -> shared_ptr<const PushForwardType>
{
    using PushForwardType1 = PushForward<Transformation::h_grad,dim,0>;
    auto grid = this->get_grid() ;
    auto push_fwd =
    PushForwardType1::create(
        IdentityMapping<dim>::create(
            make_shared<GridType>(GridType(*grid))));

    return push_fwd;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
get_degree() const -> const ComponentTable<TensorIndex<dim>> &
{
    return degree_ ;
}



template<int dim_, int dim_range_, int rank_>
void
BSplineSpace<dim_, dim_range_, rank_>::
refine_h_after_grid_refinement(
    const std::array<bool,dim> &refinement_directions,
    const GridType &grid_old)
{
    // keeping the original knots (with repetitions) before the h-refinement
    knots_with_repetitions_pre_refinement_ = knots_with_repetitions_;

    for (int direction_id = 0 ; direction_id < dim ; ++direction_id)
    {
        if (refinement_directions[direction_id])
        {
            // knots in the refined grid along the selected direction
            vector<Real> knots_new = this->get_grid()->get_knot_coordinates(direction_id);

            // knots in the original (unrefined) grid along the selected direction
            vector<Real> knots_old = grid_old.get_knot_coordinates(direction_id);

            vector<Real> knots_added(knots_new.size());

            // find the knots in the refined grid that are not present in the old grid
            auto it = std::set_difference(
                          knots_new.begin(),knots_new.end(),
                          knots_old.begin(),knots_old.end(),
                          knots_added.begin());

            knots_added.resize(it-knots_added.begin());



            for (int comp_id = 0 ; comp_id < self_t::n_components ; ++comp_id)
            {
                //--------------------------------------------------------
                // creating the new multiplicity
                const vector<int> &mult_old = mult_(comp_id).get_data_direction(direction_id);
                const int n_mult_old = mult_old.size();

                const int n_mult_to_add = n_mult_old - 1;
                const int n_mult_new = n_mult_old + n_mult_to_add;

                vector<int> mult_new(n_mult_new);
                for (int i = 0 ; i < n_mult_to_add ; ++i)
                {
                    mult_new[2*i  ] = mult_old[i];
                    mult_new[2*i+1] = 1;
                }
                mult_new[n_mult_new-1] = mult_old[n_mult_old-1];

                mult_(comp_id).copy_data_direction(direction_id,mult_new);
                //--------------------------------------------------------


                //--------------------------------------------------------
                const auto &knots_old =
                    knots_with_repetitions_pre_refinement_(comp_id).get_data_direction(direction_id);

                // knots with repetitions after refinement
                vector<Real> Ubar = knots_old;
                Ubar.insert(Ubar.end(),
                            knots_added.begin(),knots_added.end());

                sort(Ubar.begin(),Ubar.end());

                knots_with_repetitions_(comp_id).copy_data_direction(direction_id,Ubar);
                //--------------------------------------------------------
            } // end loop comp_id

        } // end if(refinement_directions[direction_id])

    } // end loop direction_id

    init_dofs();
}


template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
get_num_dofs() const -> StaticMultiArray<TensorSize<dim>,dim_range,rank>
{
    return num_dofs_;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
get_knots_with_repetitions() const -> const ComponentTable<CartesianProductArray<Real,dim>> &
{
    return knots_with_repetitions_;
}



template<int dim_, int dim_range_, int rank_>
auto
BSplineSpace<dim_, dim_range_, rank_>::
get_index_space() const -> const ComponentTable<DynamicMultiArray<Index,dim>> &
{
    return index_space_;
}



template<int dim_, int dim_range_, int rank_>
void
BSplineSpace<dim_, dim_range_, rank_>::
print_info(LogStream &out) const
{
    out << "BSplineSpace<" << dim_ << "," << dim_range_ << ">" << endl ;

    out.push("\t");

    //----------------------------------------------------------------------------------------------
    out << "Reference Patch: " ;
    this->get_grid()->print_info(out) ;
    out << endl ;
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    for (int iComp = 0 ; iComp < n_components ; iComp++)
    {
        out.push("\t") ;

        out << "Space component[" << iComp << "]: " << endl ;


        //----------------------------------------------------------------------------------------------
        out.push("\t") ;
        for (int jDim = 0 ; jDim < dim ; jDim++)
        {
            out << "Direction[" << jDim << "]:" << endl ;

            out.push("\t") ;


            out << "Degree = " << degree_(iComp)[jDim] << endl ;


            //------------------------------------------------------------------------------------------
            out << "Knot multiplicities: [ " ;
            const auto &mult_iComp_jDim = mult_(iComp).get_data_direction(jDim);
            for (const int & m : mult_iComp_jDim)
            {
                out << m << " " ;
            }
            out << "]" << endl ;
            //------------------------------------------------------------------------------------------

            //------------------------------------------------------------------------------------------
            out << "Knots vectors (with repetitions): [ " ;
            const auto &knots_iComp_jDim = knots_with_repetitions_(iComp).get_data_direction(jDim);
            for (const Real & knt : knots_iComp_jDim)
            {
                out << knt << " ";
            }
            out << "]" << endl ;
            //------------------------------------------------------------------------------------------


            const int num_intervals = bezier_op_data_(iComp).get_data_direction(jDim).size() ;

            const auto &bezier_iComp_jDim = bezier_op_(iComp).get_data_direction(jDim);
            for (int iInterv = 0 ; iInterv < num_intervals ; iInterv++)
            {
                out << "Interval[" << iInterv << "]" << endl ;

                out.push("\t") ;


                out << "Bezier extraction operator:" << endl ;

                out.push("\t") ;

                auto M = *bezier_iComp_jDim[iInterv] ;

                const int num_rows = M.size1() ;
                const int num_cols = M.size2() ;

                for (int row = 0 ; row < num_rows ; row++)
                {
                    for (int col = 0 ; col < num_cols ; col++)
                    {
                        out << M(row, col) << " " ;
                    }
                    out << endl ;
                }
                out.pop() ;

                out.pop() ;
            } // end loop iInterv

            out.pop() ;
        } // end loop iDim
        out.pop() ;
        //----------------------------------------------------------------------------------------------
        out.pop();
    } // end loop iComp


    const int num_dofs = get_num_basis() ;
    out << "Num dofs: " << num_dofs << endl ;

    //TODO: Do we need to call external functions from this output operator?
    out << "Dofs: " << dof_tools::get_dofs(*this)  << endl ;

    const SparsityPattern &sparsity_pattern = dof_tools::get_sparsity_pattern(*this) ;
    out << "Num overlapping funcs: " << sparsity_pattern.get_num_overlapping_funcs() << endl ;

    out.pop();

    out.pop();
}

IGA_NAMESPACE_CLOSE

#include <igatools/basis_functions/bspline_space.inst>

