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
/*
 *  Test for developin new BsplineSpace
 *  author: pauletti
 *  date:
 *
 */
#include "../tests.h"

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/space_spec.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/basis_functions/bernstein_extraction.h>

#include<igatools/geometry/mapping.h>
#include<igatools/geometry/push_forward.h>


#include <igatools/linear_algebra/dof_tools.h>
#include <igatools/geometry/push_forward.h>
#include <igatools/geometry/identity_mapping.h>

using std::endl;
using std::array;
using std::vector;
using std::shared_ptr;
using std::make_shared;

//IGA_NAMESPACE_OPEN

//Forward declaration to avoid including the header
template < int, int, int> class BSplineElementAccessor;

/**
 * Multivariate (tensor product) scalar, vector or k-tensor
 * valued B-spline space.
 * This object can be thought of as providing the
 * B-spline basis functions of a spline space.
 * The space is determined by:
 * - the knot vectors (implemented in the class CartesianGrid)
 * - the multiplicity vectors
 * - and the degree
 *
 * BSplineSpace allows the use of different multiplicity vectors
 * and degrees for each direction and for each component.
 *
 * \section const Constructing a BSplineSpace
 * Similar to the mechanism use in CartesianGrid we use
 * the create technique to create a smartpointer for each constructor
 * the class provides.
 * @todo enter a glossary for create idiom technique and refer from here
 * \code
 * auto space = BSplineSpace<dim>::create();
 * \endcode
 *
 * \section eval Evaluating basis function
 * Basis function are evaluated iterating on the elements.
 * Similarly we can obtain other information such as the local to global.
 *
 * \section dofs Degrees of freedom (dofs)
 * Each basis function is assigned a degree of freedom (an integer),
 * the order of this assignment is done following some prescribed
 * recipe or policy whose knowledge is not Really needed by most users.
 * They are internally stored in a grid-like multiarray container,
 * called the index space.
 * It works together with the index_space_offset which for each element
 * provides a view of the index
 * space to know which
 * are the non-zero basis function on each element, what we generally
 * refer to as the local to global mapping.
 *
 * \section bezier Storage of the basis functions (Bezier Extraction)
 * The basis functions on each element are stored in the Bspline space
 * as the 1D Bezier extraction operator.
 * When they need to be evaluated the operator applied to the Berenstein polynomials,
 * allows to compute the values at the quadrature points.
 *
 * @todo write a module about cache optimization and handling.
 *
 * \section hom_range Optimizing homogeneous range type vector spaces
 *
 * \author martinelli, 2012, 2013, 2014
 * \author pauletti, 2012, 2013
 *
 * @tparam dim Dimensionality of the parametric space (must be equal to the dimensionality
 * of the grid used top build the space
 *
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class BSplineSpace :
    public std::enable_shared_from_this<BSplineSpace<dim_,range_,rank_> >,
    public SplineSpace<dim_, range_, rank_>
{
private:
    using BaseSpace = SplineSpace<dim_, range_, rank_>;

    /** Type for current class. */
    using self_t = BSplineSpace<dim_,range_,rank_>;

public:
    /** see documentation in \ref FunctionSpaceOnGrid */
    using PushForwardType = PushForward<Transformation::h_grad,dim_,0>;

    using RefSpace = self_t;

    using GridType = typename PushForwardType::GridType;

    static const int dim = PushForwardType::dim;

    static const int codim = PushForwardType::codim;

    static const int space_dim = PushForwardType::space_dim;

    static const int range = range_;

    static const int rank = rank_;

    using BaseSpace::n_components;

    static const bool has_weights = false;

public:
    /** Type for the reference face space.*/
    using RefFaceSpace = BSplineSpace<dim-1,range,rank>;

    /** Type for the element accessor. */
    using ElementAccessor = BSplineElementAccessor<dim,range,rank>;

    /** Type for iterator over the elements.  */
    using ElementIterator = GridForwardIterator<ElementAccessor>;

    using typename BaseSpace::DegreeTable;
    using typename BaseSpace::MultiplicityTable;
    using typename BaseSpace::KnotsTable;

    enum class EndBehaviour
    {
        interpolatory, periodic, end_knots
    };
    using EndBehaviourTable = typename BaseSpace::template
            ComponentContainer<std::array<EndBehaviour, dim> >;

public:
#if 0
    /** @name Constructor and destructor */
    ///@{
    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous
     * in all components.
     */
    explicit BSplineSpace(std::shared_ptr<GridType> knots, const int degree);

    /**
     * Smart pointer create construction technique, see more detail
     * in the corresponding wrapped constructor before.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<GridType> knots, const int degree);

    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    explicit BSplineSpace(std::shared_ptr<GridType> knots,
                          const TensorIndex<dim> &degree);

    /**
     * Smart pointer create construction technique, see more detail
     * in the corresponding wrapped constructor before.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<GridType> knots, const TensorIndex<dim> &degree);

    /**
     * Constructs a maximum regularity BSpline space over CartesianGrid
     * @p knots for the given @p degree for each direction and for each
     * component.
     */
    explicit BSplineSpace(std::shared_ptr<GridType> knots,
                          const DegreeTable &degree,
                          const bool homogeneous_range = false);

    /**
     * Smart pointer create construction technique, see more detail
     * in the corresponding wrapped constructor before.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<GridType> knots,
           const DegreeTable &degree);
#endif
    /**
     * Constructs a BSpline space over the CartesianGrid
     * @p knots with the given multiplicity vector @p mult_vectors
     * for each component
     * and the given @p degree for each direction and for each
     * component.
     */
    explicit BSplineSpace(const DegreeTable &deg,
                          std::shared_ptr<GridType> knots,
                          std::shared_ptr<const MultiplicityTable> interior_mult,
                          const EndBehaviourTable & ends
                          );


    /**
     * Smart pointer create construction technique, see more detail
     * in the corresponding wrapped constructor before.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<GridType> knots,
           std::shared_ptr<const MultiplicityTable> interior_mult,
           const DegreeTable &deg);

    /** Destructor */
    ~BSplineSpace() = default;
    ///@}


    /** @name Assignment operators */
    ///@{

    /** Copy assignment. Not allowed to be used. */
    self_t &
    operator=(const self_t &space) = delete;
    ///@}

    /** @name Getting information about the space */
    ///@{
    /**
     * Returns true if all component belong to the same scalar valued
     * space.
     */
 //   bool is_range_homogeneous() const;


    std::shared_ptr<const self_t >
    get_reference_space() const;

    /**
     * Returns a reference to the dense multi array storing the global dofs.
     * Each element has a statically defined zone to read their dofs from,
     * independent of the distribution policy in use.
     */
  //  const IndexSpaceTable &get_index_space() const;
    ///@}

    /** @name Functions involving the element iterator */
    ///@{
    /**
     * Returns a element iterator to the first element of the patch
     */
    ElementIterator begin() const;

    /**
     * Returns a element iterator to the last element of the patch
     */
    ElementIterator last() const;


    /**
     * Returns a element iterator to one-pass the end of patch.
     */
    ElementIterator end() const;


    ///@}

    /**
     * Prints internal information about the space.
     * @note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;









    /** Return the push forward (non-const version). */
    std::shared_ptr<PushForwardType> get_push_forward();


    /** Return the push forward (const version). */
    std::shared_ptr<const PushForwardType> get_push_forward() const;



private:

    /** Container with the local to global basis indices */
    DofDistribution<dim, range, rank> basis_indices_;

    /** @name Bezier extraction operator. */
    BernsteinExtraction<dim, range, rank> operators_;

    ///@{
protected:
    /**
     * True if each component of the vector valued space belongs
     * to the same scalar valued space.
     */
  //  const bool homogeneous_range_;

    //TODO(pauletti, Apr 27, 2014): make this private w/getter


private:
    ///@}


    friend class BSplineElementAccessor<dim, range, rank>;


//    /**
//     * Refines the function space after a grid uniform refinement.
//     *
//     * @param[in] refinement_directions Directions along which the refinement is performed.
//     * @param[in] grid_old Grid before the refinement.
//     *
//     * @pre Before invoking this function, must be invoked the function grid_->refine().
//     * @note This function is connected to the CartesianGrid's signal for the refinement, and
//     * it is necessary in order to avoid infinite loops in the refine() function calls.
//     *
//     * @ingroup h_refinement
//     */
//    void refine_h_after_grid_refinement(
//        const std::array<bool,dim> &refinement_directions,
//        const GridType &grid_old) ;


    //TODO(pauletti, Apr 27, 2014): make this private and use a getter
//public:
//    /** Knots with repetitions before refinement */
//    KnotsTable knots_with_repetitions_pre_refinement_;

public:
    DeclException1(ExcScalarRange, int,
                   << "Range " << arg1 << "should be 0 for a scalar valued"
                   << " space.");
};


#if 0
template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(shared_ptr<GridType> grid, const int degree)
    :
    BSplineSpace(grid, TensorIndex<dim>(degree))
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(shared_ptr< GridType > knots, const int degree) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(knots, degree));
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(shared_ptr<GridType> knots, const TensorIndex<dim> &degree)
    :
    BSplineSpace(knots, DegreeTable(degree), true)
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(shared_ptr<GridType> knots,
       const TensorIndex<dim> &degree) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(knots, degree));
}



template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(shared_ptr<GridType> knots,
             const DegreeTable &degree,
             const bool homogeneous_range)
    :
    BSplineSpace(knots,
                 MultiplicityTable(knots, degree, true),
                 homogeneous_range)
{}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(shared_ptr<GridType> knots,
       const DegreeTable &degree) -> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(knots, degree));
}

#endif

template<int dim_, int range_, int rank_>
BSplineSpace<dim_, range_, rank_>::
BSplineSpace(const DegreeTable &deg,
             std::shared_ptr<GridType> knots,
             std::shared_ptr<const MultiplicityTable> interior_mult,
             const EndBehaviourTable & ends = EndBehaviourTable())
    :
    BaseSpace(deg, knots, interior_mult),
    basis_indices_(knots,BaseSpace::accumulated_interior_multiplicities(),
                   BaseSpace::get_num_basis_table(),BaseSpace::get_num_basis_per_element_table()),
    operators_(knots, BaseSpace::compute_knots_with_repetition(BaseSpace::EndBehaviour::interpolatory),
               BaseSpace::accumulated_interior_multiplicities(), deg)
{
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
create(std::shared_ptr<GridType> knots,
       std::shared_ptr<const MultiplicityTable> interior_mult,
       const DegreeTable &deg)
-> shared_ptr<self_t>
{
    return shared_ptr<self_t>(new self_t(deg, knots, interior_mult));
}








//template<int dim_, int range_, int rank_>
//void
//BSplineSpace<dim_, range_, rank_>::
//init()
//{
//    init_dofs();
//
//    // create a signal and a connection for the grid refinement
//    this->connect_refinement_h_function(
//        std::bind(&self_t::refine_h_after_grid_refinement, this,
//                  std::placeholders::_1,std::placeholders::_2));
//}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_reference_space() const -> shared_ptr<const self_t>
{
    return this->shared_from_this();
}






//template<int dim_, int range_, int rank_>
//bool
//BSplineSpace<dim_, range_, rank_>::
//is_range_homogeneous() const
//{
//    return homogeneous_range_;
//}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::begin() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(), 0);
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::last() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           this->get_grid()->get_num_elements() - 1);
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::end() const -> ElementIterator
{
    return ElementIterator(this->shared_from_this(),
                           IteratorState::pass_the_end);
}


#if 0

template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
flat_to_tensor(const Index index, const Index comp) const -> TensorIndex<dim>
{
    return index_space_(comp).flat_to_tensor(index);
}



template<int dim_, int range_, int rank_>
Index
BSplineSpace<dim_, range_, rank_>::
tensor_to_flat(const TensorIndex<dim> &tensor_index,
               const Index comp) const
{
    return index_space_(comp).tensor_to_flat(tensor_index);
}

#endif



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_push_forward() -> shared_ptr<PushForwardType>
{

    return
    PushForwardType::create(IdentityMapping<dim>::create(this->get_grid()));
}



template<int dim_, int range_, int rank_>
auto
BSplineSpace<dim_, range_, rank_>::
get_push_forward() const -> shared_ptr<const PushForwardType>
{
    using PushForwardType1 = PushForward<Transformation::h_grad,dim,0>;
    auto grid = this->get_grid();
    auto push_fwd =
    PushForwardType1::create(
        IdentityMapping<dim>::create(
            make_shared<GridType>(GridType(*grid))));

    return push_fwd;
}


#if 0
template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
refine_h_after_grid_refinement(
    const std::array<bool,dim> &refinement_directions,
    const GridType &grid_old)
{
    // keeping the original knots (with repetitions) before the h-refinement
    knots_with_repetitions_pre_refinement_ = knots_with_repetitions_;

    for (int direction_id = 0; direction_id < dim; ++direction_id)
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



            for (int comp_id = 0; comp_id < self_t::n_components; ++comp_id)
            {
                //--------------------------------------------------------
                // creating the new multiplicity
                const vector<int> &mult_old = mult_(comp_id).get_data_direction(direction_id);
                const int n_mult_old = mult_old.size();

                const int n_mult_to_add = n_mult_old - 1;
                const int n_mult_new = n_mult_old + n_mult_to_add;

                vector<int> mult_new(n_mult_new);
                for (int i = 0; i < n_mult_to_add; ++i)
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
#endif




template<int dim_, int range_, int rank_>
void
BSplineSpace<dim_, range_, rank_>::
print_info(LogStream &out) const
{
    BaseSpace::print_info(out);

    basis_indices_.print_info(out);

    operators_.print_info(out);

#if 0

    //TODO: Do we need to call external functions from this output operator?
    out << "Dofs: " << dof_tools::get_dofs(this->shared_from_this())  << endl;

    const SparsityPattern &sparsity_pattern =
        dof_tools::get_sparsity_pattern(this->shared_from_this());
    out << "Num overlapping funcs: ";
    out << sparsity_pattern.get_num_overlapping_funcs() << endl;
#endif
}
//IGA_NAMESPACE_CLOSE



int main()
{
    out.depth_console(10);

    const int dim=1;
    using BSplineSpace = BSplineSpace<dim>;
    using MultiplicityTable = typename BSplineSpace::MultiplicityTable;

    typename BSplineSpace::DegreeTable deg{{2}};
    auto grid = CartesianGrid<dim>::create(4);
    auto int_mult = shared_ptr<MultiplicityTable>(new MultiplicityTable ({ {{1,3}} }));
    auto space = BSplineSpace::create(grid, int_mult, deg);

    space->print_info(out);

    return 0;
}
