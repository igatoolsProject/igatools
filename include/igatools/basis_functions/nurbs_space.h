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

#ifndef __NURBS_SPACE_H_
#define __NURBS_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/basis_functions/bspline_space.h>

IGA_NAMESPACE_OPEN

//Forward declaration to avoid including the header
template < int, int, int > class NURBSElementAccessor;

/**
 * This class represent a function space in which the basis functions are NURBS.
 *
 */
template <int dim_, int range_ = 1, int rank_ = 1>
class NURBSSpace :
    public std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_> >,
    public FunctionSpaceOnGrid<CartesianGrid<dim_> >

{
private:
    using self_t = NURBSSpace<dim_, range_, rank_>;
    using base_t = BSplineSpace<dim_, range_, rank_>;
    using BaseSpace = FunctionSpaceOnGrid<CartesianGrid<dim_> >;

public:
    /** see documentation in \ref FunctionSpaceOnGrid */
    using PushForwardType = PushForward<Transformation::h_grad,dim_,0>;

    using RefSpace = NURBSSpace<dim_, range_, rank_>;

    using GridType = typename PushForwardType::GridType;

    static const int dim = PushForwardType::dim;

    static const int codim = PushForwardType::codim;

    static const int space_dim = PushForwardType::space_dim;

    static const int range = range_;

    static const int rank = rank_;

    static constexpr int n_components = constexpr_pow(range, rank);

    /** Type for the reference face space.*/
    using RefFaceSpace = NURBSSpace<dim-1,range,rank>;

    using DegreeTable = typename base_t::DegreeTable;

    using MultiplicityTable = typename base_t::MultiplicityTable;

    using WeightsTable = typename base_t::template ComponentTable<DynamicMultiArray<Real,dim> >;


public:
    static const bool has_weights = true;
    /**
     * Type for element accessor.
     */
    typedef NURBSElementAccessor<dim, range, rank> ElementAccessor;

    /**
     * Type for iterator over the elements.
     */
    typedef GridForwardIterator<ElementAccessor> ElementIterator;

    /**
     * Type for the face space.
     */
    //TODO rename FaceSpace_t to face_space_t
    typedef NURBSSpace<dim-1, range, rank> FaceSpace_t;


public :
    /** @name Constructor and destructor */
    ///@{
    /**
     * Constructs a maximum regularity NURBSSpace over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous in all components.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    NURBSSpace(std::shared_ptr<GridType> knots, const int &degree);

    /**
     * Returns a shared_ptr wrapping a maximum regularity NURBSSpace over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous in all components.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    static std::shared_ptr< self_t >
    create(std::shared_ptr< GridType > knots, const int &degree);

    /**
     * Constructs a maximum regularity NURBSSpace over CartesianGrid
     * @p knots for the given @p degree for each direction and for each component.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    NURBSSpace(
        std::shared_ptr<GridType> knots,
        const DegreeTable &degree);

    /**
     * Returns a shared_ptr wrapping a maximum regularity NURBSSpace over CartesianGrid
     * @p knots for the given @p degree for each direction and for each component.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    static std::shared_ptr< self_t >
    create(std::shared_ptr<GridType> knots,
           const DegreeTable &degree);

    /**
     * Construct a NURBSSpace over the CartesianGrid @p knots with
     * the given multiplicity vector @p mult_vector for each component
     * and for the given @p degree for each direction and for each component.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    NURBSSpace(
        std::shared_ptr< GridType > knots,
        const MultiplicityTable &mult_vector,
        const DegreeTable &degree);

    /**
     * Returns shared_ptr wrapping a NURBSSpace over the CartesianGrid @p knots with
     * the given multiplicity vector @p mult_vector for each component
     * and for the given @p degree for each direction and for each component.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    static std::shared_ptr< self_t >
    create(std::shared_ptr< GridType > knots,
           const MultiplicityTable &mult_vector,
           const DegreeTable &degree);

    /**
     * Construct a NURBSSpace over the CartesianGrid @p knots with
     * the given multiplicity vector @p mult_vector for each component
     * and for the given @p degree for each direction and for each component.
     */
    NURBSSpace(
        std::shared_ptr< GridType > knots,
        const MultiplicityTable &mult_vector,
        const DegreeTable &degree,
        const WeightsTable &weights);

    /**
     * Returns a shared_ptr wrapping a NURBSSpace over the CartesianGrid @p knots with
     * the given multiplicity vector @p mult_vector for each component
     * and for the given @p degree for each direction and for each component.
     */
    static std::shared_ptr< self_t >
    create(std::shared_ptr< GridType > knots,
           const MultiplicityTable &mult_vector,
           const DegreeTable &degree,
           const WeightsTable &weights);

    /** Destructor */
    ~NURBSSpace() = default;

    ///@}


    /** @name Getting information about the space */
    ///@{
    /**
     * Returns true if all component belong to the same scalar valued
     * space.
     */
    bool is_range_homogeneous() const
    {
        return sp_space_->is_range_homogeneous();
    }
    /**
     * Total number of dofs (i.e number of basis functions aka dimensionality)
     * of the space.
     */
    Size get_num_basis() const
    {
        return sp_space_->get_num_basis();
    }

    /** Return the total number of dofs in each space component. */
    std::array<Size,n_components> get_component_num_basis() const
    {
        return sp_space_->get_component_num_basis();
    }

    /** Return the total number of dofs for the i-th space component. */
    Size get_component_num_basis(int i) const
    {
        return sp_space_->get_component_num_basis(i);
    }
    /**
     *  Return the total number of dofs for the i-th space component
     *  and the j-th direction.
     */
    Size get_component_dir_num_basis(int comp, int dir) const
    {
        return sp_space_->get_component_dir_num_basis(comp, dir);
    }
    /**
     * Returns the number of dofs per element.
     */
    Size get_num_basis_per_element() const
    {
        return sp_space_->get_num_basis_per_element();
    }
    /**
     *  Return the number of dofs per element for the i-th space component.
     */
    Size get_component_num_basis_per_element(int i) const
    {
        return sp_space_->get_component_num_basis_per_element(i);
    }
    /**
     * Returns the degree of the BSpline space for each component and for each coordinate direction.
     * \return The degree of the BSpline space for each component and for each coordinate direction.
     * The first index of the returned object is the component id, the second index is the direction id.
     */
    const DegreeTable &get_degree() const
    {
        return sp_space_->get_degree();
    }
    ///@}

    /** @name Getting the space data */
    ///@{
    /**
     * Return the knots with repetitions, in each direction, for each component of the space.
     */
    const typename base_t::template ComponentTable<CartesianProductArray<Real,dim> > &
    get_knots_with_repetitions() const
    {
        return sp_space_->get_knots_with_repetitions();
    }
    ///@}


    const std::shared_ptr<base_t> get_spline_space() const
    {
        return sp_space_;
    }

    /**
     * Returns a reference to the dense multi array storing the global dofs.
     * Each element has a statically defined zone to read their dofs from,
     * independent of the distribution policy in use.
     */
    const typename base_t::template ComponentTable<DynamicMultiArray<Index,dim>> &get_index_space() const
    {
        return sp_space_->get_index_space();
    }

    /**
     * Transforms basis flat index of the component comp to a basis
     * tensor index.
     */
    TensorIndex<dim>
    flat_to_tensor(const Index index, const Index comp = 0) const
    {
        return sp_space_->flat_to_tensor(index, comp);
    }
    /**
     * Transforms a basis tensor index of the component comp to the
     * corresponding basis flat index.
     */
    Index
    tensor_to_flat(const TensorIndex<dim> &tensor_index,
                   const Index comp = 0) const
    {
        return sp_space_->tensor_to_flat(tensor_index, comp);
    }

    /** Return the push forward (non-const version). */
    std::shared_ptr<PushForwardType> get_push_forward()
    {
        return sp_space_->get_push_forward();
    }

    /** Return the push forward (const version). */
    std::shared_ptr<const PushForwardType> get_push_forward() const
    {
        return sp_space_->get_push_forward();
    }

    std::shared_ptr<const self_t >
    get_reference_space() const
    {
        return this->shared_from_this();
    }

    /**
     * Return the knot multiplicities for each component of the space.
     */
    const typename base_t::template ComponentTable<Multiplicity<dim> > &
    get_multiplicities() const
    {
        return sp_space_->get_multiplicities();
    }


    typename base_t::template ComponentTable<TensorSize<dim>> get_num_dofs() const
	{
    	return sp_space_->get_num_dofs();
	}

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

    /**
     * todo document me
     *
     */
    void print_info(LogStream &out) const;

    /**
     * Get the weights of the NURBSSpace.
     */
    const StaticMultiArray< DynamicMultiArray<Real,dim>,range,rank>
    get_weights() const;

    /**
     * Reset the weights of the NURBSSpace.
     */
    void reset_weights(const StaticMultiArray<DynamicMultiArray<iga::Real,dim>,range,rank> &weights);

private:
    std::shared_ptr<BSplineSpace<dim_, range_, rank_> > sp_space_;
    /**
     * Weights associated to the basis functions.
     */
    StaticMultiArray<DynamicMultiArray<iga::Real,dim>,range,rank> weights_;

    /**
     * Refines the NURBSSpace after the uniform refinement of the BSplineSpace.
     *
     * @param[in] refinement_directions Directions along which the refinement is performed.
     * @param[in] grid_old Grid before the refinement.
     *
     * @pre Before invoking this function, must be invoked the function BSplineSpace::refine().
     * @note This function is connected to the CartesianGrid's signal for the refinement, and
     * it is necessary in order to avoid infinite loops in the CartesianGrid::refine() function calls.
     * @note The implementation of this function is based on "The NURBS Book" Algorithm A5.4.
     *
     * @ingroup h_refinement
     */
    void refine_h_weights(const std::array<bool,dim> &refinement_directions,
                          const GridType &grid_old);

    /**
     * Create a signal and a connection for the refinement.
     */
    void create_refinement_connection();

    friend ElementAccessor;

    /**
     * Performs checks after the construction of the object.
     * In debug mode, if something is going wrong, an assertion will be raised.
     *
     * @warning This function should be used as last line in the implementation of each constructor.
     */
    void perform_post_construction_checks() const;
};



IGA_NAMESPACE_CLOSE


#endif /* __NURBS_SPACE_H_ */


