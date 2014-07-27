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
 * Multivariate (tensor product) scalar, vector or k-tensor
 * valued NURBS space.
 *
 * @ingroup containers
 */
template <int dim_, int range_ = 1, int rank_ = 1>
class NURBSSpace :
    public std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_> >,
    public FunctionSpaceOnGrid<CartesianGrid<dim_> >

{
private:
    using BaseSpace = FunctionSpaceOnGrid<CartesianGrid<dim_>>;
    using self_t = NURBSSpace<dim_, range_, rank_>;
    using spline_space_t = BSplineSpace<dim_, range_, rank_>;

public:
    /** see documentation in \ref FunctionSpaceOnGrid */
    using PushForwardType = PushForward<Transformation::h_grad,dim_,0>;

    using RefSpace = self_t;

    using GridType = typename PushForwardType::GridType;

    static const int dim = PushForwardType::dim;

    static const int codim = PushForwardType::codim;

    static const int space_dim = PushForwardType::space_dim;

    static const int range = spline_space_t::range;

    static const int rank = spline_space_t::rank;

    static constexpr int n_components = spline_space_t::n_components;

    static const std::array<int, n_components> components;

    static const bool has_weights = true;

    static const std::array<int, dim> dims;



public:
    using Func = typename spline_space_t::Func;
    template <int order>
    using Derivative = typename spline_space_t::template Derivative<order>;
    using Point = typename spline_space_t::Point;
    using Value = typename spline_space_t::Value;
    using Div   = typename spline_space_t::Div;

public:
    /** Type for the reference face space.*/
    using RefFaceSpace = Conditional<(dim>0),
          NURBSSpace<dim-1, range, rank>,
          NURBSSpace<0, range, rank> >;

    using FaceSpace = PhysicalSpace<RefFaceSpace, typename PushForwardType::FacePushForward>;

    /** Type for the element accessor. */
    using ElementAccessor = NURBSElementAccessor<dim, range, rank> ;

    /** Type for iterator over the elements.  */
    using ElementIterator = GridForwardIterator<ElementAccessor>;

public:
//    /** Container indexed by the components of the space */
    template< class T>
    using ComponentContainer = typename spline_space_t::template ComponentContainer<T>;

    using DegreeTable = typename spline_space_t::DegreeTable;
    using MultiplicityTable = typename spline_space_t::MultiplicityTable;

    using EndBehaviour = typename spline_space_t::EndBehaviour;
    using EndBehaviourTable = typename spline_space_t::EndBehaviourTable;
    using InteriorReg= typename spline_space_t::InteriorReg;
    using SpaceDimensionTable = typename spline_space_t::SpaceDimensionTable;

    using Weights = DynamicMultiArray<Real,dim>;
    using WeightsTable = ComponentContainer<Weights>;

public:
    /** @name Constructor and destructor */
    ///@{

    /**
     * Constructs a maximum regularity NURBSSpace over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous in all components.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    explicit NURBSSpace(const int degree,
                        std::shared_ptr< GridType > knots,
                        const WeightsTable &weights);

    /**
     * Returns a shared_ptr wrapping a maximum regularity NURBSSpace over CartesianGrid
     * @p knots for the given @p degree in all directions and homogeneous in all components.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    static std::shared_ptr< self_t >
    create(const int degree,
           std::shared_ptr< GridType > knots,
           const WeightsTable &weights);

    /**
     * Constructs a maximum regularity NURBSSpace over CartesianGrid
     * @p knots for the given @p degree for each direction and for each component.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    explicit NURBSSpace(const DegreeTable &degree,
                        std::shared_ptr<GridType> knots,
                        const WeightsTable &weights);

    /**
     * Returns a shared_ptr wrapping a maximum regularity NURBSSpace over CartesianGrid
     * @p knots for the given @p degree for each direction and for each component.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    static std::shared_ptr< self_t >
    create(const DegreeTable &degree,
           std::shared_ptr<GridType> knots,
           const WeightsTable &weights);

    /**
     * Construct a NURBSSpace over the CartesianGrid @p knots with
     * the given multiplicity vector @p mult_vector for each component
     * and for the given @p degree for each direction and for each component.
     */
    explicit  NURBSSpace(const DegreeTable &deg,
                         std::shared_ptr<GridType> knots,
                         std::shared_ptr<const MultiplicityTable> interior_mult,
                         const EndBehaviourTable &ends,
                         const WeightsTable &weights);


    /**
     * Smart pointer create construction technique, see more detail
     * in the corresponding wrapped constructor before.
     */
    static std::shared_ptr<self_t>
    create(const DegreeTable &deg,
           std::shared_ptr<GridType> knots,
           std::shared_ptr<const MultiplicityTable> interior_mult,
           const EndBehaviourTable &ends = EndBehaviourTable(),
           const WeightsTable &weights = WeightsTable());

    explicit  NURBSSpace(std::shared_ptr<spline_space_t> bs_space,
                         const WeightsTable &weights);

    static std::shared_ptr<self_t>
    create(std::shared_ptr<spline_space_t> bs_space,
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
#if 0
    bool is_range_homogeneous() const
    {
        return sp_space_->is_range_homogeneous();
    }
#endif
    /**
     * Total number of dofs (i.e number of basis functions aka dimensionality)
     * of the space.
     */
    Size get_num_basis() const
    {
        return sp_space_->get_num_basis();
    }

    /** Return the total number of dofs for the i-th space component. */
    Size get_num_basis(const int i) const
    {
        return sp_space_->get_num_basis(i);
    }
    /**
     *  Return the total number of dofs for the i-th space component
     *  and the j-th direction.
     */
    Size get_num_basis(const int comp, const int dir) const
    {
        return sp_space_->get_num_basis(comp, dir);
    }
    /**
     * Returns the number of dofs per element.
     */
    Size get_num_basis_per_element() const
    {
        return sp_space_->get_num_basis_per_element();
    }

    const SpaceDimensionTable get_num_basis_per_element_table() const
    {
        return sp_space_->get_num_basis_per_element_table();
    }
    /**
     *  Return the number of dofs per element for the i-th space component.
     */
    Size get_num_basis_per_element(int i) const
    {
        return sp_space_->get_num_basis_per_element(i);
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

    std::vector<Index> get_loc_to_global(const TensorIndex<dim> &j) const
    {
        return sp_space_->get_loc_to_global(j);
    }
    ///@}

    /** @name Getting the space data */
    ///@{
    /**
     * Return the knots with repetitions, in each direction, for each component of the space.
     */
#if 0
    const typename spline_space_t::template ComponentTable<CartesianProductArray<Real,dim> > &
    get_knots_with_repetitions() const
    {
        return sp_space_->get_knots_with_repetitions();
    }
#endif
    ///@}


    const std::shared_ptr<spline_space_t> get_spline_space() const
    {
        return sp_space_;
    }

    const DofDistribution<dim, range, rank> &get_basis_indices() const
    {
        return sp_space_->get_basis_indices();
    }


    DofDistribution<dim, range, rank> &get_basis_indices()
    {
        return sp_space_->get_basis_indices();
    }

    std::shared_ptr<DofsManager> get_dofs_manager() const;


#if 0
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
#endif
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

#if 0
    /**
     * Return the knot multiplicities for each component of the space.
     */
    const MultiplicityTable &
    get_multiplicities() const
    {
        return sp_space_->get_multiplicities();
    }
#endif
#if 0
    const SpaceDimensionTable &get_num_basis_table() const
    {
        return sp_space_->get_num_basis_table();
    }



    /**
     * @todo Missing documentation
     */
    const std::vector<std::vector<Index>> &get_element_global_dofs() const
    {
        return sp_space_->get_element_global_dofs();
    }
#endif

    std::shared_ptr<RefFaceSpace>
    get_ref_face_space(const Index face_id,
                       std::vector<Index> &face_to_element_dofs,
                       std::map<int, int> &elem_map) const;

    std::shared_ptr<FaceSpace>
    get_face_space(const Index face_id,
                   std::vector<Index> &face_to_element_dofs) const;



    /**
     * Adds an @p offset to the values of the dof ids.
     */
    void add_dofs_offset(const Index offset)
    {
        sp_space_->add_dofs_offset(offset);
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
     * Get the weights of the NURBSSpace.
     */
    const WeightsTable &get_weights() const;

    /**
     * Reset the weights of the NURBSSpace.
     */
    void reset_weights(const WeightsTable &weights);

    /**
     * Prints internal information about the space.
     * @note Mostly used for debugging and testing.
     */
    void print_info(LogStream &out) const;


private:
    /**
     * B-spline space
     */
    std::shared_ptr<spline_space_t> sp_space_;

    /**
     * Weights associated to the basis functions.
     */
    WeightsTable weights_;

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


