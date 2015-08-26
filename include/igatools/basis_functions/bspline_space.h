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

#ifndef __BSPLINE_SPACE_H_
#define __BSPLINE_SPACE_H_

#include <igatools/base/config.h>
#include <igatools/base/logstream.h>

#include <igatools/basis_functions/reference_space.h>
#include <igatools/basis_functions/dof_distribution.h>
#include <igatools/basis_functions/bernstein_extraction.h>
#include <igatools/geometry/mapping.h>
#include <igatools/basis_functions/physical_space.h>

IGA_NAMESPACE_OPEN


template <int, int, int> class BSplineElement;
template <int, int, int> class BSplineElementHandler;
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
 * When they need to be evaluated the operator applied to the
 * Berenstein polynomials,
 * allows to compute the values at the quadrature points.
 *
 * @todo write a module about cache optimization and handling.
 *
 * \section hom_range Optimizing homogeneous range type vector spaces
 *
 * \author martinelli, 2012, 2013, 2014
 * \author pauletti, 2012, 2013, 2014
 *
 *
 * @ingroup containers
 * @ingroup serializable
 */
template<int dim_, int range_ = 1, int rank_ = 1>
class BSplineSpace :
    public ReferenceSpace<dim_, range_, rank_>
{
private:
    using BaseSpace = ReferenceSpace<dim_, range_, rank_>;

    /** Type for current class. */
    using self_t = BSplineSpace<dim_,range_,rank_>;

public:
    /** see documentation in \ref Space */

    using GridType = CartesianGrid<dim_>;
    using ElementHandler = BSplineElementHandler<dim_, range_, rank_>;


    static const int dim       = dim_;
    static const int codim     = 0;
    static const int space_dim = dim_;
    static const int range     = range_;
    static const int rank      = rank_;
    static const bool is_physical_space = false;

public:
    using typename BaseSpace::Func;
    using typename BaseSpace::Point;
    using typename BaseSpace::Value;

    template <int order>
    using Derivative = typename BaseSpace::template Derivative<order>;

    using typename BaseSpace::Div;

    /**
     * See documentation in \ref Space
     *
     * @see Space
     */
    using PushForwardType = typename BaseSpace::PushForwardType;

    using RefSpace = typename BaseSpace::RefSpace;

    using RefPoint = Point;

public:
    /** Type for the element accessor. */
    using ElementAccessor = BSplineElement<dim,range,rank>;

    /** Type for iterator over the elements.  */
    using ElementIterator = GridIterator<ReferenceElement<dim,range,rank>>;


    using SpaceData = SplineSpace<dim_,range_,rank_>;

    using Degrees = typename SpaceData::Degrees;
    using Multiplicity = typename SpaceData::Multiplicity;
    using EndBehaviour = typename SpaceData::EndBehaviour;
    using Periodicity = typename SpaceData::Periodicity;

    using KnotsTable = typename SpaceData::KnotsTable;
    using DegreeTable = typename SpaceData::DegreeTable;
    using MultiplicityTable = typename SpaceData::MultiplicityTable;
    using TensorSizeTable = typename SpaceData::TensorSizeTable;
    using PeriodicityTable = typename SpaceData::PeriodicityTable;
    using EndBehaviourTable = typename SpaceData::EndBehaviourTable;

    using BaseSpace::ComponentContainer;

    using IndexType = TensorIndex<dim_>;
    using PropertyList = PropertiesIdContainer<IndexType>;
    using List = typename PropertyList::List;
    using ListIt = typename PropertyList::List::iterator;

public:
    /**
     * @name Creators.
     */
    ///@{
    /**
     * Builds and returns a maximum regularity (non-const) BSpline space
     * over a (non-const) CartesianGrid
     * @p grid for the given @p degree in all directions and homogeneous
     * in all components.
     */
    static std::shared_ptr<self_t>
    create_nonconst(const int degree,
                    const std::shared_ptr<GridType> &grid,
                    const InteriorReg interior_reg = InteriorReg::maximum,
                    const bool periodic = false,
                    const BasisEndBehaviour end_b = BasisEndBehaviour::interpolatory);

    /**
     * Builds and returns a maximum regularity (const) BSpline space
     * over a (const) CartesianGrid
     * @p grid for the given @p degree in all directions and homogeneous
     * in all components.
     */
    static std::shared_ptr<const self_t>
    create(const int degree,
           const std::shared_ptr<const GridType> &grid,
           const InteriorReg interior_reg = InteriorReg::maximum,
           const bool periodic = false,
           const BasisEndBehaviour end_b = BasisEndBehaviour::interpolatory);

    /**
     * Builds and returns a maximum regularity (non-const) BSpline space
     * over a (non-const) CartesianGrid
     * @p grid for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    static std::shared_ptr<self_t>
    create_nonconst(const Degrees &degree,
                    const std::shared_ptr<GridType> &grid,
                    const InteriorReg interior_reg = InteriorReg::maximum,
                    const Periodicity &periodic = Periodicity(false),
                    const EndBehaviour &end_b = EndBehaviour(BasisEndBehaviour::interpolatory));

    /**
     * Builds and returns a maximum regularity (const) BSpline space
     * over a (const) CartesianGrid
     * @p grid for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    static std::shared_ptr<const self_t>
    create(const Degrees &degree,
           const std::shared_ptr<const GridType> &grid,
           const InteriorReg interior_reg = InteriorReg::maximum,
           const Periodicity &periodic = Periodicity(false),
           const EndBehaviour &end_b = EndBehaviour(BasisEndBehaviour::interpolatory));

    /**
     * Builds and returns a (non-const) BSpline space
     * over a (non-const) CartesianGrid
     * @p grid with the given multiplicity vector @p mult_vectors
     * for each component
     * and the given @p degree for each direction and for each
     * component.
     */
    static std::shared_ptr<self_t>
    create_nonconst(const DegreeTable &deg,
                    const std::shared_ptr<GridType> &grid,
                    const MultiplicityTable &interior_mult,
                    const PeriodicityTable &periodic,
                    const EndBehaviourTable &end_b);

    /**
     * Builds and returns a (const) BSpline space
     * over a (const) CartesianGrid
     * @p grid with the given multiplicity vector @p mult_vectors
     * for each component
     * and the given @p degree for each direction and for each
     * component.
     */
    static std::shared_ptr<const self_t>
    create(const DegreeTable &deg,
           const std::shared_ptr<const GridType> &grid,
           const MultiplicityTable &interior_mult,
           const PeriodicityTable &periodic,
           const EndBehaviourTable &end_b);

    static std::shared_ptr<self_t>
    create_nonconst(const std::shared_ptr<SpaceData> &space_data,
                    const EndBehaviourTable &end_b);

    static std::shared_ptr<const self_t>
    create(const std::shared_ptr<const SpaceData> &space_data,
           const EndBehaviourTable &end_b);
    ///@}

    /**
     * Create an element (defined on this space) with a given flat_index.
     */
    virtual std::shared_ptr<SpaceElement<dim_,0,range_,rank_,Transformation::h_grad> >
    create_element(const ListIt &index, const PropId &property) const override final;


    /** Destructor. */
    virtual ~BSplineSpace() = default;

protected:
    /** @name Constructors */
    ///@{
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    BSplineSpace() = default;

    /**
     * Constructs a maximum regularity BSpline space
     * over a (non-const) CartesianGrid
     * @p grid for the given @p degree in all directions and homogeneous
     * in all components.
     */
    explicit BSplineSpace(const int degree,
                          const std::shared_ptr<GridType> &grid,
                          const InteriorReg interior_reg,
                          const bool periodic,
                          const BasisEndBehaviour end_b);

    /**
     * Constructs a maximum regularity BSpline space
     * over a (const) CartesianGrid
     * @p grid for the given @p degree in all directions and homogeneous
     * in all components.
     */
    explicit BSplineSpace(const int degree,
                          const std::shared_ptr<const GridType> &grid,
                          const InteriorReg interior_reg,
                          const bool periodic,
                          const BasisEndBehaviour end_b);

    /**
     * Constructs a maximum regularity BSpline space over
     * over a (non-const) CartesianGrid
     * @p grid for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    explicit BSplineSpace(const Degrees &degree,
                          const std::shared_ptr<GridType> &grid,
                          const InteriorReg interior_reg,
                          const Periodicity &periodic,
                          const EndBehaviour &end_b);

    /**
     * Constructs a maximum regularity BSpline space over
     * over a (const) CartesianGrid
     * @p grid for the given @p degree[i] in the i-th direction and homogeneous
     * in all components.
     */
    explicit BSplineSpace(const Degrees &degree,
                          const std::shared_ptr<const GridType> &grid,
                          const InteriorReg interior_reg,
                          const Periodicity &periodic,
                          const EndBehaviour &end_b);

    /**
     * Constructs a BSpline space over
     * over a (non-const) CartesianGrid
     * @p grid with the given multiplicity vector @p mult_vectors
     * for each component
     * and the given @p degree for each direction and for each
     * component.
     */
    explicit BSplineSpace(const DegreeTable &deg,
                          const std::shared_ptr<GridType> &grid,
                          const MultiplicityTable &interior_mult,
                          const PeriodicityTable &periodic,
                          const EndBehaviourTable &end_b);

    /**
     * Constructs a BSpline space over
     * over a (const) CartesianGrid
     * @p grid with the given multiplicity vector @p mult_vectors
     * for each component
     * and the given @p degree for each direction and for each
     * component.
     */
    explicit BSplineSpace(const DegreeTable &deg,
                          const std::shared_ptr<const GridType> &grid,
                          const MultiplicityTable &interior_mult,
                          const PeriodicityTable &periodic,
                          const EndBehaviourTable &end_b);


    explicit BSplineSpace(const std::shared_ptr<SpaceData> &space_data,
                          const EndBehaviourTable &end_b);

    explicit BSplineSpace(const std::shared_ptr<const SpaceData> &space_data,
                          const EndBehaviourTable &end_b);

    /**
     * Copy constructor. Not allowed to be used.
     */
    BSplineSpace(const self_t &space) = delete;
    ///@}

    /** @name Assignment operators */
    ///@{
    /** Copy assignment. Not allowed to be used. */
    self_t &
    operator=(const self_t &space) = delete;
    ///@}

public:



    virtual const DegreeTable &get_degree_table() const override final;


    virtual void get_element_dofs(
        const IndexType element_id,
        SafeSTLVector<Index> &dofs_global,
        SafeSTLVector<Index> &dofs_local_to_patch,
        SafeSTLVector<Index> &dofs_local_to_elem,
        const std::string &dofs_property = DofProperties::active) const override final;


    template <int sdim>
    using SubGridMap = typename GridType::template SubGridMap<sdim>;

//    using InterGridMap = std::map<Index,Index>;

    template <int k>
    using InterSpaceMap = SafeSTLVector<Index>;

    template <int k>
    using SubRefSpace = ReferenceSpace<k, range, rank>;

    template <int k>
    using SubSpace = PhysicalSpace<k,range,rank,dim-k, Transformation::h_grad>;

    /**
     * Construct a sub space of dimension k conforming to
     * the subspace sub element sub_elem_id and a map from the elements of
     * the sub_element grid to the corresponding element of the current
     * grid.
     */
    template<int k>
    std::shared_ptr<SubRefSpace<k> >
    get_ref_sub_space(const int sub_elem_id,
                      InterSpaceMap<k> &dof_map,
                      std::shared_ptr<CartesianGrid<k>> sub_grid = nullptr) const;

    template<int k>
    std::shared_ptr<SubSpace<k> >
    get_sub_space(const int s_id, InterSpaceMap<k> &dof_map,
                  std::shared_ptr<CartesianGrid<k>> sub_grid,
                  SubGridMap<k> &elem_map) const;







public:
    const PeriodicityTable &get_periodicity() const override final;

    /**
     * Returns a reference to the end behaviour table of the BSpline space.
     */
    virtual const EndBehaviourTable &get_end_behaviour_table() const override final;

    /**
     * Prints internal information about the space.
     * @note Mostly used for debugging and testing.
     */
    virtual void print_info(LogStream &out) const override final;

    std::shared_ptr<const DofDistribution<dim, range, rank> >
    get_ptr_const_dof_distribution() const override final;

    std::shared_ptr<DofDistribution<dim, range, rank> >
    get_ptr_dof_distribution() override final;




private:

    SharedPtrConstnessHandler<SpaceData > space_data_;


    EndBehaviourTable end_b_;


    /** Bezier extraction operator. */
    BernsteinExtraction<dim, range, rank> operators_;


    /** If end knots are not in the repeated knot vector */
    using EndIntervalTable = typename BaseSpace::template
                             ComponentContainer<SafeSTLArray<std::pair<Real, Real>, dim>>;
    EndIntervalTable end_interval_;

    /**
     * Container with the local to global basis indices
     * @note The concept of global indices refers to a global numbering of the
     * dofs of all the spaces.
     */
    std::shared_ptr<DofDistribution<dim,range,rank> > dof_distribution_;


    friend class BSplineElement<dim, range, rank>;
    friend class BSplineElementHandler<dim, range, rank>;

#ifdef MESH_REFINEMENT
    /**
     * Rebuild the internal state of the object after an insert_knots() function is invoked.
     *
     * @pre Before invoking this function, must be invoked the function grid_->insert_knots().
     * @note This function is connected to the CartesianGrid's signal for the refinement, and
     * it is necessary in order to avoid infinite loops in the insert_knots() function calls.
     *
     * @ingroup h_refinement
     */
    void rebuild_after_insert_knots(
        const SafeSTLArray<SafeSTLVector<Real>,dim> &knots_to_insert,
        const CartesianGrid<dim> &old_grid);

    void create_connection_for_insert_knots(std::shared_ptr<self_t> space);
#endif // MESH_REFINEMENT

    /**
     * Returns the current object wrapped by a std::shared_ptr.
     *
     * @note Internally uses the shared_from_this() function.
     */
    std::shared_ptr<const self_t > get_this_space() const;


public:
    DeclException1(ExcScalarRange, int,
                   << "Range " << arg1 << "should be 0 for a scalar valued"
                   << " space.");


    virtual bool is_bspline() const override final;

    virtual std::shared_ptr<SpaceElementHandler<dim_,0,range_,rank_,Transformation::h_grad>>
            get_elem_handler() const override final;


private:
#if 0
    /**
     * Lookup table for the local dof id in each element component
     */
    typename SpaceData::template ComponentContainer<SafeSTLVector<TensorIndex<dim> > >
    dofs_tensor_id_elem_table_;
#endif

#ifdef SERIALIZATION
    /**
     * @name Functions needed for boost::serialization
     * @see <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     */
    ///@{
    friend class boost::serialization::access;

    template<class Archive>
    void
    serialize(Archive &ar, const unsigned int version);
    ///@}
#endif // SERIALIZATION
};


IGA_NAMESPACE_CLOSE


#endif /* __BSPLINE_SPACE_H_ */

