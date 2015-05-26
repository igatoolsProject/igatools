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

#ifndef NURBS_SPACE_H_
#define NURBS_SPACE_H_

#include <igatools/base/config.h>

#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element.h>
#include <igatools/functions/ig_function.h>

#ifdef NURBS

IGA_NAMESPACE_OPEN

template <int, int, int> class NURBSElement;
template <int, int, int> class NURBSElementHandler;

/**
 * Multivariate (tensor product) scalar, vector or k-tensor
 * valued NURBS space.
 *
 * @ingroup containers
 * @ingroup serializable
 */
template <int dim_, int range_ = 1, int rank_ = 1>
class NURBSSpace :
//    public std::enable_shared_from_this<NURBSSpace<dim_,range_,rank_> >,
    public ReferenceSpace<dim_,range_,rank_>
{
private:
    using BaseSpace = ReferenceSpace<dim_,range_,rank_>;
    using self_t = NURBSSpace<dim_, range_, rank_>;

public:
    using SpSpace = BSplineSpace<dim_, range_, rank_>;

    /** see documentation in \ref Space */

    using GridType = CartesianGrid<dim_>;
    static const int dim       = dim_;
    static const int codim     = 0;
    static const int space_dim = dim_;
    static const int range     = range_;
    static const int rank      = rank_;
    static const bool is_physical_space = false;

    static const auto n_components = SplineSpace<dim_, range_, rank_>::n_components;


    /**
     * See documentation in \ref Space
     *
     * @see Space
     */
    using PushForwardType = typename BaseSpace::PushForwardType;

    using RefSpace = typename BaseSpace::RefSpace;



public:
    using Func = typename SpSpace::Func;
    template <int order>
    using Derivative = typename SpSpace::template Derivative<order>;
    using Point = typename SpSpace::Point;
    using Value = typename SpSpace::Value;
    using Div   = typename SpSpace::Div;

    using RefPoint = typename SpSpace::RefPoint;

public:

    /** Type for the element accessor. */
    using ElementAccessor = NURBSElement<dim, range, rank> ;

    /** Type for iterator over the elements.  */
    using ElementIterator = CartesianGridIterator<ReferenceElement<dim,range,rank> >;

    using ElementHandler = NURBSElementHandler<dim_, range_, rank_>;



    template <int k>
    using InterGridMap = typename GridType::template InterGridMap<k>;

    template <int k>
    using InterSpaceMap = SafeSTLVector<Index>;

    template <int k>
    using SubRefSpace = NURBSSpace<k, range, rank>;

    template <int k>
    using SubSpace = PhysicalSpace<k,range,rank, dim-k, Transformation::h_grad>;

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
                  std::shared_ptr<typename GridType::template InterGridMap<k>> elem_map) const;

public:
//    /** Container indexed by the components of the space */
    template< class T>
    using ComponentContainer = typename SpSpace::template ComponentContainer<T>;

    using Degrees = typename SpSpace::Degrees;
    using Multiplicity = typename SpSpace::Multiplicity;
    using EndBehaviour = typename SpSpace::EndBehaviour;
    using Periodicity = typename SpSpace::Periodicity;

    using KnotsTable = typename SpSpace::KnotsTable;
    using DegreeTable = typename SpSpace::DegreeTable;
    using MultiplicityTable = typename SpSpace::MultiplicityTable;
    using TensorSizeTable = typename SpSpace::TensorSizeTable;
    using PeriodicityTable = typename SpSpace::PeriodicityTable;
    using EndBehaviourTable = typename SpSpace::EndBehaviourTable;


    using WeightSpace = BSplineSpace<dim_,1,1>;
    using WeightFunction = IgFunction<dim_,0,1,1>;
    using WeightFunctionPtr = std::shared_ptr<WeightFunction>;
    using WeightFunctionPtrTable = ComponentContainer<WeightFunctionPtr>;
    using Weights = DynamicMultiArray<Real,dim>;
    using WeightsTable = ComponentContainer<Weights>;

public:
    /**
     * Returns a shared_ptr wrapping a NURBSSpace from a BSplineSpace and a scalar weight function.
     */
    static std::shared_ptr<self_t>
    create(std::shared_ptr<SpSpace> bs_space,
           const WeightFunctionPtrTable &weight_func_table);

    ///@}

    /**
     * Create an element (defined on this space) with a given flat_index.
     */
    virtual std::shared_ptr<SpaceElement<dim_,0,range_,rank_> >
    create_element(const Index flat_index) const override final;

    /** Destructor */
    virtual ~NURBSSpace() = default;

protected:
    /** @name Constructor */
    ///@{
    /**
     * Default constructor. It does nothing but it is needed for the
     * <a href="http://www.boost.org/doc/libs/release/libs/serialization/">boost::serialization</a>
     * mechanism.
     */
    NURBSSpace() = default;

    /**
     * Construct a NURBSSpace from a BSplineSpace and a table of weights.
     */
    explicit NURBSSpace(std::shared_ptr<SpSpace> bs_space,
                        const WeightFunctionPtrTable &weight_func_table);

    /**
     * Copy constructor. Not allowed to be used.
     */
    NURBSSpace(const self_t &space) = delete;

    ///@}

public:
    /** @name Getting information about the space */
    ///@{

    virtual bool is_bspline() const override final;

    virtual const DegreeTable &get_degree() const override final;

    virtual void get_element_dofs(
        const Index element_id,
        SafeSTLVector<Index> &dofs_global,
        SafeSTLVector<Index> &dofs_local_to_patch,
        SafeSTLVector<Index> &dofs_local_to_elem,
        const std::string &dofs_property = DofProperties::active) const override final;


    std::set<Index> get_interior_dofs() const;

    using typename BaseSpace::topology_variant;

    std::set<Index>
    get_boundary_dofs(const int s_id, const topology_variant &k) const override final;

    ///@}

    const std::shared_ptr<SpSpace> get_spline_space() const;



    /**
     * Get the weights of the NURBSSpace.
     */
    const WeightFunctionPtrTable &get_weights() const;

    const PeriodicityTable &get_periodicity() const override final;


    /**
     * Returns a const reference to the end behaviour table of the BSpline space.
     */
    virtual const EndBehaviourTable &get_end_behaviour_table() const override final;


    /**
     * Prints internal information about the space.
     * @note Mostly used for debugging and testing.
     */
    virtual void print_info(LogStream &out) const override final;

private:
    /**
     * B-spline space
     */
    std::shared_ptr<SpSpace> sp_space_;

    /**
     * Weight function component table. Each component of the table is a shared_ptr to an IgFunction
     * built using a scalar BSplineSpace.
     */
    WeightFunctionPtrTable weight_func_table_;


    friend ElementAccessor;
    friend ElementHandler;

    /**
     * Returns the weight coefficient associated with a given basis function.
     */
    Real get_weight_coef_from_basis_id(const Index basis_id) const;


public:
    virtual std::shared_ptr<SpaceElementHandler<dim_,0,range_,rank_>>
            get_elem_handler() const override final;



private:

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
#endif /* #ifdef NURBS */

#endif /* NURBS_SPACE_H_ */


