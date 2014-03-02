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
 * @ingroup refinement
 */
template <int dim_, int dim_range_ = 1, int rank_ = 1>
class NURBSSpace : public BSplineSpace<dim_, dim_range_, rank_>
{
public:
    using base_t = BSplineSpace<dim_, dim_range_, rank_>;

    /** Type for the grid. */
    using typename base_t::GridType;
    using base_t::dim;
    using base_t::space_dim;
    using base_t::dim_range;
    using base_t::rank;
    using base_t::n_components;
    using Multiplicities = StaticMultiArray<Multiplicity<dim_>,dim_range_,rank_>;

private:
    using self_t = NURBSSpace<dim, dim_range, rank>;

public:
    static const bool has_weights = true;
    /**
     * Type for element accessor.
     */
    typedef NURBSElementAccessor<dim, dim_range, rank> ElementAccessor;

    /**
     * Type for iterator over the elements.
     */
    typedef PatchIterator<ElementAccessor> ElementIterator;

    /**
     * Type for the face space.
     */
    //TODO rename FaceSpace_t to face_space_t
    typedef NURBSSpace<dim-1, dim_range, rank> FaceSpace_t;


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
        const StaticMultiArray<std::array<int, dim>, dim_range, rank> &degree);

    /**
     * Returns a shared_ptr wrapping a maximum regularity NURBSSpace over CartesianGrid
     * @p knots for the given @p degree for each direction and for each component.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    static std::shared_ptr< self_t >
    create(std::shared_ptr<GridType> knots,
           const StaticMultiArray<std::array<int,dim>,dim_range,rank> &degree);

    /**
     * Construct a NURBSSpace over the CartesianGrid @p knots with
     * the given multiplicity vector @p mult_vector for each component
     * and for the given @p degree for each direction and for each component.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    NURBSSpace(
        std::shared_ptr< GridType > knots,
        const Multiplicities &mult_vector,
        const StaticMultiArray<std::array<int,dim>,dim_range,rank> &degree);

    /**
     * Returns shared_ptr wrapping a NURBSSpace over the CartesianGrid @p knots with
     * the given multiplicity vector @p mult_vector for each component
     * and for the given @p degree for each direction and for each component.
     * @note All weights are set to 1.0, so the resulting space has the same structure of a BSpline space.
     */
    static std::shared_ptr< self_t >
    create(std::shared_ptr< GridType > knots,
           const Multiplicities &mult_vector,
           const StaticMultiArray<std::array<int,dim>,dim_range,rank> &degree);

    /**
     * Construct a NURBSSpace over the CartesianGrid @p knots with
     * the given multiplicity vector @p mult_vector for each component
     * and for the given @p degree for each direction and for each component.
     */
    NURBSSpace(
        std::shared_ptr< GridType > knots,
        const Multiplicities &mult_vector,
        const StaticMultiArray<std::array<int,dim>,dim_range,rank> &degree,
        const StaticMultiArray<DynamicMultiArray<Real,dim>,dim_range,rank> &weights);

    /**
     * Returns a shared_ptr wrapping a NURBSSpace over the CartesianGrid @p knots with
     * the given multiplicity vector @p mult_vector for each component
     * and for the given @p degree for each direction and for each component.
     */
    static std::shared_ptr< self_t >
    create(std::shared_ptr< GridType > knots,
           const Multiplicities &mult_vector,
           const StaticMultiArray<std::array<int,dim>,dim_range,rank> &degree,
           const StaticMultiArray<DynamicMultiArray<Real,dim>,dim_range,rank> &weights);

    /** Destructor */
    ~NURBSSpace() = default;

    ///@}

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
    const StaticMultiArray< DynamicMultiArray<Real,dim>,dim_range,rank>
    get_weights() const;

    /**
     * Reset the weights of the NURBSSpace.
     */
    void reset_weights(const StaticMultiArray<DynamicMultiArray<iga::Real,dim>,dim_range,rank> &weights);

private:
    /**
     * Weights associated to the basis functions.
     */
    StaticMultiArray<DynamicMultiArray<iga::Real,dim>,dim_range,rank> weights_;

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
     */
    void refine_h_weights(
        const std::array<bool,dim> &refinement_directions,
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


