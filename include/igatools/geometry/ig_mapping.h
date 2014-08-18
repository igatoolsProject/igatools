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

#ifndef IG_MAPPING_H_
#define IG_MAPPING_H_

#include <igatools/base/config.h>
#include <igatools/geometry/spline_mapping.h>
#include <igatools/utils/dynamic_multi_array.h>
#include <igatools/utils/static_multi_array.h>
#include <igatools/basis_functions/nurbs_space.h>

IGA_NAMESPACE_OPEN

template <class RefSpace>
class IgMapping
    : public SplineMapping<RefSpace>
{
private:
    using base_t = SplineMapping<RefSpace>;

    using base_t::dim;
    using base_t::codim;
    using base_t::space_dim;

    using typename base_t::Point;
    using typename base_t::Value;
    using typename base_t::Gradient;
    using typename base_t::Hessian;
    using typename base_t::GridType;
    using typename base_t::GridIterator;

    using typename base_t::ElementIterator;

    using self_t = IgMapping<RefSpace>;

public:
    /**
     * Default constructor.
     */
    IgMapping() = delete;

    /**
     * \brief Constructor. It builds a mapping from a function space and a vector of
     * control points.
     * \param[in] space The function space (i.e. a set of basis function) used to represents the mapping.
     * \param[in] control_points The coefficients of the linear combination of the basis function
     * of the function space used to represents the mapping.
     */
    IgMapping(const std::shared_ptr<RefSpace> space,
              const std::vector<Real> &control_points);


    /**
     * It builds a Mapping object wrapped in a std::shared_ptr,
     * from a function space and a vector of control points.
     */
    static std::shared_ptr<Mapping<dim,codim>>
                                            create(const std::shared_ptr<RefSpace> space, const std::vector<Real> &control_points);

    /**
     * Copy constructor. Performs a deep copy of the object.
     */
    IgMapping(const self_t &map);

    /**
     * Copy assignment operator.
     */
    self_t &operator=(const self_t &map) = delete;

    void init_element(const ValueFlags flag,
                      const Quadrature<dim> &quad) const override;

    void set_element(const GridIterator &elem) const override;

    void set_face_element(const Index face_id,
                          const GridIterator &elem) const override;

    /** @name Mapping as a standard function */
    ///@{
    virtual void evaluate(std::vector<Value> &values) const override;

    virtual void evaluate_gradients
    (std::vector<Gradient> &gradients) const override;

    virtual void evaluate_hessians
    (std::vector<Hessian> &hessians) const override;

    virtual void evaluate_face
    (const Index face_id, std::vector<Value> &values) const override;

    virtual void evaluate_face_gradients
    (const Index face_id, std::vector<Gradient> &gradients) const override;

    virtual void evaluate_face_hessians
    (const Index face_id, std::vector<Hessian> &hessians) const override;

    ///@}

    /** @name Function used to modify the position of the control points */
    ///@{
    /**
     * Sets the control points defining the map.
     * @param[in] control_points - Coordinates of the control points in the Euclidean space.
     */
    void set_control_points(const std::vector<Real> &control_points) override final;
    ///@}

    std::shared_ptr<RefSpace> get_iga_space() override final
    {
        return data_->ref_space_;
    }

    /**
     * Prints internal information about the mapping.
     * @note Mostly used for debugging and testing.
     * Try to call the same function on a derived class.
     */
    void print_info(LogStream &out) const override;


    /** @name Dealing with the element-based iterator. */
    ///@{
    /**
     * Returns a element iterator to the first element of the patch.
     */
    ElementIterator begin() const override final ;

    /**
     * Returns a element iterator to the last element of the patch.
     */
    ElementIterator last() const override final ;

    /**
     * Returns a element iterator to one-pass the end of patch.
     */
    ElementIterator end() const override final;
    ///@}


    /** @name Evaluating the quantities related to the IgMapping without the use of the cache. */
    ///@{
    void evaluate_at_points(const std::vector<Point> &points, std::vector<Value> &values) const override final;
    void evaluate_gradients_at_points(const std::vector<Point> &points, std::vector<Gradient> &gradients) const override final;
    void evaluate_hessians_at_points(const std::vector<Point> &points, std::vector<Hessian> &hessians) const override final;
    ///@}


private:


    template< class T >
    using ComponentTable = StaticMultiArray<T,base_t::space_dim,1>;


    class IgMappingData
    {
    public:
        /** Coordinates of the control points in the Euclidean space. */
        std::vector<Real> control_points_;

        /**
         * Weights associated with the control points (if NURBSpace is used).
         *
         * @note The weights are necessary in order to perform the h-refinement,
         * because the control points are in the euclidean space, while the
         * h-refinement algorithm (based on knot insertion) require them to be in
         * the projective space.
         */
        typename NURBSSpace<dim,space_dim,1>::WeightsTable weights_pre_refinement_;


        /** Control mesh (the coordinates are in the projective space). */
        ComponentTable<DynamicMultiArray<Real,dim>> ctrl_mesh_;


        /** The function space used to represents the mapping.*/
        std::shared_ptr<RefSpace> ref_space_;
    };


    /**
     * Data that defines the mapping, i.e. control points, weights, knots etc.
     */
    std::shared_ptr<IgMappingData> data_;

    /**
     * This variable represents the cache, that for the IgMapping is an element accessor
     * of the reference space upon which the IgMapping is built upon.
     * @note This cache should be freshly created every time a MappingElementAccessor is
     * referring to the IgMapping object.
     */
    mutable typename RefSpace::ElementIterator cache_;

public:

    /**
     * Returns the pointer wrapping the IgMappingData (i.e. control points, weights, etc.
     * except the cache.
     */
    std::shared_ptr<IgMappingData> get_data() const;


    /**
     * Constructor. Builds an igMapping with new cache, providing the data.
     * @note This is done in order to have the same mapping data shared between many IgMapping
     * objects but with different caches in order to be used independently.
     * The cache is initialized with RefSpace::begin().
     */
    IgMapping(const std::shared_ptr<IgMappingData> mapping_data);


private:
    /**
     * h-refines the control mesh of the mapping space after a grid uniform refinement.
     *
     * @param[in] refinement_directions Directions along which the refinement is performed.
     * @param[in] grid_old Grid before the refinement.
     *
     * @pre Before invoking this function, must be invoked the function grid_->refine().
     * @note This function is connected to the CartesianGrid's signal for the refinement, and
     * it is necessary in order to avoid infinite loops in the refine() function calls.
     * @note The implementation of this function is based on "The NURBS Book" Algorithm A5.4.
     *
     * @ingroup h_refinement
     *
     */
    void refine_h_control_mesh(
        const std::array<bool,dim> &refinement_directions,
        const typename base_t::GridType &grid_old);

    /**
     * Returns the control points that are active on the element represented by the cache.
     */
    std::vector<Real> get_control_points_elem() const;


};

IGA_NAMESPACE_CLOSE

#endif // #ifndef IG_MAPPING_
