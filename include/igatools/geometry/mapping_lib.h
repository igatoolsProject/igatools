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

#ifndef MAPPING_LIB_H_
#define MAPPING_LIB_H_

#include <igatools/base/config.h>
#include <igatools/geometry/analytical_mapping.h>

IGA_NAMESPACE_OPEN



/**
 * @brief Affine mapping from \f$\mathbb{R}^{dim}\f$ to \f$\mathbb{R}^{dim+codim}\f$.
 *
 * \f[
 * \mathbf{x} = \mathbf{F}(\mathbf{\hat{x}}) = \mathbf{A} \mathbf{\hat{x}} + \mathbf{b}.
 * \f]
 *
 * @author Pauletti, 2012
 * @todo Document more. Missing member documentation.
 */
template<int dim_, int codim_>
class LinearMapping : public AnalyticalMapping <dim_, codim_>
{
private:
    //TODO(pauletti, Apr 26, 2014): base_t should be mapping
    using base_t = AnalyticalMapping<dim_, codim_>;

    using base_t::dim;
    using base_t::codim;
    using base_t::space_dim;

    using typename base_t::Point;
    using typename base_t::Value;
    using typename base_t::Gradient;
    using typename base_t::Hessian;
    using typename base_t::GridType;
    using typename base_t::GridIterator;

    using self_t = LinearMapping<dim, codim>;

public:
    LinearMapping() = delete;

    LinearMapping(const Gradient &A, const Value &b);

    LinearMapping(const std::shared_ptr<GridType> grid,
                  const Gradient &A,
                  const Value &b);
    /**
     * This function creates a LinearMapping object
     * and wraps it into a shared_ptr.
     */
    static std::shared_ptr<base_t> create(const Gradient &A, const Value &b);

    static std::shared_ptr<base_t>
    create(const std::shared_ptr<GridType> grid,
           const Gradient &A, const Value &b);

    ValueFlags required_flags() const;

    void set_element(const GridIterator &elem) const;

    void set_face_element(const Index face_id,
                          const GridIterator &elem) const;

    void evaluate(std::vector<Value> &values) const override;

    void evaluate_gradients(std::vector<Gradient> &gradients) const override;

    void evaluate_hessians(std::vector<Hessian> &hessians) const override;

    void evaluate_face(const Index face_id,std::vector<Value> &values) const override;

    void evaluate_face_gradients(const Index face_id,std::vector<Gradient> &gradients) const override;

    void evaluate_face_hessians(const Index face_id,std::vector<Hessian> &hessians) const override;

private:
    const Gradient A_;
    const Value b_;


    // TODO (pauletti, Apr 23, 2014): why not have this in the base class?
    mutable std::vector<Point> points_;
    mutable std::array<std::vector<Point>, UnitElement<dim_>::faces_per_element> face_points_;
};


// TODO (pauletti, Nov 13, 2013): Document this in a proper way
/**
 * Maps a hyper rectangle into a spherical ball sector using the
 * dim-dimensional spherical coordinates, maps a hyper-rectangle
 * r in [0,R], phi_1 in [0, 2 pi], and phi_2, phi_dim-1 in [0,pi]
 * such that
 * x1 = r cos (phi_1)
 * x2 = r sin (phi_1) cos (phi_2)
 * etc
 *
 *
 * @todo Missing documentation
 */
template<int dim_>
class BallMapping : public AnalyticalMapping<dim_, 0>
{
private:
    using base_t = AnalyticalMapping<dim_, 0>;

    using base_t::dim;
    using base_t::codim;
    using base_t::space_dim;

    using typename base_t::Point;
    using typename base_t::Value;
    using typename base_t::Gradient;
    using typename base_t::Hessian;
    using typename base_t::GridType;
    using typename base_t::GridIterator;

    using self_t = BallMapping<dim>;

public:
    BallMapping(const std::shared_ptr<GridType> grid);

    static std::shared_ptr<base_t>
    create(const std::shared_ptr<GridType> grid);

    ValueFlags required_flags() const;

    void set_element(const GridIterator &elem)  const;

    void set_face_element(const Index face_id, const GridIterator &elem) const;

    void evaluate(std::vector<Value> &values) const override;

    void evaluate_gradients(std::vector<Gradient> &gradients) const override;

    void evaluate_hessians(std::vector<Hessian> &hessians) const override;

    void evaluate_face(const Index face_id, std::vector<Value> &values) const override;

    void evaluate_face_gradients(const Index face_id, std::vector<Gradient> &gradients) const override;

    void evaluate_face_hessians(const Index face_id, std::vector<Hessian> &hessians) const override;


    /** @name Evaluating the quantities related to BallMapping without the use of the cache. */
    ///@{
    void evaluate_gradients_at_points(const std::vector<Point> &points, std::vector<Gradient> &gradients) const override final;
    void evaluate_hessians_at_points(const std::vector<Point> &points, std::vector<Hessian> &hessians) const override final;
    ///@}

private:
    static const int order = 3;

    mutable std::vector<Points<dim>> points_;
    mutable std::array<std::vector<Point>, UnitElement<dim>::faces_per_element> face_points_;
    mutable std::array<std::vector<std::array<double, dim> >, order> cos_val;
    mutable std::array<std::vector<std::array<double, dim> >, order> sin_val;
    mutable std::array<std::array<std::vector<std::array<double, dim> >, order>,
            UnitElement<dim>::faces_per_element> face_cos_val;
    mutable std::array<std::array<std::vector<std::array<double, dim> >, order>,
            UnitElement<dim>::faces_per_element> face_sin_val;

};



// TODO (pauletti, Nov 21, 2013): Document this in a proper way
/**
 * @todo Missing documentation
 */
template<int dim_>
class SphereMapping : public AnalyticalMapping<dim_, 1>
{
private:
    using base_t = AnalyticalMapping<dim_, 1>;

    using base_t::dim;
    using base_t::codim;
    using base_t::space_dim;

    using typename base_t::Point;
    using typename base_t::Value;
    using typename base_t::Gradient;
    using typename base_t::Hessian;
    using typename base_t::GridType;
    using typename base_t::GridIterator;

    using self_t = SphereMapping<dim>;

public:
    SphereMapping(const std::shared_ptr<GridType> grid, const Real R = 1.);

    static std::shared_ptr<base_t> create(const std::shared_ptr<GridType> grid);

    ValueFlags required_flags() const;

    void set_element(const GridIterator &elem)  const;

    void set_face_element(const Index face_id, const GridIterator &elem) const;

    void evaluate(std::vector<Value> &values) const override;

    void evaluate_gradients(std::vector<Gradient> &gradients) const override;

    void evaluate_hessians(std::vector<Hessian> &hessians) const override;

    void evaluate_face(const Index face_id, std::vector<Value> &values) const override;

    void evaluate_face_gradients(const Index face_id, std::vector<Gradient> &gradients) const override;

    void evaluate_face_hessians(const Index face_id, std::vector<Hessian> &hessians) const override;

private:
    static const int order = 3;
    mutable std::vector<Points<dim>> points_;
    mutable std::array<std::vector<Point>, UnitElement<dim>::faces_per_element> face_points_;
    mutable std::array<std::vector<std::array<double, space_dim> >, order> cos_val;
    mutable std::array<std::vector<std::array<double, space_dim> >, order> sin_val;
    mutable std::array<std::array<std::vector<std::array<double, dim> >, order>,
            UnitElement<dim>::faces_per_element> face_cos_val;
    mutable std::array<std::array<std::vector<std::array<double, dim> >, order>,
            UnitElement<dim>::faces_per_element> face_sin_val;
    const double R = 1.;
};



/**
 * \brief This class represent a cylindrical annulus mapping.
 *
 * The mapping is
 * \f{equation*}{
 *    \begin{aligned}
 *    F(\hat{\theta},\hat{r},\hat{z}) \colon [0,1] \times[0,1]\times[0,1] & \to \Omega \\
 *    (\hat{\theta},\hat{r},\hat{z}) & \mapsto
 *    F(\hat{\theta},\hat{r},\hat{z}) =
 *    \begin{pmatrix}
 *      \bigl[ (r_1-r_0) \hat{r} + r_0 \bigr] \cos\bigl[ (\theta_1-\theta_0) \hat{\theta} \bigr] \\
 *      \bigl[ (r_1-r_0) \hat{r} + r_0 \bigr] \sin\bigl[ (\theta_1-\theta_0) \hat{\theta} \bigr] \\
 *      h_0 + (h_1-h_0) \hat{z}
 *    \end{pmatrix}
 *    \end{aligned}
 * \f}
 * where \f$ \Omega \f$ is a section of cylindrical annulus with the following characteristics:
 *   section angle \f$ \theta \in [\theta_0,\theta_1] \f$,
 *   radius \f$ r \in [r_0,r_1] \f$,
 *   height \f$ z \in[h_0,h_1] \f$.
 *
 * \author M.Martinelli
 * \date 31 Jan 2013
 */
class CylindricalAnnulus : public AnalyticalMapping<3, 0>
{
private:
    using base_t = AnalyticalMapping<3, 0>;

    using base_t::dim;
    using base_t::codim;
    using base_t::space_dim;

    using typename base_t::Point;
    using typename base_t::Value;
    using typename base_t::Gradient;
    using typename base_t::Hessian;
    using typename base_t::GridType;
    using typename base_t::GridIterator;

    using self_t = CylindricalAnnulus;

public:
    CylindricalAnnulus(
        const Real r0,
        const Real r1,
        const Real h0,
        const Real h1,
        const Real theta0,
        const Real theta1);

    CylindricalAnnulus(
        const Real r0,
        const Real r1,
        const Real h1,
        const Real theta1);

    CylindricalAnnulus(
        const std::shared_ptr<GridType> grid,
        const Real r0,
        const Real r1,
        const Real h0,
        const Real h1,
        const Real theta0,
        const Real theta1);

    CylindricalAnnulus(
        const std::shared_ptr<GridType> grid,
        const Real r0,
        const Real r1,
        const Real h1,
        const Real theta1);

    /**
     * This function creates a CylindricalAnnulus object
     * and wraps it into a shared_ptr.
     */
    static std::shared_ptr<base_t> create(
        const Real r0,
        const Real r1,
        const Real h0,
        const Real h1,
        const Real theta0,
        const Real theta1);

    static std::shared_ptr<base_t> create(
        const Real r0,
        const Real r1,
        const Real h1,
        const Real theta1);

    static std::shared_ptr<base_t> create(
        const std::shared_ptr<GridType> grid,
        const Real r0,
        const Real r1,
        const Real h0,
        const Real h1,
        const Real theta0,
        const Real theta1);

    static std::shared_ptr<base_t> create(
        const std::shared_ptr<GridType> grid,
        const Real r0,
        const Real r1,
        const Real h1,
        const Real theta1);

    ValueFlags required_flags() const;

    void set_element(const GridIterator &elem) const;

    void set_face_element(const Index face_id, const GridIterator &elem) const;

    void evaluate(std::vector<Point> &values) const override;

    void evaluate_gradients(std::vector<Gradient> &gradients) const override;

    void evaluate_hessians(std::vector<Hessian> &hessians) const override;

    void evaluate_face(const Index face_id, std::vector<Point> &values) const override;

    void evaluate_face_gradients(const Index face_id, std::vector<Gradient> &gradients) const override;

    void evaluate_face_hessians(const Index face_id, std::vector<Hessian> &hessians) const override;

    /** @name Evaluating the quantities related to CylindricalAnnulus without the use of the cache. */
    ///@{
    void evaluate_at_points(const std::vector<Point> &points, std::vector<Value> &values) const override final;
    void evaluate_gradients_at_points(const std::vector<Point> &points, std::vector<Gradient> &gradients) const override final;
    void evaluate_hessians_at_points(const std::vector<Point> &points, std::vector<Hessian> &hessians) const override final;
    ///@}

private:
    const Real r0_;
    const Real r1_;
    const Real h0_;
    const Real h1_;
    const Real theta0_;
    const Real theta1_;

    const Real dR_;
    const Real dT_;
    const Real dH_;
    mutable std::vector<Point> points_;
    mutable std::array<std::vector<Point>, UnitElement<3>::faces_per_element> face_points_;
};

IGA_NAMESPACE_CLOSE


#endif /* MAPPING_LIB_H_ */
