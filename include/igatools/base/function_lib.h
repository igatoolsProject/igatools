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

#ifndef __IGA_FUNCTION_LIB_H_
#define __IGA_FUNCTION_LIB_H_

#include <igatools/base/formula_function.h>

IGA_NAMESPACE_OPEN

/**
 * Collection of useful functions derived from the Function class.
 */
namespace functions
{
template<int dim, int codim, int range, int rank = 1>
class ConstantFunction : public FormulaFunction<dim, codim, range, rank>
{
private:
    using base_t = Function<dim, codim, range, rank>;
    using parent_t = FormulaFunction<dim, codim, range, rank>;
    using self_t = ConstantFunction<dim, codim, range, rank>;
    using typename base_t::GridType;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;
    using typename parent_t::Map;

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map,
           const Value &b);

    std::shared_ptr<base_t> clone() const override final
    {
        return std::make_shared<self_t>(self_t(*this));
    }


    ConstantFunction(const self_t &) = default;

protected:
    ConstantFunction(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map,
                     const Value &b);

private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const override;

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const override;

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const override;

private:
    const Value b_;
};



//------------------------------------------------------------------------------

template<int dim, int codim, int range>
class LinearFunction :
    public FormulaFunction<dim, codim, range, 1>
//    public std::enable_shared_from_this<LinearFunction<dim,codim,range> >
{

public:
    using base_t = Function<dim, codim, range, 1>;
    using parent_t = FormulaFunction<dim, codim, range, 1>;
    using self_t = LinearFunction<dim, codim, range>;
    using typename base_t::GridType;
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;
    using typename parent_t::Map;

private:
#if 0
    std::shared_ptr<const base_t> shared_from_derived() const override final
    {
        return this->shared_from_this();
    }
#endif

public:
    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map,
           const Gradient &A,
           const Value &b);

    virtual std::shared_ptr<base_t> clone() const override final
    {
        return std::make_shared<self_t>(self_t(*this));
    }

    LinearFunction(const self_t &) = default;


protected:
    LinearFunction(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map,
                   const Gradient &A, const Value &b);



private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const override;

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const override;

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const override;

private:
    const Gradient A_;
    const Value    b_;
};

//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
/**
 * Maps a hyper rectangle into a spherical ball sector using the
 * dim-dimensional spherical coordinates, maps a hyper-rectangle
 * r in [0,R], phi_1 in [0, 2 pi], and phi_2, phi_dim-1 in [0,pi]
 * such that
 * x1 = r cos (phi_1)
 * x2 = r sin (phi_1) cos (phi_2)
 * etc
 *
 */
template<int dim>
class BallFunction : public FormulaFunction<dim, 0, dim, 1>
{
private:
    using base_t = Function<dim, 0, dim, 1>;
    using parent_t = FormulaFunction<dim, 0, dim, 1>;
    using self_t = BallFunction<dim>;
    using typename base_t::GridType;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;
    using typename parent_t::Map;

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map);

    std::shared_ptr<base_t> clone() const override final
    {
        return std::make_shared<self_t>(self_t(*this));
    }

    BallFunction(const self_t &) = default;

protected:
    BallFunction(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map);

private:
    template<int order>
    auto
    get_aux_vals(const ValueVector<Point> &points) const;

private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const;

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const;

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const;

private:
    const Value b_;
};
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
/**
 * Maps a hyper rectangle into a spherical ball sector using the
 * dim-dimensional spherical coordinates, maps a hyper-rectangle
 * r in [0,R], phi_1 in [0, 2 pi], and phi_2, phi_dim-1 in [0,pi]
 * such that
 * x1 = r cos (phi_1)
 * x2 = r sin (phi_1) cos (phi_2)
 * etc
 *
 */
template<int dim>
class SphereFunction : public FormulaFunction<dim, 0, dim+1, 1>
{
private:
    static const int space_dim = dim + 1;
    using base_t = Function<dim, 0, dim+1, 1>;
    using parent_t = FormulaFunction<dim, 0, dim+1, 1>;
    using self_t = SphereFunction<dim>;
    using typename base_t::GridType;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;
    using typename parent_t::Map;

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map);

    std::shared_ptr<base_t> clone() const override final
    {
        return std::make_shared<self_t>(self_t(*this));
    }

    SphereFunction(const self_t &) = default;

protected:
    SphereFunction(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map);

private:
    template<int order>
    auto
    get_aux_vals(const ValueVector<Point> &points) const;

private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const;

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const;

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const;

    static const int R = 1.;

};
//------------------------------------------------------------------------------





//------------------------------------------------------------------------------
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
template<int dim>
class CylindricalAnnulus : public FormulaFunction<dim, 0, dim, 1>
{
private:
    using base_t = Function<dim, 0, dim, 1>;
    using parent_t = FormulaFunction<dim, 0, dim, 1>;
    using self_t = CylindricalAnnulus;
    using typename base_t::GridType;
public:
    using typename parent_t::Point;
    using typename parent_t::Value;
    using typename parent_t::Gradient;
    using typename parent_t::ElementIterator;
    using typename parent_t::ElementAccessor;
    template <int order>
    using Derivative = typename parent_t::template Derivative<order>;
    using typename parent_t::Map;

    static std::shared_ptr<base_t>
    create(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map,
           const Real r0,
           const Real r1,
           const Real h0,
           const Real h1,
           const Real theta0,
           const Real theta1);

    std::shared_ptr<base_t> clone() const override final
    {
        return std::make_shared<self_t>(self_t(*this));
    }

    CylindricalAnnulus(const self_t &) = default;

protected:
    CylindricalAnnulus(std::shared_ptr<GridType> grid, std::shared_ptr<Map> map,
                       const Real r0,
                       const Real r1,
                       const Real h0,
                       const Real h1,
                       const Real theta0,
                       const Real theta1);

private:
    template<int order>
    auto
    get_aux_vals(const ValueVector<Point> &points) const;

private:
    void evaluate_0(const ValueVector<Point> &points,
                    ValueVector<Value> &values) const;

    void evaluate_1(const ValueVector<Point> &points,
                    ValueVector<Derivative<1>> &values) const;

    void evaluate_2(const ValueVector<Point> &points,
                    ValueVector<Derivative<2>> &values) const;

private:
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
};


} // of namespace functions.

IGA_NAMESPACE_CLOSE

#endif
