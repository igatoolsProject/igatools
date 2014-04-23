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

#ifndef __PUSH_FORWARD_ELEMENT_ACCESSOR_H_
#define __PUSH_FORWARD_ELEMENT_ACCESSOR_H_

#include <igatools/base/config.h>
#include <igatools/utils/value_table.h>
#include <igatools/utils/value_vector.h>
#include <igatools/base/quadrature.h>
#include <igatools/geometry/mapping_element_accessor.h>

IGA_NAMESPACE_OPEN

/**
 *
 * See module on @ref accessors_iterators for a general overview.
 * @ingroup accessors_iterators
 */
template< class PushForward_ >
class PushForwardElementAccessor
    : public MappingElementAccessor<PushForward_::dim, PushForward_::codim>
{
public:
    using ContainerType =  const PushForward_;
    using base_t = MappingElementAccessor<PushForward_::dim, PushForward_::codim>;
    using base_t::dim;
    using base_t::codim;
    using base_t::space_dim;

    static const Transformation type = PushForward_::transformation_type;

    static const Transformation transformation_type = PushForward_::transformation_type;

    using base_t::ValueMap;

    template <int range, int rank>
    using RefValue = typename PushForward_::template RefValue<range, rank>;

    template <int range, int rank>
    using PhysValue = typename PushForward_::template PhysValue<range, rank>;

    template <int range, int rank, int order>
    using RefDerivative = typename PushForward_::template RefDerivative<range, rank, order>;

    template <int range, int rank, int order>
    using PhysDerivative = typename PushForward_::template PhysDerivative<range, rank, order>;
    /**
     * Default constructor. Not allowed to be used.
     */
    PushForwardElementAccessor() = delete ;

    explicit PushForwardElementAccessor(const std::shared_ptr<ContainerType> push_forward,
                                        const int index) ;

    /**
     * Copy constructor.
     */
    PushForwardElementAccessor(const PushForwardElementAccessor<PushForward_> &element) = default;

    /**
     * Copy assignment operator.
     */
    PushForwardElementAccessor<PushForward_> &operator=(const PushForwardElementAccessor<PushForward_> &element) = default ;

    /**
     * Move constructor.
     */
    PushForwardElementAccessor(PushForwardElementAccessor<PushForward_> &&element) = default;

    /**
     * Move assignment operator.
     */
    PushForwardElementAccessor<PushForward_> &operator=(PushForwardElementAccessor<PushForward_> &&element) = default ;

    /**
     * Destructor.
     */
    ~PushForwardElementAccessor() = default ;

    /**
     * This function and fill_cache must be called before
     * to call the transformation functions.
     * It takes care of allocating the necessary space in the cache
     * according to the provided evaluation flags.
     * @note The cache needs to filled by fill_cache
     * before calling the transformations functions.
     *
     */
    void init_values(const ValueFlags fill_flag,
                     const Quadrature<dim> &quad) ;

    void init_face_values(const Index face_id,
                          const ValueFlags fill_flag,
                          const Quadrature<dim-1> &quad) ;

    /** @name Mapping used for transforming quantities   */
    ///@{
    /**
     * Transform values of scalar, vector or tensor fields from
     * the reference domain to the physical one.
     * The templates arguments indicate the value type and the container type
     * of object to be transformed.
    */
    template < int dim_range, int rank, template<class T> class Container, Transformation ttype=type >
    void
    transform_values(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        Container< PhysValue<dim_range, rank> > &D0v,
        const TopologyId<dim> &topology_id = ElemTopology<dim>(),
        typename std::enable_if<ttype == Transformation::h_grad>::type * = 0) const;

    template < int dim_range, int rank, template<class T> class Container, Transformation ttype=type >
    void
    transform_values(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        Container< PhysValue<dim_range, rank> > &D0v,
        const TopologyId<dim> &topology_id = ElemTopology<dim>(),
        typename std::enable_if<ttype == Transformation::h_div>::type * = 0) const;

    template <int dim_range, int rank, template<class T> class Container, Transformation ttype=type>
    void
    transform_gradients(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
        Container< PhysDerivative<dim_range, rank, 1> > &D1v,
        const TopologyId<dim> &topology_id = ElemTopology<dim>(),
        typename std::enable_if<ttype == Transformation::h_grad>::type * = 0) const;

    template <int dim_range, int rank, template<class T> class Container, Transformation ttype=type>
    void
    transform_gradients(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
        Container< PhysDerivative<dim_range, rank, 1> > &D1v,
        const TopologyId<dim> &topology_id = ElemTopology<dim>(),
        typename std::enable_if<ttype == Transformation::h_div>::type * = 0) const;

    /**
     * Transform second derivatives of scalar, vector or tensor functions from
     * the reference domain to the physical on covariantly.
     */
    template <int dim_range, int rank, template<class T> class Container, Transformation ttype=type>
    void
    transform_hessians(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
        const Container< RefDerivative<dim_range,rank,2> > &D2v_hat,
        Container< PhysDerivative<dim_range, rank, 2> > &D2v,
        const TopologyId<dim> &topology_id = ElemTopology<dim>(),
        typename std::enable_if<ttype == Transformation::h_grad>::type * = 0) const;

    template <int dim_range, int rank, template<class T> class Container, Transformation ttype=type>
    void
    transform_hessians(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
        const Container< RefDerivative<dim_range,rank,2> > &D2v_hat,
        Container< PhysDerivative<dim_range, rank, 2> > &D2v,
        const TopologyId<dim> &topology_id = ElemTopology<dim>(),
        typename std::enable_if<ttype == Transformation::h_div>::type * = 0) const;

    /**
     * Returns det(DF(qp))
     * todo: should this have another name
     */
    ValueVector<Real> transform_measure(const TopologyId<dim> &topology_id = ElemTopology<dim>()) const;

    /**
     * Returns det(DF(qp))
     * todo: should this have another name
     */
    ValueVector<Real> transform_face_measure(const Index face_id) const;


    ///@}



    void print_info(LogStream &out,const VerbosityLevel verbosity_level = VerbosityLevel::normal) const ;

    void print_memory_info(LogStream &out) const ;

private:

    std::shared_ptr<ContainerType> push_forward_ ;



    /**
     * MappingFillFlags activated by ValueFlags
     *
     * fill_value:    (h_grad && value)
     * fill_gradient: (h_grad && gradient)||(h_div && values)  ||(h_div && gradient)||(h_curl && values)||(h_curl && gradient)
     * fill_hessian:  (h_grad && hessian )||(h_div && gradient)||(h_curl && gradient)
     * fill_inv_gradient:                                     (h_grad && gradient)
     * fill_dets:     (w_measure)
     *
     */
    ValueFlags value_to_mapping_flag(const ValueFlags v_flag) const ;



#if 0
    /** @name Transformation used for the h_grad push-forward   */
    ///@{

    /**
     * Transform values of scalar, vector or tensor fields from
     * the reference domain to the physical one covariantly.
     * The templates arguments indicate the value type and the container type
     * of object to be transformed.
    */
    template < int dim_range, int rank,template<class T> class Container >
    void
    transform_values_h_grad(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        Container< typename Derivatives< space_dim,dim_range, rank, 0>::value_t > &D0v) const;


    /**
     * Transform first derivatives of scalar, vector or tensor field from
     * the reference domain to the physical on covariantly.
     * The templates arguments indicate the value type and the container type
     * of object to be transformed.
     */
    template < int dim_range, int rank,template<class T> class Container >
    void
    transform_gradients_h_grad(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        const Container< Derivatives<dim,dim_range,rank,1> > &D1v_hat,
        Container< Derivatives< space_dim,dim_range, rank, 1> > &D1v) const;

    ///@}


    /** @name Transformation used for the h_div (a.k.a. Piola) push-forward   */
    ///@{

    /**
     * Transform the values of vector basis functions from
     * the reference domain to the physical one using the h_div transformation
     * (a.k.a. Piola transformation):
     * \f[ (v \circ F) = \frac{1}{\det DF } DF \hat{v} \; . \f]
     * The templates arguments indicate the value type and the container type of
     * objects to be transformed.
     *
     * @note In order to call this function the following condition must be satisfied
     * @code
     * dim_range == space_dim && rank == 1
     * @endcode
     */
    template < int dim_range, int rank,template<class T> class Container >
    void
    transform_values_h_div(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        Container<Values< space_dim,dim_range, rank> > &D0v,
        typename std::enable_if< dim_range == space_dim && rank == 1 >::type * = 0) const;

    /**
     * This function is there only for compatibility with the template parameters.
     * If this function is called an exception will be raised.
     * \warning If this function is called an exception will be raised.
     */
    template < int dim_range, int rank,template<class T> class Container >
    void
    transform_values_h_div(
        const Container< RefValue<dim_range, rank> > &D0v_hat,
        Container<Values< space_dim,dim_range, rank>> &D0v,
        typename std::enable_if< !(dim_range == space_dim && rank == 1) >::type * = 0) const;


    /**
     * Transform the gradients of vector basis functions from
     * the reference domain to the physical one using the h_div transformation
     * (a.k.a. Piola transformation). For the gradients the transformation writes as
     * \f[ Dv DF = \frac{1}{\det DF } \biggl( DF [D\hat{v}] + D^2F [\hat{v}] - (DF [\hat{v}]) \otimes (D^2F : DF^{-T}) \biggr) \; . \f]
     * The templates arguments indicate the value type and the container type of
     * objects to be transformed.
     *
     * @note In order to call this function the following condition must be satisfied
     * @code
     * dim_range == space_dim && rank == 1
     * @endcode
     */
    template < int dim_range, int rank,template<class T> class Container >
    void
    transform_gradients_h_div(
        const Container< Values<dim, dim_range, rank> > &D0v_hat,
        const Container< Derivatives<dim,dim_range,rank,1> > &D1v_hat,
        Container< Derivatives<space_dim,dim_range,rank,1> > &D1v,
        typename std::enable_if< dim_range == space_dim && rank == 1 >::type * = 0) const;

    /**
     * This function is there only for compatibility with the template parameters.
     * If this function is called an exception will be raised.
     * \warning If this function is called an exception will be raised.
     */
    template < int dim_range, int rank,template<class T> class Container >
    void
    transform_gradients_h_div(
        const Container< Values<dim, dim_range, rank> > &D0v_hat,
        const Container< Derivatives<dim,dim_range,rank,1> > &D1v_hat,
        Container< Derivatives<space_dim,dim_range,rank,1> > &D1v,
        typename std::enable_if< !(dim_range == space_dim && rank == 1) >::type * = 0) const;



    ///@}
#endif
} ;



IGA_NAMESPACE_CLOSE

#endif // __PUSH_FORWARD_ELEMENT_ACCESSOR_H_
