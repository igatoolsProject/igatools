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

#ifndef PUSH_FORWARD_H_
#define PUSH_FORWARD_H_

#include <igatools/geometry/mapping.h>

IGA_NAMESPACE_OPEN

constexpr
int physical_range(const int ref_range, const int space_dim, const Transformation type)
{
    return type == Transformation::h_grad ? ref_range : space_dim;
}

template < class PushForward > class PushForwardElementAccessor;

/**
 * * The main use of mapping is the transformation of objects
 * in the reference domain to the physical domain.
 * Such objects include:
 * - values and derivatives of basis functions
 * - measures (length, area, volume, etc)
 * - external normals
 * - boundary normals
 *
 * Given the objects we want to transform, there are different ways
 * to transform the same object, we call this the transformation type.
 * Such transformation types are:
 * - H_GRAD also known as covariant
 * - H_DIV  also known as Piola
 * - H_CURL
 * @tparam transformation_type_ The type of transformation. @sa Transformation
 * @ingroup containers
 */
template<Transformation transformation_type_, int dim_, int codim_ = 0 >
class PushForward
{
public:
    using Map = Mapping<dim_, codim_>;

    using GridType = typename Map::GridType;
    /**
     * see Mapping<dim_, codim_>::dim
     */
    static const int dim       = Map::dim;
    static const int codim     = Map::codim;
    static const int space_dim = Map::space_dim;


    static const Transformation type = transformation_type_;

    template<int ref_range>
    struct PhysRange
    {
        static const int value = physical_range(ref_range, space_dim, type);

    };

    template <int range, int rank>
    using RefValue = Values<range, rank>;

    template <int range, int rank, int order>
    using RefDerivative = Derivatives<dim, range, rank, order>;

    template <int range, int rank, int order>
    using RefFaceDerivative = Conditional< dim == 0,
          Derivatives<0,     range, rank, order>,
          Derivatives<dim-1, range, rank, order> >;

    template <int range, int rank>
    using PhysValue = Values<PhysRange<range>::value, rank>;

    template <int range, int rank, int order>
    using PhysDerivative = Derivatives<space_dim, PhysRange<range>::value, rank, order>;


    static const Transformation transformation_type = transformation_type_ ;

    typedef PushForward<transformation_type, dim, codim > Self ;

    /**
     * Typedef for the PushForward on a face.
     */
    typedef PushForward<transformation_type,dim-1, codim + 1> FacePushForward;


    /**
     * Type for the element accessor.
     */
    using ElementAccessor = PushForwardElementAccessor<Self>;


    /**
     * Default constructor. Not allowed to be used.
     */
    PushForward() = delete ;


    /**
     * Constructor. As input argument, it takes the Mapping @p map used to define the push-forward.
     */
    PushForward(const std::shared_ptr< Map > map) ;


    static std::shared_ptr<Self>
    create(const std::shared_ptr< Map > map) ;

    /**
     * Copy constructor.
     */
    PushForward(const Self &push_forward) ;


    void reset_map(const std::shared_ptr< Map > map) ;

    std::shared_ptr<Self> clone() const ;


    /**
     * Copy assignment operator. Not allowed to be used.
     */
    Self &operator=(const Self &push_forward) = delete ;


    /**
     * Retrieve the Mapping used to define the push-forward.
     */
    std::shared_ptr< const Map > get_mapping() const ;



    void print_info(LogStream &out) const ;

    void print_memory_info(LogStream &out) const ;

private:

    /**
     * Mapping used to define the push-forward.
     */
    std::shared_ptr< Map > map_ ;


    friend class PushForwardElementAccessor<Self> ;
} ;



IGA_NAMESPACE_CLOSE

#endif // PUSH_FORWARD_H_
