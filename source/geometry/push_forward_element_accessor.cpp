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


#include <igatools/geometry/push_forward_element_accessor.h>

#include <igatools/base/exceptions.h>
#include <igatools/geometry/unit_element.h>

#include <algorithm>

using std::array ;
using std::vector ;
using std::shared_ptr ;


IGA_NAMESPACE_OPEN

template< class PushForward >
PushForwardElementAccessor<PushForward>::
PushForwardElementAccessor(PushForward &push_forward,
                           const int index)
    :
    MappingElementAccessor<dim, codim>(*push_forward.map_, index),
    push_forward_(&push_forward)
{}



template< class PushForward >
auto
PushForwardElementAccessor<PushForward>::
value_to_mapping_flag(
    const ValueFlags v_flag) const -> ValueFlags
{
    const ValueFlags common_flag =
        ValueFlags::point|ValueFlags::map_gradient|ValueFlags::map_hessian|
        ValueFlags::w_measure|ValueFlags::face_point|ValueFlags::map_face_gradient|
        ValueFlags::map_face_hessian|ValueFlags::face_w_measure|ValueFlags::face_normal;

    /*
     * For each MappingValueFlags there is an if that checks for all
     * ValueFlags that activate the given value flag.
     */
    ValueFlags fill_flag = common_flag & v_flag;

    if (contains(v_flag, ValueFlags::point))
        fill_flag |= ValueFlags::point;

    if (contains(v_flag, ValueFlags::w_measure))
        fill_flag |= (ValueFlags::measure |
                      ValueFlags::map_gradient);

    if (contains(v_flag, ValueFlags::face_point))
        fill_flag |= ValueFlags::face_point;

    if (contains(v_flag, ValueFlags::face_w_measure))
        fill_flag |= (ValueFlags::face_measure |
                      ValueFlags::map_face_gradient);

    if (contains(v_flag, ValueFlags::face_normal))
        fill_flag |= (ValueFlags::map_face_inv_gradient |
                      ValueFlags::map_face_gradient) ;




    if (transformation_type == Transformation::h_grad)
    {
        if (contains(v_flag,ValueFlags::tran_value))
            fill_flag |= (ValueFlags::point | ValueFlags::face_point) ;
        if (contains(v_flag,ValueFlags::tran_gradient))
            fill_flag |= (ValueFlags::map_gradient |
                          ValueFlags::map_inv_gradient|
                          ValueFlags::map_face_gradient |
                          ValueFlags::map_face_inv_gradient) ;
        if (contains(v_flag,ValueFlags::tran_hessian))
            fill_flag |= (ValueFlags::map_hessian|
                          ValueFlags::map_inv_hessian |
                          ValueFlags::map_inv_gradient|
                          ValueFlags::map_face_hessian|
                          ValueFlags::map_face_inv_hessian |
                          ValueFlags::map_face_inv_gradient) ;
    }
    else if (transformation_type == Transformation::h_div)
    {
        if (contains(v_flag,ValueFlags::tran_value))
            fill_flag |= (ValueFlags::map_gradient |
                          ValueFlags::map_face_gradient) ;
        if (contains(v_flag,ValueFlags::tran_gradient))
            fill_flag |= (ValueFlags::map_gradient |
                          ValueFlags::map_hessian |
                          ValueFlags::map_face_gradient |
                          ValueFlags::map_face_hessian) ;
        if (contains(v_flag,ValueFlags::tran_hessian))
            AssertThrow(false,ExcNotImplemented()) ;
    }
    else if (transformation_type == Transformation::h_curl)
    {
        AssertThrow(false,ExcNotImplemented()) ;
        if (contains(v_flag,ValueFlags::tran_value))
            fill_flag |= (ValueFlags::map_gradient |
                          ValueFlags::map_face_gradient);
        if (contains(v_flag,ValueFlags::tran_gradient))
            fill_flag |= (ValueFlags::map_gradient |
                          ValueFlags::map_hessian |
                          ValueFlags::map_face_gradient |
                          ValueFlags::map_face_hessian) ;
        if (contains(v_flag,ValueFlags::tran_hessian))
            AssertThrow(false,ExcNotImplemented()) ;
    }
    else if (transformation_type == Transformation::l_2)
    {
        AssertThrow(false,ExcNotImplemented()) ;
        if (contains(v_flag,ValueFlags::tran_value))
            AssertThrow(false,ExcNotImplemented()) ;
        if (contains(v_flag,ValueFlags::tran_gradient))
            AssertThrow(false,ExcNotImplemented()) ;
        if (contains(v_flag,ValueFlags::tran_hessian))
            AssertThrow(false,ExcNotImplemented()) ;
    }



    // We fill extra stuff as the computation is performed anyways
    if (contains(fill_flag , ValueFlags::measure))
        fill_flag |= (ValueFlags::map_gradient |
                      ValueFlags::map_face_gradient) ;

    if (contains(fill_flag , ValueFlags::map_inv_gradient))
        fill_flag |= (ValueFlags::map_gradient |
                      ValueFlags::measure |
                      ValueFlags::map_face_gradient |
                      ValueFlags::face_measure) ;

    if (contains(fill_flag , ValueFlags::map_inv_hessian))
        fill_flag |= (ValueFlags::map_hessian |
                      ValueFlags::map_face_hessian) ;

    return fill_flag;
}



template< class PushForward >
void
PushForwardElementAccessor<PushForward>::
init_values(const ValueFlags fill_flag,
            const Quadrature<dim> &quad)
{
    base_t::init_values(value_to_mapping_flag(fill_flag), quad);
}



template< class PushForward >
void
PushForwardElementAccessor<PushForward>::
init_face_values(const Index face_id,
                 const ValueFlags fill_flag,
                 const Quadrature<dim-1> &quad)
{
    AssertThrow(false,ExcNotImplemented()) ;
}



template< class PushForward >
template < int dim_range, int rank,template<class T> class Container, Transformation ttype >
void
PushForwardElementAccessor<PushForward>::
transform_values(
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    Container< PhysValue<dim_range, rank> > &D0v,
    typename std::enable_if<ttype == Transformation::h_grad>::type *) const
{
    Assert(D0v.size() == D0v_hat.size(),
           ExcDimensionMismatch(D0v.size(),D0v_hat.size()));

    //TODO(pauletti, Jan 17, 2014): why not copy?
    auto D0v_iterator = D0v.begin() ;
    auto D0v_hat_iterator     = D0v_hat.cbegin() ;
    auto D0v_hat_iterator_end = D0v_hat.cend() ;

    for (; D0v_hat_iterator != D0v_hat_iterator_end ;
         ++D0v_iterator, ++D0v_hat_iterator)
        *D0v_iterator = *D0v_hat_iterator ;
}


template< class PushForward >
template < int dim_range, int rank,template<class T> class Container, Transformation ttype >
void
PushForwardElementAccessor<PushForward>::
transform_values(
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    Container< PhysValue<dim_range, rank> > &D0v,
    typename std::enable_if<ttype == Transformation::h_div>::type *) const
{
    AssertThrow(false, ExcMessage("This function is implemented but is not tested!")) ;

    Assert(D0v_hat.size() >= 0 , ExcEmptyObject()) ;
    Assert(D0v.size() == D0v_hat.size(), ExcDimensionMismatch(D0v.size(), D0v_hat.size())) ;

    const int num_points = this->get_num_points() ;
    Assert(num_points >= 0, ExcLowerRange(num_points,0));


    // the next two lines are written to retrieve the number of basis function in the case Container is a ValueTable object.
    // if Container is ValueVector, n_func will be equal to 1.
    Assert((D0v_hat.size() % num_points) == 0, ExcMessage("The size of the container must be a multiple of num_points.")) ;
    const int n_func = D0v_hat.size() / num_points ;


    auto D0v_iterator     = D0v.begin() ;
    auto D0v_hat_iterator = D0v_hat.cbegin() ;

    const auto &gradients_map = this->get_gradients_map() ;
    const auto &dets_map = this->get_dets_map() ;

    for (int i = 0; i < n_func; ++i)
        for (Index j_pt = 0; j_pt < num_points; ++j_pt)
        {
            const auto &DF  = gradients_map[j_pt];
            const Real det = dets_map[j_pt];

            (*D0v_iterator) = action(DF, (*D0v_hat_iterator));
            (*D0v_iterator) /= det;

            ++D0v_hat_iterator ;
            ++D0v_iterator ;
        }
}



template< class PushForward >
template < int dim_range, int rank,template<class T> class Container, Transformation ttype >
void
PushForwardElementAccessor<PushForward>::
transform_face_values(
    const Index face_id,
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    Container< PhysValue<dim_range, rank> > &D0v,
    typename std::enable_if<ttype == Transformation::h_grad>::type *) const
{
    Assert(face_id < UnitElement<PushForward::dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<PushForward::dim>::faces_per_element));
    this->transform_values(D0v_hat, D0v) ;
}


template< class PushForward >
template < int dim_range, int rank,template<class T> class Container, Transformation ttype >
void
PushForwardElementAccessor<PushForward>::
transform_face_values(
    const Index face_id,
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    Container< PhysValue<dim_range, rank> > &D0v,
    typename std::enable_if<ttype == Transformation::h_div>::type *) const
{
    AssertThrow(false, ExcMessage("This function is implemented but is not tested!")) ;

    Assert(face_id < UnitElement<PushForward::dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<PushForward::dim>::faces_per_element));

    Assert(D0v_hat.size() >= 0 , ExcEmptyObject()) ;
    Assert(D0v.size() == D0v_hat.size(), ExcDimensionMismatch(D0v.size(), D0v_hat.size())) ;

    const int num_points = this->get_num_face_points(face_id) ;
    Assert(num_points >= 0, ExcLowerRange(num_points,0));

    // the next two lines are written to retrieve the number of basis function in the case Container is a ValueTable object.
    // if Container is ValueVector, n_func will be equal to 1.
    Assert((D0v_hat.size() % num_points) == 0, ExcMessage("The size of the container must be a multiple of num_points.")) ;
    const int n_func = D0v_hat.size() / num_points ;


    auto D0v_iterator     = D0v.begin() ;
    auto D0v_hat_iterator = D0v_hat.cbegin() ;

    const auto &gradients_map = this->get_face_gradients_map(face_id) ;
    const auto &dets_map = this->get_face_dets_map(face_id) ;

    for (int i = 0; i < n_func; ++i)
        for (Index j_pt = 0; j_pt < num_points; ++j_pt)
        {
            const auto &DF  = gradients_map[j_pt];
            const Real det = dets_map[j_pt];

            (*D0v_iterator) = action(DF, (*D0v_hat_iterator));
            (*D0v_iterator) /= det;

            ++D0v_hat_iterator ;
            ++D0v_iterator ;
        }
}



template< class PushForward >
template <int dim_range, int rank, template<class T> class Container, Transformation ttype>
void
PushForwardElementAccessor<PushForward>::
transform_gradients(
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
    Container< PhysDerivative<dim_range, rank, 1> > &D1v,
    typename std::enable_if<ttype == Transformation::h_grad>::type *) const
{
    Assert(D1v_hat.size() >= 0 , ExcEmptyObject()) ;
    Assert(D1v.size() == D1v_hat.size(), ExcDimensionMismatch(D1v.size(), D1v_hat.size())) ;

    const int num_points = this->get_num_points() ;
    Assert(num_points >= 0, ExcLowerRange(num_points,0));

    // the next two lines are written to retrieve the number of basis function in the case Container is a ValueTable object.
    // if Container is ValueVector, n_func will be equal to 1.
    Assert((D1v_hat.size() % num_points) == 0,
           ExcMessage("The size of the container must be a multiple of num_points.")) ;
    const int n_func = D1v_hat.size() / num_points ;

    auto D1v_iterator     = D1v.begin() ;
    auto D1v_hat_iterator = D1v_hat.cbegin() ;

    const auto &inv_gradients_map = this->get_inv_gradients_map() ;

    for (int i_fn = 0; i_fn < n_func; ++i_fn)
    {
        for (Index j_pt = 0; j_pt < num_points; ++j_pt)
        {
            (*D1v_iterator) = compose((*D1v_hat_iterator), inv_gradients_map[j_pt]);
            ++D1v_hat_iterator ;
            ++D1v_iterator ;
        }
    }
}



template< class PushForward >
template <int dim_range, int rank, template<class T> class Container, Transformation ttype >
void
PushForwardElementAccessor<PushForward>::
transform_gradients(
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
    Container< PhysDerivative<dim_range, rank, 1> > &D1v,
    typename std::enable_if<ttype == Transformation::h_div>::type *) const
{
    AssertThrow(false, ExcMessage("This function is implemented but is not tested!")) ;
    Assert(D0v_hat.size() >= 0 , ExcEmptyObject()) ;
    Assert(D1v_hat.size() >= 0 , ExcEmptyObject()) ;
    Assert(D0v_hat.size() == D1v_hat.size(), ExcDimensionMismatch(D0v_hat.size(), D1v_hat.size())) ;
    Assert(D1v.size() == D1v_hat.size(), ExcDimensionMismatch(D1v.size(), D1v_hat.size())) ;

    const int num_points = this->get_num_points() ;
    Assert(num_points >= 0, ExcLowerRange(num_points,0));

    // the next two lines are written to retrieve the number of basis function in the case Container is a ValueTable object.
    // if Container is ValueVector, n_func will be equal to 1.
    Assert((D1v_hat.size() % num_points) == 0,
           ExcMessage("The size of the container must be a multiple of num_points.")) ;
    const int n_func = D1v_hat.size() / num_points ;

    const auto &gradients_map = this->get_gradients_map() ;
    const auto &inv_gradients_map = this->get_inv_gradients_map() ;
    const auto &hessians_map = this->get_hessians_map() ;
    const auto &dets_map = this->get_dets_map() ;

    auto Dv_iterator     = D1v.begin() ;
    auto Dv_hat_iterator = D1v_hat.cbegin() ;
    auto  v_hat_iterator =  D0v_hat.cbegin() ;

    Point<dim> D2F_invDFt_tmp ;
    const int sizeof_D2F_invDFt_tmp = sizeof(D2F_invDFt_tmp) ;
    for (int i = 0; i < n_func; ++i)
        for (Index j_pt = 0; j_pt < num_points; ++j_pt)
        {
            const auto &DF     = gradients_map[j_pt];
            const auto &DF_inv = inv_gradients_map[j_pt];
            const auto &D2F    = hessians_map[j_pt];
            const Real det   = dets_map[j_pt];

            const auto DF_Dv_hat = compose(DF, (*Dv_hat_iterator)) ;
            const auto D2F_v_hat = action(D2F,(*v_hat_iterator)) ;

            const auto DF_v_hat = action(DF, (*v_hat_iterator)) ;   // this is a Point<space_dim>

            const Tensor<dim,1,tensor::covariant, Tdouble> D2F_invDFt = contract_1(D2F,co_tensor(transpose(DF_inv))) ;

            // we copy the memory of D2F_invDFt in D2F_invDFt_tmp in order to avoid aliasing
            memcpy(&D2F_invDFt_tmp, &D2F_invDFt, sizeof_D2F_invDFt_tmp) ;

            const Tensor<dim,1,tensor::covariant,Tensor<space_dim,1,tensor::contravariant, Tdouble> >
            tens_prod = tensor_product(DF_v_hat, D2F_invDFt_tmp);

            const auto DvDF = DF_Dv_hat + D2F_v_hat - tens_prod ;

            (*Dv_iterator) = compose(DvDF ,DF_inv);
            (*Dv_iterator) /= det;

            ++Dv_iterator ;
            ++Dv_hat_iterator ;
            ++v_hat_iterator ;
        }
}



template< class PushForward >
template <int dim_range, int rank, template<class T> class Container, Transformation ttype>
void
PushForwardElementAccessor<PushForward>::
transform_face_gradients(
    const Index face_id,
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
    Container< PhysDerivative<dim_range, rank, 1> > &D1v,
    typename std::enable_if<ttype == Transformation::h_grad>::type *) const
{
    AssertThrow(false, ExcMessage("This function is implemented but is not tested!")) ;
    Assert(face_id < UnitElement<PushForward::dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<PushForward::dim>::faces_per_element));

    Assert(D1v_hat.size() >= 0 , ExcEmptyObject()) ;
    Assert(D1v.size() == D1v_hat.size(), ExcDimensionMismatch(D1v.size(), D1v_hat.size())) ;

    const int num_points = this->get_num_face_points(face_id) ;
    Assert(num_points >= 0, ExcLowerRange(num_points,0));

    // the next two lines are written to retrieve the number of basis function in the case Container is a ValueTable object.
    // if Container is ValueVector, n_func will be equal to 1.
    Assert((D1v_hat.size() % num_points) == 0,
           ExcMessage("The size of the container must be a multiple of num_points.")) ;
    const int n_func = D1v_hat.size() / num_points ;

    auto D1v_iterator     = D1v.begin() ;
    auto D1v_hat_iterator = D1v_hat.cbegin() ;

    const auto &inv_gradients_map = this->get_inv_gradients_map(FaceTopology(face_id)) ;

    for (int i_fn = 0; i_fn < n_func; ++i_fn)
    {
        for (Index j_pt = 0; j_pt < num_points; ++j_pt)
        {
            (*D1v_iterator) = compose((*D1v_hat_iterator), inv_gradients_map[j_pt]);
            ++D1v_hat_iterator ;
            ++D1v_iterator ;
        }
    }
}



template< class PushForward >
template <int dim_range, int rank, template<class T> class Container, Transformation ttype >
void
PushForwardElementAccessor<PushForward>::
transform_face_gradients(
    const Index face_id,
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
    Container< PhysDerivative<dim_range, rank, 1> > &D1v,
    typename std::enable_if<ttype == Transformation::h_div>::type *) const
{
    AssertThrow(false, ExcMessage("This function is implemented but is not tested!")) ;
    Assert(face_id < UnitElement<PushForward::dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<PushForward::dim>::faces_per_element));

    Assert(D0v_hat.size() >= 0 , ExcEmptyObject()) ;
    Assert(D1v_hat.size() >= 0 , ExcEmptyObject()) ;
    Assert(D0v_hat.size() == D1v_hat.size(), ExcDimensionMismatch(D0v_hat.size(), D1v_hat.size())) ;
    Assert(D1v.size() == D1v_hat.size(), ExcDimensionMismatch(D1v.size(), D1v_hat.size())) ;

    const int num_points = this->get_num_face_points(face_id) ;
    Assert(num_points >= 0, ExcLowerRange(num_points,0));

    // the next two lines are written to retrieve the number of basis function in the case Container is a ValueTable object.
    // if Container is ValueVector, n_func will be equal to 1.
    Assert((D1v_hat.size() % num_points) == 0,
           ExcMessage("The size of the container must be a multiple of num_points.")) ;
    const int n_func = D1v_hat.size() / num_points ;

    const auto &gradients_map = this->get_face_gradients_map(face_id) ;
    const auto &inv_gradients_map = this->get_face_inv_gradients_map(face_id) ;
    const auto &hessians_map = this->get_face_hessians_map(face_id) ;
    const auto &dets_map = this->get_face_dets_map(face_id) ;

    auto Dv_iterator     = D1v.begin() ;
    auto Dv_hat_iterator = D1v_hat.cbegin() ;
    auto  v_hat_iterator =  D0v_hat.cbegin() ;

    Point<dim> D2F_invDFt_tmp ;
    const int sizeof_D2F_invDFt_tmp = sizeof(D2F_invDFt_tmp) ;
    for (int i = 0; i < n_func; ++i)
        for (Index j_pt = 0; j_pt < num_points; ++j_pt)
        {
            const auto &DF     = gradients_map[j_pt];
            const auto &DF_inv = inv_gradients_map[j_pt];
            const auto &D2F    = hessians_map[j_pt];
            const Real det   = dets_map[j_pt];

            const auto DF_Dv_hat = compose(DF, (*Dv_hat_iterator)) ;
            const auto D2F_v_hat = action(D2F,(*v_hat_iterator)) ;

            const auto DF_v_hat = action(DF, (*v_hat_iterator)) ;   // this is a Point<space_dim>

            const Tensor<dim,1,tensor::covariant, Tdouble> D2F_invDFt = contract_1(D2F,co_tensor(transpose(DF_inv))) ;

            // we copy the memory of D2F_invDFt in D2F_invDFt_tmp in order to avoid aliasing
            memcpy(&D2F_invDFt_tmp, &D2F_invDFt, sizeof_D2F_invDFt_tmp) ;

            const Tensor<dim,1,tensor::covariant,Tensor<space_dim,1,tensor::contravariant, Tdouble> >
            tens_prod = tensor_product(DF_v_hat, D2F_invDFt_tmp);

            const auto DvDF = DF_Dv_hat + D2F_v_hat - tens_prod ;

            (*Dv_iterator) = compose(DvDF ,DF_inv);
            (*Dv_iterator) /= det;

            ++Dv_iterator ;
            ++Dv_hat_iterator ;
            ++v_hat_iterator ;
        }
}



template< class PushForward >
template <int dim_range, int rank, template<class T> class Container, Transformation ttype>
void
PushForwardElementAccessor<PushForward>::
transform_hessians(
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
    const Container< RefDerivative<dim_range,rank,2> > &D2v_hat,
    Container< PhysDerivative<dim_range, rank, 2> > &D2v,
    typename std::enable_if<ttype == Transformation::h_grad>::type *) const
{
    const int n_func = D1v_hat.get_num_functions();
    const int num_points = this->get_num_points() ;
    Assert(num_points >= 0, ExcLowerRange(num_points,0));

    Assert(num_points == D1v_hat.get_num_points(),
           ExcDimensionMismatch(num_points,D1v_hat.get_num_points()));
    Assert(num_points == D2v_hat.get_num_points(),
           ExcDimensionMismatch(num_points,D2v_hat.get_num_points()));
    Assert(n_func == D2v_hat.get_num_functions(),
           ExcDimensionMismatch(n_func,D2v_hat.get_num_functions()));
    Assert(n_func == D2v.get_num_functions(),
           ExcDimensionMismatch(n_func,D2v.get_num_functions()));
    Assert(num_points == D2v.get_num_points(),
           ExcDimensionMismatch(num_points,D2v.get_num_points()));

    const auto &inv_gradients_map = this->get_inv_gradients_map() ;
    const auto &inv_hessians_map = this->get_inv_hessians_map() ;

    for (int i = 0; i < n_func; ++i)
        for (int j = 0; j < num_points; ++j)
        {
            const auto &DF_inv  = inv_gradients_map[j];
            const auto &D2F_inv = inv_hessians_map[j];

            //TODO: create a tensor compose to get rid of for loop here
            for (int u = 0; u<dim; u++)
            {
                D2v[i][j][u] =  compose(D2v_hat[i][j][u], DF_inv);
                D2v[i][j][u] += compose(D1v_hat[i][j], D2F_inv[u]);
            }
        }
}



template< class PushForward >
template <int dim_range, int rank, template<class T> class Container, Transformation ttype>
void
PushForwardElementAccessor<PushForward>::
transform_face_hessians(
    const Index face_id,
    const Container< RefValue<dim_range, rank> > &D0v_hat,
    const Container< RefDerivative<dim_range,rank,1> > &D1v_hat,
    const Container< RefDerivative<dim_range,rank,2> > &D2v_hat,
    Container< PhysDerivative<dim_range, rank, 2> > &D2v,
    typename std::enable_if<ttype == Transformation::h_grad>::type *) const
{
    AssertThrow(false, ExcMessage("This function is implemented but is not tested!")) ;
    Assert(face_id < UnitElement<PushForward::dim>::faces_per_element && face_id >= 0,
           ExcIndexRange(face_id,0,UnitElement<PushForward::dim>::faces_per_element));

    const int n_func = D1v_hat.get_num_functions();
    const int num_points = this->get_num_face_points(face_id) ;
    Assert(num_points >= 0, ExcLowerRange(num_points,0));

    Assert(num_points == D1v_hat.get_num_points(),
           ExcDimensionMismatch(num_points,D1v_hat.get_num_points()));
    Assert(num_points == D2v_hat.get_num_points(),
           ExcDimensionMismatch(num_points,D2v_hat.get_num_points()));
    Assert(n_func == D2v_hat.get_num_functions(),
           ExcDimensionMismatch(n_func,D2v_hat.get_num_functions()));
    Assert(n_func == D2v.get_num_functions(),
           ExcDimensionMismatch(n_func,D2v.get_num_functions()));
    Assert(num_points == D2v.get_num_points(),
           ExcDimensionMismatch(num_points,D2v.get_num_points()));

    const auto &inv_gradients_map = this->get_face_inv_gradients_map(face_id) ;
    const auto &inv_hessians_map = this->get_face_inv_hessians_map(face_id) ;

    for (int i = 0; i < n_func; ++i)
        for (int j = 0; j < num_points; ++j)
        {
            const auto &DF_inv  = inv_gradients_map[j];
            const auto &D2F_inv = inv_hessians_map[j];

            //TODO: create a tensor compose to get rid of for loop here
            for (int u = 0; u<dim; u++)
            {
                D2v[i][j][u] =  compose(D2v_hat[i][j][u], DF_inv);
                D2v[i][j][u] += compose(D1v_hat[i][j], D2F_inv[u]);
            }
        }
}



template< class PushForward >
ValueVector<Real>
PushForwardElementAccessor<PushForward>::
transform_measure() const
{
    return this->get_dets_map();
}



template< class PushForward >
ValueVector<Real>
PushForwardElementAccessor<PushForward>::
transform_face_measure(const Index face_id) const
{
    return this->get_dets_map(FaceTopology(face_id));
}



template< class PushForward >
void
PushForwardElementAccessor<PushForward>::
print_info(LogStream &out) const
{
    using std::endl ;
    out << "PushForwardElementAccessor info" << endl ;

    out.push("\t") ;
    out << "transformation type = " << int(transformation_type) << endl ;
    push_forward_->print_info(out) ;
    base_t::print_info(out) ;

    out.pop() ;
}



template< class PushForward >
void
PushForwardElementAccessor<PushForward>::
print_memory_info(LogStream &out) const
{
    using std::endl ;
    out << "PushForwardElementAccessor memory info" << endl ;
    out << "this address = " << this << endl;

    out.push("\t") ;
    out << "push_forward_ address = " << push_forward_ << endl ;
    base_t::print_memory_info(out) ;
    out.pop() ;
}

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/push_forward_element_accessor.inst>
