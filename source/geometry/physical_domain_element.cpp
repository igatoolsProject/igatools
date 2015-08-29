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

#include <igatools/geometry/physical_domain_element.h>
#include <igatools/functions/function_element.h>
//#include <igatools/linear_algebra/dense_matrix.h>
//#include <boost/numeric/ublas/operation.hpp>

IGA_NAMESPACE_OPEN


template<int dim_, int codim_, class ContainerType_>
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
PhysicalDomainElementBase(std::shared_ptr<ContainerType_> phys_dom,
                          const ListIt &index,
                          const PropId &prop)
    :
    phys_dom_(phys_dom),
    grid_elem_(phys_dom_->get_grid()->create_element(index,prop)),
    func_elem_(phys_dom_->get_function()->create_element(index,prop))
{}



template<int dim_, int codim_, class ContainerType_>
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
PhysicalDomainElementBase(const self_t &elem,
                          const CopyPolicy &copy_policy)
    :
    phys_dom_(elem.phys_dom_)
{
    if (copy_policy == CopyPolicy::shallow)
    {
        grid_elem_ = elem.grid_elem_;
        /// func_elem_ = elem.func_elem_;
        local_cache_ = elem.local_cache_;
    }
    else
    {
        local_cache_ =
            std::shared_ptr<CacheType>(new CacheType(*elem.local_cache_));
        grid_elem_ = std::make_shared<GridElem>(*elem.grid_elem_,CopyPolicy::deep);
        //   func_elem_ = std::make_shared<FuncElem>(*elem.func_elem_,CopyPolicy::deep);
    }
}


template<int dim_, int codim_, class ContainerType_>
void
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
deep_copy_from(const self_t &elem)
{
    Assert(false, ExcNotImplemented());
}


template<int dim_, int codim_, class ContainerType_>
void
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
shallow_copy_from(const self_t &elem)
{
    Assert(false, ExcNotImplemented());
}



template<int dim_, int codim_, class ContainerType_>
bool
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
operator ==(const self_t &elem) const
{
    Assert(phys_dom_ == elem.phys_dom_,
           ExcMessage("Cannot compare elements on different grid."));
    return (grid_elem_ == elem.grid_elem_);
}


template<int dim_, int codim_, class ContainerType_>
bool
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
operator !=(const self_t &elem) const
{
    Assert(phys_dom_ == elem.phys_dom_,
           ExcMessage("Cannot compare elements on different grid."));
    return (grid_elem_ != elem.grid_elem_);
}

template<int dim_, int codim_, class ContainerType_>
bool
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
operator <(const self_t &elem) const
{
    Assert(phys_dom_ == elem.phys_dom_,
           ExcMessage("Cannot compare elements on different grid."));
    return (grid_elem_ < elem.grid_elem_);
}

template<int dim_, int codim_, class ContainerType_>
bool
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
operator >(const self_t &elem) const
{
    Assert(phys_dom_ == elem.phys_dom_,
           ExcMessage("Cannot compare elements on different grid."));
    return (grid_elem_ > elem.grid_elem_);
}



#if 0
template<int dim_, int codim_, class ContainerType_>
auto
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
compute_inv_first_fundamental_form() const -> ValueVector<MetricTensor>
{
    ValueVector<MetricTensor> res;
    const auto &DF = this->template get_values<_Gradient, dim>(0);
    const auto n_points = DF.get_num_points();

    res.resize(n_points);
    Real det;
    for (int i = 0; i< n_points; ++i)
    {
        const auto &A = DF[i];
        const auto A_t   = co_tensor(transpose(A));
        const auto G     = compose(A_t, A);
        res[i] = inverse(G, det);
    }

    return res;
}



template<int dim_, int codim_, class ContainerType_>
auto
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
compute_second_fundamental_form() const -> ValueVector<MetricTensor>
{
    Assert(codim==1, ExcNotImplemented());

    const auto &D2_F  = this->template get_values<_Hessian, dim>(0);
    const auto normal = this->get_external_normals();

    const auto n_points = D2_F.get_num_points();

    // const auto G_inv = compute_inv_first_fundamental_form();

    ValueVector<MetricTensor> res;
    res.resize(n_points);

    MetricTensor A;
    for (int pt = 0; pt < n_points; ++pt)
    {
        const auto &D2_F_pt = D2_F[pt];
        for (int u=0; u<dim; ++u)
        {
            const auto B = co_tensor(transpose(D2_F_pt[u]));
            A[u] = action(B, normal[pt]);
        }
        res[pt] = -A;

        // res[pt] = -compose(A, G_inv[pt]);

    }

    return res;
}


template<int dim_, int codim_, class ContainerType_>
auto
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
get_principal_curvatures() const -> const ValueVector<SafeSTLVector<Real>> &
{
#if 0
    Assert(codim==1, ExcNotImplemented());

    const auto H = compute_second_fundamental_form();
    const auto G_inv = compute_inv_first_fundamental_form();

    const auto n_points = H.get_num_points();

    ValueVector<SafeSTLVector<Real>> res(n_points);

    for (int pt = 0; pt < n_points; ++pt)
    {
        const MetricTensor B = compose(H[pt], G_inv[pt]);
        const auto A = unroll_to_matrix(B);
        res[pt] = A.eigen_values();
    }
    return res;
#endif
    Assert(codim==1, ExcNotImplemented());
    return get_values_from_cache<_Curvature,dim_>(0);
}



template<int dim_, int codim_, class ContainerType_>
auto
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
get_external_normals() const -> const ValueVector<Points<space_dim> > &
{
#if 0
    Assert(codim==1, ExcNotImplemented());
    ValueVector<Points<space_dim> > res;
    const auto &DF = this->template get_values<_Gradient, dim>(0);
    const auto n_points = DF.get_num_points();

    res.resize(n_points);
    for (int pt = 0; pt < n_points; ++pt)
    {
        res[pt] = cross_product<dim, codim>(DF[pt]);
        res[pt] /= res[pt].norm();
    }

    return res;
#endif

    Assert(codim==1, ExcNotImplemented());
    return get_values_from_cache<_OuterNormal,dim_>(0);
}



template<int dim_, int codim_, class ContainerType_>
auto
PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
get_D_external_normals() const -> ValueVector< Derivative<1> >
{
    Assert(codim==1, ExcNotImplemented());

    const auto H = compute_second_fundamental_form();
    const auto &DF = this->template get_values<_Gradient, dim>(0);
    const auto G_inv = compute_inv_first_fundamental_form();

    const auto n_points = H.get_num_points();
    ValueVector< Derivative<1> > Dn(n_points);

    for (int pt = 0; pt< n_points; ++pt)
    {
        auto L = compose(DF[pt], G_inv[pt]);
        Dn[pt] = compose(L, H[pt]);
    }
    return Dn;
}
#endif


//template<int dim_, int codim_, class ContainerType_>
//auto
//PhysicalDomainElementBase<dim_, codim_, ContainerType_>::
//clone() const -> std::shared_ptr<self_t>
//{
//    auto elem = std::make_shared<self_t>(*this,CopyPolicy::deep);
//    Assert(elem != nullptr, ExcNullPtr());
//    return elem;
//}



//  const auto &D2_F  = this->template get_values<2, dim>(0);
//  const auto normal = this->get_external_normals();
//
//  const auto n_points = D2_F.get_num_points();
//  ValueVector<SafeSTLVector<Real>> res(n_points);
//  const auto G_inv = compute_inv_first_fundamental_form();
//  DenseMatrix A(dim, dim);
//  for (int pt = 0; pt < n_points; ++pt)
//  {
//      const auto B = unroll_to_matrix(G_inv[pt]);
//      for (int i = 0; i<dim; ++i)
//          for (int j = 0; j<dim; ++j)
//              A(i,j) = -scalar_product(D2_F[pt][i][j], normal[pt]);
//      DenseMatrix C(dim, dim);
//      boost::numeric::ublas::axpy_prod(B, A, C, true);
//
//      res[pt] = C.eigen_values();
//  }



//    template<int sub_dim>
//    ValueVector<Points<space_dim> >
//    get_boundary_normals(const int s_id) const
//    {
//        Assert(dim==sub_dim+1, ExcNotImplemented());
//        ValueVector<Points<space_dim>> res;
//        const auto &DF_inv = get_inverse_values<1, sub_dim>(s_id);
//        const auto n_hat  = get_grid()->template get_boundary_normals<sub_dim>(s_id)[0];
//
//        const auto n_points = DF_inv.get_num_points();
//        res.resize(n_points);
//        for (int i = 0; i< n_points; ++i)
//        {
//            const auto DF_inv_t = co_tensor(transpose(DF_inv[i]));
//            res[i] = action(DF_inv_t, n_hat);
//            res[i] /= res[i].norm();
//        }
//        return res;
//    }

IGA_NAMESPACE_CLOSE

#include <igatools/geometry/physical_domain_element.inst>

