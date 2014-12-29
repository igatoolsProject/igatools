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
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/base/exceptions.h>

using std::make_shared;

#ifdef USE_TRILINOS
using Teuchos::rcp;
#endif

#ifdef USE_PETSC
#include <petscvec.h>
#endif

IGA_NAMESPACE_OPEN

namespace
{
DeclException3(ExcVectorAccessToNonLocalElement,
               Index, Index, Index,
               << "You tried to access element (" << arg1 << ")"
               << " of a distributed vector, but only rows "
               << arg2 << " through " << arg2
               << " are stored locally and can be accessed.");

};


#ifdef USE_TRILINOS

template <TrilinosImpl trilinos_impl>
VectorTrilinos<trilinos_impl>::
VectorTrilinos(const MapPtr map)
    :
    vector_(Tools::build_vector(map))
{}


template <TrilinosImpl trilinos_impl>
auto
VectorTrilinos<trilinos_impl>::
get_trilinos_vector() const -> Teuchos::RCP<const WrappedVector>
{
    return vector_ ;
}


template <TrilinosImpl trilinos_impl>
auto
VectorTrilinos<trilinos_impl>::
get_trilinos_vector() -> Teuchos::RCP<WrappedVector>
{
    return vector_ ;
}




Vector<LAPack::trilinos_tpetra>::
Vector(const Index num_global_dofs, CommPtr comm)
    :
    VectorTrilinos<TrilinosImpl::tpetra>(
        Tpetra::createUniformContigMap<LO,GO>(num_global_dofs,comm))
{}



Vector<LAPack::trilinos_tpetra>::
Vector(const vector<Index> &dofs_id, CommPtr comm)
    :
    VectorTrilinos<TrilinosImpl::tpetra>(
        Tpetra::createNonContigMap<LO,GO>(dofs_id,comm))
{}

Vector<LAPack::trilinos_tpetra>::
Vector(MapPtr map)
    :
    VectorTrilinos<TrilinosImpl::tpetra>(map)
{}



auto
Vector<LAPack::trilinos_tpetra>::
create(const Index size) -> std::shared_ptr<self_t>
{
    return make_shared<self_t>(self_t(size));
}



auto
Vector<LAPack::trilinos_tpetra>::
create(const vector<Index> &dof_ids) -> std::shared_ptr<self_t>
{
    return make_shared<self_t>(self_t(dof_ids));
}



void
Vector<LAPack::trilinos_tpetra>::
add_entry(const Index i, const Real value)
{
    Assert(!std::isnan(value),ExcNotANumber());
    Assert(!std::isinf(value),ExcNumberNotFinite());

    vector_->sumIntoGlobalValue(i,0,value);
}



auto
Vector<LAPack::trilinos_tpetra>::
operator+=(const self_t &vec) -> self_t &
{
    vector_->update(1., *(vec.vector_), 1.);
    return *this;
}


auto
Vector<LAPack::trilinos_tpetra>::
norm2() const -> Real
{
    return vector_->getVector(0)->norm2();
}

Real
Vector<LAPack::trilinos_tpetra>::
dot(const self_t &A) const
{
    return vector_->getVector(0)->dot(*A.vector_->getVector(0));
}



void
Vector<LAPack::trilinos_tpetra>::
clear()
{
    vector_->putScalar(0.);
}

auto
Vector<LAPack::trilinos_tpetra>::
update(const Real scalar_A, const self_t &A, const Real scalar_this) -> self_t &
{
    vector_->update(scalar_A,*A.get_trilinos_vector(),scalar_this);
    return *this;
}


const Real &
Vector<LAPack::trilinos_tpetra>::
operator()(const Index global_id) const
{
    /*
    Assert(global_id < Index(vector_->getGlobalLength()),
           ExcIndexRange(global_id,0,Index(vector_->getGlobalLength()))) ;
    //*/

    const auto map = vector_->getMap();
    const auto local_id = map->getLocalElement(global_id) ;

    Assert(local_id != Teuchos::OrdinalTraits<Index>::invalid(),
           ExcVectorAccessToNonLocalElement(
               global_id,
               map->getMinGlobalIndex(),
               map->getMaxGlobalIndex()));

    return (vector_->get2dView()[0][local_id]) ;
}



Real &
Vector<LAPack::trilinos_tpetra>::
operator()(const Index global_id)
{
    /*
    Assert(global_id < Index(vector_->getGlobalLength()),
           ExcIndexRange(global_id,0,Index(vector_->getGlobalLength()))) ;
    //*/
    const auto map = vector_->getMap();
    const auto local_id = map->getLocalElement(global_id) ;

    Assert(local_id != Teuchos::OrdinalTraits<Index>::invalid(),
           ExcVectorAccessToNonLocalElement(
               global_id,
               map->getMinGlobalIndex(),
               map->getMaxGlobalIndex()));

    return (vector_->get2dViewNonConst()[0][local_id]) ;
}


Index
Vector<LAPack::trilinos_tpetra>::
size() const
{
    return vector_->getGlobalLength() ;
}


void
Vector<LAPack::trilinos_tpetra>::
add_block(
    const vector< Index > &local_to_global,
    const DenseVector &local_vector)
{
    Assert(!local_to_global.empty(), ExcEmptyObject()) ;
    const Index num_dofs = local_to_global.size() ;

    Assert(Index(local_vector.size()) == num_dofs,
           ExcDimensionMismatch(local_vector.size(), num_dofs)) ;

    for (Index i = 0 ; i < num_dofs ; ++i)
    {
        Assert(!std::isnan(local_vector(i)),ExcNotANumber());
        Assert(!std::isinf(local_vector(i)),ExcNumberNotFinite());
        vector_->sumIntoGlobalValue(local_to_global[i],0,local_vector(i)) ;
    }
}


vector<Real>
Vector<LAPack::trilinos_tpetra>::
get_local_coefs(const vector<Index> &local_to_global_ids) const
{
    vector<Real> local_coefs;
    for (const auto &global_id : local_to_global_ids)
        local_coefs.emplace_back((*this)(global_id));

    return local_coefs;
}


vector<Real>
Vector<LAPack::trilinos_tpetra>::
get_as_vector() const
{
    Teuchos::ArrayRCP<const Real> vec = vector_->getData(0) ;

    const Index n_entries = vec.size();
    vector<Real> coefs(n_entries);

    const auto map = vector_->getMap();
    for (Index i = 0 ; i < n_entries ; ++i)
    {
        const auto global_id = map->getGlobalElement(i) ;
    	coefs[i] = (*this)(global_id);
    }
    return coefs;
}

void
Vector<LAPack::trilinos_tpetra>::
print_info(LogStream &out) const
{
    using std::endl;

    out << "-----------------------------" << endl;
    // Commented as different trilinos version show different information here
    // out << *vector_ << endl;

    Teuchos::ArrayRCP<const Real> vec = vector_->getData(0) ;

    const Index n_entries = vec.size();

    const auto map = vector_->getMap();

    out << "Global_ID        Value" << endl;
    for (Index i = 0 ; i < n_entries ; ++i)
    {
        const auto global_id = map->getGlobalElement(i) ;

        out << global_id << "        " << (*this)(global_id) <<endl ;
    }
    out << "-----------------------------" << endl;
}








Vector<LAPack::trilinos_epetra>::
Vector(const Index num_global_dofs, CommPtr comm)
    :
    VectorTrilinos<TrilinosImpl::epetra>(
        Teuchos::rcp(new Map(num_global_dofs,0,*comm)))
{}


Vector<LAPack::trilinos_epetra>::
Vector(const vector<Index> &dofs_id, CommPtr comm)
    :
    VectorTrilinos<TrilinosImpl::epetra>(
        Teuchos::rcp(new Map(-1,dofs_id.size(),dofs_id.data(),0,*comm)))
{}

auto
Vector<LAPack::trilinos_epetra>::
create(const Index num_global_dofs) -> std::shared_ptr<self_t>
{
    return make_shared<self_t>(self_t(num_global_dofs));
}

auto
Vector<LAPack::trilinos_epetra>::
create(const vector<Index> &dof_ids) -> std::shared_ptr<self_t>
{
    return make_shared<self_t>(self_t(dof_ids));
}

void
Vector<LAPack::trilinos_epetra>::
clear()
{
    vector_->PutScalar(0.);
}


Real
Vector<LAPack::trilinos_epetra>::
dot(const self_t &A) const
{
    const auto &this_vec = *(*vector_)(0);
    const auto &A_vec = *((*(A.vector_))(0));
    Real dot;
    this_vec.Dot(A_vec,&dot);
    return dot;
}


auto
Vector<LAPack::trilinos_epetra>::
update(const Real scalar_A, const self_t &A, const Real scalar_this) -> self_t &
{
    vector_->Update(scalar_A,*A.get_trilinos_vector(),scalar_this);
    return *this;
}

auto
Vector<LAPack::trilinos_epetra>::
norm2() const -> Real
{
    /*
      shared_ptr<Real> res;
      solution_->get_trilinos_vector ()->Norm2 (res.get());
      return *res;
    //*/
    Real res;
    (*vector_)(0)->Norm2(&res);

    return res;
}


void
Vector<LAPack::trilinos_epetra>::
add_entry(const Index i, const Real value)
{
    Assert(!std::isnan(value),ExcNotANumber());
    Assert(!std::isinf(value),ExcNumberNotFinite());

    vector_->SumIntoGlobalValue(i,0,value);
};



const Real &
Vector<LAPack::trilinos_epetra>::
operator()(const Index global_id) const
{
    const auto &map = vector_->Map();
    const auto local_id = map.LID(global_id) ;

    Assert(local_id != Teuchos::OrdinalTraits<Index>::invalid(),
           ExcVectorAccessToNonLocalElement(
               global_id,
               map.MinAllGID(),
               map.MaxAllGID()));

    return ((*vector_)[0][local_id]) ;
}



Real &
Vector<LAPack::trilinos_epetra>::
operator()(const Index global_id)
{
    const auto &map = vector_->Map();
    const auto local_id = map.LID(global_id) ;

    Assert(local_id != Teuchos::OrdinalTraits<Index>::invalid(),
           ExcVectorAccessToNonLocalElement(
               global_id,
               map.MinAllGID(),
               map.MaxAllGID()));

    return ((*vector_)[0][local_id]) ;
}




vector<Real>
Vector<LAPack::trilinos_epetra>::
get_local_coefs(const vector<Index> &local_to_global_ids) const
{
    vector<Real> local_coefs;
    for (const auto &global_id : local_to_global_ids)
        local_coefs.emplace_back((*this)(global_id));

    return local_coefs;
}

vector<Real>
Vector<LAPack::trilinos_epetra>::
get_as_vector() const
{
    const Index n_entries = vector_->GlobalLength();
    vector<Real> coefs(n_entries);

    const auto &map = vector_->Map();
    const auto &data = (*vector_)[0];
    for (Index i = 0 ; i < n_entries ; ++i)
    	coefs[i] = data[i];

    return coefs;
}

void
Vector<LAPack::trilinos_epetra>::
print_info(LogStream &out) const
{
    using std::endl;

    out << "-----------------------------" << endl;

    const Index n_entries = vector_->GlobalLength();
    const auto &map = vector_->Map();
    out << "Global_ID        Value" << endl;

    const auto &data = (*vector_)[0];
    for (Index i = 0 ; i < n_entries ; ++i)
    {
        const auto global_id = map.GID(i) ;

        out << global_id << "        " << data[i] <<endl ;
    }
    out << "-----------------------------" << endl;
}



void
Vector<LAPack::trilinos_epetra>::
add_block(
    const vector< Index > &local_to_global,
    const DenseVector &local_vector)
{
    Assert(!local_to_global.empty(), ExcEmptyObject()) ;
    const Index num_dofs = local_to_global.size() ;

    Assert(Index(local_vector.size()) == num_dofs,
           ExcDimensionMismatch(local_vector.size(), num_dofs)) ;

    for (Index i = 0 ; i < num_dofs ; ++i)
    {
        Assert(!std::isnan(local_vector(i)),ExcNotANumber());
        Assert(!std::isinf(local_vector(i)),ExcNumberNotFinite());
        vector_->SumIntoGlobalValue(local_to_global[i],0,local_vector(i)) ;
    }
}


template class VectorTrilinos<TrilinosImpl::tpetra>;
template class VectorTrilinos<TrilinosImpl::epetra>;

#endif // #ifdef USE_TRILINOS



#ifdef USE_PETSC

Vector<LAPack::petsc>::
Vector(const Index num_global_dofs)
{
    PetscErrorCode ierr;
    comm_ = PETSC_COMM_WORLD;
    ierr = VecCreate(comm_, &vector_);  // CHKERRQ(ierr);
    ierr = VecSetSizes(vector_, PETSC_DECIDE, num_global_dofs); // CHKERRQ(ierr);
    ierr = VecSetFromOptions(vector_); // CHKERRQ(ierr);
    ierr = VecZeroEntries(vector_); // CHKERRQ(ierr);
}



Vector<LAPack::petsc>::
Vector(const vector<Index> &dofs_id)
    :
    Vector<LinearAlgebraPackage::petsc>(dofs_id.size())
{
    //TODO: chencge this constructor in order to use the dofs_id entries and not only its size.
//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
}


auto
Vector<LAPack::petsc>::
create(const Index size) -> std::shared_ptr<self_t>
{
    return make_shared<self_t>(self_t(size));
}

auto
Vector<LAPack::petsc>::
create(const vector<Index> &dof_ids) -> std::shared_ptr<self_t>
{
    return make_shared<self_t>(self_t(dof_ids));
}



void
Vector<LAPack::petsc>::
add_entry(const Index i, const Real value)
{
    Assert(!std::isnan(value),ExcNotANumber());
    Assert(!std::isinf(value),ExcNumberNotFinite());

    PetscErrorCode ierr;
    ierr = VecSetValue(vector_, i, value, ADD_VALUES); // CHKERRQ(ierr);
};


const Real &
Vector<LAPack::petsc>::
operator()(const Index global_id) const
{
    VecGetValues(vector_,1,&global_id,const_cast<Real *>(&real_tmp_));

    return real_tmp_;
//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
}



Real &
Vector<LAPack::petsc>::
operator()(const Index global_id)
{
    VecGetValues(vector_,1,&global_id,&real_tmp_);

    return real_tmp_;
//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
}


Index
Vector<LAPack::petsc>::
size() const
{
    PetscErrorCode ierr;
    PetscInt vector_size;
    ierr = VecGetSize(vector_, &vector_size);
    CHKERRQ(ierr);
    return vector_size;
}

auto
Vector<LAPack::petsc>::
get_petsc_vector() const -> Vec
{
//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
    return vector_ ;
}

auto
Vector<LAPack::petsc>::
get_petsc_vector() -> Vec
{
//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
    return vector_ ;
}
//*/

void
Vector<LAPack::petsc>::
add_block(
    const vector< Index > &local_to_global,
    const DenseVector &local_vector)
{
    PetscErrorCode ierr;

    Assert(!local_to_global.empty(), ExcEmptyObject()) ;
    const Index num_dofs = local_to_global.size() ;

    Assert(Index(local_vector.size()) == num_dofs,
           ExcDimensionMismatch(local_vector.size(), num_dofs)) ;

    vector<PetscScalar> values;

    for (Index i = 0 ; i < num_dofs ; ++i)
    {
        Assert(!std::isnan(local_vector(i)),ExcNotANumber());
        Assert(!std::isinf(local_vector(i)),ExcNumberNotFinite());
        values.push_back(local_vector(i));
    }
    ierr = VecSetValues(vector_, num_dofs, local_to_global.data(), values.data(), ADD_VALUES); //CHKERRQ(ierr);

}


vector<Real>
Vector<LAPack::petsc>::
get_local_coefs(const vector<Index> &local_to_global_ids) const
{
    vector<Real> local_coefs;

    int num_local_dofs = local_to_global_ids.size();
    PetscScalar values[num_local_dofs];

    VecGetValues(vector_, num_local_dofs, local_to_global_ids.data(), values);

    local_coefs.assign(values, values+num_local_dofs);

    return local_coefs;

}


void
Vector<LAPack::petsc>::
print_info(LogStream &out) const
{
    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
    /*
        using std::endl;

        out << "-----------------------------" << endl;
        // Commented as different petsc version show different information here
        // out << *vector_ << endl;

        Teuchos::ArrayRCP<const Real> vec = vector_->getData(0) ;

        const Index n_entries = vec.size();

        out << "Global_ID        Value" << std::endl;
        for (Index i = 0 ; i < n_entries ; ++i)
            out << i << "        " << vec[i] << endl ;
        out << "-----------------------------" << endl;
        //*/
}

#endif // #ifdef USE_PETSC

IGA_NAMESPACE_CLOSE

