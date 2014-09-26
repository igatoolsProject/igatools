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

Vector<LAPack::trilinos>::
Vector(const Index num_global_dofs, CommPtr comm)
    :
    vector_(Tpetra::createMultiVector<Real,LO,GO>(
                Tpetra::createUniformContigMap<LO,GO>(
                    num_global_dofs,
                    comm),
                1))
{}



Vector<LAPack::trilinos>::
Vector(const vector<Index> &dofs_id, CommPtr comm)
    :
    vector_(Tpetra::createMultiVector<Real,LO,GO>(
                Tpetra::createNonContigMap<LO,GO>(
                    dofs_id,
                    comm),
                1))
{}

Vector<LAPack::trilinos>::
Vector(DofsMapPtr map)
    :
    vector_(Tpetra::createMultiVector<Real,LO,GO>(map,1))
{}


auto
Vector<LAPack::trilinos>::
create(const Index size) -> std::shared_ptr<self_t>
{
    return make_shared<self_t>(self_t(size));
}

auto
Vector<LAPack::trilinos>::
create(const vector<Index> &dof_ids) -> std::shared_ptr<self_t>
{
    return make_shared<self_t>(self_t(dof_ids));
}



void
Vector<LAPack::trilinos>::
add_entry(const Index i, const Real value)
{
    Assert(!std::isnan(value),ExcNotANumber());
    Assert(!std::isinf(value),ExcNumberNotFinite());

    vector_->sumIntoGlobalValue(i,0,value);
};


const Real &
Vector<LAPack::trilinos>::
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
Vector<LAPack::trilinos>::
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
Vector<LAPack::trilinos>::
size() const
{
    return vector_->getGlobalLength() ;
}

auto
Vector<LAPack::trilinos>::
get_trilinos_vector() const -> Teuchos::RCP<const WrappedVectorType>
{
    return vector_ ;
}

auto
Vector<LAPack::trilinos>::
get_trilinos_vector() -> Teuchos::RCP<WrappedVectorType>
{
    return vector_ ;
}

void
Vector<LAPack::trilinos>::
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
Vector<LAPack::trilinos>::
get_local_coefs(const vector<Index> &local_to_global_ids) const
{
    vector<Real> local_coefs;
    for (const auto &global_id : local_to_global_ids)
        local_coefs.emplace_back((*this)(global_id));

    return local_coefs;
}


void
Vector<LAPack::trilinos>::
print(LogStream &out) const
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


#if 0
    out << "-----------------------------" << endl;
    out << "Local_ID        Value" << endl;
    auto vector_data = vector_->get1dView();
    const int local_length = vector_->getLocalLength();
    for (int local_id = 0; local_id < local_length; ++local_id)
    {
        out << local_id << "        " << vector_data[local_id] <<endl ;
    }
    out << "-----------------------------" << endl;
#endif
}

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
print(LogStream &out) const
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

