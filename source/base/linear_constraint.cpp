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


#include <igatools/base/linear_constraint.h>
#include <igatools/base/exceptions.h>


IGA_NAMESPACE_OPEN

LinearConstraint::
LinearConstraint(const vector<Index> &dofs,const vector<Real> &coeffs,const Real rhs)
{
    Assert(dofs.size() == coeffs.size(),ExcDimensionMismatch(dofs.size(),coeffs.size()));
    Assert(!dofs.empty(),ExcEmptyObject());

    const Index n_dofs = dofs.size();
    for (Index i = 0 ; i < n_dofs ; ++i)
    {
        Assert(dofs[i] >= 0,ExcLowerRange(dofs[i],0));
        this->first.emplace_back(std::make_pair(dofs[i],coeffs[i]));
    }
    this->second = rhs;
}

Real
LinearConstraint::
get_rhs() const
{
    return this->second;
}

void
LinearConstraint::
set_rhs(const Real rhs)
{
    this->second = rhs;
}

Index
LinearConstraint::
get_num_lhs_terms() const
{
    Assert(!this->first.empty(),ExcEmptyObject());
    return this->first.size();
}


const std::pair<Index,Real> &
LinearConstraint::
get_lhs_term(const int i) const
{
    Assert(i >= 0 && i < this->get_num_lhs_terms(),ExcIndexRange(i,0,this->get_num_lhs_terms()));
    return this->first[i];
}


Index
LinearConstraint::
get_dof_index(const int i) const
{
    return this->get_lhs_term(i).first;
}

Real
LinearConstraint::
get_dof_coeff(const int i) const
{
    return this->get_lhs_term(i).second;
}


IGA_NAMESPACE_CLOSE




