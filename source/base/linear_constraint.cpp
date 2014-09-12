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
LinearConstraint(
    const LinearConstraintType &type,
    const vector<Index> &dofs,
    const vector<Real> &coeffs,
    const Real rhs)
{
    Assert(dofs.size() == coeffs.size(),ExcDimensionMismatch(dofs.size(),coeffs.size()));
    Assert(!dofs.empty(),ExcEmptyObject());

    const Index n_dofs = dofs.size();
    for (Index i = 0 ; i < n_dofs ; ++i)
    {
        Assert(dofs[i] >= 0,ExcLowerRange(dofs[i],0));
        lhs_.emplace_back(std::make_pair(dofs[i],coeffs[i]));
    }
    rhs_ = rhs;

    type_ = type;
}

std::shared_ptr<LinearConstraint>
LinearConstraint::
create(const LinearConstraintType &type,
       const vector<Index> &dofs,const vector<Real> &coeffs,const Real rhs)
{
    return std::shared_ptr<LinearConstraint>(new LinearConstraint(type,dofs,coeffs,rhs));
}

Real
LinearConstraint::
get_rhs() const
{
    return rhs_;
}

void
LinearConstraint::
set_rhs(const Real rhs)
{
    rhs_ = rhs;
}

Index
LinearConstraint::
get_num_lhs_terms() const
{
    Assert(!lhs_.empty(),ExcEmptyObject());
    return lhs_.size();
}


const std::pair<Index,Real> &
LinearConstraint::
get_lhs_term(const int i) const
{
    return lhs_[i];
}


Index
LinearConstraint::
get_dof_index(const int i) const
{
    return lhs_[i].first;
}

Real
LinearConstraint::
get_coeff(const int i) const
{
    return lhs_[i].second;
}

Index
LinearConstraint::
get_num_dofs() const
{
    return this->get_num_lhs_terms();
}

Index
LinearConstraint::
get_num_coeffs() const
{
    return this->get_num_lhs_terms();
}

vector<Index>
LinearConstraint::
get_dofs_id() const
{
    vector<Index> dofs_id;
    for (const auto &lhs_term : lhs_)
        dofs_id.push_back(lhs_term.first);

    return dofs_id;
}

vector<Real>
LinearConstraint::
get_coefficients() const
{
    vector<Real> coeffs;
    for (const auto &lhs_term : lhs_)
        coeffs.push_back(lhs_term.second);

    return coeffs;
}

LinearConstraintType
LinearConstraint::
get_type() const
{
    return type_;
}

bool
LinearConstraint::
is_dof_present(const Index dof) const
{
    bool is_dof_present = false;

    for (const auto &lhs_term : lhs_)
        if (lhs_term.first == dof)
        {
            is_dof_present = true;
            break;
        }

    return is_dof_present;
}


Real
LinearConstraint::
eval_absolute_error(const vector<Real> &dof_coeffs) const
{
    const auto n_terms = this->get_num_coeffs();

    Real error = -rhs_;
    for (Index i = 0 ; i < n_terms ; ++i)
        error += this->get_coeff(i) * dof_coeffs[this->get_dof_index(i)];

    return fabs(error);
}

void
LinearConstraint::
print_info(LogStream &out) const
{
    const int n_dofs = this->get_num_dofs();

    for (int i = 0 ; i < n_dofs ; ++i)
    {
        out << "Dof[" << i << "] = " << this->get_dof_index(i)
            << "      Coef = " << this->get_coeff(i) << std::endl;
    }
    out << "Rhs value = " << this->get_rhs() << std::endl;
}


IGA_NAMESPACE_CLOSE




