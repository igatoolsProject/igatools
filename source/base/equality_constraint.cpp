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


#include <igatools/base/equality_constraint.h>

IGA_NAMESPACE_OPEN

EqualityConstraint::
EqualityConstraint(const Index dof_id_master,const Index dof_id_slave)
    :
    dof_id_master_(dof_id_master),
    dof_id_slave_(dof_id_slave)
{
    Assert(dof_id_master_ != dof_id_slave_,
           ExcMessage("The dofs used to specify the equality constraint cannot be equal."));
    Assert(dof_id_master_ >= 0,ExcLowerRange(dof_id_master_,0));
    Assert(dof_id_slave_  >= 0,ExcLowerRange(dof_id_slave_,0));
}


Index
EqualityConstraint::
get_dof_id_master() const
{
    return dof_id_master_;
}


Index
EqualityConstraint::
get_dof_id_slave() const
{
    return dof_id_slave_;
}


void
EqualityConstraint::
print_info(LogStream &out) const
{
    out << "[ master = "<< dof_id_master_ << " , slave = " << dof_id_slave_ << " ]";
}

IGA_NAMESPACE_CLOSE
