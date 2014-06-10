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


#include <igatools/basis_functions/dofs_manager.h>
#include <igatools/base/exceptions.h>


IGA_NAMESPACE_OPEN




DofsManager::
DofsManager()
    :
    is_dofs_arrangement_open_(false),
    dofs_view_(nullptr)
{}


void
DofsManager::
dofs_arrangement_open()
{
    is_dofs_arrangement_open_ = true;
}

DofsManager::
SpaceInfo::
SpaceInfo(const Index n_dofs, const SpaceDofsView &dofs_view)
    :
    n_dofs_(n_dofs),
    offset_(0),
    dofs_view_(dofs_view)
{
    Assert(n_dofs > 0,ExcEmptyObject());
}

void
DofsManager::
add_dofs_space_view(
    const int space_id,
    const Index num_dofs_space,
    const SpaceDofsView &dofs_space_view)
{
    Assert(space_id >= 0,ExcLowerRange(space_id,0));

    spaces_info_.emplace(
        space_id,
        SpaceInfo(num_dofs_space,dofs_space_view));
}



void
DofsManager::
dofs_arrangement_close()
{
    Assert(is_dofs_arrangement_open_ == true,ExcInvalidState());

    Assert(!spaces_info_.empty(),ExcEmptyObject());

    Index offset = 0;
    for (auto &space : spaces_info_)
    {
//      auto space_id = space.first;
        auto num_dofs = space.second.n_dofs_;
        auto &dofs_view = space.second.dofs_view_;
        space.second.offset_ = offset;

        for (Index &dof : dofs_view)
            dof += offset;

        offset += num_dofs;

        auto dofs_view_begin = dofs_view.begin();
        auto view_ranges = dofs_view_begin.get_ranges();
        dofs_components_view_.insert(dofs_components_view_.end(),view_ranges.begin(),view_ranges.end());
    }

    DofsIterator dofs_begin(dofs_components_view_,0);;
    DofsIterator dofs_end(dofs_components_view_,IteratorState::pass_the_end);

    Assert(dofs_view_ == nullptr, ExcInvalidState())
    dofs_view_ = std::unique_ptr<DofsView>(new DofsView(dofs_begin,dofs_end));

    is_dofs_arrangement_open_ = false;
}

auto
DofsManager::
get_dofs_view() const -> const DofsView &
{
    Assert(is_dofs_arrangement_open_ == false,ExcInvalidState());

    Assert(dofs_view_ != nullptr, ExcNullPtr())
    return *dofs_view_;
}


int
DofsManager::
get_num_linear_constraints() const
{
    return linear_constraints_.size();
}

int
DofsManager::
get_num_equality_constraints() const
{
    return equality_constraints_.size();
}

Index
DofsManager::
get_global_dof(const int space_id, const Index local_dof) const
{
    Assert(is_dofs_arrangement_open_ == false,ExcInvalidState());

    Assert(space_id >= 0,ExcLowerRange(space_id,0));

//    const auto &space = spaces_info_.at(space_id);

    return spaces_info_.at(space_id).dofs_view_[local_dof];
}

std::vector<Index>
DofsManager::
get_global_dofs(const int space_id, const std::vector<Index> &local_dofs) const
{
    Assert(!local_dofs.empty(),ExcEmptyObject());

    std::vector<Index> global_dofs;

    for (const Index local_dof : local_dofs)
        global_dofs.emplace_back(this->get_global_dof(space_id,local_dof));

    return global_dofs;
}

void
DofsManager::
print_info(LogStream &out) const
{
    using std::endl;

    std::string tab("    ");

    out << "DofsManager infos:" << endl;

    out.push(tab);


    Assert(is_dofs_arrangement_open_ == false,ExcInvalidState());

    Assert(dofs_view_ != nullptr, ExcNullPtr())
    out << "DOFs = [ ";
    for (Index &dof : *dofs_view_)
        out << dof << " ";
    out << "]" << endl;


    Assert(!spaces_info_.empty(),ExcEmptyObject());
    int i = 0;
    for (auto &space_info : spaces_info_)
    {

        out << "Space["<< i <<"]:   ID=" << space_info.first
            << "   n_dofs=" << space_info.second.n_dofs_
            << "   DOFs=[ ";

        SpaceDofsView &dofs_space_view = const_cast<SpaceDofsView &>(space_info.second.dofs_view_);
        for (Index &dof : dofs_space_view)
            out << dof << " ";
        out << "]" << endl;

        i++;
        //*/
    }




    out << "Num. linear   constraints = " << this->get_num_linear_constraints() << endl;
    out << "Num. equality constraints = " << this->get_num_equality_constraints() << endl;

    out.pop();
}





IGA_NAMESPACE_CLOSE
