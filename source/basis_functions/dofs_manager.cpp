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

#include <map>
#include <set>

using std::vector;
using std::map;
using std::set;
using std::pair;

IGA_NAMESPACE_OPEN




DofsManager::
DofsManager()
    :
    is_dofs_arrangement_open_(false),
    dofs_view_(nullptr),
    num_unique_dofs_(0)
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

    Assert(num_dofs_space == dofs_space_view.get_num_entries(),
           ExcDimensionMismatch(num_dofs_space,dofs_space_view.get_num_entries()));

    spaces_info_.emplace(
        space_id,
        SpaceInfo(num_dofs_space,dofs_space_view));
}



void
DofsManager::
dofs_arrangement_close(const bool automatic_dofs_renumbering)
{
    Assert(is_dofs_arrangement_open_ == true,ExcInvalidState());

    Assert(!spaces_info_.empty(),ExcEmptyObject());

    Index offset = 0;
    for (auto &space : spaces_info_)
    {
        auto &dofs_view = space.second.dofs_view_;

        if (automatic_dofs_renumbering)
        {
            auto num_dofs = space.second.n_dofs_;
            space.second.offset_ = offset;

            for (Index &dof : dofs_view)
                dof += offset;

            offset += num_dofs;
        }

        auto view_ranges = dofs_view.begin().get_ranges();
        dofs_components_view_.insert(dofs_components_view_.end(),view_ranges.begin(),view_ranges.end());
    }

    DofsIterator dofs_begin(dofs_components_view_,0);;
    DofsIterator dofs_end(dofs_components_view_,IteratorState::pass_the_end);

    Assert(dofs_view_ == nullptr, ExcInvalidState())
    dofs_view_ = std::unique_ptr<DofsView>(new DofsView(dofs_begin,dofs_end));

    is_dofs_arrangement_open_ = false;


    num_unique_dofs_ = this->count_unique_dofs();
}

auto
DofsManager::
get_dofs_view() -> DofsView &
{
    Assert(is_dofs_arrangement_open_ == false,ExcInvalidState());

    Assert(dofs_view_ != nullptr, ExcNullPtr())
    return *dofs_view_;
}

auto
DofsManager::
get_dofs_view() const -> DofsConstView
{
    Assert(is_dofs_arrangement_open_ == false,ExcInvalidState());

    Assert(dofs_view_ != nullptr, ExcNullPtr())

    return DofsConstView(*dofs_view_);
}


Index
DofsManager::
get_num_dofs() const
{
    return num_unique_dofs_;
}


Index
DofsManager::
get_num_linear_constraints() const
{
    return linear_constraints_.size();
}

Index
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
equality_constraints_open()
{
    Assert(are_equality_constraints_open_ == false,
           ExcMessage("Equality constraints already opened."));
    are_equality_constraints_open_ = true;
}

void
DofsManager::
equality_constraints_close()
{
    Assert(are_equality_constraints_open_ == true,
           ExcMessage("Equality constraints already closed."));
    are_equality_constraints_open_ = false;

//    Assert(false,ExcNotImplemented());
//    AssertThrow(false,ExcNotImplemented());
}


void
DofsManager::
add_equality_constraint(const Index dof_id_master,const Index dof_id_slave)
{
    Assert(are_equality_constraints_open_ == true,
           ExcMessage("Equality constraints already closed."));

    equality_constraints_.emplace_back(EqualityConstraint(dof_id_master,dof_id_slave));
}



void
DofsManager::
linear_constraints_open()
{
    Assert(are_linear_constraints_open_ == false,
           ExcMessage("Linear constraints already opened."));
    are_linear_constraints_open_ = true;

    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}

void
DofsManager::
linear_constraints_close()
{
    Assert(are_linear_constraints_open_ == true,
           ExcMessage("Linear constraints already closed."));
    are_linear_constraints_open_ = false;

    Assert(false,ExcNotImplemented());
    AssertThrow(false,ExcNotImplemented());
}


void
DofsManager::
remove_equality_constraints_redundancies()
{
    map<Index,set<Index>> upper_sparsity_pattern_pre;

    for (const auto &eq_constr : equality_constraints_)
    {
        // retrieving the set of slaves dofs corresponding to the master dof
        auto &set_slaves_id = upper_sparsity_pattern_pre[eq_constr.get_dof_id_master()];

        // inserting the slave dof
        set_slaves_id.insert(eq_constr.get_dof_id_slave());
    }
    equality_constraints_.clear();

    map<Index,set<Index>> upper_sparsity_pattern_post;

//  LogStream out;
//  out <<"PRE" << std::endl;
    for (const auto &row_m : upper_sparsity_pattern_pre)
    {
        const Index master_id = row_m.first;
//      out << "Master = " << master_id << " ---- ";
//      out << "Slaves = [ " ;
        for (const auto &slave_id : row_m.second)
        {
//          out << slave_id << " ";

            auto &set_slaves_id_post = upper_sparsity_pattern_post[master_id];
            set_slaves_id_post.insert(slave_id);
            if (upper_sparsity_pattern_pre.count(slave_id) > 0)
            {
                auto &row_s = upper_sparsity_pattern_pre[slave_id];

                for (const auto &new_slave_id: row_s)
                    set_slaves_id_post.insert(new_slave_id);

                row_s.clear();
            }
        }
//      out << "]" << std::endl;
    }

//  out << std::endl;
//  out <<"POST" << std::endl;
    for (const auto &row_m : upper_sparsity_pattern_post)
    {
        const Index master_id = row_m.first;
//      out << "Master = " << master_id << " ---- ";
//      out << "Slaves = [ " ;
        for (const auto &slave_id : row_m.second)
        {
            equality_constraints_.emplace_back(EqualityConstraint(master_id,slave_id));
//          out << slave_id << " ";
        }
//      out << "]" << std::endl;
    }
}



Index
DofsManager::
count_unique_dofs() const
{
    Assert(is_dofs_arrangement_open_ == false,ExcInvalidState());

//    const auto &dofs = this->get_dofs_view();

    set<Index> unique_dofs(dofs_view_->begin(),dofs_view_->end());

    return unique_dofs.size();
}


SparsityPattern
DofsManager::
get_sparsity_pattern() const
{
    Assert(is_dofs_arrangement_open_ == false,ExcInvalidState());


    // build the dofs graph
    const auto & dofs_view = this->get_dofs_view();

    vector<Index> dofs_copy;
    for (const auto &dof : dofs_view)
    	dofs_copy.push_back(dof);

    Assert(!dofs_copy.empty(),ExcEmptyObject());

    SparsityPattern sparsity_pattern(dofs_copy, dofs_copy);

    using DofsInRow = set<Index>;
    DofsInRow empty_set;

    // adding the global dof keys to the map representing the dof connectivity
    for (const auto &dof : dofs_copy)
        sparsity_pattern.insert(pair<Index,DofsInRow>(dof,empty_set));

    for (const auto element_dofs : elements_dofs_view_)
    	for (const auto &dof : element_dofs)
    		sparsity_pattern[dof].insert(element_dofs.begin(),element_dofs.end());

    return (sparsity_pattern);

	Assert(false,ExcNotImplemented());
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
    Assert(are_equality_constraints_open_ == false,ExcInvalidState());
    Assert(are_linear_constraints_open_ == false,ExcInvalidState());

    Assert(dofs_view_ != nullptr, ExcNullPtr())
    out << "DOFs = [ ";
    for (Index &dof : *dofs_view_)
        out << dof << " ";
    out << "]" << endl;


    Assert(!spaces_info_.empty(),ExcEmptyObject());
    Index i = 0;
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

    out << "Num. unique dofs          = " << this->get_num_dofs() << endl;
    out << "Num. linear   constraints = " << this->get_num_linear_constraints() << endl;



    out << "Num. equality constraints = " << this->get_num_equality_constraints() << endl;
    out.push(tab);
    i = 0;
    for (const auto &eq_constr : equality_constraints_)
    {
        out << "Eq. constraint[" << i++ << "]: ";
        eq_constr.print_info(out);
        out << endl;
    }
    out.pop();


    out.pop();
}





IGA_NAMESPACE_CLOSE
