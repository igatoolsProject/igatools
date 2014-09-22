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


#include <igatools/basis_functions/space_manager.h>
#include <igatools/base/exceptions.h>

#include <map>
#include <set>


using std::map;
using std::set;
using std::pair;
using std::shared_ptr;

IGA_NAMESPACE_OPEN





SpaceManager::
SpaceManager()
    :
    is_spaces_insertion_open_(false),
    is_spaces_connectivity_open_(false),
    num_unique_dofs_(0)
{}



void
SpaceManager::
spaces_insertion_open()
{
    is_spaces_insertion_open_ = true;
}


void
SpaceManager::
spaces_insertion_close(const bool automatic_dofs_renumbering)
{
    Assert(is_spaces_insertion_open_ == true,ExcInvalidState());

    Assert(!spaces_info_.empty(),ExcEmptyObject());


    //--------------------------------------------------------------------------
    vector<DofsComponentView> dofs_components_view;

    Index offset = 0;
    for (auto &space_info_map_entry : spaces_info_)
    {
        auto &space_info = space_info_map_entry.second;

        if (automatic_dofs_renumbering)
        {
            space_info->add_dofs_offset(offset);

            offset += space_info->get_num_dofs();
        }

        auto view_ranges = space_info->get_dofs_view().begin().get_ranges();
        dofs_components_view.insert(dofs_components_view.end(),view_ranges.begin(),view_ranges.end());
    }
    //--------------------------------------------------------------------------

    dofs_view_ = DofsView(
                     DofsIterator(dofs_components_view,0),
                     DofsIterator(dofs_components_view,IteratorState::pass_the_end));

    is_spaces_insertion_open_ = false;


    num_unique_dofs_ = this->count_unique_dofs();
}


void
SpaceManager::
spaces_connectivity_open()
{
    Assert(is_spaces_connectivity_open_ == false,ExcInvalidState());
    is_spaces_connectivity_open_ = true;
}

void
SpaceManager::
spaces_connectivity_close()
{
    Assert(is_spaces_connectivity_open_ == true,ExcInvalidState());
    is_spaces_connectivity_open_ = false;
}


bool
SpaceManager::
is_spaces_connectivity_open() const
{
    return is_spaces_connectivity_open_;
}

/*
SpaceManager::
SpaceInfo::
SpaceInfo()
    :
    num_dofs_(0),
    min_dofs_id_(-1),
    max_dofs_id_(-1)
{}
//*/


SpaceManager::
SpaceInfo::
SpaceInfo(const SpacePtrVariant &space,
          const Index id,
          const int dim,
          const int codim,
          const int space_dim,
          const int range,
          const int rank,
          const Index num_dofs,
          const Index min_dofs_id,
          const Index max_dofs_id,
          const DofsView &dofs_view,
          const std::shared_ptr<const std::map<Index,DofsConstView>> elements_dofs_view)
    :
    space_(space),
    id_(id),
    dim_(dim),
    codim_(codim),
    space_dim_(space_dim),
    range_(range),
    rank_(rank),
    num_dofs_(num_dofs),
    min_dofs_id_(min_dofs_id),
    max_dofs_id_(max_dofs_id),
    dofs_view_(dofs_view),
    elements_dofs_view_(elements_dofs_view)
{
    Assert(dim_ >= 0,ExcLowerRange(dim_,0));
    Assert(codim_ >= 0,ExcLowerRange(codim_,0));
    Assert(space_dim_ > 0,ExcLowerRange(space_dim_,1));
    Assert(range_ > 0,ExcLowerRange(range_,1));
    Assert(rank_ > 0,ExcLowerRange(rank_,1));

    Assert(num_dofs_ > 0,ExcEmptyObject());
    Assert(elements_dofs_view_ != nullptr,ExcNullPtr());
    Assert(!elements_dofs_view_->empty(), ExcEmptyObject());
}

void
SpaceManager::
SpaceInfo::
add_dofs_offset(const Index offset)
{
    Assert(offset >= 0,ExcLowerRange(offset,0));
    min_dofs_id_ += offset;
    max_dofs_id_ += offset;

    for (Index &dof : dofs_view_)
        dof += offset;
}


Index
SpaceManager::
SpaceInfo::
get_num_dofs() const
{
    return num_dofs_;
}

Index
SpaceManager::
SpaceInfo::
get_min_dofs_id() const
{
    return min_dofs_id_;
}

Index
SpaceManager::
SpaceInfo::
get_max_dofs_id() const
{
    return max_dofs_id_;
}

auto
SpaceManager::
SpaceInfo::
get_space_variant() -> SpacePtrVariant &
{
    return space_;
}

auto
SpaceManager::
SpaceInfo::
get_dofs_view() -> DofsView &
{
    return dofs_view_;
}

auto
SpaceManager::
SpaceInfo::
get_dofs_view() const -> const DofsView &
{
    return dofs_view_;
}

Index
SpaceManager::
SpaceInfo::
get_id() const
{
    return id_;
}

bool
SpaceManager::
SpaceInfo::
operator==(const SpaceInfo &sp) const
{
    return (id_ == sp.id_);
}





SpaceManager::
SpacesConnection::
SpacesConnection(const SpaceInfoPtr &space_row,const SpaceInfoPtr &space_col)
    :
    space_row_(space_row),
    space_col_(space_col)
{}


bool
SpaceManager::
SpacesConnection::
operator==(const SpacesConnection &conn) const
{
    return (space_row_ == conn.space_row_) && (space_col_ == conn.space_col_);
}














auto
SpaceManager::
get_dofs_view() -> DofsView &
{
    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());

//    Assert(dofs_view_ != nullptr, ExcNullPtr())
    return dofs_view_;
}


auto
SpaceManager::
get_dofs_view() const -> DofsConstView
{
    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());

//    Assert(dofs_view_ != nullptr, ExcNullPtr())
    return DofsConstView(dofs_view_);
}



Index
SpaceManager::
get_num_dofs() const
{
    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());
    return num_unique_dofs_;
}



Index
SpaceManager::
get_num_linear_constraints() const
{
    return linear_constraints_.size();
}


Index
SpaceManager::
get_num_equality_constraints() const
{
    return equality_constraints_.size();
}


Index
SpaceManager::
get_global_dof(const int space_id, const Index local_dof) const
{
    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());

    Assert(space_id >= 0,ExcLowerRange(space_id,0));

    return spaces_info_.at(space_id)->get_dofs_view()[local_dof];
}


vector<Index>
SpaceManager::
get_global_dofs(const int space_id, const vector<Index> &local_dofs) const
{
    Assert(!local_dofs.empty(),ExcEmptyObject());

    vector<Index> global_dofs;

    for (const Index local_dof : local_dofs)
        global_dofs.emplace_back(this->get_global_dof(space_id,local_dof));

    return global_dofs;
}



bool
SpaceManager::
is_spaces_insertion_open() const
{
    return is_spaces_insertion_open_;
}


auto
SpaceManager::
SpaceInfo::
get_elements_dofs_view() const -> const std::map<Index,DofsConstView> &
{
    return *elements_dofs_view_;
}


auto
SpaceManager::
get_spaces_info() const -> const std::map<int,shared_ptr<SpaceInfo>> &
{
    return spaces_info_;
}

void
SpaceManager::
equality_constraints_open()
{
    Assert(are_equality_constraints_open_ == false,
           ExcMessage("Equality constraints already opened."));
    are_equality_constraints_open_ = true;
}


void
SpaceManager::
equality_constraints_close()
{
    Assert(are_equality_constraints_open_ == true,
           ExcMessage("Equality constraints already closed."));
    are_equality_constraints_open_ = false;
}



void
SpaceManager::
add_equality_constraint(const Index dof_id_master,const Index dof_id_slave)
{
    Assert(are_equality_constraints_open_ == true,
           ExcMessage("Equality constraints already closed."));

    equality_constraints_.emplace_back(EqualityConstraint(dof_id_master,dof_id_slave));
}




void
SpaceManager::
linear_constraints_open()
{
    Assert(are_linear_constraints_open_ == false,
           ExcMessage("Linear constraints already opened."));
    are_linear_constraints_open_ = true;
}


void
SpaceManager::
linear_constraints_close()
{
    Assert(are_linear_constraints_open_ == true,
           ExcMessage("Linear constraints already closed."));
    are_linear_constraints_open_ = false;
}



void
SpaceManager::
add_linear_constraint(const Index global_dof_id,
                      const LinearConstraintType &type,
                      const vector<Index> &dofs, const vector<Real> &coeffs, const Real rhs)
{
    Assert(are_linear_constraints_open_ == true,
           ExcMessage("Linear constraints already closed."));

    linear_constraints_.emplace(type,LC::create(global_dof_id,type,dofs,coeffs,rhs));
}


void
SpaceManager::
add_linear_constraint(std::shared_ptr<LinearConstraint> linear_constraint)
{
    Assert(are_linear_constraints_open_ == true,
           ExcMessage("Linear constraints already closed."));

    Assert(linear_constraint != nullptr, ExcNullPtr());
    linear_constraints_.emplace(linear_constraint->get_type(),linear_constraint);
}


vector<std::shared_ptr<LinearConstraint> >
SpaceManager::
get_linear_constraints(const LinearConstraintType &type_in) const
{
    //--------------------------------------------------
    // create a vector of "pure" LCType from the input argument type_in,
    // in order to use the equal_range_function
    vector<LinearConstraintType> pure_types;
    if (contains(type_in, LinearConstraintType::lagrange))
        pure_types.push_back(LinearConstraintType::lagrange);
    if (contains(type_in, LinearConstraintType::augmented_lagrange))
        pure_types.push_back(LinearConstraintType::augmented_lagrange);
    if (contains(type_in, LinearConstraintType::penalty))
        pure_types.push_back(LinearConstraintType::penalty);
    //--------------------------------------------------



    //--------------------------------------------------
    // copying the LinearConstraints of appropriate type in the return variable
    vector<std::shared_ptr<LC>> lcs;
    for (const auto type : pure_types)
    {
        auto range = linear_constraints_.equal_range(type);
        std::for_each(
            range.first,
            range.second,
            [&lcs](const auto &lc_pair)
        {
            lcs.push_back(lc_pair.second);
        });
    }
    //--------------------------------------------------

    return lcs;
}


void
SpaceManager::
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

    for (const auto &row_m : upper_sparsity_pattern_pre)
    {
        const Index master_id = row_m.first;
        for (const auto &slave_id : row_m.second)
        {
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
    }

    for (const auto &row_m : upper_sparsity_pattern_post)
    {
        const Index master_id = row_m.first;
        for (const auto &slave_id : row_m.second)
            equality_constraints_.emplace_back(EqualityConstraint(master_id,slave_id));
    }
}




Index
SpaceManager::
count_unique_dofs() const
{
    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());

//    const auto &dofs = this->get_dofs_view();

    set<Index> unique_dofs(dofs_view_.begin(),dofs_view_.end());

    return unique_dofs.size();
}

vector<std::shared_ptr<LinearConstraint> >
SpaceManager::
verify_linear_constraints(const vector<Real> &dof_values, const Real tol) const
{
    Assert(dof_values.size() == this->get_num_dofs(),
           ExcDimensionMismatch(dof_values.size(),this->get_num_dofs()));

    vector<shared_ptr<LinearConstraint> > failed_linear_constraints;
    for (const auto &lc_pair : linear_constraints_)
    {
        const auto &lc = lc_pair.second;
        if (lc->eval_absolute_error(dof_values) >= tol)
            failed_linear_constraints.push_back(lc);
    }

    return failed_linear_constraints;
}




void
SpaceManager::
print_info(LogStream &out) const
{
    using std::endl;

    std::string tab("    ");

    out << "SpaceManager infos:" << endl;

    out.push(tab);


    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());
    Assert(are_equality_constraints_open_ == false,ExcInvalidState());
    Assert(are_linear_constraints_open_ == false,ExcInvalidState());

//    Assert(dofs_view_ != nullptr, ExcNullPtr())
    out << "DOFs = [ ";
    for (const Index &dof : dofs_view_)
        out << dof << " ";
    out << "]" << endl;


    Assert(!spaces_info_.empty(),ExcEmptyObject());
    Index i = 0;
    for (const auto &space_info : spaces_info_)
    {

        out << "Space["<< i <<"]:   ID=" << space_info.first
            << "   n_dofs=" << space_info.second->get_num_dofs()
            << "   DOFs=[ ";

        const DofsView &dofs_space_view = space_info.second->get_dofs_view();
        for (const Index &dof : dofs_space_view)
            out << dof << " ";
        out << "]" << endl;

        i++;
        //*/
    }

    out << "Num. unique dofs          = " << this->get_num_dofs() << endl;


    out << "Num. linear   constraints = " << this->get_num_linear_constraints() << endl;
    out.push("   ");
    int lc_counter = 0 ;
    for (const auto &lc : linear_constraints_)
    {
        out << "Linear constraint[" << lc_counter++ << "] :" << endl;
        out.push("   ");
        lc.second->print_info(out);
        out.pop();
    }
    out.pop();


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
