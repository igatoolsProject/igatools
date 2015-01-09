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
    space_dofs_offset_(0)
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

    is_spaces_insertion_open_ = false;

    Assert(!spaces_info_.empty(),ExcEmptyObject());


    if (automatic_dofs_renumbering)
        this->perform_space_dofs_renumbering();

    this->update_dofs_view();
}

void
SpaceManager::
perform_space_dofs_renumbering()
{
    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());

    auto space_info_it = spaces_with_original_dofs_.begin();
    while (space_info_it != spaces_with_original_dofs_.end())
    {
        auto space_info = (*space_info_it);

        space_info->add_dofs_offset(space_dofs_offset_);

        space_dofs_offset_ = space_info->get_max_dofs_id() + 1;

        spaces_with_renumbered_dofs_.push_back(space_info);

        //erase the renumbered space and move the iterator to the next space to renumber
        space_info_it = spaces_with_original_dofs_.erase(space_info_it);
    }
}


void
SpaceManager::
update_dofs_view()
{
    Assert(is_spaces_insertion_open_ == false,ExcInvalidState());

    vector<DofsComponentView> dofs_components_view;
    for (auto &space_info_map_entry : spaces_info_)
    {
        auto view_ranges = space_info_map_entry.second->get_dofs_view().begin().get_ranges();
        dofs_components_view.insert(
            dofs_components_view.end(),
            view_ranges.begin(),
            view_ranges.end());
    }
    dofs_view_ = DofsView(
                     DofsIterator(dofs_components_view,0),
                     DofsIterator(dofs_components_view,IteratorState::pass_the_end));
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



SpaceManager::
SpaceInfo::
SpaceInfo(const SpacePtrVariant &space,
          const Index space_id,
          const int dim,
          const int codim,
          const int space_dim,
          const int range,
          const int rank,
          const Transformation transf_type,
          const Index num_dofs,
          const Index min_dofs_id,
          const Index max_dofs_id,
          const DofsView &dofs_view,
          const std::shared_ptr<const std::map<Index,DofsConstView>> elements_dofs_view)
    :
    space_(space),
    space_id_(space_id),
    dim_(dim),
    codim_(codim),
    space_dim_(space_dim),
    range_(range),
    rank_(rank),
    transf_type_(transf_type),
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
get_space_id() const
{
    return space_id_;
}

int
SpaceManager::
SpaceInfo::
get_dim() const
{
    return dim_;
}

int
SpaceManager::
SpaceInfo::
get_codim() const
{
    return codim_;
}

int
SpaceManager::
SpaceInfo::
get_space_dim() const
{
    return space_dim_;
}

int
SpaceManager::
SpaceInfo::
get_range() const
{
    return range_;
}

int
SpaceManager::
SpaceInfo::
get_rank() const
{
    return rank_;
}

Transformation
SpaceManager::
SpaceInfo::
get_transformation_type() const
{
    return transf_type_;
}


bool
SpaceManager::
SpaceInfo::
check_parameters(const int dim,
                 const int codim,
                 const int space_dim,
                 const int range,
                 const int rank,
                 const Transformation transf_type) const
{
    return ((dim == dim_) &&
            (codim == codim_) &&
            (space_dim == space_dim_) &&
            (range == range_) &&
            (rank == rank_) &&
            (transf_type == transf_type_));
}

bool
SpaceManager::
SpaceInfo::
operator==(const SpaceInfo &sp) const
{
    return (space_id_ == sp.space_id_);
}




SpaceManager::
SpacesConnection::
SpacesConnection(const SpaceInfoPtr &space,const bool use_dofs_connectivity_from_space)
    :
    SpacesConnection(space,space)
{
    use_dofs_connectivity_from_space_ = use_dofs_connectivity_from_space;
}

SpaceManager::
SpacesConnection::
SpacesConnection(const SpaceInfoPtr &space_row,const SpaceInfoPtr &space_col)
    :
    space_row_(space_row),
    space_col_(space_col),
    use_dofs_connectivity_from_space_(false)
{}


bool
SpaceManager::
SpacesConnection::
operator==(const SpacesConnection &conn) const
{
    return (space_row_ == conn.space_row_) && (space_col_ == conn.space_col_);
}


bool
SpaceManager::
SpacesConnection::
is_unique_space() const
{
    return (space_row_ == space_col_);
}



void
SpaceManager::
SpacesConnection::
add_dofs_connectivity(const DofsConnectivity &dofs_connectivity)
{
    extra_dofs_connectivity_.merge(dofs_connectivity);
#if 0
    for (const auto &dofs_connectivity_map_entry : dofs_connectivity)
    {
        const auto row_dof = dofs_connectivity_map_entry.first;
        const auto &col_dofs = dofs_connectivity_map_entry.second;

        extra_dofs_connectivity_[row_dof].insert(col_dofs.begin(),col_dofs.end());
    }
#endif
}


void
SpaceManager::
add_dofs_connectivity(const DofsConnectivity &dofs_connectivity)
{
    extra_dofs_connectivity_.merge(dofs_connectivity);
#if 0
    for (const auto &dofs_connectivity_map_entry : dofs_connectivity)
    {
        const auto row_dof = dofs_connectivity_map_entry.first;
        const auto &col_dofs = dofs_connectivity_map_entry.second;

        extra_dofs_connectivity_[row_dof].insert(col_dofs.begin(),col_dofs.end());
    }
#endif
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

    return DofsConstView(dofs_view_);
}



Index
SpaceManager::
get_num_linear_constraints() const
{
    Index num_linear_constraints = 0;
    for (const auto &lcs_same_type : linear_constraints_)
        num_linear_constraints += lcs_same_type.second.size();

    return num_linear_constraints;
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


auto
SpaceManager::
get_spaces_info() -> std::map<int,shared_ptr<SpaceInfo>> &
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

    auto &lcs_same_type = linear_constraints_[type];

    Assert(lcs_same_type[global_dof_id] == nullptr,
           ExcMessage("Linear constraint already added."));
    AssertThrow(lcs_same_type[global_dof_id] == nullptr,
                ExcMessage("Linear constraint already added."));
    lcs_same_type[global_dof_id] = LC::create(global_dof_id,type,dofs,coeffs,rhs);
}


void
SpaceManager::
add_linear_constraint(std::shared_ptr<LinearConstraint> linear_constraint)
{
    Assert(are_linear_constraints_open_ == true,
           ExcMessage("Linear constraints already closed."));

    Assert(linear_constraint != nullptr, ExcNullPtr());
    auto &lcs_same_type = linear_constraints_[linear_constraint->get_type()];

    Assert(lcs_same_type[linear_constraint->get_global_dof_id()] == nullptr,
           ExcMessage("Linear constraint already added."));
    AssertThrow(lcs_same_type[linear_constraint->get_global_dof_id()] == nullptr,
                ExcMessage("Linear constraint already added."));
    lcs_same_type[linear_constraint->get_global_dof_id()] = linear_constraint;
}


std::map<Index,std::shared_ptr<const LinearConstraint> >
SpaceManager::
get_linear_constraints(const LinearConstraintType &type_in) const
{
    //--------------------------------------------------
    // copying the LinearConstraints of appropriate type in the return variable
    std::map<Index,std::shared_ptr<const LC>> lcs;
    for (const auto &linear_constr_same_type : linear_constraints_)
        for (const auto &lc_pair : linear_constr_same_type.second)
            lcs[lc_pair.first] = lc_pair.second;
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



vector<std::shared_ptr<const LinearConstraint> >
SpaceManager::
verify_linear_constraints(const vector<Real> &dof_values, const Real tol) const
{
    Assert(dof_values.size() == this->get_num_col_dofs(),
           ExcDimensionMismatch(dof_values.size(),this->get_num_col_dofs()));

    vector<shared_ptr<const LinearConstraint> > failed_linear_constraints;
    for (const auto &lcs_same_type : linear_constraints_)
    {
        for (const auto &lc_pair : lcs_same_type.second)
        {
            const auto &lc = lc_pair.second;
            if (lc->eval_absolute_error(dof_values) >= tol)
                failed_linear_constraints.push_back(lc);
        }
    }


    return failed_linear_constraints;
}



std::set<Index>
SpaceManager::
get_row_dofs() const
{
    std::set<Index> row_dofs;

    for (const auto &sp_conn : spaces_connections_)
    {
        const auto row_dofs_current_space = sp_conn.get_row_dofs();
        row_dofs.insert(row_dofs_current_space.begin(),row_dofs_current_space.end());
    }

    for (const auto extra_row : extra_dofs_connectivity_)
        row_dofs.insert(extra_row.first);

    return row_dofs;
}



std::set<Index>
SpaceManager::
get_col_dofs() const
{
    std::set<Index> col_dofs;

    for (const auto &sp_conn : spaces_connections_)
    {
        const auto col_dofs_current_space = sp_conn.get_col_dofs();
        col_dofs.insert(col_dofs_current_space.begin(),col_dofs_current_space.end());
    }

    for (const auto extra_row : extra_dofs_connectivity_)
        col_dofs.insert(extra_row.second.begin(),extra_row.second.end());


    return col_dofs;
}

Index
SpaceManager::
get_num_row_dofs() const
{
    return this->get_row_dofs().size();
}

Index
SpaceManager::
get_num_col_dofs() const
{
    return this->get_col_dofs().size();
}



auto
SpaceManager::
get_sparsity_pattern() const -> shared_ptr<const DofsConnectivity>
{
    auto sparsity_pattern = shared_ptr<DofsConnectivity>(new DofsConnectivity);

    Assert(!spaces_connections_.empty(),ExcEmptyObject());
    for (const auto &sp_conn : spaces_connections_)
    {
        if (sp_conn.is_unique_space() && sp_conn.use_dofs_connectivity_from_space())
        {
            // adding the contribution of the dofs defined within the space itself
            const auto &space = sp_conn.get_space_row();
            for (const auto element_dofs : space.get_elements_dofs_view())
                for (const auto &dof : element_dofs.second)
                    (*sparsity_pattern)[dof].insert(element_dofs.second.begin(),element_dofs.second.end());
        }

        // adding the extra contribution to the connectivity defined within the spaces connection
        sparsity_pattern->merge(sp_conn.get_extra_dofs_connectivity());
    }

    // adding the extra contribution to the remaining extra connectivity
    // (i.e. the connectivity not declared within a SpacesConnection obj)
    sparsity_pattern->merge(extra_dofs_connectivity_);

    return sparsity_pattern;
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

    out << "Num. row   dofs = " << this->get_num_row_dofs() << endl;
    auto row_dofs_set = this->get_row_dofs();
    vector<Index> row_dofs(row_dofs_set.begin(),row_dofs_set.end());
    out.begin_item("Row dofs:");
    row_dofs.print_info(out);
    out.end_item();


    out << "Num. colum dofs = " << this->get_num_col_dofs() << endl;
    auto col_dofs_set = this->get_col_dofs();
    vector<Index> col_dofs(col_dofs_set.begin(),col_dofs_set.end());
    out.begin_item("Col dofs:");
    col_dofs.print_info(out);
    out.end_item();


    out << "Num. linear   constraints = " << this->get_num_linear_constraints() << endl;
    out.push("   ");
    int lc_counter = 0 ;
    const auto lcs = this->get_linear_constraints();
    for (const auto &lc : lcs)
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
