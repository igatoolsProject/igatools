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
#include "../tests.h"


#include <igatools/base/function_lib.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element.h>
#include <igatools/geometry/topology.h>

#include <igatools/basis_functions/space_tools.h>
#include <igatools/basis_functions/bspline_space.h>
#include <igatools/basis_functions/bspline_element_accessor.h>
#include <igatools/basis_functions/physical_space.h>
#include <igatools/basis_functions/physical_space_element_accessor.h>
#include <igatools/basis_functions/multi_patch_space.h>
#include <igatools/io/writer.h>
#include <iostream>

#include <igatools/base/quadrature_lib.h>
#include <igatools/linear_algebra/dense_matrix.h>
#include <igatools/linear_algebra/dense_vector.h>
#include <igatools/linear_algebra/distributed_matrix.h>
#include <igatools/linear_algebra/distributed_vector.h>
#include <igatools/linear_algebra/linear_solver.h>
#include <igatools/linear_algebra/dof_tools.h>

#include <math.h>

using namespace iga;
using namespace std;





//--------------------------------------------------------------------
template<int dim, int dim_field>
class Mortar_Interface
{
public:
    using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
    using RefSpaceField_border  = BSplineSpace<dim-1,dim_field,1>;
    using RefSpaceField_ptr     = shared_ptr<RefSpaceField>;
    using PushFw                = PushForward<Transformation::h_grad, dim,0>;
    using PushFw_border         = PushForward<Transformation::h_grad, dim-1,1>;
    using PhySpace              = PhysicalSpace<RefSpaceField, PushFw>;
    using PhySpace_border       = PhysicalSpace<RefSpaceField_border, PushFw_border>;
    using PhySpace_ptr          = shared_ptr<PhySpace>;

    using FuncMap = Function<dim, dim>;
    using GradientMap = typename FuncMap::Gradient;

    using FuncMap_border = Function<dim-1, dim>;
    using GradientMap_border = typename FuncMap_border::Gradient;

private:
    boundary_id mortar_id_;
    int degree_multiplier_;// has to change
    //array<int,dim> degree_multiplier_;

    int slave_patch_nb_;
    int slave_mortar_face_nb_;
    shared_ptr<PhySpace> slave_space_;

    int master_patch_nb_;
    int master_mortar_face_nb_;
    shared_ptr<PhySpace> master_space_;

    shared_ptr<PhySpace_border> slave_face_space_;
    vector<Index> slave_face_dof_num_;
    shared_ptr<PhySpace_border> master_face_space_;
    vector<Index> master_face_dof_num_;
    vector<Index> face_all_dof_num_;

    shared_ptr<PhySpace> multiplier_space_;
    shared_ptr<PhySpace_border> multiplier_face_space_;
    vector<Index> multiplier_face_dof_num_;



    bool active_crosspoint_;

    using pair_coef_value  = pair<Index, Real>;
    using LinearConstraint = pair<vector<pair_coef_value>, Real>;

    const int hdirch_boundary_id_val=2;
    const int dirch_boundary_id_val=3;
    const int mortar_boundary_id_val=7;

    //check_quad_pt_nb() const {
    //}

    vector<LinearConstraint> LCs_;

    void integration();

    void get_loc_face_crosspoint_modif(shared_ptr<CartesianGrid<dim>> slave_grid);
    vector<int> loc_face_crosspoint_;


public:
    /** @name Constructors and destructor */
    ///@{
    /** Default constructor. */
    Mortar_Interface() = delete;

    /** Constructor. */
    Mortar_Interface(const boundary_id &mortar_id, const int &degree_multiplier, const int &slave_patch_nb, const int &slave_mortar_face_nb, shared_ptr<PhySpace> slave_space,
                     const int &master_patch_nb, const int &master_mortar_face_nb, shared_ptr<PhySpace> master_space, shared_ptr<DofsManager> dof_manager, bool active_crosspoint)
        :
        mortar_id_(mortar_id),
        degree_multiplier_(degree_multiplier),
        slave_patch_nb_(slave_patch_nb),
        slave_mortar_face_nb_(slave_mortar_face_nb),
        slave_space_(slave_space),
        master_patch_nb_(master_patch_nb),
        master_mortar_face_nb_(master_mortar_face_nb),
        master_space_(master_space),
        multiplier_space_(
            PhySpace::create(
                RefSpaceField::create(degree_multiplier_,slave_space_->get_grid()),
                const_pointer_cast<PushFw>(slave_space_->get_push_forward()),
                1)
        ),
        active_crosspoint_(active_crosspoint)
    {
        cout<<"A LC Interface built"<<endl;
        slave_face_space_=slave_space_->get_face_space(slave_mortar_face_nb_, slave_face_dof_num_);
        slave_face_dof_num_=dof_manager->get_global_dofs(slave_patch_nb_, slave_face_dof_num_);


        master_face_space_=master_space_->get_face_space(master_mortar_face_nb_, master_face_dof_num_);
        master_face_dof_num_=dof_manager->get_global_dofs(master_patch_nb_, master_face_dof_num_);

        face_all_dof_num_=slave_face_dof_num_;
        for (int i=slave_face_dof_num_.size(); i<slave_face_dof_num_.size()+master_face_dof_num_.size(); ++i)
        {
            face_all_dof_num_.push_back(master_face_dof_num_[i-(slave_face_dof_num_.size())]);
        }

        multiplier_face_space_=multiplier_space_->get_face_space(slave_mortar_face_nb_, multiplier_face_dof_num_);


        //auto au=slave_face_space_->get_reference_space();
        //auto a=au->get_degree();
        //out<<a(1)[0]<<endl;
        //a.print_info(out);
    }


    /** Copy constructor. */
    Mortar_Interface(const  Mortar_Interface &lc_interface) = delete;


    /** Move constructor. */
    Mortar_Interface(Mortar_Interface &&lc_interface) = default;


    /** Destructor. */
    ~Mortar_Interface() = default;
    ///@}


    /** Getters. */
    const int get_degree_multiplier() const
    {
        return degree_multiplier_;
    }

    const  vector<LinearConstraint> get_LCs()
    {
        this->integration();
        return LCs_;
    }


    int is_member_vec(vector<Index>   &vec, int curr_val);
    vector<int> get_local_remove(vector<int> &remove_vec, vector<int> &vec);

};


template<int dim, int dim_field>
int Mortar_Interface<dim, dim_field>::is_member_vec(vector<Index>   &vec, int curr_val)
{
    bool find(false);
    int i(0);
    while (!((find==true) | (i>vec.size())))
    {
        if (curr_val==vec[i])
        {
            find=true;
        }
        ++i;
    }
    if (find==true)
    {
        return i-1;
    }
    else
    {
        return -1;
    }
}

template<int dim, int dim_field>
vector<int> Mortar_Interface<dim, dim_field>::get_local_remove(vector<int> &remove_vec, vector<int> &vec)
{
    vector<int> temp_vec;

    for (int i=0; i<vec.size(); ++i)
    {
        auto tmp=is_member_vec(remove_vec, vec[i]);
        if (tmp!=-1)
        {
            temp_vec.push_back(tmp);
        }
    }
    return temp_vec;
}




template<int dim, int dim_field>
void Mortar_Interface<dim, dim_field>::get_loc_face_crosspoint_modif(shared_ptr<CartesianGrid<dim>> slave_grid)
{

    vector<int> temp_list;
    if (active_crosspoint_)
    {
        if (dim==2)
        {
            if ((slave_mortar_face_nb_==0) | (slave_mortar_face_nb_==1))
                temp_list=vector<int> {2,3};
            else
                temp_list=vector<int> {0,1};
        }
        else if (dim==3)
        {
            if ((slave_mortar_face_nb_==0) | (slave_mortar_face_nb_==1))
                temp_list=vector<int> {2,3,4,5};
            else if ((slave_mortar_face_nb_==2) | (slave_mortar_face_nb_==3))
                temp_list=vector<int> {0,1,4,5};
            else
                temp_list=vector<int> {0,1,2,3};

        }
        else
        {
            Assert(false,ExcNotImplemented());
        }


        for (int i=0; i<temp_list.size(); ++i)
        {
            auto temp_bnd_id=slave_grid->get_boundary_id(temp_list[i]);
            if ((temp_bnd_id==hdirch_boundary_id_val) | (temp_bnd_id==dirch_boundary_id_val) | (temp_bnd_id==mortar_boundary_id_val))
                loc_face_crosspoint_.push_back(i) ;
        }
    }
}



template<int dim, int dim_field>
void Mortar_Interface<dim, dim_field>::integration()
{


    // Get the function number living on one face element of each face space
    const int n_slave_basis      = slave_face_space_->get_num_basis_per_element();
    const int n_master_basis     = master_face_space_->get_num_basis_per_element();
    const int n_multiplier_basis = multiplier_face_space_->get_num_basis_per_element();


    // Face quadrature
    const int int_pts(degree_multiplier_+1);
    const Quadrature<dim-1> joint_face_quad(QGauss<dim-1>(int_pts+1));
    const int n_qp                         = joint_face_quad.get_num_points();
    auto  face_w_unit_domain               = joint_face_quad.get_weights().get_flat_tensor_product();
    auto  face_quad_pts_unit_domain        = joint_face_quad.get_points().get_flat_cartesian_product();




    // Face joint grid construction -> parametric joint mesh
    auto slave_grid  = slave_space_->get_grid();

    this->get_loc_face_crosspoint_modif(slave_grid);
    //out<<"slave face"<<slave_mortar_face_nb_<<endl;
    //out<<"loc_face_crosspoint_"<<loc_face_crosspoint_<<endl;
    auto master_grid = master_space_->get_grid();

    vector< Index >  joint_face_grid_to_joint_grid_slave;
    vector< Index >  joint_face_grid_to_joint_grid_master;
    auto slave_grid_to_face  = make_shared<map<int,int> >();
    auto master_grid_to_face = make_shared<map<int,int> >();
    auto face_slave_grid     = slave_grid->get_face_grid(slave_mortar_face_nb_, *slave_grid_to_face);
    auto face_master_grid    = master_grid->get_face_grid(master_mortar_face_nb_, *master_grid_to_face);

    shared_ptr<CartesianGrid<dim-1>> face_joint_grid;
    face_joint_grid=make_shared<CartesianGrid<dim-1>>(build_cartesian_grid_union(*face_slave_grid, *face_master_grid, joint_face_grid_to_joint_grid_slave,
                                                      joint_face_grid_to_joint_grid_master));
    out<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    out<<"Face joint grid info"<<endl;
    face_joint_grid->print_info(out);
    out<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;


    // Local mortar matrix
    DenseMatrix loc_e(multiplier_face_dof_num_.size(),face_all_dof_num_.size());
    loc_e.clear();



    //Local mortar correction matrix
    DenseMatrix loc_ce(multiplier_face_dof_num_.size(),face_all_dof_num_.size());
    loc_ce.clear();
    vector<vector<Index>> loc_face_crosspoint_dof;
    for (int i=0; i<loc_face_crosspoint_.size(); ++i)
    {
        vector<Index> temp_num;
        auto temp=multiplier_face_space_->get_face_space(loc_face_crosspoint_[i],temp_num);
        loc_face_crosspoint_dof.push_back(temp_num);
    }

    vector<vector<int>> remove_multiplier_dof;
    vector<vector<int>> modify_multiplier_dof;
    for (int i=0; i<loc_face_crosspoint_dof.size(); ++i)
    {
        vector<int> temp_dof;
        vector<int> temp_dof2;
        for (int j=0; j<loc_face_crosspoint_dof[i].size(); ++j)
        {
            temp_dof.push_back(multiplier_face_dof_num_[loc_face_crosspoint_dof[i][j]]);
        }
        remove_multiplier_dof.push_back(temp_dof);
        for (int j=0; j<multiplier_face_dof_num_.size(); ++j)
        {
            if (is_member_vec(temp_dof, multiplier_face_dof_num_[j])==-1)
            {
                temp_dof2.push_back(multiplier_face_dof_num_[j]);
            }
        }
        modify_multiplier_dof.push_back(temp_dof2);
    }
    //out<<"LISTE DOF"<<loc_face_crosspoint_dof[0]<<endl;
    //using     Derivative = Derivatives< dim, 1, 1, 1>;


    for (int i=0; i<loc_face_crosspoint_dof.size(); ++i)
    {
        for (int j=0; j<loc_face_crosspoint_dof[i].size(); ++j)
        {
            out<<"LISTE DOF "<<"loc face nb: "<<i<<", face nb:"<<loc_face_crosspoint_[i]<<", loc dof nb: "<<j<<", face dof nb: "<<loc_face_crosspoint_dof[i][j]<<endl;
        }
        out<<"remove"<<remove_multiplier_dof[i]<<endl;
        out<<"modify"<<modify_multiplier_dof[i]<<endl;
    }


    // Loop on the face joint grid to lead the integration on the slave side
    for (auto felem_j: *face_joint_grid)
    {
        auto felem_j_grid = felem_j.get_grid();
        auto felem_j_id   = felem_j.get_flat_index();

        auto face_quad_pts_ref = felem_j.transform_points_unit_to_reference(face_quad_pts_unit_domain);


        out<<"face joint grid element nb:"<<felem_j_id<<endl;
        auto face_meas=static_cast<CartesianGridElement<dim-1>&>(felem_j).get_measure();
        out<<"face_meas"<<face_meas<<endl;




        auto felem_slave        = slave_face_space_->get_element(joint_face_grid_to_joint_grid_slave[felem_j_id]);
        auto dofs_face_slave    = felem_slave.get_local_to_global();

        auto felem_master        = master_face_space_->get_element(joint_face_grid_to_joint_grid_master[felem_j_id]);
        auto dofs_face_master    = felem_master.get_local_to_global();

        auto felem_multiplier     = multiplier_face_space_->get_element(joint_face_grid_to_joint_grid_slave[felem_j_id]);
        auto dofs_face_multiplier = felem_multiplier.get_local_to_global();

        out<<"All slave dofs: "<<slave_face_dof_num_<<", current face slave dofs"<<dofs_face_slave<<endl;
        out<<"All master dofs: "<<master_face_dof_num_<<", current face master dofs"<<dofs_face_master<<endl;
        out<<"All multiplier dofs: "<<multiplier_face_dof_num_<<", current face multiplier dofs"<<dofs_face_multiplier<<endl;


        //get the slave mapping from the slave physical space
        vector<GradientMap_border> face_slave_map_grad(n_qp);
        auto temp_face_slave_map=slave_face_space_->get_push_forward()->get_mapping();
        temp_face_slave_map->evaluate_gradients_at_points(face_quad_pts_ref, face_slave_map_grad);


        // get the location of the current points on each unit domain
        auto cgea_felem_slave                = felem_slave.as_cartesian_grid_element_accessor();
        auto curr_slave_face_quad_pts_unit   = cgea_felem_slave.transform_points_reference_to_unit(face_quad_pts_ref);

        auto cgea_felem_master               = felem_master.as_cartesian_grid_element_accessor();
        auto curr_master_face_quad_pts_unit  = cgea_felem_master.transform_points_reference_to_unit(face_quad_pts_ref);



        // get the function values at the current unit points
        auto basis_slave      = felem_slave.evaluate_basis_values_at_points(curr_slave_face_quad_pts_unit);
        auto basis_multiplier = felem_multiplier.evaluate_basis_values_at_points(curr_slave_face_quad_pts_unit);
        auto basis_master     = felem_master.evaluate_basis_values_at_points(curr_master_face_quad_pts_unit);


        //out<<felem_j.is_boundary(0)<<felem_j.is_boundary(1)<<endl;
        //out<<"is bound"<<felem_j.is_boundary()<<endl;
        //out<<"is bound"<<felem_j.is_boundary(0)<<endl;

        auto felem_multiplier_ref=felem_multiplier.get_ref_space_accessor();
        ////auto basis_slave_ref=felem_slave_ref.evaluate_basis_values_at_points(curr_slave_face_quad_pts_unit);
        ////out<<basis_slave.get_function_view(0)[0]<<basis_slave_ref.get_function_view(0)[0];


        //vector<int> loc_der_var;
        int face_el_type(0);
        int curr_face_crosspoint(-1);
        //if (!loc_face_crosspoint_.empty()){
        //vector<vector<int>> curr_elem_dof;
        for (int i=0; i<loc_face_crosspoint_.size(); ++i)
        {
            if ((felem_j.is_boundary(loc_face_crosspoint_[i]) && face_el_type==0))
            {
                face_el_type=1;
                curr_face_crosspoint=loc_face_crosspoint_[i];
            }
            else if ((felem_j.is_boundary(loc_face_crosspoint_[i]) && face_el_type!=0))
            {
                face_el_type=2;
                curr_face_crosspoint=curr_face_crosspoint+2*loc_face_crosspoint_[i];
            }
        }

        /*if (dim==2){
            loc_der_var.push_back(0);
        }
        else {
            if ((loc_face_crosspoint_[i]==0)|(loc_face_crosspoint_[i]==1))
                       loc_der_var.push_back(0); //to check 1 ?
            else
                       loc_der_var.push_back(1);
        }*/

        //vector<int> temp_curr_elem_dof;
        //for (int s=0; s<dofs_face_multiplier.size(); ++s){ // pas trjs member c est ici que l on classe cx a enlever et cx a garder
        //  temp_curr_elem_dof.push_back(is_member_vec(loc_face_crosspoint_dof[i], dofs_face_multiplier[s]));
        //}
        //curr_elem_dof.push_back(temp_curr_elem_dof);  // ce que l on enleve !!!

        if (face_el_type!=0)
        {
            if (degree_multiplier_==1)
                auto deri_multiplier=felem_multiplier_ref. template evaluate_basis_derivatives_at_points<1>(curr_slave_face_quad_pts_unit);
            else if (degree_multiplier_==2)
                auto deri_multiplier=felem_multiplier_ref. template evaluate_basis_derivatives_at_points<2>(curr_slave_face_quad_pts_unit);
            //auto deri_phi_m=deri_multiplier.get_function_view(0);}
            //out<<"deri"<<deri_phi_m[0]<<endl;
            //out<<"deri"<<deri_phi_m[2]<<endl;}
            //else if (degree_multiplier_==3)
            //  auto deri_multiplier=felem_multiplier_ref. template evaluate_basis_derivatives_at_points<3>(curr_slave_face_quad_pts_unit);
            //else if (degree_multiplier_==4)
            //  auto deri_multiplier=felem_multiplier_ref. template evaluate_basis_derivatives_at_points<4>(curr_slave_face_quad_pts_unit);
            //else if (degree_multiplier_==5)
            //  auto deri_multiplier=felem_multiplier_ref. template evaluate_basis_derivatives_at_points<5>(curr_slave_face_quad_pts_unit);
            else
                Assert(false,ExcNotImplemented());
        } //if (face_el_type!=0)


        //}



        /*auto phi_m = basis_multiplier.get_function_view(i);
        if (!loc_face_crosspoint_.empty()){
            if (felem_j.is_boundary(loc_face_crosspoint_[0]){
                for (qp = 0; qp < n_qp; ++qp)
                    if (is_member_vec(loc_face_crosspoint_dof[0], dofs_face_multiplier[i])~=-1)
                        phi_m[qp]=phi_m[qp]*0.;
                    else {
                            int loc_l=dofs_face_multiplier[i]%curr_elem_dof[0].size();
                            phi_l=basis_multiplier.get_function_view(loc_l);
                            deri_phi_l=deri_multiplier.get_function_view(loc_l);

                            deri_phi_m=deri_multiplier.get_function_view(i);
                            //phi_m[qp]=phi_m[qp]-deri_phi_m[qp][loc_der_var[0]/deri_phi_l[qp][loc_der_var[0]*phi_l[qp];
                            if (curr_elem_dof.size()>1)
                            auto phi_l_hat= phi_m





                            //a

                    }
                }
            }
        }*/


        auto multiplier_degree_dir=multiplier_face_space_->get_degree();


        //add the contributions of the current element
        for (int i = 0; i < n_multiplier_basis; ++i)
        {
            auto phi_m = basis_multiplier.get_function_view(i);
            for (int j = 0; j < n_slave_basis; ++j)
            {
                auto phj = basis_slave.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                {


                    ///////////////////
                    if (face_el_type==0)
                        auto phi_m_qp=phi_m[qp];
                    else if (face_el_type==1)
                    {
                        if (dim==2)
                            auto loc_l=curr_face_crosspoint*(n_multiplier_basis-1);
                        else if (dim==3)
                        {
                            //vector<int> curr_glob_dof;
                            //for (int nb=0; nb<dofs_face_multiplier.size(); ++nb){
                            //  curr_glob_dof.push_back(multiplier_face_dof_num_[dofs_face_multiplier[nb]]);
                            //}
                            if ((curr_face_crosspoint==0) || (curr_face_crosspoint==1))
                            {
                                //loc_n2=curr_glob_dof.size();
                                int loc_l2=multiplier_degree_dir(0)[1];
                                int loc_l1=n_multiplier_basis/loc_l2;
                                if (curr_face_crosspoint==0)
                                    int loc_l=(i/loc_l1)*loc_l1;
                                else
                                    int loc_l=(i/loc_l1)*loc_l1+(loc_l1-1);
                            }
                            else
                            {
                                //loc_n1=curr_glob_dof.size();
                                int loc_l1=multiplier_degree_dir(0)[0];
                                int loc_l2=n_multiplier_basis/loc_l1;
                                if (curr_face_crosspoint==2)
                                    int loc_l=i-(i/loc_l1)*loc_l1;
                                else
                                    int loc_l=i+(loc_l2-((i/loc_l1)+1))*loc_l1;
                            }
                        } // if (dim==3)


                        auto phi_m=basis_multiplier.get_function_view(i);
                        auto deri_phi_m=deri_multiplier.get_function_view(i);
                        if (is_member_vec(remove_multiplier_dof,multiplier_face_dof_num_[dofs_face_multiplier[i]])!=-1)
                            phi_m_qp=0.;
                        else
                        {
                            auto phi_l=basis_multiplier.get_function_view(loc_l);
                            auto phi_l_qp=phi_l[qp];
                            auto deri_phi_l=deri_multiplier.get_function_view(loc_l);
                            if (dim==1)
                            {
                                auto deri_phi_l_qp=deri_phi_l[qp];
                            }
                            else if (dim==2)
                            {
                                if ((curr_face_crosspoint==0) || (curr_face_crosspoint==1))
                                    int dir_var(0);
                                else
                                    int dir_var(1);
                                auto deri_phi_l_qp=deri_phi_l[qp][dir_var];
                                auto deri_phi_m_qp=deri_phi_m[qp][dir_var];
                            }


                            phi_m_qp=phi_m[qp]-(deri_phi_m_qp)/(deri_phi_l_qp)*phi_l_qp;
                        } // is_member



                    }// if face_el_type==1
                    /*  else {

                            loc_k2=multiplier_degree_dir(0)[1];
                            loc_k1=n_multiplier_basis/loc_k2;
                            if ((curr_face_crosspoint==4) || (curr_face_crosspoint==6))
                                auto loc_k=(i/loc_k1)*loc_k1;
                            else
                                auto loc_k=(i/loc_k1)*loc_k1+(loc_k1-1);


                            loc_m1=multiplier_degree_dir(0)[0];
                            loc_m2=n_multiplier_basis/loc_m1;
                            if ((curr_face_crosspoint==4) || (curr_face_crosspoint==5))
                                auto loc_m=i-(i/loc_m1)*loc_m1;
                            else
                                auto loc_m=i+(loc_m2-((i/loc_m1)+1))*loc_m1;



                            auto phi_m=basis_multiplier.get_function_view(i);
                            auto deri_phi_m=deri_multiplier.get_function_view(i);
                            if (is_member_vec(remove_multiplier_dof,multiplier_face_dof_num_[dofs_face_multiplier[i]])!=-1)
                                phi_m_qp=0.;
                            else {
                                auto phi_k=basis_multiplier.get_function_view(loc_k);
                                auto phi_k_qp=phi_k[qp];
                                auto deri_phi_k=deri_multiplier.get_function_view(loc_k);

                                auto phi_m=basis_multiplier.get_function_view(loc_m);
                                auto phi_m_qp=phi_m[qp];
                                auto deri_phi_m=deri_multiplier.get_function_view(loc_m);

                                auto deri_phi_l0_qp=deri_phi_l[qp][0];
                                auto deri_phi_m0_qp=deri_phi_m[qp][0];
                                auto deri_phi_l1_qp=deri_phi_l[qp][1];
                                auto deri_phi_m1_qp=deri_phi_m[qp][1];

                                //phi_m_qp=phi_m[qp]-(deri_phi_m_qp)/(deri_phi_l_qp)*phi_l_qp;
                            }





                        }
                        //
                        */


                    loc_e(dofs_face_multiplier[i],dofs_face_slave[j])=loc_e(dofs_face_multiplier[i],dofs_face_slave[j])+
                                                                      scalar_product(phi_m_qp,phj[qp])*face_meas*face_w_unit_domain[qp]*determinant<dim-1,dim>(face_slave_map_grad[qp]);
                }
            }
            for (int j = 0; j < n_master_basis; ++j)
            {
                auto phj = basis_master.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                {

                    //modif here
                    //if
                    auto phi_m_qp=phi_m[qp];
                    //else
                    //end

                    loc_e(dofs_face_multiplier[i],dofs_face_master[j]+slave_face_dof_num_.size())=loc_e(dofs_face_multiplier[i],dofs_face_master[j]+slave_face_dof_num_.size())+
                            + scalar_product(phi_m_qp,phj[qp])*face_meas*face_w_unit_domain[qp]*determinant<dim-1,dim>(face_slave_map_grad[qp]);
                }
            }// for j
        }// for i
        out<<"--------------------"<<endl;
    } // for elem_j



    // fill the linear constraint data LCs_
    // as LCs[curr multiplier dof nb].first[curr global dof nb].first=global dof nb (i.e. slave or master dof nb)
    // as LCs[curr multiplier dof nb].first[curr global dof nb].second=coefficient related to the corresponding global dof nb
    // as LCs[curr multiplier dof nb].second=value in the rhs
    for (int i=0; i<multiplier_face_dof_num_.size(); ++i)
    {
        vector<pair_coef_value> curr_line;
        for (int k=0; k!=face_all_dof_num_.size(); ++k)
        {
            curr_line.push_back(pair<Index, Real>(face_all_dof_num_[k],loc_e(i,k)));
        }
        LCs_.push_back(pair<vector<pair_coef_value>, Real>(curr_line, 0.0));
    }


}

//--------------------------------------------------------------------







template <int dim>
shared_ptr<Mapping<dim,0>> create_map(const vector<Real> &control_pts)
{
    auto knots = CartesianGrid<dim>::create(2);

    using RefSpaceMap = BSplineSpace<dim,dim,1>;

    const int p = 1;
    auto ref_space_map = RefSpaceMap::create(p,knots);

    Assert(control_pts.size() == ref_space_map->get_num_basis(),
           ExcDimensionMismatch(control_pts.size(),ref_space_map->get_num_basis()));

    return IgMapping<RefSpaceMap>::create(ref_space_map, control_pts);
}




template <int dim>
shared_ptr<Mapping<dim,0>> create_geo_0()
{
    const int n_basis = dim * std::pow(2,dim);
    vector<Real> control_pts(n_basis);

    if (dim == 1)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;
    }
    else if (dim == 2) // x coordinates of all the points, then y
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;
    }
    else if (dim == 3)
    {
        int id = 0 ;
        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 0.0 ;
        control_pts[id++] = 0.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

        control_pts[id++] = 1.0 ;
        control_pts[id++] = 1.0 ;

    }
    return create_map<dim>(control_pts);
}



template <int dim>
shared_ptr<Mapping<dim,0>> create_geo_1()
{
    const int n_basis = dim * std::pow(2,dim);
    vector<Real> control_pts(n_basis);

    if (dim == 1)
    {
        control_pts[0] = 1.0 ;
        control_pts[1] = 2.0 ;
    }
    else if (dim == 2)
    {
        control_pts[0] = 1.0 ;
        control_pts[1] = 2.0 ;

        control_pts[2] = 1.0 ;
        control_pts[3] = 2.0 ;

        control_pts[4] = 0.0 ;
        control_pts[5] = 0.0 ;

        control_pts[6] = 1.0 ;
        control_pts[7] = 1.0 ;
    }
    else if (dim == 3)
    {
        control_pts[0] = 1.0 ;
        control_pts[1] = 2.0 ;
        control_pts[2] = 1.0 ;

        control_pts[3] = 2.0 ;
        control_pts[4] = 1.0 ;
        control_pts[5] = 2.0 ;

        control_pts[6] = 1.0 ;
        control_pts[7] = 2.0 ;
        control_pts[8] = 0.0 ;

        control_pts[9] = 0.0 ;
        control_pts[10] = 1.0 ;
        control_pts[11] = 1.0 ;

        control_pts[12] = 0.0 ;
        control_pts[13] = 0.0 ;
        control_pts[14] = 1.0 ;

        control_pts[15] = 1.0 ;
        control_pts[16] = 0.0 ;
        control_pts[17] = 0.0 ;

        control_pts[18] = 0.0 ;
        control_pts[19] = 0.0 ;
        control_pts[20] = 1.0 ;

        control_pts[21] = 1.0 ;
        control_pts[22] = 1.0 ;
        control_pts[23] = 1.0 ;
    }

    return create_map<dim>(control_pts);
}







int main()
{


    //LogStream out;
    const int dim(2);
    const int dim_field(1);
    vector<int> degrees(2,2);
    const int deg(2);


    // Mapping $\R^2 \rightarrow \R^2$, scalar field
    vector<shared_ptr<Mapping<dim,0>>> maps_d21;
    maps_d21.push_back(create_geo_0<dim>());
    maps_d21.push_back(create_geo_1<dim>());

    using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
    using RefSpaceField_ptr     = shared_ptr<RefSpaceField>;
    using PushFw                = PushForward<Transformation::h_grad, dim,0>;
    using PhySpace              = PhysicalSpace<RefSpaceField, PushFw>;
    using PhySpace_ptr          = shared_ptr<PhySpace>;





    vector<RefSpaceField_ptr>   ref_spaces_field;
    vector<PhySpace_ptr>        spaces;

    MultiPatchSpace<PhySpace> domain_multip;
    domain_multip.patch_insertion_open();
    ////////
    const boundary_id mortar_id = 7;
    ////////
    vector<set<boundary_id>> mortar_faces;
    //2D
    //mortar_faces.push_back({1});
    //mortar_faces.push_back({0});
    //3D
    mortar_faces.push_back({1});
    mortar_faces.push_back({0});

    for (uint i=0; i!=maps_d21.size(); ++i)
    {
        ref_spaces_field.push_back(RefSpaceField::create(degrees[i], maps_d21[i]->get_grid()));
        spaces.push_back(PhySpace::create(ref_spaces_field[i], PushFw::create(maps_d21[i]),i));
        //      spaces[i]->print_info(out);
        if (i==1)
        {
            //spaces[i]->refine_h(2);
        }// Problem here with a value different than 2
        if (i==0)
        {
            spaces[i]->refine_h(2);
        }// Problem here with a value different than 2
        if (i==0)
        {
            spaces[i]->get_grid()->set_boundary_id(2,3);
            spaces[i]->get_grid()->set_boundary_id(4,3);
        }
        //spaces[i]->print_info(out);
        auto tmp_mortar_faces=mortar_faces[i];
        auto crr_mortar=tmp_mortar_faces.begin(), end_mortar=tmp_mortar_faces.end();
        for (; crr_mortar!=end_mortar; ++crr_mortar)
        {
            spaces[i]->get_grid()->set_boundary_id(*crr_mortar,mortar_id);
        }

        domain_multip.add_patch(spaces[i]);
    }

    domain_multip.patch_insertion_close();
    auto dof_manager=domain_multip.get_dofs_manager();

    //spaces[0]->get_grid()->set_boundary_id(2,2)


    auto my_lc=Mortar_Interface<dim,dim_field>(mortar_id, deg, 0, 1, spaces[0], 1, 0, spaces[1], dof_manager,true);
    //my_lc.integration();
    auto my_lc_sol=my_lc.get_LCs();
    out<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    out<<"Number of multiplier dof:"<<my_lc_sol.size()<<endl;
    out<<"-----------"<<endl;
    for (uint i=0; i!=my_lc_sol.size(); ++i)
    {
        auto tpa=my_lc_sol[i].first;
        out<<"Line nb "<<i<<", Number of primal dof :"<<tpa.size()<<endl;
        for (uint k=0; k!=tpa.size(); ++k)
        {
            out<<"mult dof nb: "<<i<<", primal dof index: "<<k<<", primal dof nb: "<<get<0>(tpa[k])<<", primal dof coef val: "<<get<1>(tpa[k])<<endl;
        }
        out<<", in the rhs: "<<my_lc_sol[i].second<<endl;
        out<<"-----------"<<endl;
    }
    out<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    //cout<<"la"<<get<0>(my_lc_sol[0].first[0])<<get<1>(my_lc_sol[0].first[0])<<my_lc_sol[0].second<<endl;

    return  0;
}
