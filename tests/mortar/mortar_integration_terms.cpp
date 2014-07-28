#include "../tests.h"


#include <igatools/base/function_lib.h>
#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/geometry/topology.h>

#include <igatools/basis_functions/space_tools.h>
#include <igatools/basis_functions/bspline_space.h>
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
    using RefSpaceField_ptr     = shared_ptr<RefSpaceField>;
    using PushFw                = PushForward<Transformation::h_grad, dim,0>;
    using PushFwt                = PushForward<Transformation::h_grad, dim-1,0>;
    using PhySpace              = PhysicalSpace<RefSpaceField, PushFw>;
    using PhySpace_ptr          = shared_ptr<PhySpace>;

    using FuncMap = Function<dim, dim>;
    using GradientMap = typename FuncMap::Gradient;

    using FuncMapt = Function<dim-1, dim>;
    using GradientMapt = typename FuncMapt::Gradient;
private:
    boundary_id mortar_id_;
    int degree_multiplier_;

    int slave_patch_nb_;
    int slave_mortar_face_nb_;
    shared_ptr<PhySpace> slave_space_;

    int master_patch_nb_;
    int master_mortar_face_nb_;
    shared_ptr<PhySpace> master_space_;


    vector<Index> face_dof_num_;
    vector<Index> multiplier_face_dof_num_;

    shared_ptr<PhySpace> multiplier_space_;

    using pair_coef_value  = pair<Index, Real>;
    using LinearConstraint = pair<vector<pair_coef_value>, Real>;



    //check_quad_pt_nb() const {
    //}

    vector<LinearConstraint> LCs_;
    void integration();


public:
    /** @name Constructors and destructor */
    ///@{
    /** Default constructor. */
    Mortar_Interface() = delete;

    /** Constructor. */
    Mortar_Interface(const boundary_id &mortar_id, const int &degree_multiplier, const int &slave_patch_nb, const int &slave_mortar_face_nb, shared_ptr<PhySpace> slave_space,
                     const int &master_patch_nb, const int &master_mortar_face_nb, shared_ptr<PhySpace> master_space, shared_ptr<DofsManager> dof_manager)
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
        )
    {
        cout<<"A LC Interface built"<<endl;
        vector<Index> slave_face_dof_num_;
        auto temp=slave_space_->get_face_space(slave_mortar_face_nb_, slave_face_dof_num_);
        slave_face_dof_num_=dof_manager->get_global_dofs(slave_patch_nb_, slave_face_dof_num_);

        vector<Index> master_face_dof_num_;
        auto tempm=master_space_->get_face_space(master_mortar_face_nb_, master_face_dof_num_);
        master_face_dof_num_=dof_manager->get_global_dofs(master_patch_nb_, master_face_dof_num_);

        face_dof_num_=slave_face_dof_num_;
        for (int i=slave_face_dof_num_.size(); i<slave_face_dof_num_.size()+master_face_dof_num_.size(); ++i)
        {
            face_dof_num_.push_back(master_face_dof_num_[i-(slave_face_dof_num_.size())]);
        };

        auto tempmm=multiplier_space_->get_face_space(slave_mortar_face_nb_, multiplier_face_dof_num_);
        multiplier_face_dof_num_=dof_manager->get_global_dofs(slave_patch_nb_, multiplier_face_dof_num_);



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


};




template<int dim, int dim_field>
void Mortar_Interface<dim, dim_field>::integration()
{


    //LogStream out;
    const int int_pts(degree_multiplier_+1);



    // Get the function number living on one element
    const int n_slave_basis      = slave_space_->get_num_basis_per_element();
    const int n_master_basis     = master_space_->get_num_basis_per_element();
    const int n_multiplier_basis = multiplier_space_->get_num_basis_per_element();



    // Get the total function number for each space to build a local matrix
    const int nt_basis0 = slave_space_->get_num_basis();
    const int nt_basis1 = master_space_->get_num_basis();
    const int nt_basism = multiplier_space_->get_num_basis();
    DenseMatrix loc_mat(nt_basism, nt_basis0+nt_basis1);
    loc_mat.clear();



    // Joint grid construction -> parametric joint mesh
    vector< Index >  joint_grid_to_slave;
    vector< Index >  joint_grid_to_master;
    auto slave_grid  = slave_space_->get_grid();
    auto master_grid = master_space_->get_grid();

    vector<shared_ptr<CartesianGrid<dim>>> joint_grid;
    joint_grid.push_back(make_shared<CartesianGrid<dim>>(build_cartesian_grid_union(*slave_grid, *master_grid, joint_grid_to_slave, joint_grid_to_master)));
    joint_grid.push_back(make_shared<CartesianGrid<dim>>(build_cartesian_grid_union(*slave_grid, *master_grid, joint_grid_to_slave, joint_grid_to_master)));

    //joint_grid[0]->print_info(out);
    //joint_grid[1]->print_info(out);

    // Apply the boundary id on the slave and master version of the joint grid
    for (uint j=0; j!=2*dim; ++j)
    {
        joint_grid[0]->set_boundary_id(j,slave_grid->get_boundary_id(j));
        joint_grid[1]->set_boundary_id(j,master_grid->get_boundary_id(j));
    }




    const Quadrature<dim-1> joint_face_quad(QGauss<dim-1>(int_pts+1));
    const int n_qp = joint_face_quad.get_num_points();
    auto  slave_joint_face_quad_gbl   = extend_face_quad(joint_face_quad, slave_mortar_face_nb_);
    auto  master_joint_face_quad_gbl  = extend_face_quad(joint_face_quad, master_mortar_face_nb_);

    auto  slave_quad_pts_unit_domain        = slave_joint_face_quad_gbl.get_points().get_flat_cartesian_product();
    auto  master_quad_pts_unit_domain       = master_joint_face_quad_gbl.get_points().get_flat_cartesian_product();
    auto  w_unit_domain                     = slave_joint_face_quad_gbl.get_weights().get_flat_tensor_product();



    vector<Index>      elem_slave_ids;
    vector<double>     elem_slave_w_meas;
    vector<Index>      elem_master_ids;
    vector<vector<Points<dim>>>  vec_slave_quad_pts_ref_domain;
    vector<vector<Points<dim>>>  vec_master_quad_pts_ref_domain;



    // Loop on the joint grid elements
    auto elem_jm=*joint_grid[1]->begin();
    for (auto elem_js: *joint_grid[0])
    {
        if (elem_js.is_boundary())
        {
            for (Index face_id = 0; face_id < UnitElement<dim>::faces_per_element; ++face_id)
            {
                if (elem_js.is_boundary(face_id))
                {
                    auto elem_js_grid = elem_js.get_grid();
                    if (elem_js_grid->get_boundary_id(face_id)==mortar_id_)
                    {
                        auto elem_js_id     = elem_js.get_flat_index();
                        elem_slave_ids.push_back(joint_grid_to_slave[elem_js_id]);
                        //out<<"SLAVE ID"<<elem_js_id<<joint_grid_to_slave[elem_js_id]<<endl;
                        vec_slave_quad_pts_ref_domain.push_back(elem_js.transform_points_unit_to_reference(slave_quad_pts_unit_domain));
                        elem_slave_w_meas.push_back(static_cast<CartesianGridElement<dim>&>(elem_js).get_measure(FaceTopology<dim>(slave_mortar_face_nb_)));
                    }
                }
                if (elem_jm.is_boundary(face_id))
                {
                    auto elem_jm_grid=elem_jm.get_grid();
                    if (elem_jm_grid->get_boundary_id(face_id)==mortar_id_)
                    {
                        auto elem_jm_id     = elem_jm.get_flat_index();
                        elem_master_ids.push_back(joint_grid_to_master[elem_jm_id]);
                        //out<<"MASTER ID"<<elem_jm_id<<joint_grid_to_master[elem_jm_id]<<endl;
                        vec_master_quad_pts_ref_domain.push_back(elem_jm.transform_points_unit_to_reference(master_quad_pts_unit_domain));
                    }
                }
            } // for face_id
        } // elem_js.is_boundary()
        ++elem_jm;
    } // for elem_js










    // Face joint grid construction
    vector< Index >  joint_face_grid_to_slave;
    vector< Index >  joint_face_grid_to_master;
    auto slave_grid_to_face  = make_shared<map<int,int> >();
    auto master_grid_to_face = make_shared<map<int,int> >();
    auto face_slave_grid     = slave_grid->get_face_grid(slave_mortar_face_nb_,*slave_grid_to_face);
    auto face_master_grid    = master_grid->get_face_grid(master_mortar_face_nb_,*master_grid_to_face);

    vector<shared_ptr<CartesianGrid<dim-1>>> face_joint_grid;
    face_joint_grid.push_back(make_shared<CartesianGrid<dim-1>>(build_cartesian_grid_union(*face_slave_grid, *face_master_grid, joint_face_grid_to_slave,
                                                                joint_face_grid_to_master)));
    face_joint_grid.push_back(make_shared<CartesianGrid<dim-1>>(build_cartesian_grid_union(*face_slave_grid, *face_master_grid, joint_face_grid_to_slave,
                                                                joint_face_grid_to_master)));

    face_joint_grid[0]->print_info(out);
    face_joint_grid[1]->print_info(out);



    vector< Index >     slave_space_to_face_dof;
    vector< Index >     master_space_to_face_dof;
    auto slave_space_face  = slave_space_->get_face_space(slave_mortar_face_nb_, slave_space_to_face_dof);
    auto master_space_face = master_space_->get_face_space(master_mortar_face_nb_, master_space_to_face_dof);

    auto  face_slave_quad_pts_unit_domain        = joint_face_quad.get_points().get_flat_cartesian_product();
    auto  face_master_quad_pts_unit_domain       = joint_face_quad.get_points().get_flat_cartesian_product();
    // no extention to have all the unit coordinates
    //auto  face_slave_quad_pts_unit_domain        = slave_joint_face_quad_gbl.get_points().get_flat_cartesian_product();
    //auto  face_master_quad_pts_unit_domain       = master_joint_face_quad_gbl.get_points().get_flat_cartesian_product();



    // Optimal loop:
    auto felem_jm=*face_joint_grid[1]->begin();
    for (auto felem_js: *face_joint_grid[0])
    {
        auto felem_js_grid = felem_js.get_grid();
        auto felem_js_id     = felem_js.get_flat_index();
        auto face_slave_quad_pts_ref = felem_js.transform_points_unit_to_reference(face_slave_quad_pts_unit_domain);

        auto face_meas=static_cast<CartesianGridElement<dim-1>&>(felem_js).get_measure(FaceTopology<dim-1>(slave_mortar_face_nb_));
        out<<"face_meas"<<face_meas<<endl;


        auto felem_slave  = slave_space_face->get_element(joint_face_grid_to_slave[felem_js_id]);
        auto elem_slave   = slave_space_->get_element(joint_grid_to_slave[slave_space_to_face_dof[felem_js_id]]);
        auto dofs_face_slave    = felem_slave.get_local_to_global();


        vector<GradientMapt> face_slave_map_grad(n_qp);
        //get the slave mapping from the slave physical space
        auto temp_face_slave_map=slave_space_face->get_push_forward()->get_mapping();
        //NON FUNZIONA
        temp_face_slave_map->evaluate_gradients_at_points(face_slave_quad_pts_ref, face_slave_map_grad);


        auto tmp_face_elem_slave          = felem_slave.as_cartesian_grid_element_accessor();
        auto temp_slave_face_pts_unit   = tmp_face_elem_slave.transform_points_reference_to_unit(face_slave_quad_pts_ref);

        //NON FUNZIONA
        auto basis_slave      = felem_slave.evaluate_basis_values_at_points(temp_slave_face_pts_unit);


        ++felem_jm;
    } // for elem_js
    // Pb: it is not possible to get the gradients at points, as the basis function value




    for (int elem_nb = 0; elem_nb!=elem_slave_ids.size(); ++elem_nb)
    {

        auto slave_quad_pts_ref_domain  = vec_slave_quad_pts_ref_domain[elem_nb];
        auto master_quad_pts_ref_domain = vec_master_quad_pts_ref_domain[elem_nb];


        auto elem_slave      = slave_space_->get_element(elem_slave_ids[elem_nb]);
        auto dofs_slave      = elem_slave.get_local_to_global();
        auto elem_multiplier = multiplier_space_->get_element(elem_slave_ids[elem_nb]);
        auto dofs_multiplier = elem_multiplier.get_local_to_global();
        auto elem_master     = master_space_->get_element(elem_master_ids[elem_nb]);
        auto dofs_master     = elem_master.get_local_to_global();

        //out<<"dof_master"<<elem_master_ids[elem_nb]<<dofs_master<<endl;
        //out<<"dof_slave"<<elem_slave_ids[elem_nb]<<dofs_slave<<endl;

        vector<GradientMap> slave_map_grad;
        //get the slave mapping from the slave physical space
        auto temp_slave_map=const_pointer_cast<PushFw>(slave_space_->get_push_forward())->get_mapping();
        temp_slave_map->evaluate_gradients_at_points(slave_quad_pts_ref_domain, slave_map_grad);
        auto w_meas = elem_slave_w_meas[elem_nb];

        auto tmp_elem_slave       = elem_slave.as_cartesian_grid_element_accessor();
        auto temp_slave_pts_unit  = tmp_elem_slave.transform_points_reference_to_unit(slave_quad_pts_ref_domain);
        auto tmp_elem_master      = elem_master.as_cartesian_grid_element_accessor();
        auto temp_master_pts_unit = tmp_elem_master.transform_points_reference_to_unit(master_quad_pts_ref_domain);


        auto basis_slave      = elem_slave.evaluate_basis_values_at_points(temp_slave_pts_unit);
        auto basis_multiplier = elem_multiplier.evaluate_basis_values_at_points(temp_slave_pts_unit);
        auto basis_master     = elem_master.evaluate_basis_values_at_points(temp_master_pts_unit);

        out<<"ss"<<basis_slave.size()<<"ms"<<basis_multiplier.size()<<"ms"<<basis_master.size()<<endl;
        out<<"ns"<<n_slave_basis<<"nm"<<n_multiplier_basis<<"nm"<<n_master_basis<<endl;

        for (int i = 0; i < n_multiplier_basis; ++i)
        {
            auto phi_m = basis_multiplier.get_function_view(i);
            for (int j = 0; j < n_slave_basis; ++j)
            {
                auto phj = basis_slave.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                {
                    loc_mat(dofs_multiplier[i],dofs_slave[j])=loc_mat(dofs_multiplier[i],dofs_slave[j])
                                                              +scalar_product(phi_m[qp],phj[qp])*w_meas*w_unit_domain[qp]*determinant<dim,dim>(slave_map_grad[qp]);
                }
            }
            for (int j = 0; j < n_master_basis; ++j)
            {
                auto phj = basis_master.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                {
                    loc_mat(dofs_multiplier[i],dofs_master[j])=loc_mat(dofs_multiplier[i],dofs_master[j])
                                                               +scalar_product(phi_m[qp],phj[qp])*w_meas*w_unit_domain[qp]*determinant<dim,dim>(slave_map_grad[qp]);
                }
            }
        }// for i
    } // for nb_elem



    for (int i=0; i<multiplier_face_dof_num_.size(); ++i)
    {
        vector<pair_coef_value> ttp;
        for (int k=0; k!=face_dof_num_.size(); ++k)
        {
            ttp.push_back(pair<Index, Real>(face_dof_num_[k],loc_mat(multiplier_face_dof_num_[i],face_dof_num_[k])));
        }
        LCs_.push_back(pair<vector<pair_coef_value>, Real>(ttp,0.0));
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
    vector<int> degrees(2,1);
    const int deg(1);


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
    domain_multip.arrangement_open();
    ////////
    const boundary_id mortar_id = 7;
    ////////
    vector<set<boundary_id>> mortar_faces;
    mortar_faces.push_back({1});
    mortar_faces.push_back({0});


    for (uint i=0; i!=maps_d21.size(); ++i)
    {
        ref_spaces_field.push_back(RefSpaceField::create(degrees[i], maps_d21[i]->get_grid()));
        spaces.push_back(PhySpace::create(ref_spaces_field[i], PushFw::create(maps_d21[i]),i));
        //      spaces[i]->print_info(out);
        if (i==1)
        {
            spaces[i]->refine_h(2);
        }// Problem here with a value different than 2
        if (i==0)
        {
            spaces[i]->get_grid()->set_boundary_id(2,3);
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

    domain_multip.arrangement_close();
    auto dof_manager=domain_multip.get_dofs_manager();




    auto my_lc=Mortar_Interface<dim,dim_field>(mortar_id, deg, 0, 1, spaces[0], 1, 0, spaces[1], dof_manager);
    auto my_lc_sol=my_lc.get_LCs();
    out<<"Number of multiplier dof"<<my_lc_sol.size()<<endl;
    for (uint i=0; i!=my_lc_sol.size(); ++i)
    {
        auto tpa=my_lc_sol[i].first;
        out<<"Number of primal dof"<<tpa.size()<<endl;
        for (uint k=0; k!=tpa.size(); ++k)
        {
            out<<"mult dof nb: "<<i<<", primal dof index: "<<k<<", primal dof nb: "<<get<0>(tpa[k])<<", primal dof coef val: "<<get<1>(tpa[k])<<endl;
        }
        out<<", in the rhs: "<<my_lc_sol[i].second<<endl;
        out<<"-----------"<<endl;
    }
    cout<<"la"<<get<0>(my_lc_sol[0].first[0])<<get<1>(my_lc_sol[0].first[0])<<my_lc_sol[0].second<<endl;
    return  0;
}
