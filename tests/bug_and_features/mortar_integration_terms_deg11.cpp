#include "../tests.h"


#include <igatools/base/function_lib.h>
#include <igatools/base/linear_constraint.h>

#include <igatools/geometry/mapping.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/ig_mapping.h>
#include <igatools/geometry/cartesian_grid.h>
#include <igatools/geometry/cartesian_grid_element_accessor.h>
#include <igatools/geometry/grid_tools.h>
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
//using namespace std;


using std::shared_ptr;
using std::pair;
using std::map;
using std::set;
using std::cout;
using std::endl;

using iga::Vector;

//--------------------------------------------------------------------
template<int dim, int dim_field>
class Mortar_Interface
{

public:
    using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
    using RefSpaceFieldt        = BSplineSpace<dim-1,dim_field,1>;
    using RefSpaceField_ptr     = shared_ptr<RefSpaceField>;
    using PushFw                = PushForward<Transformation::h_grad, dim,0>;
    using PushFwt               = PushForward<Transformation::h_grad, dim-1,1>;
    using PhySpace              = PhysicalSpace<RefSpaceField, PushFw>;
    using PhySpacet             = PhysicalSpace<RefSpaceFieldt, PushFwt>;
    using PhySpace_ptr          = shared_ptr<PhySpace>;

    using FuncMap = Function<dim, dim>;
    using GradientMap = typename FuncMap::Gradient;

    using FuncMapt = Function<dim-1, dim>;
    using GradientMapt = typename FuncMapt::Gradient;
private:

    boundary_id mortar_id_;
    int degree_multiplier_;

    int slave_mortar_face_nb_;
    shared_ptr<PhySpace> slave_space_;

    int master_mortar_face_nb_;
    shared_ptr<PhySpace> master_space_;

    shared_ptr<PhySpacet> slave_face_space_;
    iga::vector<Index> slave_face_dof_num_;
    shared_ptr<PhySpacet> master_face_space_;
    iga::vector<Index> master_face_dof_num_;
    iga::vector<Index> face_dof_num_;
    shared_ptr<PhySpacet> multiplier_face_space_;
    iga::vector<Index> multiplier_face_dof_num_;

    shared_ptr<PhySpace> multiplier_space_;

    using pair_coef_value  = std::pair<Index, Real>;



    //check_quad_pt_nb() const {
    //}

    iga::vector<shared_ptr<LinearConstraint>> LCs_;

    void integration();


public:
    /** @name Constructors and destructor */
    ///@{
    /** Default constructor. */
    Mortar_Interface() = delete;

    /** Constructor. */
    Mortar_Interface(const boundary_id &mortar_id,
                     const int &degree_multiplier,
                     const int &slave_mortar_face_nb,
                     shared_ptr<PhySpace> slave_space,
                     const int &master_mortar_face_nb,
                     shared_ptr<PhySpace> master_space,
                     shared_ptr<SpaceManager> space_manager)
        :
        mortar_id_(mortar_id),
        degree_multiplier_(degree_multiplier),
        slave_mortar_face_nb_(slave_mortar_face_nb),
        slave_space_(slave_space),
        master_mortar_face_nb_(master_mortar_face_nb),
        master_space_(master_space),
        multiplier_space_(
            PhySpace::create(
                RefSpaceField::create(degree_multiplier_,slave_space_->get_grid()),
                std::const_pointer_cast<PushFw>(slave_space_->get_push_forward())))
    {
        cout<<"A LC Interface built"<<endl;
        slave_face_space_=slave_space_->get_face_space(slave_mortar_face_nb_, slave_face_dof_num_);
        //slave_face_dof_num_=space_manager->get_global_dofs(slave_space_->get_id(), slave_face_dof_num_);


        master_face_space_=master_space_->get_face_space(master_mortar_face_nb_, master_face_dof_num_);
        //master_face_dof_num_=space_manager->get_global_dofs(master_space_->get_id(), master_face_dof_num_);

        face_dof_num_=slave_face_dof_num_;
        for (int i=slave_face_dof_num_.size(); i<slave_face_dof_num_.size()+master_face_dof_num_.size(); ++i)
        {
            face_dof_num_.push_back(master_face_dof_num_[i-(slave_face_dof_num_.size())]);
        }

        multiplier_face_space_=multiplier_space_->get_face_space(slave_mortar_face_nb_, multiplier_face_dof_num_);
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

    const iga::vector<shared_ptr<LinearConstraint>> get_LCs()
    {
        this->integration();
        return LCs_;
    }


};




template<int dim, int dim_field>
void Mortar_Interface<dim, dim_field>::integration()
{
    // Face quadrature
    const int int_pts(degree_multiplier_+1);
    const Quadrature<dim-1> joint_face_quad(QGauss<dim-1>(int_pts+1));
    const int n_qp                         = joint_face_quad.get_num_points();
    auto  face_w_unit_domain               = joint_face_quad.get_weights().get_flat_tensor_product();
    auto  face_quad_pts_unit_domain        = joint_face_quad.get_points().get_flat_cartesian_product();




    // Face joint grid construction -> parametric joint mesh
    auto slave_grid  = slave_space_->get_grid();
    auto master_grid = master_space_->get_grid();

    using FaceGridMap = typename CartesianGrid<dim>::FaceGridMap;
    auto slave_grid_to_face  = make_shared<FaceGridMap>();
    auto master_grid_to_face = make_shared<FaceGridMap>();
    auto face_slave_grid     = slave_grid->get_face_grid(slave_mortar_face_nb_, *slave_grid_to_face);
    auto face_master_grid    = master_grid->get_face_grid(master_mortar_face_nb_, *master_grid_to_face);

    using FaceGrid = CartesianGrid<dim-1>;
    shared_ptr<FaceGrid> face_joint_grid;
    typename grid_tools::InterGridMap<dim-1> joint_face_grid_to_joint_grid_slave;
    typename grid_tools::InterGridMap<dim-1> joint_face_grid_to_joint_grid_master;
    face_joint_grid = grid_tools::build_cartesian_grid_union(
                          *face_slave_grid,
                          *face_master_grid,
                          joint_face_grid_to_joint_grid_slave,
                          joint_face_grid_to_joint_grid_master);

    out<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    out<<"Face joint grid info"<<endl;
    face_joint_grid->print_info(out);
    out<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;


    // Local mortar matrix
    DenseMatrix loc_e(multiplier_face_dof_num_.size(),face_dof_num_.size());
    loc_e.clear();




    // Loop on the face joint grid to lead the integration on the slave side
    auto felem_j     = face_joint_grid->begin();
    auto felem_j_end = face_joint_grid->end();

//    for (auto felem_j: *face_joint_grid)
    for (; felem_j != felem_j_end ; ++felem_j)
    {
        auto felem_j_grid = felem_j->get_grid();
        auto felem_j_id   = felem_j->get_flat_index();

        auto face_quad_pts_ref = felem_j->transform_points_unit_to_reference(face_quad_pts_unit_domain);


        out<<"face joint grid element nb:"<<felem_j_id<<endl;
        auto face_meas=static_cast<CartesianGridElement<dim-1>&>(*felem_j).get_measure();
        out<<"face_meas"<<face_meas<<endl;

        auto felem_slave     = slave_face_space_->get_element(
                                   joint_face_grid_to_joint_grid_slave[felem_j]->get_flat_index());
        auto dofs_face_slave = felem_slave.get_local_to_global();
        auto n_slave_basis   = dofs_face_slave.size();

        auto felem_master     = master_face_space_->get_element(
                                    joint_face_grid_to_joint_grid_master[felem_j]->get_flat_index());
        auto dofs_face_master = felem_master.get_local_to_global();
        auto n_master_basis   = dofs_face_master.size();

        auto felem_multiplier     = multiplier_face_space_->get_element(
                                        joint_face_grid_to_joint_grid_slave[felem_j]->get_flat_index());
        auto dofs_face_multiplier = felem_multiplier.get_local_to_global();
        auto n_multiplier_basis   = dofs_face_multiplier.size();

        out<<"All slave dofs: ";
        slave_face_dof_num_.print_info(out);
        out<<", current face slave dofs";
        dofs_face_slave.print_info(out);
        out<<endl;

        out<<"All master dofs: ";
        master_face_dof_num_.print_info(out);
        out<<", current face master dofs";
        dofs_face_master.print_info(out);
        out<<endl;

        out<<"All multiplier dofs: ";
        multiplier_face_dof_num_.print_info(out);
        out<<", current face multiplier dofs";
        dofs_face_multiplier.print_info(out);
        out<<endl;


        //get the slave mapping from the slave physical space
        ValueVector<GradientMapt> face_slave_map_grad(n_qp);
        auto temp_face_slave_map = slave_face_space_->get_push_forward()->get_mapping();
        temp_face_slave_map->evaluate_gradients_at_points(face_quad_pts_ref, face_slave_map_grad);


        // get the location of the current points on each unit domain
        auto cgea_felem_slave               = felem_slave.as_cartesian_grid_element_accessor();
        auto curr_slave_face_quad_pts_unit  = cgea_felem_slave.transform_points_reference_to_unit(face_quad_pts_ref);

        auto cgea_felem_master              = felem_master.as_cartesian_grid_element_accessor();
        auto curr_master_face_quad_pts_unit = cgea_felem_master.transform_points_reference_to_unit(face_quad_pts_ref);



        // get the function values at the current unit points
        auto basis_slave      = felem_slave.evaluate_basis_values_at_points(curr_slave_face_quad_pts_unit);
        auto basis_multiplier = felem_multiplier.evaluate_basis_values_at_points(curr_slave_face_quad_pts_unit);
        auto basis_master     = felem_master.evaluate_basis_values_at_points(curr_master_face_quad_pts_unit);


        //add the contributions of the current element
        for (int i = 0; i < n_multiplier_basis; ++i)
        {
            auto phi_m = basis_multiplier.get_function_view(i);
            for (int j = 0; j < n_slave_basis; ++j)
            {
                auto phj = basis_slave.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                {
                    loc_e(dofs_face_multiplier[i], dofs_face_slave[j]) =
                        loc_e(dofs_face_multiplier[i],dofs_face_slave[j])+ scalar_product(phi_m[qp],phj[qp]) *
                        face_meas * face_w_unit_domain[qp] * determinant<dim-1,dim>(face_slave_map_grad[qp]);
                }
            }
            for (int j = 0; j < n_master_basis; ++j)
            {
                auto phj = basis_master.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                {
                    loc_e(dofs_face_multiplier[i],dofs_face_master[j]+slave_face_dof_num_.size())=loc_e(dofs_face_multiplier[i],dofs_face_master[j]+slave_face_dof_num_.size())+
                            + scalar_product(phi_m[qp],phj[qp])*face_meas*face_w_unit_domain[qp]*determinant<dim-1,dim>(face_slave_map_grad[qp]);
                }
            }// for j
        }// for i
        out<<"--------------------"<<endl;
    } // for elem_j



    // fill the linear constraint data LCs_
    // as LCs[curr multiplier dof nb].first[curr global dof nb].first=global dof nb (i.e. slave or master dof nb)
    // as LCs[curr multiplier dof nb].first[curr global dof nb].second=coefficient related to the corresponding global dof nb
    // as LCs[curr multiplier dof nb].second=value in the rhs
    const auto n_dofs_multipliers = multiplier_face_dof_num_.size();
    const auto n_dofs_face        = face_dof_num_.size();
    for (int i = 0 ; i < n_dofs_multipliers ; ++i)
    {
        vector<Index> primal_dofs_id;
        vector<Real> coefficients;

        for (int k = 0 ; k != n_dofs_face; ++k)
        {
            primal_dofs_id.push_back(face_dof_num_[k]);
            coefficients.push_back(loc_e(i,k));
        }
        LCs_.push_back(shared_ptr<LinearConstraint>(new LinearConstraint(primal_dofs_id,coefficients,0.0)));
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
    domain_multip.patch_insertion_open();
    ////////
    const boundary_id mortar_id = 7;
    ////////
    vector<set<boundary_id>> mortar_faces;
    mortar_faces.push_back( {1});
    mortar_faces.push_back( {0});


    for (int i=0; i!=maps_d21.size(); ++i)
    {
        ref_spaces_field.push_back(RefSpaceField::create(degrees[i], maps_d21[i]->get_grid()));

        spaces.push_back(
            PhySpace::create(ref_spaces_field[i],
                             PushFw::create(maps_d21[i])));
        //      spaces[i]->print_info(out);
        if (i==1)
        {
            spaces[i]->refine_h(2);
        }// Problem here with a value different than 2
        else if (i==0)
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

    domain_multip.patch_insertion_close();
    auto space_manager = domain_multip.get_space_manager();




    auto my_lc = Mortar_Interface<dim,dim_field>(mortar_id,deg,1,spaces[0],0,spaces[1], space_manager);
    //my_lc.integration();
    auto my_lc_sol = my_lc.get_LCs();

    space_manager->linear_constraints_open();

    out<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    out<<"Number of multiplier dof:"<<my_lc_sol.size()<<endl;
    out<<"-----------"<<endl;
    int lc_counter = 0;
    for (const auto &lc : my_lc_sol)
    {
        out << "Line nb " << lc_counter << ", Number of primal dof :" << lc->get_num_lhs_terms() << endl;
        lc->print_info(out);
        out<<"-----------"<<endl;


        space_manager->add_linear_constraint(lc);



        ++lc_counter;
    }
    out<<"+++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

    space_manager->linear_constraints_close();

    space_manager->print_info(out);

    //cout<<"la"<<get<0>(my_lc_sol[0].first[0])<<get<1>(my_lc_sol[0].first[0])<<my_lc_sol[0].second<<endl;
    return  0;
}
