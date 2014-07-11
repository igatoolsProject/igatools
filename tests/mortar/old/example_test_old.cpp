
#include <igatools/base/function_lib.h>
#include <igatools/geometry/mapping_lib.h>
#include <igatools/geometry/identity_mapping.h>
#include <igatools/geometry/ig_mapping.h>
//#include <igatools/io/reader.h>
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
using functions::ConstantFunction;
using space_tools::project_boundary_values;
using dof_tools::apply_boundary_values;
using dof_tools::get_sparsity_pattern;



template<const int dim>
shared_ptr<Mapping<dim>> Create_geo0(){
	using RefSpace_map         = BSplineSpace<dim,dim,1>;
	shared_ptr<RefSpace_map>   ref_space_map;
	
	const int p = 1;
	auto knots     = CartesianGrid<dim>::create(2);
	ref_space_map  = RefSpace_map::create(p,knots);
	vector<Real>   control_pts(ref_space_map->get_num_basis());
		
	if (dim == 1)
		{
			int id = 0 ;
			control_pts[id++] = 0.0 ;
			control_pts[id++] = 1.0 ;
		}
	else if (dim == 2)
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
	auto map = IgMapping<RefSpace_map>::create(ref_space_map, control_pts);
	//auto map=IdentityMapping<dim,0>::create(knots);
	return map;
		
} 
	
	







template<const int dim>
shared_ptr<Mapping<dim>> Create_geo1(){
	using RefSpace_map         = BSplineSpace<dim,dim,1>;
	shared_ptr<RefSpace_map>   ref_space_map;
	
	const int p = 1;
	auto knots     = CartesianGrid<dim>::create(2);
	ref_space_map  = RefSpace_map::create(p,knots);
	vector<Real>   control_pts(ref_space_map->get_num_basis());

		
	if (dim == 1)
		{
			int id = 0 ;
			control_pts[id++] = 1.0 ;
			control_pts[id++] = 2.0 ;
		}
	else if (dim == 2)
		{
			int id = 0 ;
			control_pts[id++] = 1.0 ;
			control_pts[id++] = 2.0 ;
			
			control_pts[id++] = 1.0 ;
			control_pts[id++] = 2.0 ;
			
			control_pts[id++] = 0.0 ;
			control_pts[id++] = 0.0 ;
			
			control_pts[id++] = 1.0 ;
			control_pts[id++] = 1.0 ;
		}
	else if (dim == 3)
		{
			int id = 0 ;
			control_pts[id++] = 1.0 ;
			control_pts[id++] = 2.0 ;
			
			control_pts[id++] = 1.0 ;
			control_pts[id++] = 2.0 ;
			
			control_pts[id++] = 1.0 ;
			control_pts[id++] = 2.0 ;
		
			control_pts[id++] = 1.0 ;
			control_pts[id++] = 2.0 ;
			
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
	auto map = IgMapping<RefSpace_map>::create(ref_space_map, control_pts);
	//map->get_grid()->refine();
	//auto map=IdentityMapping<dim,0>::create(knots);
	
	return map;
}



//template<const int dim>
//shared_ptr<Mapping<dim>> nurb_geometry_from_file0()
//{
//    string input_file = "nurb_geometry0.xml";
//    auto map = get_mapping_from_file<dim>(input_file);
//	
//    const int n_plot_points = 10;
//    Writer<dim> writer(map, n_plot_points);
//    string filename = "view_nurb_geometry0";
//    writer.save(filename);
//	return map;
//}

//template<const int dim>
//shared_ptr<Mapping<dim>> nurb_geometry_from_file1()
//{
//    string input_file = "nurb_geometry1.xml";
//    auto map = get_mapping_from_file<dim>(input_file);
//	
//   const int n_plot_points = 10;
//    Writer<dim> writer(map, n_plot_points);
//    string filename = "view_nurb_geometry1";
//    writer.save(filename);
//	return map;
//}
//



// Dirichlet_function
template <int dim, int dim_field, int rank=1>
class D_function : public Function<dim,dim_field,rank>
{
public:
	D_function(){};
	
	void evaluate(const std::vector< typename Function<dim,dim_field,rank>::Point> &points, 
				  std::vector< typename Function<dim, dim_field,rank>::Value> &values) const
	
	{
		const int num_points = points.size() ;
        Assert(num_points == values.size(), ExcDimensionMismatch(num_points, values.size())) ;
		
		for (int i =0; i<num_points; ++i){
			//auto pt_i = points[i];
			values[i]= 1.0;//sin(5*pt_i[0])*cos(6*pt_i[1]);
		}

	}
};	




int main()
{
	const int dim(2);
	const int dim_field(1);
	auto map0=Create_geo0<dim>();
	auto map1=Create_geo1<dim>();
	//auto map0=nurb_geometry_from_file0<dim>();
	//auto map1=nurb_geometry_from_file1<dim>();
	
	using RefSpace_field         = BSplineSpace<dim,dim_field,1>;
	//using RefSpace_field         = NURBSSpace<dim,dim_field,1>;
	shared_ptr<RefSpace_field>   ref_space_field0;
	shared_ptr<RefSpace_field>   ref_space_field1;
	const int deg0 = 1;
	const int deg1 = 1;
	ref_space_field0  = RefSpace_field::create(deg0, map0->get_grid());
	ref_space_field1  = RefSpace_field::create(deg1, map1->get_grid());
	
    using PushFw      = PushForward<Transformation::h_grad, dim,0>;
    using PhySpace    = PhysicalSpace<RefSpace_field, PushFw>;
    auto space0     = PhySpace::create(ref_space_field0, PushFw::create(map0),0);
	auto space1     = PhySpace::create(ref_space_field1, PushFw::create(map1),1);
	
	//map0->get_grid()->refine();
	//space0->refine_h(2);
	//space1->refine_h(2);

	
	MultiPatchSpace<PhySpace> domain_multip;
	domain_multip.arrangement_open();
	domain_multip.add_patch(space0);

	domain_multip.add_patch(space1);
	domain_multip.add_interface(InterfaceType::C0_strong,
									space0,1,
									space1,0);
	domain_multip.arrangement_close();
	LogStream out;
	domain_multip.print_info(out);
	auto dof_manager=domain_multip.get_dofs_manager();
	dof_manager->print_info(out);
	//	auto num_lin_c=dof_manager.get_num_linear_constraints();
	
	
	/////const boundary_id nh_id = 0;
	////const boundary_id n_id = 1;
	////const boundary_id dh_id = 2;
	const boundary_id d_id = 3;
	/////set<boundary_id> nh_faces = {1};
	set<boundary_id> d_faces    = {0};
	//// Change the face id 0->NH // 1->NNH // 2->DH // 3->DNH
	/////auto crr_nh=nh_faces.begin(), end_nh=nh_faces.end();
	/////for (; crr_nh!=end_nh; ++crr_nh)
	/////{
	/////	space0->get_grid()->set_boundary_id(*crr_nh,nh_id);
	/////}
	
	auto crr_d=d_faces.begin(), end_d=d_faces.end();
	for (; crr_d!=end_d; ++crr_d)
	{
		space0->get_grid()->set_boundary_id(*crr_d,d_id);
	}

	
	const Quadrature<dim>   elem_quad0(QGauss<dim>(deg0+1));
    const Quadrature<dim-1> face_quad0(QGauss<dim-1>(deg0+1));
	std::shared_ptr<Matrix<LAPack::trilinos>> matrix;
    std::shared_ptr<Vector<LAPack::trilinos>> rhs;
    std::shared_ptr<Vector<LAPack::trilinos>> solution;
	const auto n_basis = space0->get_num_basis();
    matrix   = Matrix<LAPack::trilinos>::create(get_sparsity_pattern<PhySpace>(space0));
    rhs      = Vector<LAPack::trilinos>::create(n_basis);
    solution = Vector<LAPack::trilinos>::create(n_basis);
	DenseMatrix loc_mat(n_basis, n_basis);
    DenseVector loc_rhs(n_basis);
    vector<Index> loc_dofs(n_basis);
	
    const int n_qp = elem_quad0.get_num_points();
    ConstantFunction<dim> f({0.0});
	using Value = typename Function<dim>::Value;
    vector<Value> f_values(n_qp);
	
    auto elem = space0->begin();
    const auto elem_end = space0->end();
    ValueFlags fill_flags = ValueFlags::value | ValueFlags::gradient |
	ValueFlags::w_measure | ValueFlags::point;
    elem->init_values(fill_flags, elem_quad0);
	
    for (; elem != elem_end; ++elem)
    {
        loc_mat.clear();
        loc_rhs.clear();
        elem->fill_values();
		
        auto points  = elem->get_points();
        auto phi     = elem->get_basis_values();
        auto grd_phi = elem->get_basis_gradients();
        auto w_meas  = elem->get_w_measures();
		
        f.evaluate(points, f_values);
		
        for (int i = 0; i < n_basis; ++i)
        {
            auto grd_phi_i = grd_phi.get_function_view(i);
            for (int j = 0; j < n_basis; ++j)
            {
                auto grd_phi_j = grd_phi.get_function_view(j);
                for (int qp = 0; qp < n_qp; ++qp)
                    loc_mat(i,j) +=
					scalar_product(grd_phi_i[qp], grd_phi_j[qp])
					* w_meas[qp];
            }
			
            auto phi_i = phi.get_function_view(i);
            for (int qp = 0; qp < n_qp; ++qp)
                loc_rhs(i) += scalar_product(phi_i[qp], f_values[qp])
				* w_meas[qp];
        }
		
        loc_dofs = elem->get_local_to_global();
        matrix->add_block(loc_dofs, loc_dofs, loc_mat);
        rhs->add_block(loc_dofs, loc_rhs);
    }
	
    matrix->fill_complete();
	
    //ConstantFunction<dim> g({0.0});
	D_function<dim, dim_field, 1> dfunction;
	//vector< Value > dfunc_at_points(points.size());
	//function.evaluate(points, d_func_at_points);
    std::map<Index, Real> values;
    //project_boundary_values<PhySpace,LinearAlgebraPackage::trilinos>(g,space0,face_quad0,dir_id,values);
	project_boundary_values<PhySpace,LAPack::trilinos>(dfunction,space0,face_quad0,d_id,values);
    apply_boundary_values(values, *matrix, *rhs, *solution);
	solution->print(out);

	
    using LinSolver = LinearSolver<LAPack::trilinos>;
    LinSolver solver(LinSolver::SolverType::CG);
    solver.solve(*matrix, *rhs, *solution);
	solution->print(out);

	const int n_plot_points = 2;
    Writer<dim> writer(map0, n_plot_points);
	
    writer.add_field(space0, *solution, "solution");
    string filename = "poisson_problem_P0" ;
    writer.save(filename);
	
     
	
		out<<"Essai";
	
	
	
	
    return  0;
}
