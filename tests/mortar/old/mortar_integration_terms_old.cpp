
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
shared_ptr<Mapping<dim,0>> create_geo_0(){
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
shared_ptr<Mapping<dim,0>> create_geo_1(){
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




template <int dim, int dim_field>
void do_test(vector<shared_ptr<Mapping<dim,0>>> & maps, vector<int> degrees, LogStream& out, MultiPatchSpace<PhySpace> & domain_multip){	
	
	using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
	using RefSpaceField_ptr     = shared_ptr<RefSpaceField>;
	using PushFw                = PushForward<Transformation::h_grad, dim,0>;
	//using PhySpace              = PhysicalSpace<RefSpaceField, PushFw>;
	using PhySpace_ptr          = shared_ptr<PhySpace>;	   
	
	vector<RefSpaceField_ptr>	ref_spaces_field;
	vector<PhySpace_ptr>	    spaces;
	 
	//MultiPatchSpace<PhySpace> domain_multip;
	domain_multip.arrangement_open();
	
	
	for(int i=0; i!=maps.size(); ++i){
		ref_spaces_field.push_back(RefSpaceField::create(degrees[i], maps[i]->get_grid()));
		spaces.push_back(PhySpace::create(ref_spaces_field[i], PushFw::create(maps[i]),i));
		spaces[i]->print_info(out);
		spaces[i]->refine_h(2); // Problem here with a value different than 2
		spaces[i]->print_info(out);
		domain_multip.add_patch(spaces[i]);
	}
	 
	domain_multip.arrangement_close();
}




int main()
{
	LogStream out;
	const int dim(2);
	const int dim_field(1);
	vector<int> degrees(1);
	
	using RefSpaceField         = BSplineSpace<dim,dim_field,1>;
	using PushFw                = PushForward<Transformation::h_grad, dim,0>;
	using PhySpace              = PhysicalSpace<RefSpaceField, PushFw>;
	
	// Mapping $\R^2 \rightarrow \R^2$, scalar field
	vector<shared_ptr<Mapping<dim,0>>> maps_d21;
	maps_d21.push_back(create_geo_0<dim>());
	maps_d21.push_back(create_geo_1<dim>());
	MultiPatchSpace<PhySpace> domain_multip;
    do_test<dim,dim_field>(maps_d21, degrees, out, domain_multip);
	
	return  0;
}
