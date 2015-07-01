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
/*
 *  Test for the DofDistribution class and its serialization
 *  author: pauletti, martinelli
 *  date: May 06th, 2015
 *
 */

// TODO (pauletti, Dec 26, 2014): make test dimension independent

#include "../tests.h"
#include <igatools/basis_functions/dof_distribution.h>



template <int dim>
void serialize_deserialize(const DofDistribution<dim> &dof_admin)
{
    out.begin_item("Original DofDistribution.");
    dof_admin.print_info(out);
    out.end_item();



    std::string filename = "dof_distribution_dim" + std::to_string(dim) + ".xml";
    std::string tag_name = "DofDistribution_dim" + std::to_string(dim);
    {
        // serialize the DofDistribution object to an xml file
        std::ofstream xml_ostream(filename);
        OArchive xml_out(xml_ostream);

        xml_out << boost::serialization::make_nvp(tag_name.c_str(),dof_admin);
        xml_ostream.close();
    }

    DofDistribution<dim> dof_admin_new;
    {
        // de-serialize the DofDistribution object from an xml file
        std::ifstream xml_istream(filename);
        IArchive xml_in(xml_istream);
        xml_in >> BOOST_SERIALIZATION_NVP(dof_admin_new);
        xml_istream.close();
    }
    out.begin_item("DofDistribution after serialize-deserialize.");
    dof_admin_new.print_info(out);
    out.end_item();
    //*/
}


template <int dim>
void test1()
{
    OUTSTART
    using SplineSpace = SplineSpace<dim>;
    using MultiplicityTable = typename SplineSpace::MultiplicityTable;

    typename SplineSpace::DegreeTable deg {{2}};

    auto grid = CartesianGrid<dim>::create(4);

    auto int_mult = MultiplicityTable({ {{1,3}} });
    auto sp_spec = SplineSpace::create(deg, grid, int_mult);

    CartesianProductArray<Real,2> bn_x {{-0.5, 0, 0}, {1.1, 1.2, 1.3}};
    typename SplineSpace::BoundaryKnotsTable bdry_knots { {bn_x} };
    typename SplineSpace::EndBehaviourTable end_b((typename SplineSpace::EndBehaviourTable(SafeSTLArray<BasisEndBehaviour,dim>(BasisEndBehaviour::end_knots))));

    auto rep_knots = sp_spec->compute_knots_with_repetition(end_b,bdry_knots);

    auto n_basis = sp_spec->get_num_basis_table();
    auto degree = sp_spec->get_degree_table();

    DofDistribution<dim> dof_admin(n_basis, degree, sp_spec->get_periodic_table());

    //-----------------------------------------------------------------
    const auto &dofs_view = dof_admin.get_dofs_view();
    const std::string property_active = "active";
    //dof_admin.add_dofs_property(property_active);

    for (const auto &dof : dofs_view)
    {
        if (dof % 2 == 0)
            dof_admin.set_dof_property_status(property_active, dof, true);
        else
            dof_admin.set_dof_property_status(property_active, dof, false);
    }
    //-----------------------------------------------------------------


    serialize_deserialize(dof_admin);

    OUTEND
}

template <int dim>
void test2()
{
    OUTSTART
    using SplineSpace = SplineSpace<dim>;

    typename SplineSpace::DegreeTable deg {{1,2}};

    auto grid = CartesianGrid<dim>::create({4,3});
    auto int_mult = SplineSpace::get_multiplicity_from_regularity(InteriorReg::maximum,
                    deg, grid->get_num_intervals());
    auto sp_spec = SplineSpace::create(deg, grid, int_mult);

    typename SplineSpace::EndBehaviourTable end_b((typename SplineSpace::EndBehaviourTable(SafeSTLArray<BasisEndBehaviour,dim>(BasisEndBehaviour::interpolatory))));

    auto rep_knots = sp_spec->compute_knots_with_repetition(end_b);


    auto n_basis = sp_spec->get_num_basis_table();
    auto degree = sp_spec->get_degree_table();

    DofDistribution<dim> basis_index(n_basis, degree, sp_spec->get_periodic_table());

    //-----------------------------------------------------------------
    const auto &dofs_view = basis_index.get_dofs_view();
    const std::string property_active = "active";
    // basis_index.add_dofs_property(property_active);

    for (const auto &dof : dofs_view)
    {
        if (dof % 2 == 0)
            basis_index.set_dof_property_status(property_active, dof, true);
        else
            basis_index.set_dof_property_status(property_active, dof, false);
    }
    //-----------------------------------------------------------------

    serialize_deserialize(basis_index);

    OUTEND
}

int main()
{
    out.depth_console(10);
    test1<1>();
    test2<2>();

    return 0;
}
