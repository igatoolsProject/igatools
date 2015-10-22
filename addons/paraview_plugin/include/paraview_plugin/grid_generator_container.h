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

#ifndef VTK_GRID_GENERATOR_CONTAINER_H_
#define VTK_GRID_GENERATOR_CONTAINER_H_

#include <igatools/base/config.h>

#include <igatools/base/tuple_utils.h>
#include <paraview_plugin/grid_generator.h>

class vtkMultiBlockDataSet;

IGA_NAMESPACE_OPEN

class FunctionsContainer;


class VtkIgaGridGeneratorContBase
{
private:

  /**
   * Self type.
   */
  typedef VtkIgaGridGeneratorContBase Self_;

  /**
   * Self shared pointer type.
   */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /**
   * Grid generator type for @p dim and @p codim.
   */
  template <int dim, int codim>
  using GridGen_ = VtkIgaGridGenerator<dim, codim>;

  /**
   * Shared pointer type of the grid generator.
   */
  template <int dim, int codim>
  using GridGenPtr_ = std::shared_ptr<GridGen_<dim, codim>>;

  /**
   * Type for the map containing the grid generators.
   */
  template <int dim, int codim>
  using GenMap_ = std::map<Index, GridGenPtr_<dim, codim>>;

  /**
   * Function map type.
   */
  template <int dim, int codim>
  using MapFun_ = Function<dim, 0, dim + codim, 1>;

protected:
  /**
   * Function map shared pointer type.
   */
  template <int dim, int codim>
  using MapFunPtr_ = std::shared_ptr<MapFun_<dim, codim>>;


  /**
   * Grid information shared pointer type.
   */
  typedef std::shared_ptr<VtkGridInformation> GridInfoPtr_;

  /**
   * Functions container shared pointer type.
   */
  typedef std::shared_ptr<FunctionsContainer> FunContPtr_;


  /**
   * Constructor.
   */
  VtkIgaGridGeneratorContBase(const FunContPtr_ funcs_container,
                              const GridInfoPtr_ solid_info,
                              const GridInfoPtr_ knot_info);

private:

  template <int dim>
  class VtkGridGeneratorContBaseSameDim
  {

  public:
    template <int codim>
    class VtkGridGeneratorContBaseSameDimCodim
    {


    public:

      const GenMap_<dim, codim> &get_generators() const;

      GenMap_<dim, codim> &get_generators();

    private:
      GenMap_<dim, codim> grid_generators_;

    };

  public:
    /**
     * Returns a const-reference to the data identified by the index pair <tt><dim,codim></tt>.
     *
     * @note The returned container holds the geometry parametrizations
     * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}}\f$.
     */
    template <int codim>
    const VtkGridGeneratorContBaseSameDimCodim<codim> &get_data_codim() const
    {
      return boost::fusion::at_key<Topology<codim>>(data_varying_codim_);
    };

    /**
     * Returns a const-reference to the data identified by the index pair <tt><dim,codim></tt>.
     *
     * @note The returned container holds the geometry parametrizations
     * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}}\f$.
     */
    template <int codim>
    VtkGridGeneratorContBaseSameDimCodim<codim> &get_data_codim()
    {
      return boost::fusion::at_key<Topology<codim>>(data_varying_codim_);
    };

    /**
     * Returns a const-reference to the internal data.
     *
     * @note The returned object refers to different values of the <tt>codim</tt> parameter.
     */
    const auto &get_data() const
    {
      return data_varying_codim_;
    }



    /**
     * All the data in the VtkGridGeneratorContBaseSameDim class, organized by
     * the @p codim index (starting from 0 to @p dim - 1)
     */
    DataVaryingId<VtkGridGeneratorContBaseSameDimCodim, 0, 4-dim> data_varying_codim_;

  };


protected:

  /**
   * Same dim and codim Vtk grid generator.
   */
  template <int dim, int codim>
  using VtkGridGenContSameCodim_ =
    typename VtkGridGeneratorContBaseSameDim<dim>::template VtkGridGeneratorContBaseSameDimCodim<codim>;


  /**
   * Returns a const-reference to the data identified by the index @p dim.
   */
  template <int dim>
  const VtkGridGeneratorContBaseSameDim<dim> &get_data_dim() const
  {
    return boost::fusion::at_key<Topology<dim>>(data_varying_dim_);
  };

  /**
   * Returns a reference to the data identified by the index @p dim.
   */
  template <int dim>
  VtkGridGeneratorContBaseSameDim<dim> &get_data_dim()
  {
    return boost::fusion::at_key<Topology<dim>>(data_varying_dim_);
  };


  /**
   * Returns a const-reference to the data identified by the index pair <tt><dim,codim></tt>.
   *
   * @note The returned container holds the geometry parametrizations
   * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}}\f$.
   */
  template <int dim,int codim>
  const VtkGridGenContSameCodim_<dim, codim> &
  get_data_dim_codim() const;

  /**
   * Returns a reference to the data identified by the index pair <tt><dim,codim></tt>.
   *
   * @note The returned container holds the geometry parametrizations
   * \f$ \mathbf{F}_i \colon \mathbb{R}^{\text{dim}} \to \mathbb{R}^{\text{dim}+\text{codim}}\f$.
   */
  template <int dim,int codim>
  VtkGridGenContSameCodim_<dim, codim> &
  get_data_dim_codim();

  /**
   * Functions container.
   */
  FunContPtr_ funcs_container_;

  /**
   * Grids information for the solid mesh.
   */
  const GridInfoPtr_ solid_info_;

  /**
   * Grids information for the knot mesh.
   */
  const GridInfoPtr_ knot_info_;


  /**
   * All the data in the VtkGridGeneratorContBase class, organized by the
   *  @p dim index (starting from 1 to 3)
   */
  DataVaryingId<VtkGridGeneratorContBaseSameDim, 1, 3> data_varying_dim_;

  /**
   * Container for numbering the generators included in the container.
   * Each entry of the vector is a tuple, whose components are:
   *   - 1st: global id of mapping function (global in the igatools framework).
   *   - 2nd: name associated to the mapping function.
   *   - 3rd: flag for indicating if it is active or not.
   *   - 4rd: flag for indicating if it is an ig mapping or not.
   */
  SafeSTLVector<std::tuple<Index, std::string, bool, bool>> generators_numbering_;

public:

  /**
   * Returns the number of grids.
   * */
  Size get_number_grids() const;

  /**
   * Returns the number of active grids.
   * */
  Size get_number_active_grids() const;

  /**
   * Returns the name of the @p id grid.
   * */
  const std::string &get_grid_name(const Index &id) const;

  /**
   * Returns the status (active/inactive) of the grid named
   * @p name.
   * */
  bool get_grid_status(const std::string &name) const;

  /**
   * Set the @p status (active/inactive) of the grid named
   * @p name.
   * */
  void set_grid_status(const std::string &name, const bool status);

  /**
   * TODO: to document.
   * */
  void set_solid_grids(vtkMultiBlockDataSet *const mb);

  /**
   * TODO: to document.
   * */
  void set_knot_grids(vtkMultiBlockDataSet *const mb);

};



class VtkIgaGridGeneratorContParm : public VtkIgaGridGeneratorContBase
{


private:
  /**
   * Self type.
   */
  typedef VtkIgaGridGeneratorContParm Self_;

  /**
   * Self shared pointer type.
   */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /**
   * Base class.
   */
  typedef VtkIgaGridGeneratorContBase Base_;

  /**
   * Functions container shared pointer type.
   */
  typedef typename Base_::FunContPtr_ FunContPtr_;


  /**
   * Grid information shared pointer type.
   */
  typedef typename Base_::GridInfoPtr_ GridInfoPtr_;

  /**
   * Constructor.
   */
  VtkIgaGridGeneratorContParm(const FunContPtr_ funcs_container,
                              const GridInfoPtr_ solid_info,
                              const GridInfoPtr_ knot_info);

  void fill_generators();

  template <int dim, int codim>
  void
  insert_generator(const MapFunPtr_<dim, codim> map_fun,
                   const std::string &map_name);

public:
  /**
   * TODO:
   */
  static SelfPtr_ create(const FunContPtr_ funcs_container,
                         const GridInfoPtr_ solid_info,
                         const GridInfoPtr_ knot_info);

  /**
   * TODO:
   */
  void update(const GridInfoPtr_ solid_info,
              const GridInfoPtr_ knot_info);


};



class VtkIgaGridGeneratorContPhys : public VtkIgaGridGeneratorContBase
{
private:

  /**
   * Self type.
   */
  typedef VtkIgaGridGeneratorContPhys Self_;

  /**
   * Self shared pointer type.
   */
  typedef std::shared_ptr<Self_> SelfPtr_;

  /**
   * Base class.
   */
  typedef VtkIgaGridGeneratorContBase Base_;

  /**
   * Functions container shared pointer type.
   */
  typedef typename Base_::FunContPtr_ FunContPtr_;


  /**
   * Grid information shared pointer type.
   */
  typedef typename Base_::GridInfoPtr_ GridInfoPtr_;

  /**
   * Control grid information shared pointer type.
   */
  typedef std::shared_ptr<VtkControlGridInformation> ControlGridInfoPtr_;

  /**
   * Constructor.
   */
  VtkIgaGridGeneratorContPhys(const FunContPtr_ funcs_container,
                              const GridInfoPtr_ solid_info,
                              const GridInfoPtr_ knot_info,
                              const ControlGridInfoPtr_ control_info);


public:
  /**
   * TODO:
   */
  static SelfPtr_ create(const FunContPtr_ funcs_container,
                         const GridInfoPtr_ solid_info,
                         const GridInfoPtr_ knot_info,
                         const ControlGridInfoPtr_ control_info);

  /**
   * TODO:
   */
  void update(const GridInfoPtr_ solid_info,
              const GridInfoPtr_ knot_info,
              const ControlGridInfoPtr_ control_info);

private:

  /**
   * Grids information for the control mesh.
   */
  const ControlGridInfoPtr_ control_info_;


  void fill_generators();

  template <int dim, int codim>
  void insert_generator(const MapFunPtr_<dim, codim> map_fun,
                        const std::string &map_name);

public:

  /**
   * Returns the number of active grids corresponding to ig mappings.
   * */
  Size get_number_active_grids_ig() const;


  /**
   * TODO: to document.
   * */
  void set_control_grids(vtkMultiBlockDataSet *const mb);

};

IGA_NAMESPACE_CLOSE

#endif // VTK_GRID_GENERATOR_CONTAINER_H_
