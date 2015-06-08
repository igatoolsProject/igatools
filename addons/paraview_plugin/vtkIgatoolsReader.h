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

#ifndef VTK_IGATOOLS_READER_H_
#define VTK_IGATOOLS_READER_H_

#include "vtkMultiBlockDataSetAlgorithm.h"

#include <paraview_plugin/iga_vtk.h>

class vtkStructuredGrid;

// class vtkIgatoolsReader : public vtkUnstructuredGridAlgorithm
class vtkIgatoolsReader : public vtkMultiBlockDataSetAlgorithm
{
public:

    vtkTypeMacro(vtkIgatoolsReader, vtkMultiBlockDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);

    static vtkIgatoolsReader *New();

    /*
     * Description:
     * Specify file name of the .iga file.
     */
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    /*
     * Description:
     * Set the number of VTK visualization elements per direction for each Bezier
     * element for the physical solid mesh.
     */
    vtkSetVector3Macro(NumVisualizationElementsPhysicalSolid, int);
    vtkGetVectorMacro(NumVisualizationElementsPhysicalSolid, int, 3);

    /*
     * Description:
     * Set the number of VTK visualization elements per direction for each Bezier
     * element the parametric solid mesh.
     */
    vtkSetVector3Macro(NumVisualizationElementsParametricSolid, int);
    vtkGetVectorMacro(NumVisualizationElementsParametricSolid, int, 3);

    /*
     * Description:
     * Set the number of VTK visualization elements per direction for each Bezier
     * element for the physical knot mesh.
     */
    vtkSetVector3Macro(NumVisualizationElementsPhysicalKnot, int);
    vtkGetVectorMacro(NumVisualizationElementsPhysicalKnot, int, 3);

    /*
     * Description:
     * Set the number of VTK visualization elements per direction for each Bezier
     * element the parametric knot mesh.
     */
    vtkSetVector3Macro(NumVisualizationElementsParametricKnot, int);
    vtkGetVectorMacro(NumVisualizationElementsParametricKnot, int, 3);

    /*
     * Set/Get Grid type for the physical solid mesh.
     *  - 0: unstructured grid : quadratic elements.
     *  - 1: unstructured grid : linear elements.
     *  - 2: structured grid.
     */
    vtkSetMacro(GridTypePhysicalSolid, int);
    vtkGetMacro(GridTypePhysicalSolid, int);

    /*
     * Set/Get Grid type for the parametric solid mesh
     *  - 0: unstructured grid : quadratic elements.
     *  - 1: unstructured grid : linear elements.
     *  - 2: structured grid.
     */
    vtkSetMacro(GridTypeParametricSolid, int);
    vtkGetMacro(GridTypeParametricSolid, int);

    /*
     * Set/Get Grid type for the physical knot mesh.
     *  - 0: unstructured grid : quadratic elements.
     *  - 1: unstructured grid : linear elements.
     *  - 2: structured grid.
     */
    vtkSetMacro(GridTypePhysicalKnot, int);
    vtkGetMacro(GridTypePhysicalKnot, int);

    /*
     * Set/Get Grid type for the parametric knot mesh
     *  - 0: unstructured grid : quadratic elements.
     *  - 1: unstructured grid : linear elements.
     *  - 2: structured grid.
     */
    vtkSetMacro(GridTypeParametricKnot, int);
    vtkGetMacro(GridTypeParametricKnot, int);

    /*
     * Set/Get Solid Mesh creation flag.
     *  - true:  create the mesh.
     *  - false: do not create the mesh.
     */
    vtkSetMacro(SolidMeshPhysical, bool);
    vtkGetMacro(SolidMeshPhysical, bool);

    /*
     * Set/Get Solid Mesh creation flag.
     *  - true:  create the mesh.
     *  - false: do not create the mesh.
     */
    vtkSetMacro(SolidMeshParametric, bool);
    vtkGetMacro(SolidMeshParametric, bool);

    /*
     * Set/Get Control Mesh creation flag.
     *  - true:  create the mesh.
     *  - false: do not create the mesh.
     */
    vtkSetMacro(ControlMeshPhysical, bool);
    vtkGetMacro(ControlMeshPhysical, bool);

    /*
     * Set/Get Knot Mesh creation flag.
     *  - true:  create the mesh.
     *  - false: do not create the mesh.
     */
    vtkSetMacro(KnotMeshPhysical, bool);
    vtkGetMacro(KnotMeshPhysical, bool);

    /*
     * Set/Get Knot Mesh creation flag.
     *  - true:  create the mesh.
     *  - false: do not create the mesh.
     */
    vtkSetMacro(KnotMeshParametric, bool);
    vtkGetMacro(KnotMeshParametric, bool);

    /*
     * Set/Get Parametric Mesh creation flag.
     *  - true:  create the mesh.
     *  - false: do not create the mesh.
     */
    vtkSetMacro(ParametricMesh, bool);
    vtkGetMacro(ParametricMesh, bool);

    /*
     * Set/Get Physical Mesh creation flag.
     *  - true:  create the mesh.
     *  - false: do not create the mesh.
     */
    vtkSetMacro(PhysicalMesh, bool);
    vtkGetMacro(PhysicalMesh, bool);

protected:
    /*
     * Deleter
     */
    vtkIgatoolsReader();

    /*
     * Destructor.
     */
    ~vtkIgatoolsReader() {}

    virtual int RequestData(vtkInformation *,
                            vtkInformationVector **,
                            vtkInformationVector *) override final;

    virtual int RequestInformation(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) override final;

public:
    /*
     * Test whether the file with the given name exists and can be read by this
     * reader.
     */
    int CanReadFile(const char *name);

private:
    vtkIgatoolsReader(const vtkIgatoolsReader &) = delete;
    vtkIgatoolsReader(const vtkIgatoolsReader &&) = delete;
    void operator=(const vtkIgatoolsReader &) = delete;
    void operator=(const vtkIgatoolsReader&&) = delete;


    /*
     * Retrieves the file name and the path of the file.
     */
    void get_file_and_path(std::string &file_name, std::string &file_path);

    /*
     * Check number of visualization element.
     */
    void check_number_visualization_elements();

    /*
     * Grid type for the physical solid mesh.
     */
    int GridTypePhysicalSolid;

    /*
     * Grid type for the parametric solid mesh.
     */
    int GridTypeParametricSolid;

    /*
     * Grid type for the physical knot mesh.
     */
    int GridTypePhysicalKnot;

    /*
     * Grid type for the parametric knot mesh.
     */
    int GridTypeParametricKnot;

    /*
     * Solid mesh creation flag.
     */
    bool SolidMeshPhysical;

    /*
     * Solid mesh creation flag.
     */
    bool SolidMeshParametric;

    /*
     * Control mesh creation flag.
     */
    bool ControlMeshPhysical;

    /*
     * Knot mesh creation flag.
     */
    bool KnotMeshPhysical;

    /*
     * Knot mesh creation flag.
     */
    bool KnotMeshParametric;

    /*
     * Parametric mesh creation flag.
     */
    bool ParametricMesh;

    /*
     * Physical mesh creation flag.
     */
    bool PhysicalMesh;

    /*
     * File name variable.
     */
    char *FileName = NULL;

    /*
     * Number of VTK visualization elements per direction per Bezier element
     * for the physical solid mesh.
     */
    int NumVisualizationElementsPhysicalSolid[3];

    /*
     * Number of VTK visualization elements per direction per Bezier element
     * for the parametric solid mesh.
     */
    int NumVisualizationElementsParametricSolid[3];

    /*
     * Number of VTK visualization elements per direction per Bezier element
     * for the physical knot mesh.
     */
    int NumVisualizationElementsPhysicalKnot[3];

    /*
     * Number of VTK visualization elements per direction per Bezier element
     * for the parametric knot mesh.
     */
    int NumVisualizationElementsParametricKnot[3];

    /*
     * Iga vtk object.
     */
    iga::IGAVTK iga_vtk_;
};


#endif // VTK_IGATOOLS_READER_H_