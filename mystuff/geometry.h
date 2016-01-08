IGA_NAMESPACE_OPEN

template<int dim>
class Geometry
{
public:
  TensorSize<dim>  nel;
  TensorIndex<dim> deg;
  IgCoefficients   coefs;
  IgCoefficients   weights;
  void load(const char *fname)
  {
    FILE *fp;
    int check;
    fp=fopen(fname,"r");
    if (fp==0) std::cout << "cannot open the .nurbs file!" << std::endl;
    else
    {
      // checking dimension
      fscanf(fp,"%d",&check);
      if (check!=dim) std::cout << "wrong geometry dimension!" << std::endl;
      else
      {
        // reading degrees
        for (int idim=0; idim<dim; idim++)
          fscanf(fp,"%d",&deg[idim]);
        // reading elements
        for (int idim=0; idim<dim; idim++)
          fscanf(fp,"%d",&nel[idim]);
        // reading number of control points
        int ncp;
        fscanf(fp,"%d",&ncp);
        // reading control points
        double data;
        for (int icp=0; icp<ncp*dim; icp++)
        {
          fscanf(fp,"%lf",&data);
          coefs[icp]=data;
        }
        // reading weights
        for (int icp=0; icp<ncp; icp++)
        {
          fscanf(fp,"%lf",&data);
          weights[icp]=data;
        }
      }
      fclose(fp);
    }
  }
};

IGA_NAMESPACE_CLOSE
