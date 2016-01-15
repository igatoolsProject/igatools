#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
#
# This file is part of the igatools library.
#
# The igatools library is free software: you can use it, redistribute
# it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-+--------------------------------------------------------------------



"""@package init_instantiation_data

This module is loaded from all instantiation scripts.
It will read a .txt table created by another script called
generate_instantiation_table.py

Instantiation rational

Instantiation dependencies

The user supplies at configure time the physical bases that 
the library shall be used for.

1) The user physical bases require a face physical basis

2) Each physical basis requires:
   - A ref basis
   - A domain
   - a push-foward

3) Each domain of the physical bases can be based on a IgGridFunction
   that requires a ReferenceBasis of the appropriate type
   
4) In igatools each reference basis is treated as a special physical basis,
   requiring:
   - an identity domain of codimenion  0
   - a push-foward of h_grad type
   

"""

# Removes duplicates of a list while keeping the original order
def unique(seq):
   checked = []
   for e in seq:
      if e not in checked:
         checked.append(e)
   return checked

# Class specifying the description of a physical basis.
class PhysBasisSpecs:
   # Constructor of the class.
   def __init__(self, arg_list):
      self.dim        = arg_list[0]
      self.codim      = arg_list[1]
      self.range      = arg_list[2]
      self.rank       = arg_list[3]
      self.trans_type = arg_list[4]
      self.space_dim  = self.dim +  self.codim
      self.phys_range = self.physical_range(self.range, self.space_dim, 
                                            self.trans_type)
      self.phys_rank  = self.physical_rank(self.rank)
      return None
   def __eq__(self, other):
        if isinstance(other, PhysBasisSpecs):
            return(self.dim == other.dim) and \
              (self.codim == other.codim) and \
              (self.range == other.range) and \
              (self.rank == other.rank) and \
              (self.trans_type == other.trans_type) 
        return NotImplemented
    
   def __hash__(self):
        return hash((self.dim, self.codim, self.range, self.rank, self.trans_type))  

    
   def physical_range(self, ref_range, space_dim, trans_type):
      if trans_type == 'h_grad':
         return ref_range
      if trans_type == 'h_div':
         return space_dim

   def physical_rank(self, ref_rank):
      return ref_rank


class PhysBasis:
    def __init__(self, specs, ref_basis, name):
       self.spec   = specs
       self.name   = ( '%s <' %name + 
                       '%d,%d,%d,' % (specs.dim, specs.range, specs.rank) +
                       '%d>' %(specs.codim) )

class RefBasis:
    def __init__(self, specs, ref_basis):
       self.spec   = specs
       self.name   = ref_basis + '<%d,%d,%d>' % (specs.dim, specs.range, specs.rank)     


class FunctionRow:
   #function dim, codim, range and rank
   def __init__(self, arg_list):
      self.dim        = arg_list[0]
      self.codim      = arg_list[1]
      self.range      = arg_list[2]  
      self.rank       = arg_list[3]
      return None
  
   def __eq__(self, other):
      if isinstance(other, FunctionRow):
          return(self.dim == other.dim) \
              and (self.range == other.range) \
              and (self.rank == other.rank) \
              and (self.codim == other.codim) 
      return NotImplemented
  
   def __hash__(self):
        return hash((self.dim, self.codim, self.range, self.rank))
     

class MappingRow:
   #mappings dim, codim and space_dim
   def __init__(self, arg_list):
      self.dim        = arg_list[0]
      self.codim      = arg_list[1]  
      self.space_dim  = self.dim + self.codim
      return None
  
   def __eq__(self, other):
      if isinstance(other, MappingRow):
          return(self.dim == other.dim) \
              and (self.codim == other.codim) 
      return NotImplemented
  
   def __hash__(self):
        return hash((self.dim, self.codim, self.space_dim))
    
class PForwRow:
   #mappings dim, codim and space_dim
   def __init__(self, arg_list):
      self.dim        = arg_list[0]
      self.codim      = arg_list[1]  
      self.trans_type = arg_list[2]
      return None

   def __eq__(self, other):
      if isinstance(other, MappingRow):
          return(self.dim == other.dim) \
              and (self.codim == other.codim) \
              and (self.trans_type == other.trans_type)
      return NotImplemented
  
class RefBasisRow:
    def __init__(self, arg_list):
      self.dim        = arg_list[0]
      self.range      = arg_list[1]  
      self.rank       = arg_list[2]
      return None
  
    def __eq__(self, other):
       if isinstance(other, RefBasisRow):
           return(self.dim == other.dim)  \
               and (self.range == other.range) \
               and (self.rank == other.rank)
       return NotImplemented
    def __hash__(self):
        return hash((self.dim, self.range, self.rank))
    
class InstantiationInfo:
   """ Stores "tables" with useful entries to be used for instantiations.
   
   This information is generated using a table of
   physical bases that was genererated at configure time
   by user passed options.

   """
  
   def __init__(self, filename, max_der_order, nurbs):
      """The constructor."""
      self.n_sub_element = 1 #1 for faces, 2 for faces and faces-1, max dim for all
      
      self.phy_sp_dims  =[] # Physical bases that the library is suppussed to be used on
      self.ref_sp_dims  =[]
      self.mapping_dims =[]
      self.push_fw_dims =[]
      self.function_dims=[]
      self.domain_dims = []
      
      self.all_phy_sp_dims  =[] # Physical bases that the library is suppussed to be used on
      self.all_ref_sp_dims  =[]
      self.all_mapping_dims =[]
      self.all_push_fw_dims =[]
      self.all_function_dims=[]
      self.all_domain_dims = []
      
      
      self.ig_bases = ['ReferenceBasis']
#      self.ig_bases = ['BSplineBasis'] if nurbs == 'OFF' \
#      else ['BSplineBasis', 'NURBSBasis']

      self.deriv_order = range(int(max_der_order)+1)
      self.derivatives=[]  # allderivative classes
      self.values=[]
      self.divs=[]
      
      self.iterators = ['GridIteratorBase<Accessor>',
                               'GridIterator<Accessor>',
                               'ConstGridIterator<Accessor>']
#---------------------------------------
      self.RefBases=[]     # all required reference bases
      self.UserPhysBases=[]
      self.PhysBases=[]
 #---------------------------------------
 
      
      self.read_dimensions_file(filename)
      
      
      self.create_derivatives()
      
      
      self.create_ref_bases()
      self.create_PhysBases()
      
      
      return None

   def sub_dims(self, dim):
       return range(max(0, dim - self.n_sub_element), dim + 1)

   def read_dimensions_file(self, filename):
      '''Reads a text file where each line describes a physical basis and
            genereate the main tables '''

      file_input = open(filename, 'r')
      user_bases=[]
      for i in file_input:
         row = i.strip().split()
         if (len(row) > 0) and (row[0] != '#') :
            user_bases.append( [int(x) for x in row[0:4]] + row[-1:])
      file_input.close()
      phy_sp_dims = [PhysBasisSpecs(row)  for row in user_bases]
      
                            
      #------------------------------------------
      for k in range(0, self.n_sub_element + 1):   
          bases     = [ PhysBasisSpecs([x.dim-k, x.codim+k, x.range, x.rank, x.trans_type])
                         for x in phy_sp_dims if x.dim>=k]
          
          ref_bases = unique( [RefBasisRow([x.dim, x.range, x.rank])
                                for x in bases] )
          mappings   = unique( [MappingRow([x.dim,  x.codim]) 
                                for x in bases] )
          functions  = unique( [FunctionRow([x.dim,  x.codim, x.phys_range, x.phys_rank]) 
                                for x in bases] )
      
          functions = unique( functions +
                              [FunctionRow([x.dim,  0, x.range, x.rank]) 
                               for x in ref_bases] +
                              [FunctionRow([x.dim,  0, x.dim, 1]) 
                               for x in ref_bases])
          
          map_funcs = unique([FunctionRow([x.dim,  0, x.space_dim, 1]) 
                              for x in mappings] )
          
          functions = unique(functions + map_funcs)
          
          ref_bases = unique(ref_bases + 
                              [RefBasisRow([x.dim, x.range, 1]) for x in map_funcs]) 
         
          if k == 0:
             self.phy_sp_dims   = unique(bases)
             self.ref_sp_dims   = unique(ref_bases)
             self.mapping_dims  = unique(mappings)
             self.function_dims = unique(functions)
             
             
          self.all_phy_sp_dims   = unique(self.all_phy_sp_dims + bases)
          self.all_ref_sp_dims   = unique(self.all_ref_sp_dims + ref_bases)
          self.all_mapping_dims  = unique(self.all_mapping_dims + mappings)
          self.all_function_dims = unique(self.all_function_dims + functions)


      # for the get sub bases of the reference bases    
      self.all_phy_sp_dims = unique(self.all_phy_sp_dims +
                                    [PhysBasisSpecs([x.dim, 0, x.range, x.rank, 'h_grad'])
                                      for x in self.all_ref_sp_dims])    
      self.phy_sp_dims = unique(self.phy_sp_dims +
                                    [PhysBasisSpecs([x.dim, 0, x.range, x.rank, 'h_grad'])
                                      for x in self.all_ref_sp_dims])
      for k in range(1, self.n_sub_element + 1):
         self.all_phy_sp_dims = unique(self.all_phy_sp_dims +
                                    [PhysBasisSpecs([x.dim-k, k, x.range, x.rank, 'h_grad']) 
                                     for x in self.ref_sp_dims  if x.dim>=k])

      self.all_mapping_dims = unique(self.all_mapping_dims +
                                    [MappingRow([x.dim,  x.codim])  for x in self.all_phy_sp_dims] )
      self.all_function_dims = unique(self.all_function_dims +
                                    [FunctionRow([x.dim,  x.codim, x.phys_range, x.phys_rank]) 
                                for x in self.all_phy_sp_dims] )
      self.domain_dims     = unique([x.dim for x in self.function_dims])   
      self.all_domain_dims = unique([x.dim for x in self.all_function_dims])

      self.sub_domain_dims = list(set(self.all_domain_dims) - set(self.domain_dims))
      self.sub_ref_sp_dims = list(set(self.all_ref_sp_dims) - set(self.ref_sp_dims))
      self.sub_function_dims = list(set(self.all_function_dims) - set(self.function_dims))
      self.sub_mapping_dims = list(set(self.all_mapping_dims) - set(self.mapping_dims))
      self.sub_phy_sp_dims = list(set(self.all_phy_sp_dims) - set(self.phy_sp_dims))
      return None


  



   def create_ref_bases(self):
      ''' Creates a list of Reference bases as table and as classes'''
    
      return None

  
   def create_PhysBases(self):
       
      self.PhysBases   = unique( [PhysBasis(x,sp,'PhysicalBasis')
                                 for sp in self.ig_bases
                                 for x in self.phy_sp_dims] )
      self.AllPhysBases = unique( [PhysBasis(x,sp,'PhysicalBasis')
                                 for sp in self.ig_bases
                                 for x in self.all_phy_sp_dims] )
      
      self.SubPhysBases = unique( [PhysBasis(x,sp,'PhysicalBasis')
                                 for sp in self.ig_bases
                                 for x in self.sub_phy_sp_dims] )
      
      self.RefBases = unique( [RefBasis(x,sp)
                                for sp in self.ig_bases
                                for x in self.ref_sp_dims ] )
      self.AllRefBases = unique( [RefBasis(x,sp)
                                   for sp in self.ig_bases
                                   for x in self.all_ref_sp_dims ] )


      
   def create_derivatives(self):
      '''Creates a list of the tensor types for the required values and derivatives'''
          
      inv_func = unique([FunctionRow([x.space_dim,  0, x.dim, 1]) 
                              for x in self.all_mapping_dims] )
      self.dims_list = unique (self.all_function_dims + inv_func)

      deriv_list=[]
      value_list=[]
      div_list=[]
      for order in self.deriv_order:
         for dims in self.dims_list:
            (dim, range, rank) = (dims.dim, dims.range, dims.rank)
            if order == 0:
               (dim, order) = (1,1)
            deriv ='Tensor<%d, %d, tensor::covariant, Tensor<%d, %d, tensor::contravariant, Tdouble>>' % (dim, order, range, rank)
            value ='Tensor<%d, %d, tensor::contravariant, Tdouble>' % (range, rank)
            div   ='Tensor<%d, %d, tensor::contravariant, Tdouble>' % (range, rank - 1)
            deriv_list.append(deriv)
            value_list.append(value)
            div_list.append(div)

      self.derivatives = unique(deriv_list)
      self.values = unique(value_list)
      self.divs = unique(div_list)
      return None




class Instantiation:
    """ Main function called at the beginning of all instatiation scripts."""


    def __init__(self, inc_files=[], other_inc_files=[], verbose=False):
        #Getting a dictionary or arguments.
        from sys import argv as sysargv
        from os import sep as ossep
        args = dict([arg.split('=') for arg in sysargv[1:]])    

        # Reading information from dimensions file.
        self.inst = InstantiationInfo(args['config_file'], 
                                      args['max_der_order'], 
                                      args['nurbs'])
        #  Some debug information printing
        if verbose:
            print('dim codim range rank space_dim')
            for x in self.inst.all_phy_sp_dims:
                print (x.dim, x.codim, x.range, x.rank, x.space_dim, x.trans_type)
                print (x)
                
            for x in self.inst.all_ref_sp_dims:
                print (x.dim,  x.range, x.rank)
        # Openning the output file.
        self.file_output = open(args['out_file'], 'w')
        # Writing the header.
        header = ( '// This file was automatically generated ' +
                   'from %s \n' % (sysargv[0].split(ossep)[-1]) +
                   '// DO NOT edit as it will be overwritten.\n\n')
        self.file_output.write(header)
        if inc_files:
            for file in inc_files:
                self.file_output.write('#include <igatools/%s>\n' %file)
        if other_inc_files:
            for file in other_inc_files:
                self.file_output.write('#include <%s>\n' %file)        
        self.file_output.write('IGA_NAMESPACE_OPEN\n')
       
    def __del__(self):
        self.file_output.write('IGA_NAMESPACE_CLOSE\n')
        self.file_output.close() 
        


#-------------------------------------------------
def skel_size(dim,sdim):
  res = 0
  if dim == sdim:
    res = 1
  elif ((sdim < dim) and (sdim >= 0)):
    if sdim > 0:
      res = 2 * skel_size(dim-1, sdim) + skel_size(dim-1, sdim-1)
    else:
      res = 2**dim

  return res;
#-------------------------------------------------
