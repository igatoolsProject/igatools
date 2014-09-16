#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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

The user supplies at configure time the physical spaces that 
the library shall be used for.

1) The user physical spaces require a face physical space

2) Each physical space requires:
   - A ref space
   - A mapping
   - a push-foward

3) Each mapping of the physical spaces can be an Igmapping
   that requires a refspace of the appropriate type
   
4) In igatools each reference space is treated as a special physical space,
   requiring:
   - an identity mapping of codimenion  0
   - a push-foward of h_grad type
   

"""

# Removes duplicates of a list while keeping the original order
def unique(seq):
   checked = []
   for e in seq:
      if e not in checked:
         checked.append(e)
   return checked

# Class specifying the description of a physical space.
class PhysSpaceSpecs:
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
        if isinstance(other, PhysSpaceSpecs):
            return(self.dim == other.dim) and (self.codim == other.codim) and (self.range == other.range) and (self.rank == other.rank)and (self.trans_type == other.trans_type) 
        return NotImplemented
   def physical_range(self, ref_range, space_dim, trans_type):
      if trans_type == 'h_grad':
         return ref_range
      if trans_type == 'h_div':
         return space_dim

   def physical_rank(self, ref_rank):
      return ref_rank


class PhysSpace:
    def __init__(self, specs, ref_space):
       self.spec   = specs
       self.name   = ( 'PhysicalSpace <' + ref_space + 
                       '<%d,%d,%d>' % (specs.dim, specs.range, specs.rank) +
                       ', PushForward<Transformation::%s, %d, %d> >' 
                       %(specs.trans_type, specs.dim, specs.codim) )

class RefSpace:
    def __init__(self, specs, ref_space):
       self.spec   = specs
       self.name   = ref_space + '<%d,%d,%d>' % (specs.dim, specs.range, specs.rank)     

class FunctionRow:
   #function dim, range and rank
   def __init__(self, arg_list):
      self.dim        = arg_list[0]
      self.range      = arg_list[1]  
      self.rank       = arg_list[2]
      return None



class MappingRow:
   #mappings dim, codim and space_dim
   def __init__(self, arg_list):
      self.dim        = arg_list[0]
      self.codim      = arg_list[1]  
      self.space_dim  = self.dim + self.codim
      return None


class PForwRow:
   #mappings dim, codim and space_dim
   def __init__(self, arg_list):
      self.dim        = arg_list[0]
      self.codim      = arg_list[1]  
      self.trans_type = arg_list[2]
      return None

  
class RefSpaceRow:
    def __init__(self, arg_list):
      self.dim        = arg_list[0]
      self.range      = arg_list[1]  
      self.rank       = arg_list[2]
      return None
    def __eq__(self, other):
       if isinstance(other, RefSpaceRow):
           return(self.dim == other.dim)  and (self.range == other.range) and (self.rank == other.rank)
       return NotImplemented
   
class InstantiationInfo:
   """ Stores "tables" with useful entries to be used for instantiations.
   
   This information is generated using a table of
   physical spaces that was genererated at configure time
   by user passed options.

   """
  
   def __init__(self, filename, max_der_order):
      """The constructor."""
      self.user_phy_sp_dims =[] # Physical spaces that the library is suppussed to be used on
      self.face_phy_sp_dims =[] # Physical Spaces that are faces of the user spaces
      self.all_phy_sp_dims  =[] #the physical spaces the user provides plus the one that are necesary on top
      self.igm_phy_sp_dims=[]

      self.function_dims=[] # table of dim, range, rank for functions
      
      self.user_mapping_dims =[] # table of dim codim
      self.all_mapping_dims  =[] # table of dim codim
      
      self.user_ref_sp_dims=[]
      self.face_ref_sp_dims=[]
      self.all_ref_sp_dims=[]
      self.igm_ref_sp_dims=[]
      self.really_all_ref_sp_dims=[]
       
      self.deriv_order = range(int(max_der_order)+1)
      self.derivatives=[]  # allderivative classes
      self.values=[]
      self.divs=[]
      
      self.domain_dims = [] # list all domain dimensions
      self.user_domain_dims = []
      self.face_domain_dims = []

#---------------------------------------
      self.RefSpaces=[]     # all required reference spaces
      self.UserRefSpaces=[]
      self.UserFilteredRefSpaces=[] #iga mapping required ref spaces

      
      self.UserPhysSpaces=[]
      self.PhysSpaces=[]
 #---------------------------------------
 
      
      self.read_dimensions_file(filename)
      
      self.create_mapping_dims()
      self.create_function_dims()
      self.create_derivatives()
      self.create_ref_dim()
      
      self.create_ref_spaces()
      self.create_PhysSpaces()
      return None



   def read_dimensions_file(self, filename):
      '''Reads a text file where each line describes a physical space and
            genereate the main tables '''

      file_input = open(filename, 'r')
      user_spaces=[]
      for i in file_input:
         row = i.strip().split()
         if (len(row) > 0) and (row[0] != '#') :
            user_spaces.append( [int(x) for x in row[0:4]] + row[-1:])
      file_input.close()
    
      
      for row in user_spaces:
         self.user_phy_sp_dims.append(PhysSpaceSpecs(row))
         self.all_phy_sp_dims.append(PhysSpaceSpecs(row))

      #Add the spaces for the faces     
      face_spaces = unique ([ [sp.dim-1, sp.codim+1, sp.range, sp.rank, sp.trans_type]
                    for sp in self.user_phy_sp_dims ] )

      for row in face_spaces:
         self.face_phy_sp_dims.append(PhysSpaceSpecs(row))
         self.all_phy_sp_dims.append(PhysSpaceSpecs(row))

      self.all_phy_sp_dims = unique(self.all_phy_sp_dims)
            
      self.domain_dims = unique([sp.dim for sp in self.all_phy_sp_dims])
           
      return None



   def create_mapping_dims(self):
      '''Fills all_mapping_dims with a list of all mappings '''
      dims_list = unique([ [row.dim,  row.codim] for row in self.all_phy_sp_dims])
      self.all_mapping_dims = [MappingRow(row) for row in dims_list]
               
      dims_list = unique([ [row.dim,  row.codim] for row in self.user_phy_sp_dims])
      self.user_mapping_dims = [MappingRow(row) for row in dims_list]     
      return None


   
   def create_function_dims(self):
      dims_list=[]
      for row in self.all_phy_sp_dims:
         # Add derivative for the reference space and physical space
         dims_list.append((row.dim,  row.range, row.rank))
         dims_list.append((row.space_dim, row.phys_range, row.phys_rank))
         
      for row in self.all_mapping_dims:
         dims_list.append((row.dim,  row.space_dim, 1))
         dims_list.append((row.space_dim,  row.dim, 1))
      
      for row in self.igm_ref_sp_dims:
         dims_list.append((row.dim,  row.range, 1))
                 
      for row in unique(dims_list):
          self.function_dims.append(FunctionRow(row))
     
      #print(dims_list)
      return None
 
     

   def create_ref_spaces(self):
      ''' Creates a list of Reference spaces as table and as classes'''
       
      self.user_ref_sp_dims = unique( [RefSpaceRow([x.dim, x.range, x.rank])
                                       for x in self.user_phy_sp_dims] )
      self.all_ref_sp_dims  = unique( [RefSpaceRow([x.dim, x.range, x.rank]) 
                                       for x in self.all_phy_sp_dims] )
      self.face_ref_sp_dims = unique( [RefSpaceRow([x.dim, x.range, x.rank])
                                       for x in self.face_phy_sp_dims] )

      self.igm_ref_sp_dims = unique( [RefSpaceRow([x.dim, x.space_dim, 1])
                                       for x in self.all_mapping_dims] )
      
      self.really_all_ref_sp_dims=unique(self.all_ref_sp_dims + self.igm_ref_sp_dims)
     
      #self.all_phy_sp_dims.append(  
      a= unique( [PhysSpaceSpecs([x.dim, 0, x.range, x.rank, 'h_grad'])
                                           for x in self.really_all_ref_sp_dims] +
                                         [PhysSpaceSpecs([x.dim-1, 1, x.range, x.rank, 'h_grad'])
                                          for x in self.user_ref_sp_dims])
                                    #)
      
      self.all_phy_sp_dims = unique(self.all_phy_sp_dims+a)
     
#       self.igm_phy_sp_dims = unique( [PhysSpaceRow([x.dim, 0, x.space_dim, 1, 'h_grad'])
#                                        for x in self.igm_ref_sp_dims] )
     
      RefDims = ['<%d,%d,%d>' % (x.dim, x.range, x.rank)
                                for x in self.all_ref_sp_dims ]  

      UserRefDims = ['<%d,%d,%d>' % (x.dim, x.range, x.rank)
                     for x in self.user_ref_sp_dims ] 

      IgRefDims = ['<%d,%d,%d>' % (x.dim, x.range, x.rank)
                     for x in self.igm_ref_sp_dims ]
      
   #todo : for only bspsline
      spaces = ('BSplineSpace', 'NURBSSpace')
#      spaces = ['BSplineSpace']
      self.RefSpaces = ( ['%s%s' % (sp, dims) for sp in spaces
                            for dims in RefDims] )
      self.UserRefSpaces = ( ['%s%s' % (sp, dims)
                                for sp in spaces
                                for dims in UserRefDims] )
      self.IgmRefSpaces = ( ['%s%s' % (sp, dims)
                             for sp in spaces
                             for dims in IgRefDims] )
      

      return None



   def create_PhysSpaces(self):
      
      #todo : for only bspsline
      spaces = ('BSplineSpace', 'NURBSSpace')
#      spaces = ['BSplineSpace']
      self.PhysSpaces_v2 = unique( [PhysSpace(x,sp)
                                 for sp in spaces
                                 for x in self.all_phy_sp_dims] )
      
      self.AllRefSpaces_v2 = unique( [RefSpace(x,sp)
                                      for sp in spaces
                                      for x in self.really_all_ref_sp_dims ] )
      
      
      self.AllPushForwards = unique(['PushForward<Transformation::%s, %d, %d>'
                                    %(x.trans_type, x.dim, x.codim) for x in self.all_phy_sp_dims] )
      
      self.RefPushForwards = unique(['PushForward<Transformation::%s, %d, %d>'
                                    %('h_grad', x.dim, 0) for x in self.all_ref_sp_dims + self.igm_ref_sp_dims] )
      
     
      self.PhysSpaces = unique( ['PhysicalSpace <' +
                                   '%s<%d,%d,%d>' % (sp, x.dim, x.range, x.rank) +
                                   ', PushForward<Transformation::%s, %d, %d> >'
                                   %(x.trans_type, x.dim, x.codim)
                                   for sp in spaces
                                   for x in self.all_phy_sp_dims] )
      
      
      
      
      self.UserPhysSpaces = unique( ['PhysicalSpace <' +
                                     '%s<%d,%d,%d>' % (sp, x.dim, x.range, x.rank) +
                                     ', PushForward<Transformation::%s, %d, %d> >'
                                     %(x.trans_type, x.dim, x.codim)
                                     for sp in spaces
                                     for x in self.user_phy_sp_dims] )


   def create_ref_dim(self):
      self.ref_dom_dims = unique([x.dim for x in self.all_phy_sp_dims])
      self.user_ref_dom_dims = unique([x.dim for x in self.user_phy_sp_dims])
      self.face_ref_dom_dims = unique([x.dim for x in self.face_phy_sp_dims])
      return None

  

   def create_derivatives(self):
      '''Creates a list of the tensor types for the required values and derivatives'''
      mapping_list = [FunctionRow([x.dim, x.space_dim, 1]) for x in self.all_mapping_dims]
      mapping_list.append(FunctionRow([0, 0, 1])) #todo fix approprietly
      dims_list = self.function_dims + mapping_list
      deriv ='Tensor<dim, order, tensor::covariant, Tensor<range, rank, tensor::contravariant, Tdouble>>'
      value ='Tensor<range, rank, tensor::contravariant, Tdouble>'
      div   ='Tensor<range, rank-1, tensor::contravariant, Tdouble>'

      deriv_list=[]
      value_list=[]
      div_list=[]
      for order in self.deriv_order:
         for dims in dims_list:
            (dim, range, rank) = (dims.dim, dims.range, dims.rank)
            if order == 0:
               (dim, order) = (1,1)
            replace_table = (('order', str(order)), ('dim', str(dim)),('range', str(range)),('rank', str(rank)))
            temp = deriv
            temp_v = value
            temp_d = div
            for rep in replace_table:
               temp = temp.replace(rep[0], rep[1])
               temp_v = temp_v.replace(rep[0], rep[1])
               temp_d = temp_d.replace(rep[0], rep[1])
            deriv_list.append(temp)
            value_list.append(temp_v)
            div_list.append(temp_d)
            
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
        self.inst = InstantiationInfo(args['config_file'], args['max_der_order'])
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
        