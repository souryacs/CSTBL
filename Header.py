#!/usr/bin/env python

import dendropy
from dendropy import TreeList, Tree, Taxon, TaxonSet
import time
import os
import sys
from optparse import OptionParser
import math

# this is the path of QP executable based on GNU C library
#QP_Executable = './GNU_BFGS2'

""" 
this dictionary defines the taxa pair relations
each entry is indexed by two leaves 
"""
TaxaPair_Reln_Dict = dict()

""" 
this list contains the complete set of taxa present in the input source trees 
"""
COMPLETE_INPUT_TAXA_LIST = []

# this is the debug level
# set for printing the necessary information
DEBUG_LEVEL = 0

# this is a global string containing the objective function of the quadratic 
# programming employed for branch length determination of the tree
objective_function_string = ''

# this is a global string containing the derivative of the objective function
# employed for the above mentioned quadratic programming
deriv_objective_function_string = ''

# this dictionary stores the edge information
# key of the dictionary: (terminal_node1, terminal_node2)
# value of the element: array index (which will be used to denote corresponding variable with respect to the branch)
# variables are stored as x[0], x[1], .....
# they are employed in a quadratic programming
EdgeInfoDict = dict()

""" this variable associates weight of matrix corresponding to individual source trees """
Matrix_Weight_Val = []

# this is the corrected edge length for the output weighted supertree
# when the QP computation returns a negative edge length
CORRECTED_POSITIVE_EDGE_LEN = 0.00001

##-----------------------------------------------------
""" 
this class defines a pair of taxa
key of this class --- taxa1, taxa2  
"""
class Reln_TaxaPair(object):
  def __init__(self):    
            
    """ 
    with respect to the output supertree, this array stores the branches between this couplet
    (indexed by some mechanism) which are between this couplet and their MRCA node
    """    
    self.Branch_Array_Idx_List = []
    
    """
    with respect to the input treelist, this list stores the input tree indices
    where the taxa pair is supported
    """
    self.Supported_Tree_Idx_List = []

    """ 
    this list stores the distance values between these taxa pairs
    as supported by the source trees
    """
    self.DistMatValues = []
            
  """ 
  add branch length distance value (for an input tree)
  to a list
  """
  def _AddDistMatValue(self, val):
    self.DistMatValues.append(val)

  """
  adds supporting tree
  """
  def _AddSupportTreeIndex(self, idx):
    self.Supported_Tree_Idx_List.append(idx)
  
  """
  returns the indices of trees supporting this couplet
  """
  def _GetSupportTreeList(self):
    return self.Supported_Tree_Idx_List
  
  """
  returns the number of trees supporting this couplet
  """
  def _GetNoSupportTrees(self):
    return len(self.Supported_Tree_Idx_List)
    
  """
  returns the branch indices of the supertree between this couplet
  """
  def _GetBranchArrayIdxList(self):
    return self.Branch_Array_Idx_List
    
  """
  adds index of one particular branch between this couplet
  """
  def _AddBranchArrayIdx(self, idx):
    if idx not in self.Branch_Array_Idx_List:
      self.Branch_Array_Idx_List.append(idx)
          
  """
  returns the average distance between current taxa pair
  with respect to all the supporting input trees
  this is the weighted average distance
  where weights of individual trees supporting this couplet are considered
  """
  def _GetAvgDistMatVal(self):
    num = 0
    denom = 0
    for i in range(len(self.DistMatValues)):
      num = num + self.DistMatValues[i] * Matrix_Weight_Val[self.Supported_Tree_Idx_List[i]]
      denom = denom + Matrix_Weight_Val[self.Supported_Tree_Idx_List[i]]
    return (num * 1.0) / denom
      
  """ 
  this function adds branch length distance information between a couplet
  with respect to a particular input tree and their MRCA node with respect to that tree
  """
  def _Add_Edge_Distance(self, node1_dist_from_mrca_node, node2_dist_from_mrca_node):
    self._AddDistMatValue(node1_dist_from_mrca_node + node2_dist_from_mrca_node)
        
    