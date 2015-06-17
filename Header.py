#!/usr/bin/env python

import dendropy
from dendropy import TreeList, Tree, Taxon, TaxonSet
import time
import os
import sys
from optparse import OptionParser
import math

# this is the path of QP executable based on GNU C library
QP_Executable = './GNU_BFGS2'

""" this dictionary defines the taxa pair relations
each entry is indexed by two nodes """
TaxaPair_Reln_Dict = dict()

""" this list contains the complete set of taxa present in the input source trees """
COMPLETE_INPUT_TAXA_LIST = []

# this is the debug level
# set for printing the necessary information
DEBUG_LEVEL = 0

# this text file stores all the printing output
Output_Text_File = 'complete_output_description.txt'

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
""" this class defines the connectivity relationship between a pair of taxa
initially the information are obtained from the input source trees
later the contents of these class instances are modified according to the generation of the consensus tree
key of this class --- taxa1, taxa2  
in the class, the edge type signifies the relation between a pair of taxa """
class Reln_TaxaPair(object):
  def __init__(self):    
            
    # for a particular taxa pair, when the tree topology is established, this array stires the branches 
    # those are in between this taxa pair
    self.Branch_Array_Idx_List = []
    
    # this list stores the tree index where the taxa pair is supported
    self.Supported_Tree_Idx_List = []

    # this list stores the distance values between these taxa pairs
    # as supported by the source trees
    self.DistMatValues = []
            
  # these functions are for computing the QP 
  def _AddDistMatValue(self, val):
    self.DistMatValues.append(val)
    
  def _GetNoSupportTrees(self):
    return len(self.Supported_Tree_Idx_List)
    
  def _GetBranchArrayIdxList(self):
    return self.Branch_Array_Idx_List
    
  def _AddBranchArrayIdx(self, idx):
    if idx not in self.Branch_Array_Idx_List:
      self.Branch_Array_Idx_List.append(idx)
      
  def _AddSupportTreeIndex(self, idx):
    self.Supported_Tree_Idx_List.append(idx)
        
  def _GetSupportTreeList(self):
    return self.Supported_Tree_Idx_List
    
  # this function returns the average value of distance between current taxa pair
  def _GetAvgDistMatVal(self):
    # comment - sourya
    #return (self.dist_node1_mrca_global_mean + self.dist_node2_mrca_global_mean)
    # add - sourya
    num = 0
    denom = 0
    for i in range(len(self.DistMatValues)):
      num = num + self.DistMatValues[i] * Matrix_Weight_Val[self.Supported_Tree_Idx_List[i]]
      denom = denom + Matrix_Weight_Val[self.Supported_Tree_Idx_List[i]]
    return (num * 1.0) / denom
      
  # this function adds one edge count (with a given input edge type)
  def _Add_Edge_Distance(self, node1_dist_from_mrca_node, node2_dist_from_mrca_node):
    # add the relative branch length timings (with respect to MRCA)
    self._AddDistMatValue(node1_dist_from_mrca_node + node2_dist_from_mrca_node)
        
    
