#!/usr/bin/env python


##---------------------------------------------
""" 
this program is used to generate a supertree (consensus) from a set of constituent trees
input trees contain the branch length information as well
output needs to conform to the input source trees
so as to get the minimal least square error with respect to the source trees
the input is multiple source trees
there may be conflicts among the input tree - we have to select the consensus
"""

"""
Author: Sourya Bhattacharyya
Dept of CSE, IIT Kharagpur
V1.0 - 07.04.2014 - basic code
V2.0 - 01.07.2014 - modified the clustering and conflict detection routine
V3.0 - 23.03.2015 - reduced storage complexity and cleaned the code
V4.0 - 17.06.2015 - reduced to only estimate branch lengths from a given supertree topology
		  - also reduced the code to use only GNU BFGS method
		  - no other scipy library based optimization is used
V5.0 - 06.09.2015 - added GSL based QP optimization executable path in the main command line option
		  - code clean and comment
"""

## Copyright 2014, 2015 Sourya Bhattacharyya and Jayanta Mukherjee.
## All rights reserved.
##
## See "LICENSE.txt" for terms and conditions of usage.
##
##---------------------------------------------

import Header
from Header import *
import UtilFunc
from UtilFunc import *
import Edge_Len_Adjust
from Edge_Len_Adjust import *

##-----------------------------------------------------
# this function is useful to parse various options for input data processing
def parse_options():  
  parser = OptionParser()
  
  parser.add_option("-I", "--INPFILE", \
			  type="string", \
			  action="store", \
			  dest="INP_FILENAME", \
			  default="", \
			  help="name of the file containing input phylogenetic trees")
  
  parser.add_option("-p", "--inpform", \
			  type="int", \
			  action="store", \
			  dest="inp_file_format", \
			  default=1, \
			  help="1 - format of the input file (containing the input treelist) is NEWICK (default) \
			  2 - input file format is NEXUS")
  
  parser.add_option("-T", "--topology", \
			  type="string", \
			  action="store", \
			  dest="topology_input_tree_file", \
			  default="", \
			  help="File containing custom unweighted supertree topology, such as RFS. \
			  The supertree is built using the input trees (provided with -I option above)")
			  
  parser.add_option("-t", "--topform", \
			  type="int", \
			  action="store", \
			  dest="topology_file_format", \
			  default=1, \
			  help="1 - Format of the custom unweighted supertree file is NEWICK (default) \
			  2 - file format of the topology input file is NEXUS")

  parser.add_option("-Q", "--QPExec", \
			  type="string", \
			  action="store", \
			  dest="QP_Exec_Path", \
			  default="", \
			  help="Absolute path of the executable for QP solver")
    			        
  opts, args = parser.parse_args()
  return opts, args
  
  
##-----------------------------------------------------
''' main function '''
def main():  
  opts, args = parse_options()
  
  ROOTED_TREE = True	#opts.default_rooted
  PRESERVE_UNDERSCORE = True	#opts.preserve_underscores
  if (opts.inp_file_format == 1):
    INPUT_FILE_FORMAT = 'newick'
  else:
    INPUT_FILE_FORMAT = 'nexus'
  INPUT_FILENAME = opts.INP_FILENAME
  TOPOLOGY_INPUT_TREE_FILENAME = opts.topology_input_tree_file
  if (opts.topology_file_format == 1):
    TOPOLOGY_FILE_FORMAT = 'newick'
  else:
    TOPOLOGY_FILE_FORMAT = 'nexus'
  NO_OF_LOOPS = 15
  METHOD_OF_QP = 3
  """
  abspath function converts input possibly relative path into an absoloute path
  """
  if (opts.QP_Exec_Path == ""):
    print '******** THERE IS NO PATH FOR QP SOLVER (GNU_BFGS2) IS PROVIDED - RETURN **********'
    return
  QP_EXEC_PATH = os.path.abspath(opts.QP_Exec_Path)
  
  if (INPUT_FILENAME == ""):
    print '******** THERE IS NO INPUT FILE (CONTAINING THE SOURCE TREES) SPECIFIED - RETURN **********'
    return
  print 'input filename containing the source trees: ', INPUT_FILENAME
  
  if (TOPOLOGY_INPUT_TREE_FILENAME == ""):
    print '******** THERE IS NO CUSTOM SUPERTREE TOPOLOGY FILE SPECIFIED - RETURN **********'
    return
  print 'file containing the custom supertree topology: ', TOPOLOGY_INPUT_TREE_FILENAME
  
  """
  obtain the directory of input treelist
  """
  k = INPUT_FILENAME.rfind("/")
  if (k == -1):
    dir_of_inp_file = './'
  else:
    dir_of_inp_file = INPUT_FILENAME[:(k+1)]
  if (DEBUG_LEVEL > 1):
    print 'dir_of_inp_file: ', dir_of_inp_file
  
  """
  obtain the directory of custom supertree topology file
  """
  k = TOPOLOGY_INPUT_TREE_FILENAME.rfind("/")
  if (k == -1):
    dir_of_topology_inp_file = './'
  else:
    dir_of_topology_inp_file = TOPOLOGY_INPUT_TREE_FILENAME[:(k+1)]
  if (DEBUG_LEVEL > 1):
    print 'dir_of_topology_inp_file: ', dir_of_topology_inp_file  
    
  """
  the output weighted supertree will be placed within a directory placed 
  under the directory containing custom supertree topology
  """
  dir_of_curr_exec = dir_of_topology_inp_file + 'CUSTOM_SUPERTREE_QP'   
  # create the directory
  if (os.path.isdir(dir_of_curr_exec) == False):
    mkdr_cmd = 'mkdir ' + dir_of_curr_exec
    os.system(mkdr_cmd)         
  # append the current output directory in the text file
  Output_Text_File = dir_of_curr_exec + '/' + 'Complete_Output_Description.txt'

  print 'Output_Text_File: ', Output_Text_File

  print 'QP_EXEC_PATH: ', QP_EXEC_PATH

  # note the program beginning time 
  start_timestamp = time.time()
    
  #-------------------------------------  
  """ 
  read the input treelist file
  """
  Input_Treelist = Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)  

  fp = open(Output_Text_File, 'w')      
  #-------------------------------------    
  """
  from the input trees, note the number of taxa (total)
  """
  for tr_idx in range(len(Input_Treelist)):
    taxa_labels_curr_tree = Input_Treelist[tr_idx].infer_taxa().labels()
    if (DEBUG_LEVEL > 1):
      fp.write('\n Tree no : ' + str(tr_idx+1) +  'no of leaf nodes: ' + str(len(taxa_labels_curr_tree)))
    if (DEBUG_LEVEL > 2):
      fp.write('\n taxa set belonging to current tree: ' + str(taxa_labels_curr_tree))
    for i in range(len(taxa_labels_curr_tree)):
      if taxa_labels_curr_tree[i] not in COMPLETE_INPUT_TAXA_LIST:
	COMPLETE_INPUT_TAXA_LIST.append(taxa_labels_curr_tree[i])
  
  number_of_taxa = len(COMPLETE_INPUT_TAXA_LIST)
  
  """ 
  now process individual trees to find the couplet relations within those trees
  """
  for tr_idx in range(len(Input_Treelist)):
    DeriveCoupletRelations(Input_Treelist[tr_idx], tr_idx)
      
  if (DEBUG_LEVEL >= 0):
    fp.write('\n  total no of taxa: ' + str(number_of_taxa))
  if (DEBUG_LEVEL > 1):
    fp.write('\n len COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST))
    fp.write('\n len TaxaPair_Reln_Dict : ' + str(len(TaxaPair_Reln_Dict)))
  
  # close the output file
  fp.close()
  
  """ 
  read the custom supertree topology from the specified input custom topology file
  """
  Final_Supertree = Read_Input_Tree(ROOTED_TREE, PRESERVE_UNDERSCORE, TOPOLOGY_FILE_FORMAT, TOPOLOGY_INPUT_TREE_FILENAME)    
  fp = open(Output_Text_File, 'a')
  fp.write('\n\n ---output supertree  without branch length information (in newick format): ' + Final_Supertree.as_newick_string())
  fp.close()

  # this function assigns the branch length information on the generated supertree
  AssignBranchLen(Final_Supertree, Input_Treelist, QP_EXEC_PATH, Output_Text_File)
  
  out_treefilename = dir_of_curr_exec + '/' + 'CUSTOM_SUPERTREE_with_branch_length_newick.tre'
  outfile = open(out_treefilename, 'w')
  outfile.write(Final_Supertree.as_newick_string())
  outfile.close()
  
  # note the timestamp
  # this will signify the time required for tree reading and couplet feature extraction
  end_timestamp = time.time()  
  fp = open(Output_Text_File, 'a')
  fp.write('\n \n\n ===============>>>>>>>>>>>>>>> TIME COMPLEXITY : complete method execution: ' + str(end_timestamp - start_timestamp))
  fp.close()
          
  #--------------------------------------------------------------  
  # delete the storage variables associated with the current execution 
  
  # clear the dictionaries
  TaxaPair_Reln_Dict.clear()
  EdgeInfoDict.clear()
  
  # clear the lists associated
  if (len(COMPLETE_INPUT_TAXA_LIST) > 0):
    COMPLETE_INPUT_TAXA_LIST[:] = []
  if (len(Matrix_Weight_Val) > 0):
    Matrix_Weight_Val[:] = []
      
#-----------------------------------------------------
if __name__ == "__main__":
    main() 
  
