#!/usr/bin/env python


##---------------------------------------------
''' 
this program is used to generate a supertree (consensus) from a set of constituent trees
input trees contain the branch length information as well
output needs to conform to the input source trees
so as to get the minimal least square error with respect to the source trees
the input is multiple source trees
there may be conflicts among the input tree - we have to select the consensus

Author: Sourya Bhattacharyya
Dept of CSE, IIT Kharagpur
V1.0 - 07.04.2014 - basic code
V2.0 - 01.07.2014 - modified the clustering and conflict detection routine
V3.0 - 23.03.2015 - reduced storage complexity and cleaned the code
V4.0 - 17.06.2015 - reduced to only estimate branch lengths from a given supertree topology
		  - also reduced the code to use only GNU BFGS method
		  - no other scipy library based optimization is used
'''

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
			  help="name of the input file containing input trees")
  
  parser.add_option("-p", "--inpform", \
			  type="int", \
			  action="store", \
			  dest="inp_file_format", \
			  default=1, \
			  help="1 - input file (containing the input treelist) format is NEWICK (default) \
			  2 - input file format is NEXUS")
  
  parser.add_option("-T", "--topology", \
			  type="string", \
			  action="store", \
			  dest="topology_input_tree_file", \
			  default="", \
			  help="user can provide a custom unweighted supertree topology as an input, using this option. \
			  the custom unweighted supertree is built using some other suprtree method (other than COSPEDTREE) such as RFS. \
			  The supertree is built using the source tree list provided with the following -I option \
			  from the input unweighted supertree, only the branch lengths are computed \
			  the weighted supertree is used for performance comparison \
			  here, no supertree computation is performed. \
			  Only the branch length assignment is performed")
			  
  parser.add_option("-t", "--topform", \
			  type="int", \
			  action="store", \
			  dest="topology_file_format", \
			  default=1, \
			  help="1 - if the tree topology is provided using the -T option, it is provided in NEWICK (default) formatted file \
			  2 - file format of the topology input file is NEXUS")			  
    			        
  opts, args = parser.parse_args()
  return opts, args
  
  
##-----------------------------------------------------
''' main function '''
def main():  
  opts, args = parse_options()
  
  ROOTED_TREE = False	#opts.default_rooted
  PRESERVE_UNDERSCORE = True	#opts.preserve_underscores
  if (opts.inp_file_format == 1):
    INPUT_FILE_FORMAT = 'newick'
  else:
    INPUT_FILE_FORMAT = 'nexus'
  INPUT_FILENAME = opts.INP_FILENAME
  NO_OF_LOOPS = 15
  METHOD_OF_QP = 3
  TOPOLOGY_INPUT_TREE_FILENAME = opts.topology_input_tree_file
  if (opts.topology_file_format == 1):
    TOPOLOGY_FILE_FORMAT = 'newick'
  else:
    TOPOLOGY_FILE_FORMAT = 'nexus'
  
  global Output_Text_File
  
  if (INPUT_FILENAME == ""):
    print '******** THERE IS NO INPUT FILE (CONTAINING THE SOURCE TREES) SPECIFIED - RETURN **********'
    return
  else:
    print 'input filename containing the source trees: ', INPUT_FILENAME
    
  # according to the location of input filename
  # adjust the locations of the output files as well
  k = INPUT_FILENAME.rfind("/")
  if (k == -1):
    dir_of_inp_file = './'
  else:
    dir_of_inp_file = INPUT_FILENAME[:(k+1)]
  if (DEBUG_LEVEL > 1):
    print 'dir_of_inp_file: ', dir_of_inp_file
  
  if (TOPOLOGY_INPUT_TREE_FILENAME == ""):
    print '******** NO INPUT FILE CONTAINING THE CUSTOM SUPERTREE TOPOLOGY IS PROVIDED - RETURN **********'
    return
  else:
    dir_of_curr_exec = dir_of_inp_file + 'CUSTOM_SUPERTREE_QP'   

  # create the directory
  if (os.path.isdir(dir_of_curr_exec) == False):
    mkdr_cmd = 'mkdir ' + dir_of_curr_exec
    os.system(mkdr_cmd)         
    
  # append the current output directory in the text file
  Output_Text_File = dir_of_curr_exec + '/' + Output_Text_File
    
  fp = open(Output_Text_File, 'w')    
    
  # this variable notes the count of input source trees  
  tree_count = 0
    
  # note the program beginning time 
  start_timestamp = time.time()
    
  #-------------------------------------  
  """ read the source trees collection and store it in a tree collection structure
  individual elements of this collection is thus a source tree """
  Input_Treelist = Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)  
  
  #-------------------------------------    
  # from the input source trees, note the number of taxa (total)
  # and also define the class instances corresponding to single taxa
  for tr_idx in range(len(Input_Treelist)):
    tree_count = tree_count + 1
    taxa_labels_curr_tree = Input_Treelist[tr_idx].infer_taxa().labels()
    if (DEBUG_LEVEL > 1):
      fp.write('\n Tree no : ' + str(tree_count) +  'no of leaf nodes: ' + str(len(taxa_labels_curr_tree)))
    if (DEBUG_LEVEL > 2):
      fp.write('\n taxa set belonging to current tree: ' + str(taxa_labels_curr_tree))
    for i in range(len(taxa_labels_curr_tree)):
      if taxa_labels_curr_tree[i] not in COMPLETE_INPUT_TAXA_LIST:
	COMPLETE_INPUT_TAXA_LIST.append(taxa_labels_curr_tree[i])
    
  # now process individual trees to find the couplet relations of those trees
  for tr_idx in range(len(Input_Treelist)):
    DeriveCoupletRelations(Input_Treelist[tr_idx], tr_idx)
   
  number_of_taxa = len(COMPLETE_INPUT_TAXA_LIST)
    
  if (DEBUG_LEVEL >= 0):
    fp.write('\n  total no of taxa: ' + str(number_of_taxa))
  if (DEBUG_LEVEL > 1):
    #fp.write('\n len Taxa_Info_Dict: ' + str(len(Taxa_Info_Dict)))
    fp.write('\n len COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST))
    fp.write('\n len TaxaPair_Reln_Dict : ' + str(len(TaxaPair_Reln_Dict)))
  
  fp.close()
  
  # read the custom tree topology from the specified input file
  Final_Supertree = Read_Input_Tree(ROOTED_TREE, PRESERVE_UNDERSCORE, TOPOLOGY_FILE_FORMAT, TOPOLOGY_INPUT_TREE_FILENAME)    
  suptr_taxa_labels = Final_Supertree.infer_taxa().labels()
  Final_Supertree.retain_taxa_with_labels(suptr_taxa_labels, True, True)
  fp = open(Output_Text_File, 'a')
  fp.write('\n\n --- after deleteting the internal nodes with degree 1')
  fp.write('\n\n ---output tree  without branch length information (in newick format): ' + Final_Supertree.as_newick_string())
  fp.write('\n\n ---with branch length information - plotting : ')
  fp.write(Final_Supertree.as_ascii_plot())
  fp.close()

  # this function assigns the branch length information on the generated supertree
  AssignBranchLen(Final_Supertree, Input_Treelist, Output_Text_File)
  
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
  
