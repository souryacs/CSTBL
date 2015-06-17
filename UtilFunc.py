#!/usr/bin/env python

import Header
from Header import *
                     
#--------------------------------------------------------
# this function defines relationship between a pair of nodes in a tree
# the relationship is either ancestor / descendant, or siblings, or no relationship 
# parameters: 1) node1 and node2 are two nodes
# node1_dist_from_mrca_node is the distance of node 1 from MRCA
# node2_dist_from_mrca_node is the distance of node 2 from MRCA
def DefineLeafPairReln(node1, node2, node1_dist_from_mrca_node, node2_dist_from_mrca_node, tree_idx):
  key1 = (node1.taxon.label, node2.taxon.label)
  key2 = (node2.taxon.label, node1.taxon.label)
  if key1 in TaxaPair_Reln_Dict:
    TaxaPair_Reln_Dict[key1]._AddSupportTreeIndex(tree_idx)
    TaxaPair_Reln_Dict[key1]._Add_Edge_Distance(node1_dist_from_mrca_node, node2_dist_from_mrca_node)
  elif key2 in TaxaPair_Reln_Dict:
    TaxaPair_Reln_Dict[key2]._AddSupportTreeIndex(tree_idx)
    TaxaPair_Reln_Dict[key2]._Add_Edge_Distance(node1_dist_from_mrca_node, node2_dist_from_mrca_node)
  else:
    TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
    TaxaPair_Reln_Dict[key1]._AddSupportTreeIndex(tree_idx)
    TaxaPair_Reln_Dict[key1]._Add_Edge_Distance(node1_dist_from_mrca_node, node2_dist_from_mrca_node)
      
  return      
      
#--------------------------------------------------------
# this function derives coupket relations belonging to one tree
# that is provided as an input argument to this function
def DeriveCoupletRelations(Curr_tree, tree_idx):
    
  # traverse the internal nodes of the tree in postorder fashion
  for curr_node in Curr_tree.postorder_internal_node_iter():
    # distance of this node from the root node
    curr_node_dist_from_root = curr_node.distance_from_root()

    # list the leaf and internal children of the current node
    curr_node_child_leaf_nodes = []
    curr_node_child_internal_nodes = []
    for x in curr_node.child_nodes():
      if (x.is_leaf() == True):
	curr_node_child_leaf_nodes.append(x)
      else:
	curr_node_child_internal_nodes.append(x)
    
    # pair of leaf nodes will be related by sibling relations
    if (len(curr_node_child_leaf_nodes) > 1):
      for i in range(len(curr_node_child_leaf_nodes) - 1):
	node1_dist_from_mrca_node = curr_node_child_leaf_nodes[i].distance_from_root() - curr_node_dist_from_root
	for j in range(i+1, len(curr_node_child_leaf_nodes)):
	  node2_dist_from_mrca_node = curr_node_child_leaf_nodes[j].distance_from_root() - curr_node_dist_from_root
	  DefineLeafPairReln(curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], \
	    node1_dist_from_mrca_node, node2_dist_from_mrca_node, tree_idx)
    
    # one leaf node (direct descendant) and another leaf node (under one internal node)
    # will be related by ancestor / descendant relations
    if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
      for p in curr_node_child_leaf_nodes:
	node1_dist_from_mrca_node = p.distance_from_root() - curr_node_dist_from_root
	for q in curr_node_child_internal_nodes:
	  for r in q.leaf_nodes():
	    node2_dist_from_mrca_node = r.distance_from_root() - curr_node_dist_from_root
	    DefineLeafPairReln(p, r, node1_dist_from_mrca_node, node2_dist_from_mrca_node, tree_idx)
    
    # finally a pair of leaf nodes which are descendant of internal nodes will be related by NO_EDGE relation
    if (len(curr_node_child_internal_nodes) > 1):
      for i in range(len(curr_node_child_internal_nodes) - 1):
	for j in range(i+1, len(curr_node_child_internal_nodes)):
	  for p in curr_node_child_internal_nodes[i].leaf_nodes():
	    node1_dist_from_mrca_node = p.distance_from_root() - curr_node_dist_from_root
	    for q in curr_node_child_internal_nodes[j].leaf_nodes():
	      node2_dist_from_mrca_node = q.distance_from_root() - curr_node_dist_from_root
	      DefineLeafPairReln(p, q, node1_dist_from_mrca_node, node2_dist_from_mrca_node, tree_idx)

##-----------------------------------------------------
# this function reads the input tree list file
# parameters: ROOTED_TREE - whether the treelist to be read as rooted format
# PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
# INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
# INPUT_FILENAME: file containing the input treelist

def Read_Input_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
  Inp_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, schema=INPUT_FILE_FORMAT, \
						  preserve_underscores=PRESERVE_UNDERSCORE, \
						  default_as_rooted=ROOTED_TREE)
  
  return Inp_TreeList
      
      
##-----------------------------------------------------
# this function reads an input tree from a specified file
# parameters: ROOTED_TREE - whether the treelist to be read as rooted format
# PRESERVE_UNDERSCORE: whether underscores of the taxa name will be preserved or not
# INPUT_FILE_FORMAT: data is read from the file according to NEWICK or NEXUS format
# INPUT_FILENAME: file containing the input tree

def Read_Input_Tree(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
  Inp_Tree = dendropy.Tree.get_from_path(INPUT_FILENAME, schema=INPUT_FILE_FORMAT, \
						  preserve_underscores=PRESERVE_UNDERSCORE, \
						  default_as_rooted=ROOTED_TREE)
  
  return Inp_Tree
      
      
      
