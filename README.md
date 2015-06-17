# CSTBL
Branch Length assignment to a Custom input SuperTree topology, given a set of input weighted (with branch lengths) phylogenetic trees 

Description
------------

CSTBL (Custom Super Tree with Branch Length) is a program to assign the branch length of a given input supertree topology, from a set of input weighted phylogenetic trees (carrying distance information).

Input
-----

1) A set of input phylogenetic trees with overlapping taxa set and branch length (distance) information.

2) A given supertree topology, pre-computed using any existing unweighted supertree method, such as MRP, RFS, SuperFine, MinFlip, ........  

Note: The supertree should cover all the input taxa from the source trees.

Output
------

A weighted supertree whose topology is exactly the specified supertree topology, and whose branch lengths are assigned based on the source weighted phylogenetic trees. The branch length assignment is carried out in such a fashion that the difference between the output distance matrix (generated from the derived weighted supertree) and the input distance matrices (generated from the input phylogenetic trees) is as small as possible.

Dependencies / Installation Requirements
--------------

CSTBL is developed in Linux Systems (Ubuntu 14.04), using Python 2.7.

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default)

Note: We have not tested the code on Python 3. Any user having Python 3 environment need to check the correct execution of our code, and optionally needs to upgrade it accordingly.

We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ )

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. We did not upgrade the code for Dendropy 4.0 support, so any user having this new version of Dendropy might need to check the functionalities of COSPEDBTree and possibly upgrade / replace / edit few dendrop[y related functions. So, we recommend users to use the earlier version of Dendropy, to avoid any conflict.

Support for Dendropy 4 and corresponding update of code will be done in a future release.


UBUNTU version issues
------------------

For systems having Ubuntu with lower versions (lower than 14.04), please notify in case of any errors due to OS environments.

Note: We do not support development version corresponding to Windows XP and MacOS, although that will be done in some future release.

Execution
------------

CSTBL is to be executed with the following command line options, from a terminal: (assuming the present working directory contains the source codes)

chmod +x CSTBL.py (To change its permission to make it an executable file)

./CSTBL.py [options]

Details of the options are mentioned below:

-h, --help show this help message and exit

  -I INP_FILENAME, --INPFILE=INP_FILENAME
  
                      name of the input file containing input trees
  
  -p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT
                        
                        1 - input file (containing the input treelist) format
                        is NEWICK (default)                         
                        2 - input file format is NEXUS
  
  -T TOPOLOGY_INPUT_TREE_FILE, --topology=TOPOLOGY_INPUT_TREE_FILE
                        
                        user can provide a custom unweighted supertree
                        topology as an input, using this option.
                        the custom unweighted supertree is built using some
                        other suprtree method such as RFS.
                        The supertree is built using the source tree list
                        provided with the following -I option
                        from the input unweighted supertree, only the branch
                        lengths are computed. The
                        weighted supertree is used for performance comparison
                        here, no supertree computation is performed.
                        Only the branch length assignment is performed
  
  -t TOPOLOGY_FILE_FORMAT, --topform=TOPOLOGY_FILE_FORMAT
                        
                        1 - if the tree topology is provided using the -T
                        option, it is provided in NEWICK (default) formatted file                       
                        2 - file format of the topology input file is NEXUS

Example of a command (followed for the results published in the manuscript)

./COSPEDBTree -I source_tree_input.txt -p1-T supertree_topology_file.txt -t1 > out.txt

command descriptions:

1) -I specifies the input filename

2) source_tree_input.txt : contains the input collection of trees

3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, as specified by the option (-p1) (1 stands for newick)

4) -T option is used to specify the custom supertree topology, built from the trees in the file source_tree_input.txt

5) supertree_topology_file.txt contains the supertree topology in either newick (preferable) or nexus format

6) -t option is analogous to -p option, to specify the format of supertree topology file (1 = newick, 2 = nexus)

The output texts are printed at console. User can redirect the output results to any standard text file by using standard redirection operation (>). For example, in the above command, all the detailed results (textual descriptions) are redirected to file out.txt.

In addition, one output folder 'CUSTOM_SUPERTREE_QP' is created in the current directory containing the source phylogenetic trees. The folder contains the custom input unweighted supertree, its weighted version (computed using this package), and text files as outputs of nonlinear programming functions employed in the current package.


Utilities
-----------

CSTBL requires O(MN^2) time and O(N^2) space complexity, for N input taxa and M input trees. 

For any queries, please contact
-------------------------------

Sourya Bhattacharyya 

Department of Computer Science and Engineering 

Indian Institute of Technology Kharagpur 

email: sourya.bhatta@gmail.com


