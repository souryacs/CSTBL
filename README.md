# CSTBL
-------------

Branch Length assignment to a Custom input SuperTree topology, given a set of input weighted (with branch lengths) phylogenetic trees 

Description
----------------

CSTBL (Custom Super Tree with Branch Length) is a program to assign the branch length of a 
given input supertree topology, from a set of input weighted phylogenetic trees (carrying distance information).

Input
---------

1) A set of input phylogenetic trees with overlapping taxa set and branch length (distance) information.

2) A given supertree topology, pre-computed using any existing unweighted supertree method, such as MRP, RFS, SuperFine, MinFlip, ........

Note: The supertree should cover all the input taxa from the source trees.

Output
----------

A weighted supertree whose topology is exactly the specified supertree topology, and whose branch lengths are assigned 
based on the source weighted phylogenetic trees. The branch length assignment is carried out in such a fashion that 
the difference between the output distance matrix (generated from the derived weighted supertree) and the 
input distance matrices (generated from the input phylogenetic trees) is as small as possible.

Validation
------------

For branch length prediction accuracy, use least square error between the distance matrices and the output supertree
for topological accuracy, use RF distance between the source trees and the output supertree.

Dependencies / Installation Requirements
-----------------------------------------------------

CSTBL is developed in Linux Systems (Ubuntu 14.04), using Python 2.7.

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default)

Note: We have not tested the code on Python 3. Any user having Python 3 environment need to 
check the correct execution of our code, and optionally needs to upgrade it accordingly.

********* We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ )

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. 
We did not upgrade the code for Dendropy 4.0 support, so any user having this new version of 
Dendropy might need to check the functionalities of CSTBL and possibly 
upgrade / replace / edit few dendropy related functions. So, we recommend users to use the 
earlier version of Dendropy, to avoid any conflict.

*********** Support for Dendropy 4 and corresponding update of code will be done in a future release.

3) A binary executable file GNU_BFGS2 is provided (in a zipped archieve) along with this release. 
User needs to Download, extract the archieve and save it in the location containing the source codes.

The package requires GNU scinetific library (GSL) for its execution. The package can be installed in any of the following ways:

A) User needs to go to the link http://www.gnu.org/software/gsl/ and follow the 
download and installation instructions.

B) Otherwise, for UBUNTU 14.04 system users, please check the following link: 
http://askubuntu.com/questions/490465/install-gnu-scientific-library-gsl-on-ubuntu-14-04-via-terminal

C) Or you can use the following command for CentOS or Fedora

yum install gsl-devel

Issues (Mentioned in the following (A) and (B) points)
----------------------------------------------------------------------

A) Related Ubuntu Version
----------------------------------

For systems having Ubuntu with lower versions (lower than 14.04), please notify in case of any errors due to OS environments.

Note: We do not support development version corresponding to Windows XP and MacOS, although that will be done in some future release.

B) Related GSL
--------------------

Sometimes, after installing the GSL, following error is encountered in the system:

"error while loading shared libraries: libgsl.so.0: cannot open shared object file: No such file or directory"

To resolve this error, check out the following link and perform the steps as mentioned.

https://www.gnu.org/software/gsl/manual/html_node/Shared-Libraries.html

Otherwise, locate the file '.bashrc' within your system, and append the following line at the end of it:

LD_LIBRARY_PATH=/usr/local/lib

Execution
---------------

Download the zipped archieve and extract the codes.

First, extract the zipped archieve named GNU_BFGS2.zip, to unpack the corresponding executable. 

There is also a python file named CSTBL.py, which is the main file of this package.

In terminal, go to the directory containing the source codes, and type the following commands:

chmod +x CSTBL.py (To change its permission to make it an executable file)

./CSTBL.py [options]

Details of the options are mentioned below:

-h, --help show this help message and exit

-I INP_FILENAME, --INPFILE=INP_FILENAME

                  name of the input file containing input phylogenetic trees

-p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT

                    1 - input file (containing the input treelist) format is NEWICK (default)                         
                    2 - input file format is NEXUS

-T CUSTOM_SUPERTREE_FILE, --topology=CUSTOM_SUPERTREE_FILE

                    File containing custom unweighted supertree topology, such as RFS.			  
                    The supertree is built using the input trees (provided with -I option above)

-t TOPOLOGY_FILE_FORMAT, --topform=TOPOLOGY_FILE_FORMAT

                    The value can be either 1 or 2.
                    1 - supertree topology (provided using the -T option) is a NEWICK (default) formatted file                       
                    2 - file format of the supertree topology file is NEXUS

-Q QP_Exec_Path , --QPExec QP_Exec_Path

                  Path (absolute / relative) of the executable for QP solver (for branch length assignment). 
                  User needs to provide the path of GNU_BFGS2 executable with this option.

Example of a command 
(followed for the results published in the manuscript)
--------------------------------------------------------------------------------------------------

./CSTBL.py -I source_tree_input.txt -p1-T supertree_topology_file.txt -t1 -Q './GNU_BFGS2' > out.txt

command descriptions:

1) -I specifies the input filename:
  
2)  source_tree_input.txt contains the input collection of trees

3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, 
as specified by the option (-p1) (1 stands for newick)

4) -T option is used to specify the custom supertree topology, built from the trees in the file source_tree_input.txt

5) supertree_topology_file.txt contains the supertree topology in either newick 
(preferable) or nexus format

6) -t option is analogous to -p option, to specify the format of supertree topology file 
(1 = newick, 2 = nexus)

7) GNU_BFGS2 executable location (here we assume that the executable is placed in the current directory)

The output texts are printed at console. User can redirect the output results to any standard text file by 
using standard redirection operation (>). For example, in the above command, all the detailed results 
(textual descriptions) are redirected to file out.txt.

In addition, one output folder 'CUSTOM_SUPERTREE_QP' is created in the current directory 
containing the supertree topology file. The folder contains the custom input unweighted supertree, 
its weighted version (computed using this package), and text files as outputs of nonlinear programming 
functions employed in the current package.

Utilities
-----------

CSTBL requires O(MN^2) time and O(N^2) space complexity, for N input taxa and M input trees.

For any queries, please contact
---------------------------------------

Sourya Bhattacharyya

Department of Computer Science and Engineering

Indian Institute of Technology Kharagpur

email: sourya.bhatta@gmail.com


