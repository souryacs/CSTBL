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

NOTE:

All the options except the first three, signify toggle / complement of their corresponding DEFAULT values. First option (help) displays these command line parameters.

It Is Preferable For A Beginner, To Not Use Any Option Other Than The Second And Third Options. Second option is for specifying the input filename (mandatory) Third option is for specifying the corresponding file format.

Details of the options are mentioned below:

-h, --help show this help message and exit


