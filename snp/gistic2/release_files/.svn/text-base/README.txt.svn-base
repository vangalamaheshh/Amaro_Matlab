
This GISTIC2 release has the following directory structure.

<gistic_install_directory>                        where the archive was unpacked
|-- README.txt					  this file
|-- GISTICDocumentation_110324_standalone.htm     HTML documentation
|-- GISTICDocumentation_110324_standalone_files   documentation images
|   `-- {16 image files}
|-- INSTALL.txt                                   how install the compiled program
|-- LICENSE.txt					  licensing information
|-- MATLAB_Component_Runtime {empty}              suggested location of MCR
|-- MCRInstaller.bin                              MCR v7.14 intstallation program
|-- example_results {empty}                       output directory for example             
|-- gp_gistic2_from_seg                           compiled GISTIC 2.0 executable
|-- run_gistic_example                            shell script to run the example
|-- examplefiles                                  example input files
|   |-- arraylistfile.txt
|   |-- cnvfile.txt
|   |-- markersfile.txt
|   `-- segmentationfile.txt
|-- refgenefiles                                  various reference genomes
|   |-- hg16.mat
|   |-- hg17.mat
|   `-- hg18.mat
`-- source					  MATLAB source code
    |-- {195 *.m GISTIC 2 class files}
    `-- @SegArray				  source code for the SegArray class
        `-- {95 *.m SegArray class files}


This release contains both the compiled executable and the source code from 
which it was compiled, as well as documentation and an example. 

Please read LICENSE.txt to understand the terms under which GISTIC 2.0
is licensed to end users.

The executable is the gene pattern module gp_gistic2_from_seg. This
module depends on the libraries in the MATLAB Component Runtime (MCR)
version 7.14 (from R2010b). If your Linux system does not already have
the MCR installed, you may do so by running the provided
MCRInstaller.bin from The MathWorks.  The example script,
run_gistic_example, assumes that the MCR will be installed in the
MATLAB_Component_Runtime subdirectory. If it is not, you must edit the
example. See INSTALL.txt for details.

The source code is in the 'source' subdirectory. The
<gistic_install_directory>/source directory must be added to the
matlab path before running or compiling the code. The top level module
for the executable program provided is gp_gistic2_from_seg.m.

The minimum requirement for MATLAB is R2009B (7.14).  The GISTIC
software was developed on x86-64 Linux systems and has not been tested
on other MATLAB platforms.

HTML documentation for the input parameters and important output files
is available by browsing GISTICDocumentation_110324_standalone.htm. 
The source code is extensively documented and most functions will return 
information in response to a 'help <function>' at the matlab prompt.

This release contains input files to support the run_gistic_example in
the directories examplefiles and refgenefiles.

