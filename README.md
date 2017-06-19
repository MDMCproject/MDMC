# MDMC

This software uses the XML reader/writer library called xmlf90, a copy of which
is located in src/xmlf90-1.2. For more information about this library see 
https://github.com/rscircus/xmlf90 (xmlf90 is 3-clause BSD license which is compatible with GNU GPL).

To compile
----------

* Intel Microsoft Visual studio: Open mdmc.sln in /src and compile.

* Other platforms: in /src find both a traditional makefile and SCons makefile scripts (named SConstruct and SConscript) are available (possibly a bit out of date).

To run
------

The source-code directory src contains the the sub-directory input which contains examples 
of input files in the form of MDMC job files (XML formatted files). One of these is
./input/md_control.xml which should be fast to execute and simply just run a small MD
simulation. 

From the command line type mdmc.exe (mdmc on linux). You are then asked to name a MDMC
job file. The program then runs and output all results into a folder 'output' in the directory 
of the executable.

