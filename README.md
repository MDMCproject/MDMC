# MDMC

This software was created to test a new algorithm for fitting/optimising Potential Energy (PE) parameters against dynamical information such as S(q,omega) and S(q,t).

A zip file containing Windows executable, user manual, scripts for plotting, example MDMC job files and Argon test data are available for download from https://github.com/MDMCproject/MDMC/releases.

Documentation related to and referenced to in the code can be found in source code directory doc.

Library used by this software
-----------------------------

This software uses the XML reader/writer library called xmlf90, a copy of which
is located in src/xmlf90-1.2. For more information about this library see 
https://github.com/rscircus/xmlf90 (xmlf90 is 3-clause BSD license which is compatible with GNU GPL).

To compile source code
----------------------

* Intel Microsoft Visual studio: Open mdmc.sln in /src and compile.

* Other platforms: in /src find both a traditional makefile and SCons makefile scripts (named SConstruct and SConscript) is available (possibly a bit out of date).
