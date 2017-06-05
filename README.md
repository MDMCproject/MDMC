# MDMC

This software uses the XML reader/writer library called xmlf90. The version of this library used with this software can be found in src/xmlf90-1.2. For more information about this library see https://github.com/rscircus/xmlf90.

To compile
----------

* Intel Microsoft Visual studio: Open mdmc.vfproj in /src and compile.

* Other platforms: in /src find both a traditional makefile and SCons makefile scripts (named SConstruct and SConscript) are available (possibly a bit out of date).

To run
------

A bit primitive, although this would be simple to improve. The /src directory contains the two sub-directories /input and /output. The simulation that is executed is controlled by an XML file. The particular XML file which is loaded is hard coded in mdmc.f90, as of this writing, with the line "./input/mdmc_control.xml". 

* To run a different simulation either change mdmc_control.xml or chose one of the other XML simulation files in /src/input.
* Then from Visual Studio press Crtl+F5, or from a command prompt type mdmc or mdmc.exe.
