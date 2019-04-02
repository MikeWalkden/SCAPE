This folder contains binary files compiled for SCAPE Version 1.21 together with some 
standard libraries in the "libs" folder.

These are compiled with gfortran version 4.8.1 from mingw32.org:
===============================================================

CREATE_FILES.exe - the supporting application that creates the model.omi and model.fcs files 
required to describe a SCAPE model component to FluidEarth

ENGINE.dll - the SCAPE engine

RUN_SCAPE.exe - the free-standing SCAPE application
 
*** 
Note that this has been deleted (2nd April 2019) - the source code has been updated and I don't have access to a GFortran compiler
*** 


This is compiled with Visual Studio Express 2013:
================================================

Engine_Wrapper.dll - the SCAPE engine wrapper required by FluidEarth

The libs Folder
===============

Contains a selection of libraries taken from a gfortran version 4.8.1 installation. 
These are provided for convenience and if to be used must be copied to a folder listed 
on the Windows PATH.

Reccomendation is that you do not use these libraries but instead install gfortran from the 
mingw32.org website.
