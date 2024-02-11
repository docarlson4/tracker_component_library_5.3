David Vallado's original distribution came with implementations in C++,
Fortran, Microsoft Excel, Matlab, and Pascal. Only the C++ implementation
is retained here and is called as an external library by a mex function to
run the SGP 4 propagator. The code has been slightly modified to eliminate
all warnings when compiled with the -Weverything option under Clang, except
for the Wcovered-switch-default warning. For example, a few unused
variables were eliminated. Also, the debugging code and the IO code as well
as files containing helper functions for the IO code have been removed as
they is not used in the Tracker Library. The definition of pi had to be
moved from sgp4unit.h to sgp4unit.cpp to avoid a namespace conflict with
one of Matlab's headers. The documentation files appearing in AIAA journals
have been removed. They have been superseded by the publication:
D. A. Vallado, P. Crawford, R. Hujsak, and T. S. Kelso "Revisiting
Spacetrack Report #3," AIAA/AAS Astrodynamics Specialist Conference and
Exhibit, Guidance, Navigation, and Control and Co-located Conferences,
Keystone Colorado, 21-24 August 2006.

December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
