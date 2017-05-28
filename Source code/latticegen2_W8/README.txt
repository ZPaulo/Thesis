The code has been developed with "CODEBLOCKS" IDE.

step 1 - copy .vrt and .bnd files in the folder. The geometry must be drawn in Pointwise. It is necessary to set the boundary conditions and mesh with a sufficient number of points on the curved boundaries because the connectors are straight lines between points of the geometries.

step 2 - set names of the files to open in GeomRead.cpp

step 3 - set number of nodes in x-direction in MeshGen.cpp

step 4 - Set seed point in FluidTracking.cpp. XSeed and YSeed dimensions are the real x and y of the geometry in Pointwise.

WARNING! Once happened to me that the code gave me error. In this case scaling (increasing the dimension) the geometry in Pointwise solved the problem.