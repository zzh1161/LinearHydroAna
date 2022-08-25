# Linear Hydrostatic Analysis
## Document Description
+ `include` contains all the header files.

  + `info.hpp` contains some macro definitions.
  
  + `LinSysSolver.hpp` defines an abstract base class that stores sparse matrices in CSR format for the assembly of stiffness matrices and the solution of linear systems.
  
    `EigenLibSolver.hpp` and `CHOLMODSolver.hpp` are the corresponding inherited classes, representing the use of the Eigen library or CHOLMOD library to solve.
    
  + `readTetMesh.hpp` defines the function to read the tetrahedral mesh.msh file.
  
    `readHexMesh.hpp` defines the function to read the hexahedral mesh .in files, where the .in files are parsed from the .vtk files.
  
    `readCondition.hpp` defines function to read boundary conditions and other information, with requirements for the input file format.
  
  + `shapeInterface.hpp` defines the abstract base class for the grid.
  
    `tet4nodeEle.hpp` is an inherited 4-node tetrahedron class.
  
    `hex8nodeEle.hpp` is an inherited 8-node hexahedron class.
  
  + `tet4nodeSolver.hpp` encapsulates the finite element analysis process of a tetrahedral mesh model.
  
    `hex8nodeSolver.hpp` encapsulates the finite element analysis process of a hexahedral mesh model.
  
+ `tests` contains all test examples.

  + `tet_cube` tests a simple mesh model with 8 vertices and 8 tetrahedra, and the correctness was verified by a matlab program.
  + `tet_bunny` tests the tetrahedral mesh model *bunny.msh*.
  + `tet_armadillo` tests the tetrahedral mesh model *Armadillo219K.msh*, which has 219k nodes and takes about 30 minutes to solve using the Eigen library and only 5 minutes using the CHOLMOD library.
  + `hex_smallCube` tests a small cube model consisting of 27 small hexahedra and verified the correctness using matlab.
  + `hex_cube` tests a larger cube model consisting of 512 hexahedra with a total of 729 nodes.

## Compile and Run

The code author has the following environment and software locally:

+ Ubuntu 20.04
+ g++(GCC) 11.2.0 (C++20)
+ Python 3.8.10
+ CMake 3.16.3
+ Eigen 3.3.7
+ SuiteSparse 5.7.1
+ [MshIO](https://github.com/qnzhou/MshIO)
+ [meshio 5.3.4](https://github.com/nschloe/meshio)

Compiling:

+ In the *hex* test sample, type the command

  ```bash
  mkdir build && cd build
  cmake ..
  make
  ```

  to compile. Type `./main` to run it.

+ In the *tet* sample, type `make` to compile. Type `./main` to run it.