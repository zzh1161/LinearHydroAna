#include <iostream>
#include <set>

#include <mshio/mshio.h>
#include <python3.8/Python.h>

#include "../../include/hex8nodeSolver.hpp"

using namespace Eigen;
using namespace lhfea;

int main()
{
    Py_Initialize();
    PyRun_SimpleString("import os,sys");
    PyRun_SimpleString("os.chdir('../')");

    std::string filePath = "./input/cube.in";
    std::string condPath = "./input/condition.in";

    VectorXr result_uvw;
    hex8nodeSolver(filePath, condPath, result_uvw, SolverType::cholmod);

    std::ofstream outf("./output/result_uvw.txt");
    for(int i=0; i<result_uvw.size(); i+=3){
        outf << result_uvw[i] << " "
             << result_uvw[i+1] << " "
             << result_uvw[i+2] << std::endl;
    }
    outf.close();

    Py_Finalize();

    return 0;
}