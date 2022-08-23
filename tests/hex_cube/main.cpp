#include <iostream>
#include <set>

#include <python3.8/Python.h>
#include <mshio/mshio.h>

#include "../../include/hex8nodeSolver.hpp"

using namespace Eigen;
using namespace lhfea;

void conditionProcess(std::string filePath, std::string condPath){
    realn E = 1e10, nu = 0.25;
    MatrixXr TV; MatrixXi TT;
    readHexMesh(filePath, TV, TT);
    std::ofstream out(condPath);
    out << E << " " << nu << std::endl;
    std::vector<int> bottom_numb;
    std::vector<int> top_numb;
    for(int i=0; i<TV.rows(); ++i){
        if(TV(i,1) == 0)
            bottom_numb.push_back(i);
        if(TV(i,1) == 1)
            top_numb.push_back(i);
    }
    for(auto i : bottom_numb){
        out << i << " " << 1 << " " << 1 << " " << 1 << " " 
            << 0 << " " << 0 << " " << 0 << " "
            << 0 << " " << 0 << " " << 0 << std::endl;
    }
    for(auto i : top_numb){
        out << i << " " << 0 << " " << 0 << " " << 0 << " " 
            << 0 << " " << 0 << " " << 0 << " "
            << 5e6 << " " << 0 << " " << 0 << std::endl;
    }
    out.close();
}

int main()
{
    /* PrePythonProcess */
    Py_Initialize();
    PyRun_SimpleString("import os,sys");
    PyRun_SimpleString("sys.path.append('./')");
    PyRun_SimpleString("sys.path.append('../')");
    PyRun_SimpleString("import processVTK");
    PyRun_SimpleString("processVTK.transHexFile()");
    
    std::string filePath = "./input/cube.in";
    std::string condPath = "./input/condition.in";
    conditionProcess(filePath, condPath);

    VectorXr result_uvw;
    hex8nodeSolver(filePath, condPath, result_uvw, SolverType::cholmod);

    std::ofstream outf("./output/result_uvw.txt");
    for(int i=0; i<result_uvw.size(); i+=3){
        outf << result_uvw[i] << " "
             << result_uvw[i+1] << " "
             << result_uvw[i+2] << std::endl;
    }
    outf.close();

    /* PostPythonProcess */
    PyRun_SimpleString("processVTK.outputHexVTK()");
    Py_Finalize();

    return 0;
}