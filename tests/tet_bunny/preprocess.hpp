#pragma once

#include "../../include/tet4nodeSolver.hpp"

/**************************************************************
 > @description: Pre-processing for bunny.msh
 > @param {string} filePath : .msh file path
 > @param {string} condPath : .in file path
**************************************************************/
bool pre_condition_process_bunny(std::string filePath, std::string condPath)
{
    double E = 1e10, nu = 0.25;
    lhfea::MatrixXr TV;
    Eigen::MatrixXi TT;
    if(!lhfea::readTetMesh(filePath, TV, TT))
        return false;

    std::ofstream out(condPath);
    out << E << " " << nu << std::endl;

    std::vector<int> bottom_numb;
    std::vector<int> top_numb;
    for(int i=0; i<TV.rows(); ++i){
        if(TV(i,1) < -0.24)
            bottom_numb.push_back(i);
        if(TV(i,1) > 2.50)
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
            << 2e5 << " " << 0 << " " << 0 << std::endl;
    }

    out.close();

    return true;
}