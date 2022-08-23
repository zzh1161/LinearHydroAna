#pragma once

#include "../../include/tet4nodeSolver.hpp"

/**************************************************************
 > @description: Pre-processing for Armadillo219K.msh
 > @param {string} filePath : .msh file path
 > @param {string} condPath : .in file path
**************************************************************/
bool pre_condition_process_armadillo(std::string filePath, std::string condPath)
{
    double E = 1e10, nu = 0.25;
    lhfea::MatrixXr TV;
    Eigen::MatrixXi TT;
    if(!lhfea::readTetMesh(filePath, TV, TT))
        return false;

    std::ofstream out(condPath);
    out << E << " " << nu << std::endl;

    std::vector<int> bottom_numb;
    std::vector<int> left_numb;
    std::vector<int> right_numb;
    for(int i=0; i<TV.rows(); ++i){
        if(TV(i,1) < 21.5)
            bottom_numb.push_back(i);
        if(TV(i,0) > 52 && TV(i,1) > 21.5)
            right_numb.push_back(i);
        if(TV(i,0) < -52 && TV(i,1) > 21.5)
            left_numb.push_back(i);
    }
    for(auto i : bottom_numb){
        out << i << " " << 1 << " " << 1 << " " << 1 << " " 
            << 0 << " " << 0 << " " << 0 << " "
            << 0 << " " << 0 << " " << 0 << std::endl;
    }
    for(auto i : left_numb){
        out << i << " " << 0 << " " << 0 << " " << 0 << " " 
            << 0 << " " << 0 << " " << 0 << " "
            << 2e6 << " " << 0 << " " << 0 << std::endl;
    }
    for(auto i : right_numb){
        out << i << " " << 0 << " " << 0 << " " << 0 << " " 
            << 0 << " " << 0 << " " << 0 << " "
            << -2e6 << " " << 0 << " " << 0 << std::endl;
    }

    out.close();

    return true;
}