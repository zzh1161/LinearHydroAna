#ifndef READ_HEX_MESH_HPP
#define READ_HEX_MESH_HPP

#include <fstream>
#include <vector>

#include <spdlog/spdlog.h>

#include "info.hpp"

namespace lhfea{

bool readHexMesh(const std::string filePath, MatrixXr &TV, Eigen::MatrixXi &TT)
{
    std::ifstream in(filePath);
    int rows; in >> rows;
    try{
        TV.resize(rows, 3);
        for(int i=0; i<rows; i++){
            realn a, b, c;
            in >> a >> b >> c;
            TV.row(i) << a, b, c;
        }
        in >> rows;
        TT.resize(rows, 8);
        for(int i=0; i<rows; i++){
            int a,b,c,d,e,f,g,h;
            in >> a >> b >> c >> d >> e >> f >> g >> h;
            TT.row(i) << a,b,c,d,e,f,g,h;
        }
    }
    catch(...){
        spdlog::error("The input file does not meet the format requirements.");
        exit(-1);
    }
    in.close();
    return true;
}

} // namespace lhfea



#endif