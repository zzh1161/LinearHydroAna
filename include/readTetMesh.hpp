#ifndef READ_TET_MESH_HPP
#define READ_TET_MESH_HPP

#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <vector>

#include <spdlog/spdlog.h>
#include <mshio/mshio.h>
#include <eigen3/Eigen/Eigen>
// #include <boost/filesystem.hpp>

#include "info.hpp"

namespace lhfea{

bool readTetMesh(const std::string& filePath, MatrixXr& TV, Eigen::MatrixXi& TT){
    // using namespace boost::filesystem;
    // if (!exists(path(filePath))) {
    //     return false;
    // }

    mshio::MshSpec spec;
    try {
        spec = mshio::load_msh(filePath);
    }
    catch (...) {
        spdlog::error("MshIO only supports MSH 2.2 and 4.1");
        exit(-1);
    }

    const auto& nodes = spec.nodes;
    const auto& els = spec.elements;
    const int vAmt = nodes.num_nodes;
    int elemAmt = 0;
    for (const auto& e : els.entity_blocks) {
        assert(e.entity_dim == 3);
        // assert(e.element_type == 4); // linear tet
        elemAmt += e.num_elements_in_block;
    }

    TV.resize(vAmt, 3);
    int index = 0;
    for (const auto& n : nodes.entity_blocks) {
        for (int i = 0; i < n.num_nodes_in_block * 3; i += 3) {
            TV.row(index) << n.data[i], n.data[i + 1], n.data[i + 2];
            ++index;
        }
    }

    TT.resize(elemAmt, 4);
    int elm_index = 0;
    for (const auto& e : els.entity_blocks) {
        for (int i = 0; i < e.data.size(); i += 5) {
            index = 0;
            for (int j = i + 1; j <= i + 4; ++j) {
                TT(elm_index, index++) = e.data[j] - 1;
            }
            ++elm_index;
        }
    }

    // finding the surface because $Surface is not supported by MshIO
    // spdlog::info("Finding the surface triangle mesh for {:s}", filePath);
    // find_surfTri_from_tet(TT, SF);

    // spdlog::info("tet mesh loaded with {:d} particles, {:d} tets, and {:d} surface triangles.",
    //              TV.rows(), TT.rows(), SF.rows());

    return true;
}

bool readCondition(std::string path, realn &E, realn &nu,
            std::unordered_map<int,std::vector<realn>> &bound)
{
    std::ifstream in(path);
    in >> E >> nu;
    try{
        while(!in.fail()){
            int idx; in >> idx;
            realn tmp;
            std::vector<realn> list;
            for(int i=0; i<9; i++){
                in >> tmp;
                list.push_back(tmp);
            }
            bound.emplace(idx, list);
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