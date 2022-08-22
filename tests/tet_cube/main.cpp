#include <iostream>

#include "../../include/tet4nodeSolver.hpp"

using namespace lhfea;
using namespace Eigen;

int main()
{
    std::string filePath = "./input/cube.msh";
    std::string condPath = "./input/condition.in";
    VectorXr result_uvw;

    MatrixXr TV; MatrixXi TT;
    readTetMesh(filePath, TV, TT);
    std::cout << "TV: \n" << TV << std::endl << std::endl
              << "TT: \n" << TT << std::endl << std::endl;

    tet4nodeSolver(filePath, condPath, result_uvw, SolverType::cholmod);
    
    for(int i=0; i<14; i++){
        for(int j=0; j<3; j++){
            std::cout << result_uvw(3*i+j) << "  ";
        }
        std::cout << std::endl;
    }

    mshio::MshSpec spec;
    spec = mshio::load_msh(filePath);
    int j = 0;
    for(auto &n : spec.nodes.entity_blocks){
        for (int i = 0; i < n.num_nodes_in_block*3; i += 3){
            n.data[i]   += result_uvw[j];
            n.data[i+1] += result_uvw[j+1];
            n.data[i+2] += result_uvw[j+2];
            j += 3;
        }
    }
    mshio::save_msh("./output/result_cube.msh", spec);

    return 0;
}
