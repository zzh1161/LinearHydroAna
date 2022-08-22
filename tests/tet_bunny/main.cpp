#include "preprocess.hpp"
#include "../../include/tet4nodeSolver.hpp"

using namespace lhfea;
using namespace Eigen;

int main()
{
    std::string filePath = "./input/bunny.msh";
    std::string condPath = "./input/condition.in";

    pre_condition_process_bunny(filePath, condPath);

    // MatrixXr TV;
    // MatrixXi TT;
    // readTetMesh(filePath, TV, TT);
    // std::cout << "TV:\n" << TV << std::endl << std::endl
    //           << "TT:\n" << TT << std::endl << std::endl;

    VectorXr result_uvw;
    tet4nodeSolver(filePath, condPath, result_uvw, SolverType::cholmod);

    std::ofstream out("./output/result_uvw.txt");
    for(int i=0; i<result_uvw.size(); i+=3){
        out << result_uvw(i) << " "
            << result_uvw(i+1) << " "
            << result_uvw(i+2) << std::endl;
    }
    out.close();

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
    mshio::save_msh("./output/result_bunny.msh", spec);

    return 0;
}