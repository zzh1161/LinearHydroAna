#include <iostream>
#include <set>

#include <mshio/mshio.h>
#include "../../include/readHexMesh.hpp"
#include "../../include/readCondition.hpp"
#include "../../include/hex8nodeEle.hpp"
#include "../../include/EigenLibSolver.hpp"

using namespace Eigen;
using namespace lhfea;

int main()
{
    std::string filePath = "./input/cube.in";
    std::string condPath = "./input/condition.in";

    MatrixXr TV; MatrixXi TT;
    assert(readHexMesh(filePath, TV, TT) && "Read Hex Mesh");
    realn E, nu;
    std::unordered_map<int,std::vector<realn>> conditions;
    assert(readCondition(condPath, E, nu, conditions) && "Read bounary condition file");

    /* Pre for set_pattern: set vNeighbor */
    std::vector<std::set<int>> vNeighbor(TV.rows());
    for(int i=0; i<TT.rows(); i++){
        for(int j=0; j<8; j++){
            for(int k=0; k<8; k++){
                if(k == j) continue;
                vNeighbor[TT(i,j)].insert(TT(i,k));
            }
        }
    }
    std::cout << "Set vNeigehbor success!" << std::endl;

    /* Initialize hexes */
    std::vector<hex8node> hexes;
    for(int i=0; i<TT.rows(); i++){
        auto tmp = TT.row(i);
        hex8node Hex(tmp(0),tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6),tmp(7));
        Hex.calStiffnessMat(TV, E, nu);
        hexes.push_back(Hex);
    }
    std::cout << "Initialize hexes success!" << std::endl;

    /* Assemble stiffness matrices */
    SIM::EigenLibSolver<Eigen::VectorXi, VectorXr> LHsolver;
    LHsolver.set_pattern(vNeighbor);
    LHsolver.analyze_pattern();
    for(auto &hex : hexes){
        for(int i=0; i<8; i++){
            for(int j=0; j<8; j++){
                for(int m=0; m<3; m++)
                    for(int n=0; n<3; n++)
                        LHsolver.addCoeff(3*hex.node_num[i]+m, 3*hex.node_num[j]+n, hex.ele_stiff(3*i+m, 3*j+n));
            }
        }
    }
    std::cout << "Assemble stiffness matrices success!" << std::endl;

    /* Set boundary conditions and RHS */
    VectorXr rhs(TV.rows()*3,1);
    for(auto iter=conditions.begin(); iter!=conditions.end(); ++iter){
        auto list = iter->second;
        for(int i=0; i<3; ++i){
            if(list[i] == 1){
                LHsolver.setUnit_row(3*iter->first+i);
                LHsolver.setUnit_col(3*iter->first+i, vNeighbor[iter->first]);
                rhs[iter->first*3+i] = list[i+3];
            }
            else{
                rhs[iter->first*3+i] += (list[i+3]+list[i+6]);
            }
        }
    }
    std::cout << "Set boundary conditions and RHS success!" << std::endl;

    /* Solve the linear system */
    VectorXr result_uvw;
    LHsolver.factorize();
    LHsolver.solve(rhs, result_uvw);
    std::cout << "Solve the linear system success!" << std::endl;

    std::ofstream outf("./output/result_uvw.txt");
    for(int i=0; i<result_uvw.size(); i+=3){
        outf << result_uvw[i] << " "
             << result_uvw[i+1] << " "
             << result_uvw[i+2] << std::endl;
    }
    // outf << std::endl;
    // for(auto &hex : hexes){
    //     outf << hex.ele_stiff << std::endl << std::endl;
    // }
    outf.close();

    return 0;
}