#ifndef HEX_8_NODE_SOLVER_HPP
#define HEX_8_NODE_SOLVER_HPP

#include "info.hpp"
#include "readHexMesh.hpp"
#include "readCondition.hpp"
#include "hex8nodeEle.hpp"
#include "EigenLibSolver.hpp"

namespace lhfea{

#ifdef USE_FUNCTOR_SOLVER
struct _Hex8NodeSolver{
    void solve(std::string filePath, std::string condPath, VectorXr &result_uvw) const;
    void operator()(std::string filePath, std::string condPath, VectorXr &result_uvw) const{
        solve(filePath, condPath, result_uvw);
    }
};
inline constexpr _Hex8NodeSolver hex8nodeSolver;
#endif

#ifdef USE_FUNCTOR_SOLVER
void _Hex8NodeSolver::solve(std::string filePath, std::string condPath, VectorXr &result_uvw) const
#else
void hex8nodeSolver(std::string filePath, std::string condPath, VectorXr &result_uvw)
#endif
{
    /* Get info from .in format files processed by python */
    MatrixXr TV; Eigen::MatrixXi TT;
    assert(readHexMesh(filePath, TV, TT) && "Read Hex Mesh");

    /* Get boundary conditions and other info */
    realn E, nu;
    std::unordered_map<int,std::vector<realn>> conditions;
    assert(readCondition(condPath, E, nu, conditions) && "Read bounary condition file");

    /* Preprocess for set_pattern: set vNeighbor */
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
    LHsolver.factorize();
    LHsolver.solve(rhs, result_uvw);
    std::cout << "Solve the linear system success!" << std::endl;
}

} // namespace lhfea


#endif