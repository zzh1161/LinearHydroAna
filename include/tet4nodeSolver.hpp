#ifndef TET_4_NODE_SOLVER_HPP
#define TET_4_NODE_SOLVER_HPP

#include <memory>

#include "info.hpp"
#include "readTetMesh.hpp"
#include "readCondition.hpp"
#include "EigenLibSolver.hpp"
#include "CHOLMODSolver.hpp"
#include "tet4nodeEle.hpp"

namespace lhfea{

#ifdef USE_FUNCTOR_SOLVER
struct _Tet4NodeSolver
{
    void solve(std::string meshPath, std::string condiPath, VectorXr &result_uvw, SolverType solver) const;
    void operator()(std::string meshPath, std::string condiPath, VectorXr &result_uvw, SolverType solver) const{
        solve(meshPath, condiPath, result_uvw, solver);
    }
};
inline constexpr _Tet4NodeSolver tet4nodeSolver;
#endif

#ifdef USE_FUNCTOR_SOLVER
void _Tet4NodeSolver::solve(std::string meshPath, std::string condiPath, VectorXr &result_uvw, SolverType solver) const
#else
void tet4nodeSolver(std::string meshPath, std::string condiPath, VectorXr &result_uvw, SolverType solver)
#endif
{
    /* Get info from .msh format files */
    MatrixXr TV; Eigen::MatrixXi TT;
    assert(readTetMesh(meshPath, TV, TT) && "Read Tet Mesh");

    /* Get boundary conditions and other info */
    realn E, nu;
    std::unordered_map<int,std::vector<realn>> conditions;
    assert(readCondition(condiPath, E, nu, conditions) && 
            "Read bounary condition file");

    std::vector<std::set<int>> vNeighbor(TV.rows());
    for(int i=0; i<TT.rows(); i++){
        if(TT(i,0)==TT(i,1) || TT(i,0)==TT(i,2) || TT(i,0)==TT(i,3) ||
            TT(i,1)==TT(i,2) || TT(i,1)==TT(i,3) || TT(i,2)==TT(i,3))
            continue;
       vNeighbor[TT(i,0)].insert(TT(i,1)); vNeighbor[TT(i,0)].insert(TT(i,2)); vNeighbor[TT(i,0)].insert(TT(i,3));
       vNeighbor[TT(i,1)].insert(TT(i,0)); vNeighbor[TT(i,1)].insert(TT(i,2)); vNeighbor[TT(i,1)].insert(TT(i,3));
       vNeighbor[TT(i,2)].insert(TT(i,0)); vNeighbor[TT(i,2)].insert(TT(i,1)); vNeighbor[TT(i,2)].insert(TT(i,3));
       vNeighbor[TT(i,3)].insert(TT(i,0)); vNeighbor[TT(i,3)].insert(TT(i,1)); vNeighbor[TT(i,3)].insert(TT(i,2));   
    }
    std::cout << "Set vNeigehbor successfully!" << std::endl;

    /* Initialize tethrahedrons */
    std::vector<tet4node> tethrahedrons;
    for(int i=0; i<TT.rows(); i++){
        if(TT(i,0)==TT(i,1) || TT(i,0)==TT(i,2) || TT(i,0)==TT(i,3) ||
            TT(i,1)==TT(i,2) || TT(i,1)==TT(i,3) || TT(i,2)==TT(i,3))
            continue;
        auto tmp = TT.row(i);
        tet4node Tet(tmp(0), tmp(1), tmp(2), tmp(3));
        Tet.calStiffnessMat(TV, E, nu);
        tethrahedrons.push_back(Tet);
    }
    std::cout << "Initialize tethrahedrons successfully!" << std::endl;

    /* Assemble stiffness matrices */
    std::unique_ptr<SIM::LinSysSolver<Eigen::VectorXi,VectorXr>> LHsolver;
    if(solver == SolverType::Eigen)
        LHsolver = std::make_unique<SIM::EigenLibSolver<Eigen::VectorXi,VectorXr>>();
    else if(solver == SolverType::cholmod)
        LHsolver = std::make_unique<SIM::CHOLMODSolver<Eigen::VectorXi,VectorXr>>();
    LHsolver->set_pattern(vNeighbor);
    LHsolver->analyze_pattern();
    LHsolver->setZero();
    for(auto &tet : tethrahedrons){
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                for(int m=0; m<3; m++)
                    for(int n=0; n<3; n++)
                        LHsolver->addCoeff(3*tet.node_num[i]+m, 3*tet.node_num[j]+n, tet.ele_stiff(3*i+m, 3*j+n));
            }
        }
    }
    std::cout << "Assemble stiffness matrices successfully!" << std::endl;

    /* Set boundary conditions and RHS */
    VectorXr rhs(TV.rows()*3,1);
    for(auto iter=conditions.begin(); iter!=conditions.end(); ++iter){
        auto list = iter->second;
        for(int i=0; i<3; ++i){
            if(list[i] == 1){
                LHsolver->setUnit_row(3*iter->first+i);
                LHsolver->setUnit_col(3*iter->first+i, vNeighbor[iter->first]);
                rhs[iter->first*3+i] = list[i+3];
            }
            else{
                rhs[iter->first*3+i] += (list[i+3]+list[i+6]);
            }
        }
    }
    std::cout << "Set boundary conditions and RHS successfully!" << std::endl;

    /* Solve the linear system */
    LHsolver->factorize();
    LHsolver->solve(rhs, result_uvw);
    std::cout << "Solve the linear system successfully!" << std::endl;
}


} // namespace lhfea


#endif