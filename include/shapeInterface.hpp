#ifndef SHAPE_INTERFACE_HPP
#define SHAPE_INTERFACE_HPP

#include "info.hpp"
#include "readTetMesh.hpp"

namespace lhfea{

// #define NOT_KEEP_BMat
#if __GNUC__ > 7
    #define USE_FUNCTOR_SOLVER
#endif

struct shapeInterface
{
    std::vector<int> node_num;  // number of nodes
    // MatrixXr      node_coo;  // coordinates of nodes, but not necessary
    MatrixXr         ele_stiff; // unit stiffness matrix
    // MatrixXr      ele_DMat;  // elasticity coefficient matrix, not necessary to keep it

    shapeInterface() = default;
    shapeInterface(std::vector<int> num): node_num(num){}
    virtual ~shapeInterface() = default;

    /*******************************************************************************
     > @description: Calculating Unit Stiffness Matrix
     > @param {MatrixXr} &TV : Matrix describing the edges of polyhedrons
     > @param {realn} E      : Modulus of elasticity
     > @param {realn} nu     : Poisson's ratio
    ********************************************************************************/
    virtual void calStiffnessMat(MatrixXr &TV, realn E, realn nu) = 0;
};

} // namespace lhfea


#endif