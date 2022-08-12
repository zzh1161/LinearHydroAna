#ifndef TET_4_NODE_ELE_HPP
#define TET_4_NODE_ELE_HPP

#include <cmath>

#include "shapeInterface.hpp"

namespace lhfea{

struct tet4node : public shapeInterface
{
    // realn ele_vol; // volume of the element, but not necessary
#ifndef NOT_KEEP_BMat
    MatrixXr         ele_BMat;  // strain-displacement matrix
#endif

    tet4node(){}
    tet4node(std::vector<int> num) : shapeInterface(num){}
    tet4node(int i, int j, int m, int n)
        {node_num.push_back(i); node_num.push_back(j); node_num.push_back(m); node_num.push_back(n);}
    ~tet4node(){}

    void calStiffnessMat(MatrixXr &TV, realn E, realn nu);
};

void tet4node::calStiffnessMat(MatrixXr &TV, realn E, realn nu)
{
#ifndef NOT_KEEP_BMat
    ele_BMat.resize(6,12);
#else
    MatrixXr ele_BMat(6,12);
#endif
    MatrixXr ele_DMat(6,6);
    MatrixXr coord(4,4);
    #define x(i) TV(node_num[i-1],0)
    #define y(i) TV(node_num[i-1],1)
    #define z(i) TV(node_num[i-1],2)
    coord << 1, x(1), y(1), z(1),
             1, x(2), y(2), z(2),
             1, x(3), y(3), z(3),
             1, x(4), y(4), z(4);
    realn volx6 = coord.determinant();
    realn beta[5]; realn gamma[5]; realn delta[5];
    for(int m=1; m<5; m++){
        std::vector<int> idx;
        for(int j=1; j<5; j++){ if(j != m) idx.push_back(j); }
        MatrixXr matbeta(3,3);
        MatrixXr matgam(3,3);
        MatrixXr matdel(3,3);
        matbeta << 1,y(idx[0]),z(idx[0]), 1,y(idx[1]),z(idx[1]), 1,y(idx[2]),z(idx[2]);
        matgam  << 1,x(idx[0]),z(idx[0]), 1,x(idx[1]),z(idx[1]), 1,x(idx[2]),z(idx[2]);
        matdel  << 1,x(idx[0]),y(idx[0]), 1,x(idx[1]),y(idx[1]), 1,x(idx[2]),y(idx[2]);
        beta[m]  = std::pow(-1,m)*matbeta.determinant();
        gamma[m] = std::pow(-1,m-1)*matgam.determinant();
        delta[m] = std::pow(-1,m)*matdel.determinant();
    }
    ele_BMat << beta[1],0,0,         beta[2],0,0,         beta[3],0,0,         beta[4],0,0,
                0,gamma[1],0,        0,gamma[2],0,        0,gamma[3],0,        0,gamma[4],0,
                0,0,delta[1],        0,0,delta[2],        0,0,delta[3],        0,0,delta[4],
                gamma[1],beta[1],0,  gamma[2],beta[2],0,  gamma[3],beta[3],0,  gamma[4],beta[4],0,
                0,delta[1],gamma[1], 0,delta[2],gamma[2], 0,delta[3],gamma[3], 0,delta[4],gamma[4],
                delta[1],0,beta[1],  delta[2],0,beta[2],  delta[3],0,beta[3],  delta[4],0,beta[4];
    ele_BMat = ele_BMat*(1.0/volx6);
    ele_DMat << 1-nu, nu, nu, 0, 0, 0,
                nu, 1-nu, nu, 0, 0, 0,
                nu, nu, 1-nu, 0, 0, 0,
                0, 0, 0, 0.5-nu, 0, 0,
                0, 0, 0, 0, 0.5-nu, 0,
                0, 0, 0, 0, 0, 0.5-nu;
    ele_DMat = ele_DMat*(E/((1+nu)*(1-2*nu)));
    ele_stiff = ele_BMat.transpose()*ele_DMat*ele_BMat*std::abs(volx6/6);
}

} // namespace lhfea

#endif