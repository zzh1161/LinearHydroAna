#ifndef HEX_8_NODE_ELE_HPP
#define HEX_8_NODE_ELE_HPP

#include <functional>
#include "shapeInterface.hpp"

namespace lhfea{

struct hex8node : public shapeInterface
{
    using funcMat = std::function<MatrixXr(realn,realn,realn)>;

    // realn ele_vol; // not necessary
#ifndef NOT_KEEP_BMat
    funcMat ele_BMat;
    funcMat Jacobi;
#endif

    hex8node() : shapeInterface(){}
    hex8node(std::vector<int> num) : shapeInterface(num){}
    hex8node(int a, int b, int c, int d, int e, int f, int g, int h)
            { node_num = std::vector<int>{a,b,c,d,e,f,g,h}; }
    ~hex8node(){}

    void calStiffnessMat(MatrixXr &TV, realn E, realn nu);

private:
#ifndef NOT_KEEP_BMat
    void calBMatrix(MatrixXr &TV, realn E, realn nu);
#else
    void calBMatrix(MatrixXr &TV, realn E, realn nu, funcMat &ele_BMat, funcMat &Jacobi) const;
#endif
    void calDMatrix(realn E, realn nu, MatrixXr &ele_DMat) const;
};

void hex8node::calStiffnessMat(MatrixXr &TV, realn E, realn nu)
{
// #ifndef NOT_KEEP_BMat
//     ele_BMat.resize(6,12);
// #else
//     MatrixXr ele_BMat(6,12);
// #endif
    MatrixXr ele_DMat(6,6);
    #define x(i) TV(node_num[i-1],0)
    #define y(i) TV(node_num[i-1],1)
    #define z(i) TV(node_num[i-1],2)
}

#ifndef NOT_KEEP_BMat
void hex8node::calBMatrix(MatrixXr &TV, realn E, realn nu)
#else
void hex8node::calBMatrix(MatrixXr &TV, realn E, realn nu, 
        hex8node::funcMat &ele_BMat, hex8node::funcMat &Jacobi) const
#endif
{
    #define x(i) TV(node_num[i-1],0)
    #define y(i) TV(node_num[i-1],1)
    #define z(i) TV(node_num[i-1],2)
    std::vector<int> xi    = {0,1,1,-1,-1,1,1,-1,-1};
    std::vector<int> eta   = {0,-1,1,1,-1,-1,1,1,-1};
    std::vector<int> gamma = {0,-1,-1,-1,-1,1,1,1,1};
    // std::vector<std::function<realn(realn,realn,realn)>> dNds(9), dNdt(9), dNdn(9);
    // for(int i=1; i<9; i++){
    //     dNds[i] = [&xi, &eta, &gamma, &i](realn s, realn t, realn n){
    //         return xi[i]*(1+eta[i]*t)*(1+gamma[i]*n)/8;
    //     };
    //     dNdt[i] = [&xi, &eta, &gamma, &i](realn s, realn t, realn n){
    //         return eta[i]*(1+xi[i]*s)*(1+gamma[i]*n)/8;
    //     };
    //     dNdn[i] = [&xi, &eta, &gamma, &i](realn s, realn t, realn n){
    //         return gamma[i]*(1+xi[i]*s)*(1+eta[i]*t)/8;
    //     };
    // }
    #define dNds(NUM,FIR,SEC,THR) xi[NUM]*(1+eta[NUM]*SEC)*(1+gamma[NUM]*THR)/8
    #define dNdt(NUM,FIR,SEC,THR) eta[NUM]*(1+xi[NUM]*FIR)*(1+gamma[NUM]*THR)/8
    #define dNdn(NUM,FIR,SEC,THR) gamma[NUM]*(1+xi[NUM]*FIR)*(1+eta[NUM]*SEC)/8
    std::function<realn(realn,realn,realn)> dxds, dxdt, dxdn, dyds, dydt, dydn, dzds, dzdt, dzdn;
    dxds = [&](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += x(i)*dNds(i,s,t,n); return res;};
    dyds = [&](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += y(i)*dNds(i,s,t,n); return res;};
    dzds = [&](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += z(i)*dNds(i,s,t,n); return res;};
    dxdt = [&](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += x(i)*dNdt(i,s,t,n); return res;};
    dydt = [&](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += y(i)*dNdt(i,s,t,n); return res;};
    dzdt = [&](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += z(i)*dNdt(i,s,t,n); return res;};
    dxdn = [&](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += x(i)*dNdn(i,s,t,n); return res;};
    dydn = [&](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += y(i)*dNdn(i,s,t,n); return res;};
    dzdn = [&](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += z(i)*dNdn(i,s,t,n); return res;};
    /* Calculate Jacobi matrix */
    Jacobi = [&](realn s, realn t, realn n){
        MatrixXr result(3,3);
        result << dxds(s,t,n), dyds(s,t,n), dzds(s,t,n),
                  dxdt(s,t,n), dydt(s,t,n), dzdt(s,t,n),
                  dxdn(s,t,n), dydn(s,t,n), dzdn(s,t,n);
        return result;
    };
    auto Jdet = [&](realn s, realn t, realn n){
        return Jacobi(s,t,n).determinant();  
    };
    /* Calculate B matrix */

}

void hex8node::calDMatrix(realn E, realn nu, MatrixXr &ele_DMat) const
{
    ele_DMat.resize(6,6);
    ele_DMat << 1-nu, nu, nu, 0, 0, 0,
                nu, 1-nu, nu, 0, 0, 0,
                nu, nu, 1-nu, 0, 0, 0,
                0, 0, 0, 0.5-nu, 0, 0,
                0, 0, 0, 0, 0.5-nu, 0,
                0, 0, 0, 0, 0, 0.5-nu;
    ele_DMat = ele_DMat*(E/((1+nu)*(1-2*nu)));
}

// #ifndef NOT_KEEP_BMat
// void hex8node::calJacobiMat(MatrixXr &TV)
// #else
// void hex8node::calJacobiMat(MatrixXr &TV, hex8node::funcMat &Jacobi) const
// #endif
// {
//     #define x(i) TV(node_num[i-1],0)
//     #define y(i) TV(node_num[i-1],1)
//     #define z(i) TV(node_num[i-1],2)
//     std::vector<int> xi    = {0,1,1,-1,-1,1,1,-1,-1};
//     std::vector<int> eta   = {0,-1,1,1,-1,-1,1,1,-1};
//     std::vector<int> gamma = {0,-1,-1,-1,-1,1,1,1,1};
//     std::vector<std::function<realn(realn,realn)>> dNds(9), dNdt(9), dNdn(9);
//     for(int i=1; i<9; i++){
//         dNds[i] = [&xi, &eta, &gamma, i](realn t, realn n){
//             return xi[i]*(1+eta[i]*t)*(1+gamma[i]*n)/8;
//         };
//         dNdt[i] = [&xi, &eta, &gamma, i](realn s, realn n){
//             return eta[i]*(1+xi[i]*s)*(1+gamma[i]*n)/8;
//         };
//         dNdn[i] = [&xi, &eta, &gamma, i](realn s, realn t){
//             return gamma[i]*(1+xi[i]*s)*(1+eta[i]*t)/8;
//         };
//     }
//     Jacobi = [&](realn s, realn t, realn n){
//         realn dxds=0, dxdt=0, dxdn=0;
//         realn dyds=0, dydt=0, dydn=0;
//         realn dzds=0, dzdt=0, dzdn=0;
//         for(int i=1; i<9; i++){
//             dxds += x(i)*dNds[i](t,n); dxdt += x(i)*dNdt[i](s,n); dxdn += x(i)*dNdn[i](s,t);
//             dyds += y(i)*dNds[i](t,n); dydt += y(i)*dNdt[i](s,n); dydn += y(i)*dNdn[i](s,t);
//             dzds += z(i)*dNds[i](t,n); dzdt += z(i)*dNdt[i](s,n); dzdn += z(i)*dNdt[i](s,t);
//         }
//         MatrixXr result(3,3);
//         result << dxds, dyds, dzds, dxdt, dydt, dzdt, dxdn, dydn, dzdn;
//         return result;
//     };
// }

} // namespace lhfea


#endif