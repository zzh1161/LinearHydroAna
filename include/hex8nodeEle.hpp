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
    void calBMatrix(MatrixXr &TV);
#else
    void calBMatrix(MatrixXr &TV, funcMat &ele_BMat, funcMat &Jacobi) const;
#endif
    void calDMatrix(realn E, realn nu, MatrixXr &ele_DMat) const;
};

void hex8node::calStiffnessMat(MatrixXr &TV, realn E, realn nu)
{
    MatrixXr ele_DMat(6,6);
#ifdef NOT_KEEP_BMat
    hex8node::funcMat ele_BMat;
    hex8node::funcMat Jacobi;
    calBMatrix(TV, ele_BMat, Jacobi);
#else
    calBMatrix(TV);
#endif
    calDMatrix(E, nu, ele_DMat);
    std::vector<realn> intpont = {-0.774596669241483, 0, 0.774596669241483};
    std::vector<realn> intwght = {0.555555555555556, 0.888888888888889, 0.555555555555556};
    ele_stiff.resize(24,24); ele_stiff.setZero();
    for(int m=0; m<3; m++){
        for(int j=0; j<3; j++){
            for(int i=0; i<3; i++){
                auto B = ele_BMat(intpont[i],intpont[j],intpont[m]);
                ele_stiff += B.transpose()*ele_DMat*B*
                    Jacobi(intpont[i],intpont[j],intpont[m]).determinant()*(intwght[m]*intwght[j]*intwght[i]);
            }
        }
    }
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

#ifndef NOT_KEEP_BMat
void hex8node::calBMatrix(MatrixXr &TV)
#else
void hex8node::calBMatrix(MatrixXr &TV, hex8node::funcMat &ele_BMat, hex8node::funcMat &Jacobi) const
#endif
{
    #define x(i) TV(node_num[i-1],0)
    #define y(i) TV(node_num[i-1],1)
    #define z(i) TV(node_num[i-1],2)
    std::vector<int> xi    = {0,1,1,-1,-1,1,1,-1,-1};
    std::vector<int> eta   = {0,-1,1,1,-1,-1,1,1,-1};
    std::vector<int> gamma = {0,-1,-1,-1,-1,1,1,1,1};
    #define dNds(NUM,FIR,SEC,THR) (xi[NUM]*(1+eta[NUM]*SEC)*(1+gamma[NUM]*THR)/8)
    #define dNdt(NUM,FIR,SEC,THR) (eta[NUM]*(1+xi[NUM]*FIR)*(1+gamma[NUM]*THR)/8)
    #define dNdn(NUM,FIR,SEC,THR) (gamma[NUM]*(1+xi[NUM]*FIR)*(1+eta[NUM]*SEC)/8)
    //////////!!!!! TODO: Maybe using macro is better than using lambda here !!!!!//////////
    // std::function<realn(realn,realn,realn)> dxds, dxdt, dxdn, dyds, dydt, dydn, dzds, dzdt, dzdn;
    // #define Cap_List_1  xi,eta,gamma,&TV,this
    // dxds = [Cap_List_1](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += x(i)*dNds(i,s,t,n); return res;};
    // dyds = [Cap_List_1](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += y(i)*dNds(i,s,t,n); return res;};
    // dzds = [Cap_List_1](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += z(i)*dNds(i,s,t,n); return res;};
    // dxdt = [Cap_List_1](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += x(i)*dNdt(i,s,t,n); return res;};
    // dydt = [Cap_List_1](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += y(i)*dNdt(i,s,t,n); return res;};
    // dzdt = [Cap_List_1](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += z(i)*dNdt(i,s,t,n); return res;};
    // dxdn = [Cap_List_1](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += x(i)*dNdn(i,s,t,n); return res;};
    // dydn = [Cap_List_1](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += y(i)*dNdn(i,s,t,n); return res;};
    // dzdn = [Cap_List_1](realn s, realn t, realn n){realn res=0; for(int i=1; i<9; i++) res += z(i)*dNdn(i,s,t,n); return res;};
    #define dxds(FIR,SEC,THR) (x(1)*dNds(1,FIR,SEC,THR)+x(2)*dNds(2,FIR,SEC,THR)+x(3)*dNds(3,FIR,SEC,THR)+x(4)*dNds(4,FIR,SEC,THR)\
                              +x(5)*dNds(5,FIR,SEC,THR)+x(6)*dNds(6,FIR,SEC,THR)+x(7)*dNds(7,FIR,SEC,THR)+x(8)*dNds(8,FIR,SEC,THR))
    #define dyds(FIR,SEC,THR) (y(1)*dNds(1,FIR,SEC,THR)+y(2)*dNds(2,FIR,SEC,THR)+y(3)*dNds(3,FIR,SEC,THR)+y(4)*dNds(4,FIR,SEC,THR)\
                              +y(5)*dNds(5,FIR,SEC,THR)+y(6)*dNds(6,FIR,SEC,THR)+y(7)*dNds(7,FIR,SEC,THR)+y(8)*dNds(8,FIR,SEC,THR))
    #define dzds(FIR,SEC,THR) (z(1)*dNds(1,FIR,SEC,THR)+z(2)*dNds(2,FIR,SEC,THR)+z(3)*dNds(3,FIR,SEC,THR)+z(4)*dNds(4,FIR,SEC,THR)\
                              +z(5)*dNds(5,FIR,SEC,THR)+z(6)*dNds(6,FIR,SEC,THR)+z(7)*dNds(7,FIR,SEC,THR)+z(8)*dNds(8,FIR,SEC,THR))
    #define dxdt(FIR,SEC,THR) (x(1)*dNdt(1,FIR,SEC,THR)+x(2)*dNdt(2,FIR,SEC,THR)+x(3)*dNdt(3,FIR,SEC,THR)+x(4)*dNdt(4,FIR,SEC,THR)\
                              +x(5)*dNdt(5,FIR,SEC,THR)+x(6)*dNdt(6,FIR,SEC,THR)+x(7)*dNdt(7,FIR,SEC,THR)+x(8)*dNdt(8,FIR,SEC,THR))
    #define dydt(FIR,SEC,THR) (y(1)*dNdt(1,FIR,SEC,THR)+y(2)*dNdt(2,FIR,SEC,THR)+y(3)*dNdt(3,FIR,SEC,THR)+y(4)*dNdt(4,FIR,SEC,THR)\
                              +y(5)*dNdt(5,FIR,SEC,THR)+y(6)*dNdt(6,FIR,SEC,THR)+y(7)*dNdt(7,FIR,SEC,THR)+y(8)*dNdt(8,FIR,SEC,THR))
    #define dzdt(FIR,SEC,THR) (z(1)*dNdt(1,FIR,SEC,THR)+z(2)*dNdt(2,FIR,SEC,THR)+z(3)*dNdt(3,FIR,SEC,THR)+z(4)*dNdt(4,FIR,SEC,THR)\
                              +z(5)*dNdt(5,FIR,SEC,THR)+z(6)*dNdt(6,FIR,SEC,THR)+z(7)*dNdt(7,FIR,SEC,THR)+z(8)*dNdt(8,FIR,SEC,THR))
    #define dxdn(FIR,SEC,THR) (x(1)*dNdn(1,FIR,SEC,THR)+x(2)*dNdn(2,FIR,SEC,THR)+x(3)*dNdn(3,FIR,SEC,THR)+x(4)*dNdn(4,FIR,SEC,THR)\
                              +x(5)*dNdn(5,FIR,SEC,THR)+x(6)*dNdn(6,FIR,SEC,THR)+x(7)*dNdn(7,FIR,SEC,THR)+x(8)*dNdn(8,FIR,SEC,THR))
    #define dydn(FIR,SEC,THR) (y(1)*dNdn(1,FIR,SEC,THR)+y(2)*dNdn(2,FIR,SEC,THR)+y(3)*dNdn(3,FIR,SEC,THR)+y(4)*dNdn(4,FIR,SEC,THR)\
                              +y(5)*dNdn(5,FIR,SEC,THR)+y(6)*dNdn(6,FIR,SEC,THR)+y(7)*dNdn(7,FIR,SEC,THR)+y(8)*dNdn(8,FIR,SEC,THR))
    #define dzdn(FIR,SEC,THR) (z(1)*dNdn(1,FIR,SEC,THR)+z(2)*dNdn(2,FIR,SEC,THR)+z(3)*dNdn(3,FIR,SEC,THR)+z(4)*dNdn(4,FIR,SEC,THR)\
                              +z(5)*dNdn(5,FIR,SEC,THR)+z(6)*dNdn(6,FIR,SEC,THR)+z(7)*dNdn(7,FIR,SEC,THR)+z(8)*dNdn(8,FIR,SEC,THR))
    /* Calculate Jacobi matrix */
    #define Cap_list_2  xi,eta,gamma,&TV,this/*,dxds,dxdt,dxdn,dyds,dydt,dydn,dzds,dzdt,dzdn*/
    Jacobi = [Cap_list_2](realn s, realn t, realn n){
        MatrixXr result(3,3);
        result << dxds(s,t,n), dyds(s,t,n), dzds(s,t,n),
                  dxdt(s,t,n), dydt(s,t,n), dzdt(s,t,n),
                  dxdn(s,t,n), dydn(s,t,n), dzdn(s,t,n);
        return result;
    };
    /* Calculate B matrix */
    #define coef_a(FIR,SEC,THR) (dydt(FIR,SEC,THR)*dzdn(FIR,SEC,THR)-dzdt(FIR,SEC,THR)*dydn(FIR,SEC,THR))
    #define coef_b(FIR,SEC,THR) (dyds(FIR,SEC,THR)*dzdn(FIR,SEC,THR)-dzds(FIR,SEC,THR)*dydn(FIR,SEC,THR))
    #define coef_c(FIR,SEC,THR) (dyds(FIR,SEC,THR)*dzdt(FIR,SEC,THR)-dzds(FIR,SEC,THR)*dydt(FIR,SEC,THR))
    #define coef_d(FIR,SEC,THR) (dxdt(FIR,SEC,THR)*dzdn(FIR,SEC,THR)-dzdt(FIR,SEC,THR)*dxdn(FIR,SEC,THR))
    #define coef_e(FIR,SEC,THR) (dxds(FIR,SEC,THR)*dzdn(FIR,SEC,THR)-dzds(FIR,SEC,THR)*dxdn(FIR,SEC,THR))
    #define coef_f(FIR,SEC,THR) (dxds(FIR,SEC,THR)*dzdt(FIR,SEC,THR)-dzds(FIR,SEC,THR)*dxdt(FIR,SEC,THR))
    #define coef_g(FIR,SEC,THR) (dxdt(FIR,SEC,THR)*dydn(FIR,SEC,THR)-dydt(FIR,SEC,THR)*dxdn(FIR,SEC,THR))
    #define coef_h(FIR,SEC,THR) (dxds(FIR,SEC,THR)*dydn(FIR,SEC,THR)-dyds(FIR,SEC,THR)*dxdn(FIR,SEC,THR))
    #define coef_i(FIR,SEC,THR) (dxds(FIR,SEC,THR)*dydt(FIR,SEC,THR)-dyds(FIR,SEC,THR)*dxdt(FIR,SEC,THR))
#ifndef NOT_KEEP_BMat
    #define Cap_list_3  xi,eta,gamma,&TV,this/*,dxds,dxdt,dxdn,dyds,dydt,dydn,dzds,dzdt,dzdn*/
#else
    #define Cap_list_3  xi,eta,gamma,&TV,this,&Jacobi/*,dxds,dxdt,dxdn,dyds,dydt,dydn,dzds,dzdt,dzdn*/
#endif
    ele_BMat = [Cap_list_3](realn s, realn t, realn n){
        MatrixXr res(6,24);
        for(int i=0; i<8; i++){
            res(0, 3*i+0) =  coef_a(s,t,n)*dNds(i+1,s,t,n) - coef_b(s,t,n)*dNdt(i+1,s,t,n) + coef_c(s,t,n)*dNdn(i+1,s,t,n);
            res(0, 3*i+1) =  0;
            res(0, 3*i+2) =  0;
            res(1, 3*i+0) =  0;
            res(1, 3*i+1) = -coef_d(s,t,n)*dNds(i+1,s,t,n) + coef_e(s,t,n)*dNdt(i+1,s,t,n) - coef_f(s,t,n)*dNdn(i+1,s,t,n);
            res(1, 3*i+2) =  0;
            res(2, 3*i+0) =  0;
            res(2, 3*i+1) =  0;
            res(2, 3*i+2) =  coef_g(s,t,n)*dNds(i+1,s,t,n) - coef_h(s,t,n)*dNdt(i+1,s,t,n) + coef_i(s,t,n)*dNdn(i+1,s,t,n);
            res(3, 3*i+0) = -coef_d(s,t,n)*dNds(i+1,s,t,n) + coef_e(s,t,n)*dNdt(i+1,s,t,n) - coef_f(s,t,n)*dNdn(i+1,s,t,n);
            res(3, 3*i+1) =  coef_a(s,t,n)*dNds(i+1,s,t,n) - coef_b(s,t,n)*dNdt(i+1,s,t,n) + coef_c(s,t,n)*dNdn(i+1,s,t,n);
            res(3, 3*i+2) =  0;
            res(4, 3*i+0) =  0;
            res(4, 3*i+1) =  coef_g(s,t,n)*dNds(i+1,s,t,n) - coef_h(s,t,n)*dNdt(i+1,s,t,n) + coef_i(s,t,n)*dNdn(i+1,s,t,n);
            res(4, 3*i+2) = -coef_d(s,t,n)*dNds(i+1,s,t,n) + coef_e(s,t,n)*dNdt(i+1,s,t,n) - coef_f(s,t,n)*dNdn(i+1,s,t,n);
            res(5, 3*i+0) =  coef_g(s,t,n)*dNds(i+1,s,t,n) - coef_h(s,t,n)*dNdt(i+1,s,t,n) + coef_i(s,t,n)*dNdn(i+1,s,t,n);
            res(5, 3*i+1) =  0;
            res(5, 3*i+2) =  coef_a(s,t,n)*dNds(i+1,s,t,n) - coef_b(s,t,n)*dNdt(i+1,s,t,n) + coef_c(s,t,n)*dNdn(i+1,s,t,n);
        }
        res = res/Jacobi(s,t,n).determinant();
        return res;
    };
    #undef Cap_list_3
    #undef coef_i
    #undef coef_h
    #undef coef_g
    #undef coef_f
    #undef coef_e
    #undef coef_d
    #undef coef_c
    #undef coef_b
    #undef coef_a
    #undef Cap_list_2
    #undef dzdn
    #undef dydn
    #undef dxdn
    #undef dzdt
    #undef dydt
    #undef dxdt
    #undef dzds
    #undef dyds
    #undef dxds
    // #undef Cap_list_1
    #undef dNdn
    #undef dNdt
    #undef dNds
    #undef z
    #undef y
    #undef x
}

} // namespace lhfea


#endif