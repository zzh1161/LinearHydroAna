#ifndef HEX_8_NODE_ELE_HPP
#define HEX_8_NODE_ELE_HPP

#include "shapeInterface.hpp"

namespace lhfea{

struct hex8node : public shapeInterface
{
    // realn ele_vol; // not necessary

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
    void calBMatrix(MatrixXr &TV, realn E, realn nu, MatrixXr &ele_BMat);
#endif
    void calDMatrix(MatrixXr &TV, realn E, realn nu, MatrixXr &ele_DMat);
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
void hex8node::calBMatrix(MatrixXr &TV, realn E, realn nu, MatrixXr &ele_BMat)
#endif
{

}

void hex8node::calDMatrix(MatrixXr &TV, realn E, realn nu, MatrixXr &ele_DMat)
{

}

} // namespace lhfea


#endif