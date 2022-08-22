#define USE_CHOLMOD
#ifdef USE_CHOLMOD

#ifndef CHOLMOD_SOLVER_HPP
#define CHOLMOD_SOLVER_HPP

#include "LinSysSolver.hpp"

#include <suitesparse/cholmod.h>
#include <eigen3/Eigen/Eigen>

#include <vector>
#include <set>

namespace SIM {

template <typename vectorTypeI, typename vectorTypeS>
class CHOLMODSolver : public LinSysSolver<vectorTypeI, vectorTypeS> {
    typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:
    cholmod_common cm;
    cholmod_sparse* A;
    cholmod_factor* L;
    cholmod_dense *b, *solution;
    cholmod_dense *x_cd, *y_cd; // for multiply

    void *Ai, *Ap, *Ax, *bx, *solutionx, *x_cdx, *y_cdx;

public:
    CHOLMODSolver(void);
    ~CHOLMODSolver(void);

    void set_pattern(const std::vector<std::set<int>>& vNeighbor);
    void set_pattern(const Eigen::SparseMatrix<double>& mtr); //NOTE: mtr must be SPD
    void load(const char* filePath, Eigen::VectorXd& rhs);
    void write(const char *filePath, const Eigen::VectorXd& rhs);

    void analyze_pattern(void);

    void setZero(){Base::setZero();}
    void setCoeff(int rowI, int colI, double val)/*{ Base::setCoeff(rowI, colI, val); }*/;
    void addCoeff(int rowI, int colI, double val)/*{ Base::addCoeff(rowI, colI, val); }*/;
    void setUnit_row(int rowI)/*{ Base::setUnit_row(rowI); }*/;
    void setUnit_col(int colI, const std::set<int>& rowVIs)/*{ Base::setUnit_col(colI, rowVIs); }*/;

    bool factorize(void);

    void solve(Eigen::VectorXd& rhs, Eigen::VectorXd& result);

    virtual void multiply(const Eigen::VectorXd& x, Eigen::VectorXd& Ax);

    virtual void outputFactorization(const std::string& filePath);
};

template <typename vectorTypeI, typename vectorTypeS>
CHOLMODSolver<vectorTypeI, vectorTypeS>::CHOLMODSolver(void)
{
    cholmod_start(&cm);
    A = NULL;
    L = NULL;
    b = NULL;
    x_cd = y_cd = NULL;

    Ai = Ap = Ax = NULL;
    bx = NULL;
    solutionx = x_cdx = y_cdx = NULL;
}

template <typename vectorTypeI, typename vectorTypeS>
CHOLMODSolver<vectorTypeI, vectorTypeS>::~CHOLMODSolver(void)
{
    if (A) {
        A->i = Ai;
        A->p = Ap;
        A->x = Ax;
        cholmod_free_sparse(&A, &cm);
    }
    cholmod_free_factor(&L, &cm);
    if (b) {
        b->x = bx;
        cholmod_free_dense(&b, &cm);
    }
    if (x_cd) {
        x_cd->x = x_cdx;
        cholmod_free_dense(&x_cd, &cm);
    }
    if (y_cd) {
        y_cd->x = y_cdx;
        cholmod_free_dense(&y_cd, &cm);
    }
    cholmod_finish(&cm);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::set_pattern(const std::vector<std::set<int>>& vNeighbor)
{
    Base::set_pattern(vNeighbor);

    //TODO: directly save into A
    if (!A) {
        A = cholmod_allocate_sparse(Base::numRows, Base::numRows, Base::ja.size(),
            true, true, -1, CHOLMOD_REAL, &cm);
        Ax = A->x;
        Ap = A->p;
        Ai = A->i;
        // -1: upper right part will be ignored during computation
    }
    Base::ia.array() -= 1;
    Base::ja.array() -= 1; // CHOLMOD's index starts from 0
    A->i = Base::ja.data();
    A->p = Base::ia.data();
    A->x = Base::a.data();
}
template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::set_pattern(const Eigen::SparseMatrix<double>& mtr)
{
    Base::set_pattern(mtr);

    if (!A) {
        A = cholmod_allocate_sparse(Base::numRows, Base::numRows, mtr.nonZeros(),
            true, true, -1, CHOLMOD_REAL, &cm);
        Ax = A->x;
        Ap = A->p;
        Ai = A->i;
        // -1: upper right part will be ignored during computation

        A->i = Base::ja.data();
        A->p = Base::ia.data();
        A->x = Base::a.data();
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::load(const char* filePath, Eigen::VectorXd& rhs)
{
    Base::load(filePath, rhs);

    //TODO: directly save into A
    if (!A) {
        A = cholmod_allocate_sparse(Base::numRows, Base::numRows, Base::ja.size(),
            true, true, -1, CHOLMOD_REAL, &cm);
        Ax = A->x;
        Ap = A->p;
        Ai = A->i;
        // -1: upper right part will be ignored during computation
    }
    Base::ia.array() -= 1;
    Base::ja.array() -= 1; // CHOLMOD's index starts from 0
    A->i = Base::ja.data();
    A->p = Base::ia.data();
    A->x = Base::a.data();
    cm.supernodal = CHOLMOD_SIMPLICIAL;
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::write(const char *filePath, const Eigen::VectorXd& rhs)
{
    FILE *out = fopen(filePath, "w");
    cholmod_write_sparse(out, A, NULL, NULL, &cm);
    fclose(out);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::analyze_pattern(void)
{
    // std::cout << getCurrentRSS() << std::endl;
    cholmod_free_factor(&L, &cm);
    L = cholmod_analyze(A, &cm);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::setCoeff(int rowI, int colI, double val)
{
    if(rowI <= colI){
        assert(rowI < Base::IJ2aI.size());
        const auto finder = Base::IJ2aI[rowI].find(colI);
        assert(finder != Base::IJ2aI[rowI].end());
        Base::a[finder->second] = val;
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::addCoeff(int rowI, int colI, double val)
{
    if(rowI <= colI){
        assert(rowI < Base::IJ2aI.size());
        const auto finder = Base::IJ2aI[rowI].find(colI);
        assert(finder != Base::IJ2aI[rowI].end());
        Base::a[finder->second] += val;
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::setUnit_row(int rowI)
{
    assert(Base::numRows == Base::IJ2aI.size());
    assert(rowI < Base::numRows);
    for (const auto& colIter : Base::IJ2aI[rowI]) {
        Base::a[colIter.second] = (colIter.first == rowI);
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::setUnit_col(int colI, const std::set<int>& rowVIs)
{
    assert(Base::numRows == Base::IJ2aI.size());
    assert(colI < Base::numRows);
    for (const auto& rowVI : rowVIs) {
        for (int dimI = 0; dimI < DIM_; ++dimI) {
            int rowI = rowVI * DIM_ + dimI;
            if(rowI <= colI){
                const auto finder = Base::IJ2aI[rowI].find(colI);
                if (finder != Base::IJ2aI[rowI].end()) {
                    Base::a[finder->second] = (rowI == colI);
                }
            }
        }
    }
}

template <typename vectorTypeI, typename vectorTypeS>
bool CHOLMODSolver<vectorTypeI, vectorTypeS>::factorize(void)
{
    cholmod_factorize(A, L, &cm);
    // std::cout << getCurrentRSS() << std::endl;
    // exit(0);
    return cm.status != CHOLMOD_NOT_POSDEF;
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::solve(Eigen::VectorXd& rhs,
    Eigen::VectorXd& result)
{
    //TODO: directly point to rhs?
    if (!b) {
        b = cholmod_allocate_dense(Base::numRows, 1, Base::numRows, CHOLMOD_REAL, &cm);
        bx = b->x;
    }
    b->x = rhs.data();
    cholmod_dense* x;
    x = cholmod_solve(CHOLMOD_A, L, b, &cm);
    result.conservativeResize(rhs.size());
    memcpy(result.data(), x->x, result.size() * sizeof(result[0]));
    cholmod_free_dense(&x, &cm);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::multiply(const Eigen::VectorXd& x,
    Eigen::VectorXd& Ax)
{
    assert(x.size() == Base::numRows);

    if (!x_cd) {
        x_cd = cholmod_allocate_dense(Base::numRows, 1, Base::numRows,
            CHOLMOD_REAL, &cm);
        x_cdx = x_cd->x;
    }
    x_cd->x = (void*)x.data();

    Ax.conservativeResize(Base::numRows);
    if (!y_cd) {
        y_cd = cholmod_allocate_dense(Base::numRows, 1, Base::numRows,
            CHOLMOD_REAL, &cm);
        y_cdx = y_cd->x;
    }
    y_cd->x = (void*)Ax.data();

    double alpha[2] = { 1.0, 1.0 }, beta[2] = { 0.0, 0.0 };

    cholmod_sdmult(A, 0, alpha, beta, x_cd, y_cd, &cm);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::outputFactorization(const std::string& filePath)
{
    cholmod_sparse* spm = cholmod_factor_to_sparse(L, &cm);

    FILE* out = fopen(filePath.c_str(), "w");
    assert(out);

    cholmod_write_sparse(out, spm, NULL, "", &cm);

    fclose(out);
}

} // namespace SIM

#endif /* CHOLMODSolver_hpp */

#endif /* USE_CHOLMOD */
