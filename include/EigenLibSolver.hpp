#ifndef EIGEN_LIB_SOLVER_HPP
#define EIGEN_LIB_SOLVER_HPP

#include "LinSysSolver.hpp"

#include <eigen3/Eigen/Eigen>

#include <vector>
#include <set>

namespace SIM {

template <typename vectorTypeI, typename vectorTypeS>
class EigenLibSolver : public LinSysSolver<vectorTypeI, vectorTypeS> 
{
    typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:
    Eigen::SparseMatrix<double> coefMtr;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> simplicialLDLT;

public:
    void set_pattern(const std::vector<std::set<int>>& vNeighbor);
    void set_pattern(const Eigen::SparseMatrix<double>& mtr); //NOTE: mtr must be SPD
    void write(const char *filePath, const Eigen::VectorXd& rhs);

    void analyze_pattern(void);

    bool factorize(void);

    void solve(Eigen::VectorXd& rhs, Eigen::VectorXd& result);

    double coeffMtr(int rowI, int colI) const;

    void setZero(void);

    virtual void setCoeff(int rowI, int colI, double val);

    virtual void addCoeff(int rowI, int colI, double val);

    virtual void setUnit_row(int rowI);

    virtual void setUnit_col(int colI, const std::set<int>& rowVIs);
};

template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::set_pattern(const std::vector<std::set<int>>& vNeighbor)
{
    Base::set_pattern(vNeighbor);

    //TODO: directly save into mtr
    coefMtr.resize(Base::numRows, Base::numRows);
    coefMtr.reserve(Base::ja.size());
    Base::ia.array() -= 1.0;
    Base::ja.array() -= 1.0;
    memcpy(coefMtr.innerIndexPtr(), Base::ja.data(), Base::ja.size() * sizeof(Base::ja[0]));
    memcpy(coefMtr.outerIndexPtr(), Base::ia.data(), Base::ia.size() * sizeof(Base::ia[0]));
}
template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::set_pattern(const Eigen::SparseMatrix<double>& mtr) //NOTE: mtr must be SPD
{
    Base::numRows = static_cast<int>(mtr.rows());
    coefMtr = mtr;
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::write(const char *filePath, const Eigen::VectorXd& rhs)
{
    std::ofstream out(filePath);
    out << coefMtr;
    out.close();
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::analyze_pattern(void)
{
    simplicialLDLT.analyzePattern(coefMtr);
    assert(simplicialLDLT.info() == Eigen::Success);
}

template <typename vectorTypeI, typename vectorTypeS>
bool EigenLibSolver<vectorTypeI, vectorTypeS>::factorize(void)
{
    bool succeeded = false;
    simplicialLDLT.factorize(coefMtr);
    succeeded = (simplicialLDLT.info() == Eigen::Success);
    assert(succeeded);
    return succeeded;
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::solve(Eigen::VectorXd& rhs, Eigen::VectorXd& result)
{
    result = simplicialLDLT.solve(rhs);
    assert(simplicialLDLT.info() == Eigen::Success);
}

template <typename vectorTypeI, typename vectorTypeS>
double EigenLibSolver<vectorTypeI, vectorTypeS>::coeffMtr(int rowI, int colI) const
{
    return Base::coeffMtr(rowI, colI);
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::setZero(void)
{
    //TODO: directly manipulate valuePtr without a
    Base::setZero();
    memcpy(coefMtr.valuePtr(), Base::a.data(), Base::a.size() * sizeof(Base::a[0]));
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::setCoeff(int rowI, int colI, double val)
{
    //TODO: directly manipulate valuePtr without a

    if (rowI <= colI) {
        assert(rowI < Base::IJ2aI.size());
        const auto finder = Base::IJ2aI[rowI].find(colI);
        assert(finder != Base::IJ2aI[rowI].end());
        Base::a[finder->second] = val;
        coefMtr.valuePtr()[finder->second] = val;
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::addCoeff(int rowI, int colI, double val)
{
    //TODO: directly manipulate valuePtr without a

    if (rowI <= colI) {
        assert(rowI < Base::IJ2aI.size());
        const auto finder = Base::IJ2aI[rowI].find(colI);
        assert(finder != Base::IJ2aI[rowI].end());
        Base::a[finder->second] += val;
        coefMtr.valuePtr()[finder->second] += val;
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::setUnit_row(int rowI)
{
    for (const auto& colIter : Base::IJ2aI[rowI]) {
        coefMtr.valuePtr()[colIter.second] = (colIter.first == rowI);
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenLibSolver<vectorTypeI, vectorTypeS>::setUnit_col(int colI, const std::set<int>& rowVIs)
{
    for (const auto& rowVI : rowVIs) {
        for (int dimI = 0; dimI < DIM_; ++dimI) {
            int rowI = rowVI * DIM_ + dimI;
            if (rowI <= colI) {
                const auto finder = Base::IJ2aI[rowI].find(colI);
                if (finder != Base::IJ2aI[rowI].end()) {
                    coefMtr.valuePtr()[finder->second] = (rowI == colI);
                }
            }
        }
    }
}

} // namespace SIM

#endif