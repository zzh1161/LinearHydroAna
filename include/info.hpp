#ifndef FEA_INFO_HPP
#define FEA_INFO_HPP

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>

namespace lhfea{

using realn = double;

#define MAKE_TYPEDEFS(Type, TypeSuffix, Size, SizeSuffix)                  \
typedef Eigen::Matrix<Type, Size, Size> Matrix##SizeSuffix##TypeSuffix;    \
typedef Eigen::Matrix<Type, Size, 1>    Vector##SizeSuffix##TypeSuffix;    \
typedef Eigen::Matrix<Type, 1, Size>    RowVector##SizeSuffix##TypeSuffix;

// MAKE_TYPEDEFS(real, r, 2, 2)
// MAKE_TYPEDEFS(real, r, 3, 3)
// MAKE_TYPEDEFS(real, r, 4, 4)
MAKE_TYPEDEFS(realn, r, -1, X)

#undef MAKE_TYPEDEFS

using SpMatXr = Eigen::SparseMatrix<realn>;
using TripXr  = Eigen::Triplet<realn>;

enum class SolverType{
    Eigen = 1,
    cholmod = 2
};

} // namespace lhfea


#endif