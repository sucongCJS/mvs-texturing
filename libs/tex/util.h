#ifndef TEX_UTIL_HEADER
#define TEX_UTIL_HEADER

#include <string>
#include <vector>
#include <fstream>
#include <cstring>
#include <cerrno>
#include <cassert>

#include <Eigen/Core>

#include "util/exception.h"
#include "util/file_system.h"

#include "math/vector.h"
#include "math/matrix.h"
#include "math/functions.h"

// 将MVE向量转换为Eigen行向量
template <typename T, int N>
Eigen::Matrix<T, 1, N> mve_to_eigen(math::Vector<T, N> const & vec) {
    Eigen::Matrix<T, 1, N> ret;
    for (int n = 0; n < N; ++n)
        ret(0, n) = vec(n);
    return ret;
}

/**
 * X 1*N, mu 1*N, covariance_inv N*N
 * 返回 f$\exp(-\frac{1}{2} (X-\mbox{mu})^T \cdot \mbox{covariance\_inv} \cdot (X-\mbox{mu}))\f$
 */
template <typename T, int N> T const
multi_gauss_unnormalized(Eigen::Matrix<T, 1, N> const & X, Eigen::Matrix<T, 1, N> const & mu,
    Eigen::Matrix<T, N, N> const & covariance_inv) {

    Eigen::Matrix<T, 1, N> mean_removed = X - mu;
    return std::exp(T(-0.5) * mean_removed * covariance_inv * mean_removed.adjoint());
}

#endif /* TEX_UTIL_HEADER */
