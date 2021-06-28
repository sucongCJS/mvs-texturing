#ifndef TEX_UTIL_HEADER
#define TEX_UTIL_HEADER

#include <string>
#include <vector>
#include <fstream>
#include <cstring>
#include <cerrno>
#include <cassert>

#include <Eigen/Core>
#include <util/system.h>
#include <util/exception.h>
#include <util/file_system.h>

#include "math/vector.h"
#include "math/matrix.h"
#include "math/functions.h"

inline void remove_tmp(std::string const tmp_dir){
    for (util::fs::File const & file : util::fs::Directory(tmp_dir)) {
        util::fs::unlink(util::fs::join_path(file.path, file.name).c_str());
    }
    util::fs::rmdir(tmp_dir.c_str());
}

/**
  * Write vector to binary file.
  * @throws util::FileException
  */
template <typename T> void
vector_to_file(std::string const & filename, std::vector<T> const & vector) {
    std::ofstream out(filename.c_str(), std::ios::binary);
    if (!out.good())
        throw util::FileException(filename, std::strerror(errno));

    out.write(reinterpret_cast<const char*>(&vector[0]), vector.size()*sizeof(T));
    out.close();
}

/**
  * Loads vector from binary file.
  * @throws util::FileException
  */
template <typename T> std::vector<T>
vector_from_file(std::string const & filename, std::string const tmp_dir) {
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in.good()){
        remove_tmp(tmp_dir);
        throw util::FileException(filename, std::strerror(errno));
    }
    in.seekg (0, in.end);
    const size_t filesize = in.tellg();
    in.seekg (0, in.beg);
    const size_t num_elements = filesize / sizeof(T);
    std::vector<T> vector(num_elements);
    in.read(reinterpret_cast<char*>(&vector[0]), num_elements*sizeof(T));
    in.close();
    return vector;
}

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
 * @return f$\exp(-\frac{1}{2} (X-\mbox{mu})^T \cdot \mbox{covariance\_inv} \cdot (X-\mbox{mu}))\f$
 */
template <typename T, int N> T const
multi_gauss_unnormalized(Eigen::Matrix<T, 1, N> const & X, Eigen::Matrix<T, 1, N> const & mu,
    Eigen::Matrix<T, N, N> const & covariance_inv) {

    Eigen::Matrix<T, 1, N> mean_removed = X - mu;
    return std::exp(T(-0.5) * mean_removed * covariance_inv * mean_removed.adjoint());
}

#endif /* TEX_UTIL_HEADER */
