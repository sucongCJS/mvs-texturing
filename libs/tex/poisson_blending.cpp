/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <cstdint>
#include <iostream>

#include <math/vector.h>
#include <Eigen/SparseCore>  // SparseMatrix and SparseVector classes
#include <Eigen/SparseLU>  // factorization to solve general square sparse systems

#include "poisson_blending.h"

typedef Eigen::SparseMatrix<float> SpMat;

// i是像素位置
math::Vec3f simple_laplacian(int i, mve::FloatImage::ConstPtr img){
    const int width = img->width();
    assert(i > width + 1 && i < img->get_pixel_amount() - width -1);  // i不能是img最外那一圈像素

    return -4.0f * math::Vec3f(&img->at(i, 0))
        + math::Vec3f(&img->at(i - width, 0))
        + math::Vec3f(&img->at(i - 1, 0))
        + math::Vec3f(&img->at(i + 1, 0))
        + math::Vec3f(&img->at(i + width, 0));
}

bool valid_mask(mve::ByteImage::ConstPtr mask){
    const int width = mask->width();
    const int height = mask->height();

    // mask的最外一圈像素不能为255??
    for (int x = 0; x < width; ++x)
        if (mask->at(x, 0, 0) == 255 || mask->at(x, height - 1, 0) == 255)
            return false;
    for (int y = 0; y < height; ++y)
        if (mask->at(0, y, 0) == 255 || mask->at(width - 1, y, 0) == 255)
            return false;

    //TODO check for sane boundary conditions...

    return true;
}

// src贴到dest上
// mask是ByteImage, 能存0~255的数
void poisson_blend(mve::FloatImage::ConstPtr src, mve::ByteImage::ConstPtr mask,
                    mve::FloatImage::Ptr dest, float alpha) {

    assert(src->width() == mask->width() && mask->width() == dest->width());
    assert(src->height() == mask->height() && mask->height() == dest->height());
    assert(src->channels() == 3 && dest->channels() == 3);
    assert(mask->channels() == 1);
    assert(valid_mask(mask));

    const int n = dest->get_pixel_amount();
    const int width = dest->width();
    const int height = dest->height();
    const int channels = dest->channels();

    mve::Image<int>::Ptr indices = mve::Image<int>::create(width, height, 1);
    indices->fill(-1);  // 如果在mask内就附上一个值(从0开始赋), 否则就为-1
    int index = 0;  // mask内的像素
    for (int i = 0; i < n; ++i) {
        if (mask->at(i) != 0) {
            indices->at(i) = index;
            index += 1;
        }
    }
    const int nnz = index;  // mask内像素数量

    std::vector<math::Vec3f> coefficients_b;
    coefficients_b.resize(nnz);  // 存储混合后的mask那部分像素的Laplacian散度

    std::vector<Eigen::Triplet<float, int> > coefficients_A;  // triplets是Eigen中定义的一种用于储存稀疏矩阵的非零值数据结构 <flaot, int>是给谁的类型?: <存的数的类型, 索引的类型>
    coefficients_A.reserve(nnz);// 预留空间, 但不创建元素对象, size也不会改变 //TODO better estimate...

    for (int i = 0; i < n; ++i) {  // 遍历每个像素
        const int row = indices->at(i);  // 图片中的像素在mask中的位置

        // 什么情况mask会变128, 64?? 边界
        if (mask->at(i) == 128 || mask->at(i) == 64) {  // mask会变??
            Eigen::Triplet<float, int> t(row, row, 1.0f);
            coefficients_A.push_back(t);

            coefficients_b[row] = math::Vec3f(&dest->at(i, 0));  // 边界的颜色是原来图片的颜色
        }

        // 只有是白色的才会被插入到coefficients_A
        if (mask->at(i) == 255) {
            const int i01 = indices->at(i - width);
            const int i10 = indices->at(i - 1);
            const int i11 = indices->at(i);
            const int i12 = indices->at(i + 1);
            const int i21 = indices->at(i + width);

            /* All neighbours should be eighter border conditions or part of the optimization. */
            assert(i01 != -1 && i10 != -1 && i11 != -1 && i12 != -1 && i21 != -1);  // 必须都在mask内

            Eigen::Triplet<float, int> t01(row, i01, 1.0f);  // row不是mask内第几个像素吗怎么是行呢?: 矩阵的行和列数等于mask内所有像素的数量
            Eigen::Triplet<float, int> t10(row, i10, 1.0f);  // i**也是maks内第几个像素
            Eigen::Triplet<float, int> t11(row, i11, -4.0f);
            Eigen::Triplet<float, int> t12(row, i12, 1.0f);
            Eigen::Triplet<float, int> t21(row, i21, 1.0f);
            Eigen::Triplet<float, int> triplets[] = {t01, t10, t11, t12, t21};

            coefficients_A.insert(coefficients_A.end(), triplets, triplets + 5);  // 在结尾插入triplets列表中的5个元素

            math::Vec3f l_d = simple_laplacian(i, dest);
            math::Vec3f l_s = simple_laplacian(i, src);

            coefficients_b[row] = (alpha * l_s + (1.0f - alpha) * l_d);  // 将dest和src的Laplacian做一个composition. alpha是啥??
        }
    }

    // ??不用使用pyramid吗

    SpMat A(nnz, nnz);  // 可以看 https://blog.csdn.net/hjimce/article/details/45716603 的图, 每一行要么存有一个4, 四个1, 其余都是0, 要么存的是一个1, 其余都是0
    A.setFromTriplets(coefficients_A.begin(), coefficients_A.end());  // 给稀疏矩阵赋值

    Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> > solver;
    solver.compute(A);  // 求解x  A

    for (int channel = 0; channel < channels; ++channel) {
        Eigen::VectorXf b(nnz);  // 保留一个通道的coefficients_b
        for (std::size_t i = 0; i < coefficients_b.size(); ++i)
            b[i] = coefficients_b[i][channel];

        Eigen::VectorXf x(n);  // 保留一个通道的合成后的mask部分图片颜色
        x = solver.solve(b);

        // 修改dest颜色, 将src的颜色贴上去
        for (int i = 0; i < n; ++i) {
            int index = indices->at(i);
            if (index != -1) 
                dest->at(i, channel) = x[index];
        }
    }
}
