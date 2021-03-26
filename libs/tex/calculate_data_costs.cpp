/*
 * Copyright (C) 2015, Nils Moehrle, Michael Waechter
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <numeric>

#include <mve/image_color.h>
#include <acc/bvh_tree.h>
#include <Eigen/Core>
#include <Eigen/LU>

#include "util.h"
#include "histogram.h"
#include "texturing.h"
#include "sparse_table.h"
#include "progress_counter.h"

typedef acc::BVHTree<unsigned int, math::Vec3f> BVHTree;

TEX_NAMESPACE_BEGIN

/**
 * Dampens the quality of all views in which the face's projection
 * has a much different color than in the majority of views.
 *
 * @param infos contains information about one face seen from several views
 * @param settings runtime configuration.
 * @return whether the outlier removal was successfull.
 */
bool photometric_outlier_detection(std::vector<FaceProjectionInfo> * infos, Settings const & settings) {
    if (infos->size() == 0) return true;  // 这个面没有被任何view看到

    /* Configuration variables. */

    double const gauss_rejection_threshold = 6e-3;

    /* If all covariances drop below this we stop outlier detection. */
    double const minimal_covariance = 5e-4;  // 

    int const outlier_detection_iterations = 10;
    int const minimal_num_inliers = 4;

    float outlier_removal_factor = std::numeric_limits<float>::signaling_NaN();
    switch (settings.outlier_removal) {
        case OUTLIER_REMOVAL_NONE: return true;
        case OUTLIER_REMOVAL_GAUSS_CLAMPING: outlier_removal_factor = 1.0f; break;  // ??
        case OUTLIER_REMOVAL_GAUSS_DAMPING:  outlier_removal_factor = 0.2f; break;
    }

    Eigen::MatrixX3d inliers(infos->size(), 3);  // view的个数 * 3, 这个面在每个view的mean color
    std::vector<std::uint32_t> is_inlier(infos->size(), 1);  // 初始化为1, size为view的个数
    for (std::size_t row = 0; row < infos->size(); ++row) {
        inliers.row(row) = mve_to_eigen(infos->at(row).mean_color).cast<double>();
    }

    Eigen::RowVector3d var_mean;
    Eigen::Matrix3d covariance;
    Eigen::Matrix3d covariance_inv;

    for (int i = 0; i < outlier_detection_iterations; ++i) {

        if (inliers.rows() < minimal_num_inliers)  // inlier小于一定数目就不剔除了
            return false;

        /* Calculate the inliers' mean color and color covariance. */
        var_mean = inliers.colwise().mean();  // 计算三个通道的mean值
        Eigen::MatrixX3d centered = inliers.rowwise() - var_mean;  // 和mean作差
        covariance = (centered.adjoint() * centered) / double(inliers.rows() - 1);  // 协方差公式, 越小越不相关 ?为什么要减var_mean:因为公式

        /* If all covariances are very small we stop outlier detection 
         * and only keep the inliers (set quality of outliers to zero). 如果最大值都小于某个阈值 */
        if (covariance.array().abs().maxCoeff() < minimal_covariance) {  // 给矩阵的每个元素加绝对值 
            for (std::size_t row = 0; row < infos->size(); ++row) {
                if (!is_inlier[row])  // outliers 的 quality为0
                    infos->at(row).quality = 0.0f;
            }
            return true;
        }

        /* Invert the covariance. FullPivLU is not the fastest way but
         * it gives feedback about numerical stability during inversion. */
        Eigen::FullPivLU<Eigen::Matrix3d> lu(covariance);
        if (!lu.isInvertible()) return false;

        covariance_inv = lu.inverse();

        /* Compute new number of inliers (all views with a gauss value above a threshold). */
        for (std::size_t row = 0; row < infos->size(); ++row) {
            Eigen::RowVector3d color = mve_to_eigen(infos->at(row).mean_color).cast<double>();
            double gauss_value = multi_gauss_unnormalized(color, var_mean, covariance_inv);
            is_inlier[row] = (gauss_value >= gauss_rejection_threshold ? 1 : 0);  // 越大越好
        }
        
        /* Resize Eigen matrix accordingly and fill with new inliers. */
        inliers.resize(std::accumulate(is_inlier.begin(), is_inlier.end(), 0), Eigen::NoChange);  // 缩小inliers列表, 有多少个is_inlier, inlier列表就多大
        // 把筛剩下的inlier放入inliers列表中
        for (std::size_t row = 0, inlier_row = 0; row < infos->size(); ++row) {
            if (is_inlier[row]) {
                inliers.row(inlier_row++) = mve_to_eigen(infos->at(row).mean_color).cast<double>();
            }
        }
    }

    covariance_inv *= outlier_removal_factor;
    for (FaceProjectionInfo & info : *infos) {
        Eigen::RowVector3d color = mve_to_eigen(info.mean_color).cast<double>();
        double gauss_value = multi_gauss_unnormalized(color, var_mean, covariance_inv);
        assert(0.0 <= gauss_value && gauss_value <= 1.0);
        switch(settings.outlier_removal) {
            case OUTLIER_REMOVAL_NONE:
                return true;
            case OUTLIER_REMOVAL_GAUSS_DAMPING:
                info.quality *= gauss_value; break;
            case OUTLIER_REMOVAL_GAUSS_CLAMPING:
                if (gauss_value < gauss_rejection_threshold) info.quality = 0.0f; break;
        }
    }
    return true;
}

/*
    texture_views 是图片
    face_projection_infos 保存view能看到的face(s)的组合
*/
void calculate_face_projection_infos(mve::TriangleMesh::ConstPtr mesh,
    std::vector<TextureView> * texture_views, Settings const & settings,
    FaceProjectionInfos * face_projection_infos) {

    std::vector<unsigned int> const & faces = mesh->get_faces();  // 三角形顶点索引
    std::vector<math::Vec3f> const & vertices = mesh->get_vertices();
    mve::TriangleMesh::NormalList const & face_normals = mesh->get_face_normals();

    std::size_t const num_views = texture_views->size();

    util::WallTimer timer;
    std::cout << "\tBuilding BVH from " << faces.size() / 3 << " faces... " << std::flush;
    BVHTree bvh_tree(faces, vertices);
    std::cout << "done. (Took: " << timer.get_elapsed() << " ms)" << std::endl;

    ProgressCounter view_counter("\tCalculating face qualities", num_views);
    #pragma omp parallel
    {
        std::vector<std::pair<std::size_t, FaceProjectionInfo> > projected_face_view_infos;

        #pragma omp for schedule(dynamic)
        for (std::uint16_t j = 0; j < static_cast<std::uint16_t>(num_views); ++j) {  // 遍历每个view static_cast 类型转换
            view_counter.progress<SIMPLE>();

            TextureView * texture_view = &texture_views->at(j);  // view类
            texture_view->load_image();
            texture_view->generate_validity_mask();

            if (settings.data_term == DATA_TERM_GMI) {
                texture_view->generate_gradient_magnitude();  // 每个view计算梯度图
                texture_view->erode_validity_mask();
            }

            math::Vec3f const & view_pos = texture_view->get_pos();
            math::Vec3f const & viewing_direction = texture_view->get_viewing_direction();

            for (std::size_t i = 0; i < faces.size(); i += 3) {  // 遍历每个三角形
                std::size_t face_id = i / 3;

                math::Vec3f const & v1 = vertices[faces[i]];
                math::Vec3f const & v2 = vertices[faces[i + 1]];
                math::Vec3f const & v3 = vertices[faces[i + 2]];
                math::Vec3f const & face_normal = face_normals[face_id];
                math::Vec3f const face_center = (v1 + v2 + v3) / 3.0f;

                /* Check visibility and compute quality */

                math::Vec3f view_to_face_vec = (face_center - view_pos).normalized();
                math::Vec3f face_to_view_vec = (view_pos - face_center).normalized();

                /* Backface and basic frustum culling */
                float viewing_angle = face_to_view_vec.dot(face_normal);
                if (viewing_angle < 0.0f || viewing_direction.dot(view_to_face_vec) < 0.0f)  // 面背对视角 或者 视角背对面 不管这个面
                    continue;

                if (std::acos(viewing_angle) > MATH_DEG2RAD(75.0f))  // ?? face在view的正中间也可以角度大啊
                    continue;

                /* Projects into the valid part of the TextureView? */
                if (!texture_view->inside(v1, v2, v3))  // ??将点投到view上看是否有效 mask?
                    continue;

                // frustum culling  ??
                if (settings.geometric_visibility_test) {
                    /* Viewing rays do not collide? */
                    bool visible = true;
                    math::Vec3f const * samples[] = {&v1, &v2, &v3};
                    // TODO: random monte carlo samples...

                    for (std::size_t k = 0; k < sizeof(samples) / sizeof(samples[0]); ++k) {
                        BVHTree::Ray ray;
                        ray.origin = *samples[k];
                        ray.dir = view_pos - ray.origin;
                        ray.tmax = ray.dir.norm();
                        ray.tmin = ray.tmax * 0.0001f;
                        ray.dir.normalize();

                        BVHTree::Hit hit;
                        if (bvh_tree.intersect(ray, &hit)) {
                            visible = false;
                            break;
                        }
                    }
                    if (!visible) continue;
                }

                FaceProjectionInfo info = {j, 0.0f, math::Vec3f(0.0f, 0.0f, 0.0f)};  // {哪一个view, ?, ?}

                /* Calculate quality. */
                texture_view->get_face_info(v1, v2, v3, &info, settings);

                if (info.quality == 0.0) continue;

                /* Change color space. */
                mve::image::color_rgb_to_ycbcr(*(info.mean_color));

                std::pair<std::size_t, FaceProjectionInfo> pair(face_id, info);  // view和face的信息保存到一起, label?
                projected_face_view_infos.push_back(pair);
            }

            texture_view->release_image();
            texture_view->release_validity_mask();
            if (settings.data_term == DATA_TERM_GMI) {
                texture_view->release_gradient_magnitude();
            }
            view_counter.inc();
        }

        //std::sort(projected_face_view_infos.begin(), projected_face_view_infos.end());

        #pragma omp critical
        {
            for (std::size_t i = projected_face_view_infos.size(); 0 < i; --i) {
                std::size_t face_id = projected_face_view_infos[i - 1].first;  // 面
                FaceProjectionInfo const & info = projected_face_view_infos[i - 1].second;  // 
                face_projection_infos->at(face_id).push_back(info);
            }
            projected_face_view_infos.clear();
        }
    }
}

void postprocess_face_infos(Settings const & settings, FaceProjectionInfos * face_projection_infos,
        DataCosts * data_costs) {

    ProgressCounter face_counter("\tPostprocessing face infos", face_projection_infos->size());

    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < face_projection_infos->size(); ++i) {  // 遍历每个面
        face_counter.progress<SIMPLE>();

        std::vector<FaceProjectionInfo> & infos = face_projection_infos->at(i);
        if (settings.outlier_removal != OUTLIER_REMOVAL_NONE) {
            photometric_outlier_detection(&infos, settings);

            infos.erase(std::remove_if(infos.begin(), infos.end(),
                [](FaceProjectionInfo const & info) -> bool {return info.quality == 0.0f;}),
                infos.end());  // 将infos(当前面)里quality是0的视角移除
        }
        std::sort(infos.begin(), infos.end());

        face_counter.inc();
    }

    /* Determine the function for the normlization. */
    float max_quality = 0.0f;  // 所有面的所有view中quality最大的
    for (std::size_t i = 0; i < face_projection_infos->size(); ++i) {  // 遍历每个面的信息(把outlier给剔除了)
        for (FaceProjectionInfo const & info : face_projection_infos->at(i))
            max_quality = std::max(max_quality, info.quality);
    }

    Histogram hist_qualities(0.0f, max_quality, 10000);  // 创建一个10000个bin的直方图, 最小值是0, 最大值是max_quality
    for (std::size_t i = 0; i < face_projection_infos->size(); ++i) {  // 将所有面的所有view的quality加到直方图中, 从小到大排
        for (FaceProjectionInfo const & info : face_projection_infos->at(i))
            hist_qualities.add_value(info.quality);
    }

    float percentile = hist_qualities.get_approx_percentile(0.995f);

    /* Calculate the costs. */
    for (std::uint32_t i = 0; i < face_projection_infos->size(); ++i) {
        for (FaceProjectionInfo const & info : face_projection_infos->at(i)) {

            /* Clamp to percentile and normalize. */
            float normalized_quality = std::min(1.0f, info.quality / percentile);  // quality过高的都限制在1及以内
            float data_cost = (1.0f - normalized_quality);  // quality越低, cost越高
            data_costs->set_value(i, info.view_id, data_cost);
        }

        /* Ensure that all memory is freeed. */
        face_projection_infos->at(i) = std::vector<FaceProjectionInfo>();
    }

    std::cout << "\tMaximum quality of a face within an image: " << max_quality << std::endl;
    std::cout << "\tClamping qualities to " << percentile << " within normalization." << std::endl;
}

void calculate_data_costs(mve::TriangleMesh::ConstPtr mesh, std::vector<TextureView> * texture_views,
    Settings const & settings, DataCosts * data_costs) {

    std::size_t const num_faces = mesh->get_faces().size() / 3;  // 三角形面片的个数
    std::size_t const num_views = texture_views->size();  // 视角的个数

    if (num_faces > std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("Exeeded maximal number of faces");
    if (num_views > std::numeric_limits<std::uint16_t>::max())
        throw std::runtime_error("Exeeded maximal number of views");

    FaceProjectionInfos face_projection_infos(num_faces);  // 默认值初始化num_faces个vector 一个面一个vector, 这个vector里再存能看到这个面的所有views, 以及这个面在当前view的颜色质量信息
    calculate_face_projection_infos(mesh, texture_views, settings, &face_projection_infos);
    postprocess_face_infos(settings, &face_projection_infos, data_costs);
}

TEX_NAMESPACE_END
