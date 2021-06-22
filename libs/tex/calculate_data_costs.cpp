#include <numeric>

#include <mve/image_color.h>
#include <acc/bvh_tree.h>
#include <Eigen/Core>
#include <Eigen/LU>

#include "texturing.h"
#include "sparse_table.h"
#include "progress_counter.h"
#include "util.h"
#include "histogram.h"

typedef acc::BVHTree<unsigned int, math::Vec3f> BVHTree;

TEX_NAMESPACE_BEGIN

/*
* texture_views: 图片
* face_projection_infos: 保存view能看到的face(s)的组合
* */
void calculate_face_projection_infos(mve::TriangleMesh::ConstPtr mesh, TextureViews * texture_views, Settings const & settings, FaceProjectionInfos * face_projection_infos){
    std::vector<unsigned int> const & faces = mesh->get_faces();  // 三角形顶点索引
    std::vector<math::Vec3f> const & vertices = mesh->get_vertices();
    mve::TriangleMesh::NormalList const & face_normals = mesh->get_face_normals();

    std::size_t const num_views = texture_views->size();

    std::cout<<"\t 使用"<<faces.size()<<"个面构建BVH树..."<<std::flush;
    BVHTree bvh_tree(faces, vertices);
    std::cout<<"完成. "<<std::endl;

    ProgressCounter view_counter("\t 计算每个面片的代价", num_views);
    #pragma omp parallel
    {
        std::vector<std::pair<size_t, FaceProjectionInfo>> projected_face_view_infos;

        #pragma omp for schedule(dynamic)
        for(std::uint16_t j=0; j<static_cast<std::uint16_t>(num_views); ++j){  // 遍历每个view
            view_counter.progress<SIMPLE>();

            TextureView * texture_view = &texture_views->at(j);
            texture_view->load_image();
            texture_view->generate_validity_mask();

            if(settings.data_term == DATA_TERM_GMI){
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

                // 检查可见性 和 计算quality
                // 1. 检查可见性
                math::Vec3f view_to_face_vec = (face_center - view_pos).normalized();
                math::Vec3f face_to_view_vec = (view_pos - face_center).normalized();

                // 背面剔除
                float viewing_angle = face_to_view_vec.dot(face_normal);
                if (viewing_angle < 0.0f || viewing_direction.dot(view_to_face_vec) < 0.0f)  // 面背对视角 或者 视角背对面 不管这个面
                    continue;

                if (!texture_view->inside(v1, v2, v3))  // 将点投到view上看是否有效, 以及看mask
                    continue;

                // 视域剔除
                if (settings.geometric_visibility_test) {
                    /* Viewing rays do not collide? */
                    bool visible = true;
                    math::Vec3f const * samples[] = {&v1, &v2, &v3};

                    for (std::size_t k=0; k<sizeof(samples) / sizeof(samples[0]); ++k) {  // 从三角形的三个顶点连线到view_pos看有没有被挡住, 如果其中一个被挡住, 则三角形设为不可见
                        BVHTree::Ray ray;
                        ray.origin = *samples[k];
                        ray.dir = view_pos - ray.origin;
                        ray.tmax = ray.dir.norm();  // BVH最大盒子的最大值 https://zhuanlan.zhihu.com/p/114307697
                        ray.tmin = ray.tmax * 0.0001f;  // BVH最大盒子的最小值
                        ray.dir.normalize();

                        BVHTree::Hit hit;
                        if (bvh_tree.intersect(ray, &hit)) {
                            visible = false;
                            break;
                        }
                    }
                    if (!visible) continue;
                }
                FaceProjectionInfo info = {j, 0.0f, math::Vec3f(0.0f, 0.0f, 0.0f)};  // {view_id, quality, mean_color}

                // 计算quality
                texture_view->get_face_info(v1, v2, v3, &info, settings);

                if(info.quality == 0.0) continue;

                // 转ycbcr
                mve::image::color_rgb_to_ycbcr(*(info.mean_color));  // 

                std::pair<std::size_t, FaceProjectionInfo> pair(face_id, info);  // view和face的信息保存到一起
                projected_face_view_infos.push_back(pair);
            }

            texture_view->release_image();
            texture_view->release_validity_mask();
            if (settings.data_term == DATA_TERM_GMI) {
                texture_view->release_gradient_magnitude();
            }

            view_counter.inc();
        }

        #pragma omp critical
        {
            for (std::size_t i = projected_face_view_infos.size(); 0 < i; --i){
                std::size_t face_id = projected_face_view_infos[i-1].first;  // 面
                FaceProjectionInfo const & info = projected_face_view_infos[i-1].second;  // 视角信息
                face_projection_infos->at(face_id).push_back(info);
            }
            projected_face_view_infos.clear();
        }
    }

}

/*
* 消除面片投影在的所有视角的中质量与其他大多数视角中的颜色有很大不同的面片
* param infos 包含从多个视角看到的一张脸的信息
* param settings 运行时配置。
* return 异常值删除是否成功。
* */
bool photometric_outlier_detection(std::vector<FaceProjectionInfo> * infos, Settings const & settings) {
    if (infos->size() == 0) return true;  // 这个面没有被任何view看到

    // 参数配置
    double const gauss_rejection_threshold = 6e-3;
    double const minimal_covariance = 5e-4;  // 如果所有协方差都低于这个值，就停止异常检测

    int const outlier_detection_iterations = 10;  // 迭代次数
    int const minimal_num_inliers = 4;

    float outlier_removal_factor = std::numeric_limits<float>::signaling_NaN();
    switch(settings.outlier_removal){
        case OUTLIER_REMOVAL_NONE: return true;
        case OUTLIER_REMOVAL_GAUSS_CLAMPING: outlier_removal_factor = 1.0f; break;
        case OUTLIER_REMOVAL_GAUSS_DAMPING:  outlier_removal_factor = 0.2f; break;
    }

    Eigen::MatrixX3d inliers(infos->size(), 3);  // shape=view的个数 * 3, 存这个面在每个view的mean color
    for(std::size_t row=0; row < infos->size(); ++row) {
        inliers.row(row) = mve_to_eigen(infos->at(row).mean_color).cast<double>();
    }

    std::vector<std::uint32_t> is_inlier(infos->size(), 1);  // 把所有的view都先当成inlier, 1为inlier, 0不是

    Eigen::RowVector3d var_mean;
    Eigen::Matrix3d covariance;
    Eigen::Matrix3d covariance_inv;

    for(int i = 0; i < outlier_detection_iterations; ++i) {  // mean shift 迭代一定次数
        if(inliers.rows() < minimal_num_inliers)  // inlier小于一定数目就不剔除了
            return false;

        // 计算inliers(每个view的平均颜色)的平均颜色和协方差
        var_mean = inliers.colwise().mean();  // 计算三个通道的mean值
        Eigen::MatrixX3d centered = inliers.rowwise() - var_mean;  // 和mean作差
        covariance = (centered.adjoint() * centered) / double(inliers.rows() - 1);

        // 如果所有协方差都很小，就停止，保留已有inlier, 并将outlier的质量设置为零
        if(covariance.array().abs().maxCoeff() < minimal_covariance) {  // 给矩阵的每个元素加绝对值 
            for(std::size_t row = 0; row < infos->size(); ++row) {  // 遍历每个view
                if(!is_inlier[row]){
                    infos->at(row).quality = 0.0f;
                }
            }
            return true;
        }

        Eigen::FullPivLU<Eigen::Matrix3d> lu(covariance);
        if (!lu.isInvertible()) return false;  // 如果协方差变得不可导, 退出

        covariance_inv = lu.inverse();

        // 计算新的内联线数（高斯值高于阈值的所有视图）
        for (std::size_t row = 0; row < infos->size(); ++row) {
            Eigen::RowVector3d color = mve_to_eigen(infos->at(row).mean_color).cast<double>();
            double gauss_value = multi_gauss_unnormalized(color, var_mean, covariance_inv);
            is_inlier[row] = (gauss_value >= gauss_rejection_threshold ? 1 : 0);  // 越大越好
        }

        // 根据还有多少linlier调整矩阵
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

void postprocess_face_infos(Settings const & settings, FaceProjectionInfos * face_projection_infos, DataCosts * data_costs){
    ProgressCounter face_counter("\t 处理面片信息", face_projection_infos->size());

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

    // 确定标准化函数
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

    /* 计算 costs. */
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

    std::cout << "\t 图片中质量最大的面片: " << max_quality << std::endl;
    std::cout << "\t 把质量范围控制在" << percentile << "以内" << std::endl;
}

void calculate_data_costs(mve::TriangleMesh::ConstPtr mesh, TextureViews * texture_views, Settings const & settings, DataCosts * data_costs){
    std::size_t const num_faces = mesh->get_faces().size()/3;  // 三角形面片的个数
    std::size_t const num_views = texture_views->size();  // 视角的个数

    if (num_faces > std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("超出最多面片数");
    if (num_views > std::numeric_limits<std::uint16_t>::max())
        throw std::runtime_error("低于最少面片数");

    FaceProjectionInfos face_projection_infos(num_faces);  // 默认值初始化num_faces个vector 一个面片一个vector, 这个vector里再存能看到这个面的所有views, 以及这个面在当前view的颜色质量信息. 如果面片质量为0, 则是空的. 
    calculate_face_projection_infos(mesh, texture_views, settings, &face_projection_infos);  // 计算face_projection_infos

    postprocess_face_infos(settings, &face_projection_infos, data_costs);  // 预先计算data_cost, 后面会用到
}

TEX_NAMESPACE_END
