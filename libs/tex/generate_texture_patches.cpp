#include <set>
#include <list>

#include <util/timer.h>
#include <mve/image_tools.h>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include "texturing.h"
#include "rect.h"

TEX_NAMESPACE_BEGIN

/** 
 * 候选纹理块, 最终的纹理块由多个候选纹理块构成
 */
struct TexturePatchCandidate {
    Rect<int> bounding_box;
    TexturePatch::Ptr texture_patch;
};


/** 
 * @param label: 候选纹理块的标签
 * @param texture_view: 使用的视角
 * @param subgraph: 使用了label标签的连通的面片的向量
 * @param mesh: 网格
 */
TexturePatchCandidate generate_candidate(int label, TextureView const & texture_view,
    std::vector<std::size_t> const & subgraph, mve::TriangleMesh::ConstPtr mesh,
    Settings const & settings){
    mve::ByteImage::Ptr view_image = texture_view.get_image();
    int min_x = view_image->width(), min_y = view_image->height();
    int max_x = 0, max_y = 0;
    
    mve::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();
    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();

    std::vector<math::Vec2f> texcoords;  // 所有面片的三个顶点投影到视角的二维图片上的二维坐标
    for (std::size_t i = 0; i < subgraph.size(); ++i){  // 遍历用了相同标签的连通的面片
        for (std::size_t j = 0; j < 3; ++j){
            math::Vec3f vertex = vertices[mesh_faces[subgraph[i] * 3 + j]];
            math::Vec2f pixel = texture_view.get_pixel_coords(vertex);

            texcoords.push_back(pixel);  // 面片的顶点投影到视角的二维图片上的二维坐标

            min_x = std::min(static_cast<int>(std::floor(pixel[0])), min_x);
            min_y = std::min(static_cast<int>(std::floor(pixel[1])), min_y);
            max_x = std::max(static_cast<int>(std::ceil(pixel[0])), max_x);
            max_y = std::max(static_cast<int>(std::ceil(pixel[1])), max_y);
        }
    }

    // 检查是否超过视角的成像图片的宽高
    assert(min_x >= 0);
    assert(min_y >= 0);
    assert(max_x < view_image->width());
    assert(max_y < view_image->height());

    // 
    int width  = max_x - min_x + 1;  width += 2 * texture_patch_border;
    int height = max_y - min_y + 1; height += 2 * texture_patch_border;
    min_x -= texture_patch_border;
    min_y -= texture_patch_border;

    // 将所有顶点 在视角的图片的位置 转为 在aabb中的相对位置, aabb被视角的图片所包含
    math::Vec2f min(min_x, min_y);
    for (std::size_t i = 0; i < texcoords.size(); ++i) {
        texcoords[i] = texcoords[i] - min;
    }

    mve::ByteImage::Ptr byte_image;
    byte_image = mve::image::crop(view_image, width, height, min_x, min_y, *math::Vec3uc(255, 0, 255));  // 裁剪图片, Region may exceed the input image dimensions, new pixel values are initialized with the given color(255,0,255).
    mve::FloatImage::Ptr image = mve::image::byte_to_float_image(byte_image);  // scaling from [0, 255] to [0, 1]

    // gamma矫正
    if (settings.tone_mapping == TONE_MAPPING_GAMMA) {
        mve::image::gamma_correct(image, 2.2f);
    }

    TexturePatchCandidate texture_patch_candidate = {
        Rect<int>(min_x, min_y, max_x, max_y), TexturePatch::create(label, subgraph, texcoords, image)
    };

    return texture_patch_candidate;
}

/**
 * @brief merge用了相同纹理块的顶点的面片信息, 顶点的投影信息不需要合并
 */
void merge_vertex_projection_infos(std::vector<std::vector<VertexProjectionInfo> > * vertex_projection_infos) {
    #pragma omp parallel for
    for (std::size_t i = 0; i < vertex_projection_infos->size(); ++i) {  // 遍历每个有投影信息的顶点
        std::vector<VertexProjectionInfo> & infos = vertex_projection_infos->at(i);
        
        std::map<std::size_t, VertexProjectionInfo> info_map;  // 纹理块的id:面片的信息
        std::map<std::size_t, VertexProjectionInfo>::iterator it;

        for (VertexProjectionInfo const & info : infos) {  // 遍历一个顶点的所有投影信息
            std::size_t texture_patch_id = info.texture_patch_id;
            if((it = info_map.find(texture_patch_id)) == info_map.end()) {
                info_map[texture_patch_id] = info;
            }
            else{
                it->second.faces.insert(it->second.faces.end(), info.faces.begin(), info.faces.end());  // 合并顶点的面片信息
            }
        }

        infos.clear();
        infos.reserve(info_map.size());
        for(it=info_map.begin(); it!=info_map.end(); ++it){
            infos.push_back(it->second);
        }
    }
}

/**
 * @param texture_patches: 
 */
void generate_texture_patches(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    mve::MeshInfo const & mesh_info,
    std::vector<TextureView> * texture_views, Settings const & settings,
    std::vector<std::vector<VertexProjectionInfo> > * vertex_projection_infos,
    std::vector<TexturePatch::Ptr> * texture_patches) {

    mve::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();
    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();
    vertex_projection_infos->resize(vertices.size());

    std::size_t num_patches = 0;

    std::cout << "\t执行中... " << std::flush;

    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < texture_views->size(); ++i) {  // 遍历每一个视角
        std::vector<std::vector<std::size_t> > subgraphs;  // 保存有相同标签的面片, 每一片连通的面片一个集合

        int const label = i + 1;
        graph.get_subgraphs(label, &subgraphs);

        TextureView * texture_view = &texture_views->at(i);
        texture_view->load_image();
        std::list<TexturePatchCandidate> candidates;  // 每个subgraph为一个候选纹理块, 所有的subgraph构成最终的纹理块

        for (std::size_t j = 0; j < subgraphs.size(); ++j) {  // 遍历每个连通的subgraph
            candidates.push_back(generate_candidate(label, *texture_view, subgraphs[j], mesh, settings));
        }
        texture_view->release_image();

        // 合并包含相同图像的candidate, 也就是这个视角中被用于纹理映射的图片部分
        std::list<TexturePatchCandidate>::iterator it, sit;
        for (it = candidates.begin(); it != candidates.end(); ++it) {
            for (sit = candidates.begin(); sit != candidates.end();) {
                Rect<int> bounding_box = sit->bounding_box;
                if (it != sit && bounding_box.is_inside(&it->bounding_box)) {  // 如果一个candicate的aabb包含另一个的aabb(但是两个candicate肯定是不连通的), 合并
                    // 合并面片
                    TexturePatch::Faces & faces = it->texture_patch->get_faces();
                    TexturePatch::Faces & ofaces = sit->texture_patch->get_faces();
                    faces.insert(faces.end(), ofaces.begin(), ofaces.end());

                    // 合并顶点在二维图片上的相对坐标
                    TexturePatch::Texcoords & texcoords = it->texture_patch->get_texcoords();
                    TexturePatch::Texcoords & otexcoords = sit->texture_patch->get_texcoords();
                    math::Vec2f offset;
                    offset[0] = sit->bounding_box.min_x - it->bounding_box.min_x;
                    offset[1] = sit->bounding_box.min_y - it->bounding_box.min_y;
                    for (std::size_t i = 0; i < otexcoords.size(); ++i) {
                        texcoords.push_back(otexcoords[i] + offset);
                    }

                    sit = candidates.erase(sit);
                }
                else{
                    ++sit;
                }
            }
        }

        for (it = candidates.begin(); it != candidates.end(); ++it) {
            std::size_t texture_patch_id;

            #pragma omp critical
            {
                texture_patches->push_back(it->texture_patch);
                texture_patch_id = num_patches++;
            }

            std::vector<std::size_t> const & faces = it->texture_patch->get_faces();
            std::vector<math::Vec2f> const & texcoords = it->texture_patch->get_texcoords();
            for (std::size_t i = 0; i < faces.size(); ++i) {  // 遍历patch里的每个面
                std::size_t const face_id = faces[i];
                std::size_t const face_pos = face_id * 3;
                for (std::size_t j = 0; j < 3; ++j) {  // 遍历每个顶点
                    std::size_t const vertex_id = mesh_faces[face_pos  + j];
                    math::Vec2f const projection = texcoords[i * 3 + j];

                    VertexProjectionInfo info = {texture_patch_id, projection, {face_id}};

                    #pragma omp critical
                    vertex_projection_infos->at(vertex_id).push_back(info);  // 一个顶点会被多个面用到
                }
            }
        }
    }

    merge_vertex_projection_infos(vertex_projection_infos);

    // hole_filling

    std::cout<<" 生成"<<num_patches<<"块纹理块."<<std::endl;
}

TEX_NAMESPACE_END
