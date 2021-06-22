#ifndef TEX_SEAMLEVELING_HEADER
#define TEX_SEAMLEVELING_HEADER

#include <vector>

#include <mve/mesh.h>

TEX_NAMESPACE_BEGIN

struct VertexProjectionInfo {
    std::size_t texture_patch_id;  // 纹理块id
    math::Vec2f projection;  // 顶点投影到texture_patch_id纹理块中的二维坐标
    std::vector<std::size_t> faces;  // 包含这个顶点的面片

    bool operator<(VertexProjectionInfo const & other) const {
        return texture_patch_id < other.texture_patch_id;
    }
};

struct EdgeProjectionInfo {
    std::size_t texture_patch_id;  // 纹理块id
    math::Vec2f p1;  // 顶点
    math::Vec2f p2;

    bool operator<(EdgeProjectionInfo const & other) const {
        return texture_patch_id < other.texture_patch_id;
    }
};

struct MeshEdge {
    std::size_t v1;  // 顶点
    std::size_t v2;
};

/*
* 
* */
// void find_mesh_edge_projections(std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
//     MeshEdge mesh_edge, std::vector<EdgeProjectionInfo> * edge_projection_infos){
//      // 获取边界上两个点的投影信息
//     std::vector<VertexProjectionInfo> const & v1_projection_infos = vertex_projection_infos[mesh_edge.v1];
//     std::vector<VertexProjectionInfo> const & v2_projection_infos = vertex_projection_infos[mesh_edge.v2];

//     std::set<EdgeProjectionInfo> edge_projection_infos_set;

//     for (VertexProjectionInfo v1_projection_info : v1_projection_infos) {
//         for (VertexProjectionInfo v2_projection_info : v2_projection_infos) {
//             if (v1_projection_info.texture_patch_id != v2_projection_info.texture_patch_id) continue;

//             ...
//         }
//     }
    
// }

TEX_NAMESPACE_END

#endif /* TEX_SEAMLEVELING_HEADER */
