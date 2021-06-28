#include "texturing.h"
#include "seam_leveling.h"
// #include "seam_leveling.cpp"

TEX_NAMESPACE_BEGIN

void local_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    mve::MeshInfo const & mesh_info,
    VertexProjectionInfos const & vertex_projection_infos,
    TexturePatches * texture_patches){

    std::size_t const num_vertices = vertex_projection_infos.size();
    std::vector<math::Vec3f> vertex_colors(num_vertices);
    std::vector<std::vector<math::Vec3f> > edge_colors;

    {
        std::vector<MeshEdge> seam_edges;
        find_seam_edges(graph, mesh, &seam_edges);
    }
}

TEX_NAMESPACE_END