#include "seam_leveling.h"

TEX_NAMESPACE_BEGIN

/**
  * 找缝隙
  * @return seam_edges 
  */
void find_seam_edges(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh, std::vector<MeshEdge> * seam_edges){
    mve::TriangleMesh::FaceList const & faces = mesh->get_faces();

    for(std::size_t node=0; node<graph.num_nodes(); ++node){
        std::vector<std::size_t> adj_nodes = graph.get_adj_nodes(node);
        for(std::size_t adj_nodes : adj_nodes){
            if(node > adj_nodes) continue;  // 防止重复添加啊

            int label1 = graph.get_label(node);
            int label2 = graph.get_label(adj_nodes);

            if(label1 == label2) continue;

            std::vector<std::size_t> shared_edge;
            for(int i=0; i<3; ++i){
                std::size_t v1 = faces[node*3 + i];
                for(int j=0; j<3; ++j){
                    std::size_t v2 = faces[adj_nodes*3 + j];

                    if(v1 == v2) shared_edge.push_back(v1);
                }
            }

            assert(shared_edge.size() == 2);
            std::size_t v1 = shared_edge[0];
            std::size_t v2 = shared_edge[1];

            assert(v1 != v2);
            if(v1 > v2) std::swap(v1, v2);

            seam_edges->push_back({v1, v2});
        }
    }
}

TEX_NAMESPACE_END
