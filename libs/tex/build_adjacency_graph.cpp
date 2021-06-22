#include "texturing.h"

#include "progress_counter.h"

TEX_NAMESPACE_BEGIN

/*
* 有共用的边就算相邻
* */
void build_adjacency_graph(mve::TriangleMesh::ConstPtr mesh, mve::MeshInfo const & mesh_info, Graph * graph){
    mve::TriangleMesh::FaceList const & faces = mesh->get_faces();
    std::size_t const num_faces = faces.size() / 3;

    ProgressCounter face_counter("\t 添加边到邻接图: ", num_faces);
    for(std::size_t i=0; i<faces.size(); i+=3){
        face_counter.progress<SIMPLE>();

        // 当前三角形的三个顶点
        std::size_t v1 = faces[i];
        std::size_t v2 = faces[i+1];
        std::size_t v3 = faces[i+2];

        std::vector<std::size_t> adj_faces;  // 当前三角形的相邻面
        mesh_info.get_faces_for_edge(v1, v2, &adj_faces);  // 使用v1, v2作为边的面加入到adj_faces中
        mesh_info.get_faces_for_edge(v2, v3, &adj_faces);  // 使用v2, v3作为边的面加入到adj_faces中
        mesh_info.get_faces_for_edge(v3, v1, &adj_faces);  // 使用v3, v1作为边的面加入到adj_faces中

        // 去除adj_faces中重复的面以及当前面
        for(std::size_t j=0; j<adj_faces.size(); ++j){
            std::size_t face_id = i/3;
            std::size_t adj_face_id = adj_faces[j];

            if(face_id != adj_face_id){  // 相邻不包括自己
                if(!graph->has_edge(face_id, adj_face_id)){
                    graph->add_edge(face_id, adj_face_id);
                }
            }
        }
        face_counter.inc();
    }
    std::cout<<"\t 共"<<graph->num_edges()<<"条边"<<std::endl;
}

TEX_NAMESPACE_END
