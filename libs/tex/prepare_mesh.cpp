#include "texturing.h"

TEX_NAMESPACE_BEGIN

std::size_t remove_redudant_faces(mve::MeshInfo const & mesh_info, mve::TriangleMesh::Ptr mesh){
    mve::TriangleMesh::FaceList & faces = mesh->get_faces();
    mve::TriangleMesh::FaceList new_faces;
    new_faces.reserve(faces.size());

    std::size_t num_redundant = 0;
    for(std::size_t i=0; i<faces.size(); i+=3){  // 遍历每个三角形
        std::size_t face_id = i/3;
        bool is_redundant = false;

        for(std::size_t j=0; !is_redundant && j<3; ++j){
            mve::MeshInfo::AdjacentFaces const & adj_faces = mesh_info[faces[i+j]].faces;  // 当前三角形三个顶点相邻的面
            for(std::size_t k=0; !is_redundant && k<adj_faces.size(); ++k){
                std::size_t adj_face_id = adj_faces[k];

                // 只移除有相邻id的面片
                if(face_id < adj_face_id){
                    bool identical = true;
                    // 如果由相同的顶点构成, 那么这两个三角形判为一样
                    for(std::size_t l = 0; l < 3; ++l) {
                        std::size_t vertex = faces[adj_face_id * 3 + l];
                        if (std::find(&faces[i], &faces[i + 3], vertex) == &faces[i + 3]) {
                            identical = false;
                            break;
                        }
                    }

                    is_redundant = identical;
                }
            }
        }

        if(is_redundant){
            ++num_redundant;
        }
        else{
            new_faces.insert(new_faces.end(), faces.cbegin()+i, faces.cbegin()+i+3);
        }
    }

    faces.swap(new_faces);

    return num_redundant;
}

void prepare_mesh(mve::MeshInfo *mesh_info, mve::TriangleMesh::Ptr mesh){
    std::size_t num_redundant = remove_redudant_faces(*mesh_info, mesh);
    if(num_redundant > 0){
        std::cout<<"\t 移除了"<<num_redundant<<"个冗余面"<<std::endl;
    }

    // 确保面和顶点数能够对应起来
    mesh->ensure_normals(true, true);

    // 更新顶点信息
    mesh_info->clear();
    mesh_info->initialize(mesh);
}

TEX_NAMESPACE_END
