#include <map>
#include <set>

#include <util/timer.h>
#include <math/accum.h>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include "texturing.h"
#include "seam_leveling.h"
#include "progress_counter.h"

TEX_NAMESPACE_BEGIN

/*
* 寻找以vertex为边的所有边界(一般有两条), 边界的左右两边的图片分别是label1, label2
* */
void find_seam_edges_for_vertex_label_combination(UniGraph const & graph, mve::TriangleMesh::ConstPtr & mesh,
    mve::MeshInfo const & mesh_info, std::size_t vertex, std::size_t label1, std::size_t label2,
    std::vector<MeshEdge> * seam_edges){
    
    assert(label1 != 0 && label2 != 0 && label1 < label2);

    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();

    std::vector<std::size_t> const & adj_verts = mesh_info[vertex].verts;  // 顶点相邻的顶点
    for (std::size_t i = 0; i < adj_verts.size(); ++i){
        std::size_t adj_vertex = adj_verts[i];
        if (vertex == adj_vertex) continue;
        
        std::vector<std::size_t> edge_faces;
        mesh_info.get_faces_for_edge(vertex, adj_vertex, &edge_faces);  // 返回两个顶点相邻的面片

        for (std::size_t j = 0; j < edge_faces.size(); ++j) {
            for(std::size_t k = j + 1; k < edge_faces.size(); ++k) {
                std::size_t face_label1 = graph.get_label(edge_faces[j]);
                std::size_t face_label2 = graph.get_label(edge_faces[k]);
                if (!(face_label1 < face_label2)) std::swap(face_label1, face_label2);
                
                if (face_label1 != label1 || face_label2 != label2) continue;

                math::Vec3f v1 = vertices[vertex];
                math::Vec3f v2 = vertices[adj_vertex];

                if ((v2 - v1).norm() == 0.0f) continue;  // 边长为0, 退出

                MeshEdge seam_edge = {vertex, adj_vertex};
                seam_edges->push_back(seam_edge);
            }
        }
    }
}



math::Vec3f calculate_difference(VertexProjectionInfos const & vertex_projection_infos,
    mve::TriangleMesh::ConstPtr & mesh, std::vector<TexturePatch::Ptr> const &  texture_patches,
    std::vector<MeshEdge> const & seam_edges, int label1, int label2) {

    assert(label1 != 0 && label2 != 0 && label1 < label2);
    assert(!seam_edges.empty());
    
    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();

    math::Accum<math::Vec3f> color1_accum(math::Vec3f(0.0f));  // 一个累加器(支持基本整型的数据相加)
    math::Accum<math::Vec3f> color2_accum(math::Vec3f(0.0f));

    for (MeshEdge const & seam_edge : seam_edges) {
        math::Vec3f v1 = vertices[seam_edge.v1];
        math::Vec3f v2 = vertices[seam_edge.v2];
        
        if ((v2 - v1).norm() == 0.0f) continue;  // 边长为0, 退出

        std::vector<EdgeProjectionInfo> edge_projection_infos;
        // find_mesh_edge_projections(vertex_projection_infos, seam_edge, &edge_projection_infos);

        
    }
}

void global_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    mve::MeshInfo const & mesh_info,
    std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    std::vector<TexturePatch::Ptr> * texture_patches) {
    
    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();
    std::size_t const num_vertices = vertices.size();

    std::cout << "\t 创建用于优化的矩阵... " << std::flush;

    std::vector<std::map<std::size_t, std::size_t> > vertlabel2row;  // 把顶点和view对应起来, 一个顶点可能被多个面共用, 所以一个顶点可能对多个view, 相同点不同的view,点的序号也不同  vertlabel2row[顶点_id(即在mesh中所有顶点的索引)][view_id] = 点的新序号(即x_row)
    vertlabel2row.resize(num_vertices);

    std::vector<std::vector<std::size_t> > labels;  // 索引和顶点在mesh中的索引一样 {{view的id}}
    labels.resize(num_vertices);

    // 为每个标签的每个顶点 指定一个新索引（row）
    std::size_t x_row = 0;  // 顶点的序号
    for (std::size_t i = 0; i < num_vertices; ++i) {  // 遍历每个顶点
        std::set<std::size_t> label_set;  // 集合 如果使用当前顶点的多个面来自同一个view, 只保留一个
    
        std::vector<std::size_t> faces = mesh_info[i].faces;  // 相邻面
        std::set<std::size_t>::iterator it = label_set.begin();
        for (std::size_t j = 0; j < faces.size(); ++j) {  // 遍历每个相邻面, 找出这个顶点相邻的面片用到的view有哪些
            std::size_t label = graph.get_label(faces[j]);  // 面用的是哪个view
            label_set.insert(it, label);
        }

        for (it = label_set.begin(); it != label_set.end(); ++it) {  // 遍历顶点相邻的面用到的view
            std::size_t label = *it;
            if (label == 0) continue;

            vertlabel2row[i][label] = x_row;  // 
            labels[i].push_back(label);
            ++x_row;
        }
    }

    std::size_t x_rows = x_row;
    assert(x_rows < static_cast<std::size_t>(std::numeric_limits<int>::max()));

    float const lambda = 0.1f;

    /* Fill the Tikhonov matrix Gamma(regularization constraints). */
    std::size_t Gamma_row = 0;  // 能表示有多少对在同一个patch又相邻的顶点
    std::vector<Eigen::Triplet<float, int> > coefficients_Gamma;  // A small structure to hold a non zero as a triplet (i,j,value).  用来选中在同一个纹理块的顶点, 给他们的g加个正则项防止过大  [从0开始的索引][顶点的序号]:+/-λ
    coefficients_Gamma.reserve(2 * num_vertices);
    for (std::size_t i = 0; i < num_vertices; ++i) {  // 遍历每个顶点
        for (std::size_t j = 0; j < labels[i].size(); ++j) {  // 遍历每个顶点连接的view
            std::vector<std::size_t> const & adj_verts = mesh_info[i].verts;  // 相邻顶点
            for (std::size_t k = 0; k < adj_verts.size(); ++k) {  // 遍历相邻顶点
                std::size_t adj_vertex = adj_verts[k];
                for (std::size_t l = 0; l < labels[adj_vertex].size(); ++l) {  // 遍历相邻顶点连接的view
                    std::size_t label = labels[i][j];  // 当前顶点的一个view
                    std::size_t adj_vertex_label = labels[adj_vertex][l];  // 相邻顶点的一个view
                    if (i < adj_vertex && label == adj_vertex_label) {  // 两个顶点用的是同一个view
                        Eigen::Triplet<float, int> t1(Gamma_row, vertlabel2row[i][label], lambda);
                        Eigen::Triplet<float, int> t2(Gamma_row, vertlabel2row[adj_vertex][adj_vertex_label], -lambda);
                        coefficients_Gamma.push_back(t1);
                        coefficients_Gamma.push_back(t2);
                        Gamma_row++;
                    }
                }
            }
        }
    }

    std::size_t Gamma_rows = Gamma_row;
    assert(Gamma_rows < static_cast<std::size_t>(std::numeric_limits<int>::max()));

    Eigen::SparseMatrix<float> Gamma(Gamma_rows, x_rows);  // (相邻顶点且相同view的个数, 顶点的view的组合个数) 行为同一个patch又相邻的顶点对的对数, 列为顶点(如果是不同view的顶点, 要区别对待)
    Gamma.setFromTriplets(coefficients_Gamma.begin(), coefficients_Gamma.end());

    /* Fill the matrix A and the coefficients for the Vector b of the linear equation system. */
    std::vector<Eigen::Triplet<float, int> > coefficients_A;  // [从0开始的索引][view的序号]:+/-1
    std::vector<math::Vec3f> coefficients_b;
    std::size_t A_row = 0;
    for (std::size_t i = 0; i < num_vertices; ++i) {
        for (std::size_t j = 0; j < labels[i].size(); ++j) {  // 遍历顶点连接的view
            for (std::size_t k = 0; k < labels[i].size(); ++k) {
                std::size_t label1 = labels[i][j];
                std::size_t label2 = labels[i][k];
                if (label1 < label2) {  // 不同的patch才有seam

                    std::vector<MeshEdge> seam_edges;  // 纹理的边缘
                    find_seam_edges_for_vertex_label_combination(graph, mesh, mesh_info, i, label1, label2, &seam_edges);

                    if (seam_edges.empty()) continue;

                    Eigen::Triplet<float, int> t1(A_row, vertlabel2row[i][label1], 1.0f);
                    Eigen::Triplet<float, int> t2(A_row, vertlabel2row[i][label2], -1.0f);
                    coefficients_A.push_back(t1);
                    coefficients_A.push_back(t2);

                    coefficients_b.push_back(calculate_difference(vertex_projection_infos, mesh, *texture_patches, seam_edges, label1, label2));

                    ++A_row;
                }
            }
        }
    }

    std::size_t A_rows = A_row;
    assert(A_rows < static_cast<std::size_t>(std::numeric_limits<int>::max()));

    Eigen::SparseMatrix<float> A(A_rows, x_rows);  // (接缝的数量, 顶点的view的组合个数), 每一行中有一个-1, 1, 且这两个元素位于接缝处
    A.setFromTriplets(coefficients_A.begin(), coefficients_A.end());

    Eigen::SparseMatrix<float> Lhs = A.transpose() * A + Gamma.transpose() * Gamma;
    /* Only keep lower triangle (CG only uses the lower), prune the rest and compress matrix. */
    Lhs.prune([](const int& row, const int& col, const float& value) -> bool {  // left hand side
            return col <= row && value != 0.0f;  // 下三角
        }); // value != 0.0f is only to suppress a compiler warning

    std::vector<std::map<std::size_t, math::Vec3f> > adjust_values(num_vertices);  // 每个顶点的调整值1
    std::cout << "完成." << std::endl;
    std::cout << "\tLhs的维度: " << Lhs.rows() << " x " << Lhs.cols() << std::endl;

    util::WallTimer timer;
    std::cout << "\t 计算颜色调整:"<< std::endl;
    #pragma omp parallel for
    for (std::size_t channel = 0; channel < 3; ++channel) {
        /* Prepare solver. */
        Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower> cg;
        cg.setMaxIterations(1000);
        cg.setTolerance(0.0001);   // The tolerance corresponds to the relative residual error: |Ax-b|/|b|
        cg.compute(Lhs);

        /* Prepare right hand side. */
        Eigen::VectorXf b(A_rows);
        for (std::size_t i = 0; i < coefficients_b.size(); ++i) {
            b[i] = coefficients_b[i][channel];
        }
        Eigen::VectorXf Rhs = Eigen::SparseMatrix<float>(A.transpose()) * b;

        /* Solve for x. */
        Eigen::VectorXf x(x_rows);
        x = cg.solve(Rhs);

        /* 减去平均值，因为系统是欠约束的，我们寻求最小调整的解。 */
        x = x.array() - x.mean();

        #pragma omp critical
        std::cout << "\t\t 颜色通道" << channel << ": CG 花了 "
            << cg.iterations() << "次迭代. 残差是" << cg.error() << std::endl;

        #pragma omp critical
        for (std::size_t i = 0; i < num_vertices; ++i) {
            for (std::size_t j = 0; j < labels[i].size(); ++j) {
                std::size_t label = labels[i][j];
                adjust_values[i][label][channel] = x[vertlabel2row[i][label]];
            }
        }
    }
    // std::cout << "\t\t 花费了" << timer.get_elapsed_sec() << "秒" << std::endl;

    // mve::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();

    // ProgressCounter texture_patch_counter("\t 调整纹理块", texture_patches->size());
    // #pragma omp parallel for schedule(dynamic)
    // for (std::size_t i = 0; i < texture_patches->size(); ++i) {
    //     texture_patch_counter.progress<SIMPLE>();

    //     TexturePatch::Ptr texture_patch = texture_patches->at(i);

    //     int label = texture_patch->get_label();
    //     std::vector<std::size_t> const & faces = texture_patch->get_faces();
    //     std::vector<math::Vec3f> patch_adjust_values(faces.size() * 3, math::Vec3f(0.0f));  // 三通道

    //     /* Only adjust texture_patches originating form input images. */
    //     if (label == 0) {
    //         texture_patch->adjust_colors(patch_adjust_values);
    //         texture_patch_counter.inc();
    //         continue;
    //     };

    //     for (std::size_t j = 0; j < faces.size(); ++j) {
    //         for (std::size_t k = 0; k < 3; ++k) {
    //             std::size_t face_pos = faces[j] * 3 + k;
    //             std::size_t vertex = mesh_faces[face_pos];
    //             patch_adjust_values[j * 3 + k] = adjust_values[vertex].find(label)->second;
    //         }
    //     }

    //     texture_patch->adjust_colors(patch_adjust_values);
    //     texture_patch_counter.inc();
    // }
}

TEX_NAMESPACE_END
