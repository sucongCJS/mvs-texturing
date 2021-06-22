#ifndef TEX_TEXTURING_HEADER
#define TEX_TEXTURING_HEADER

#include <vector>
#include <iostream>
#include <cassert>

#include "mve/mesh.h"
#include "mve/mesh_info.h"

#include "defines.h"
#include "texture_view.h"
#include "uni_graph.h"
#include "sparse_table.h"
#include "texture_atlas.h"
#include "texture_patch.h"
#include "seam_leveling.h"
#include "obj_model.h"

TEX_NAMESPACE_BEGIN

typedef std::vector<TextureView> TextureViews;
typedef UniGraph Graph;  // 固定节点数的单向图
typedef ObjModel Model;
typedef SparseTable<std::uint32_t, std::uint16_t, float> DataCosts;  // 每个面, 和对应视角之间代价
typedef std::vector<std::vector<FaceProjectionInfo>> FaceProjectionInfos;  // 
typedef std::vector<std::vector<VertexProjectionInfo>> VertexProjectionInfos;  // 顶点的投影信息, 一个顶点可能有多个投影信息, 
typedef std::vector<TextureAtlas::Ptr> TextureAtlases;
typedef std::vector<TexturePatch::Ptr> TexturePatches;

/*
* 准备网格以进行纹理处理
* - 删除冗余的面
* - 确保法线（面和顶点）
* */
void prepare_mesh(mve::MeshInfo * mesh_info, mve::TriangleMesh::Ptr mesh);

/*
* 从scene目录中生成纹理视角
* */
void generate_texture_views(std::string const & in_scene, TextureViews * texture_views, std::string const & tmp_dir);

/*
* 创建网格面片的邻接图
* */
void build_adjacency_graph(mve::TriangleMesh::ConstPtr mesh, mve::MeshInfo const & mesh_info, Graph * graph);

/*
* 如果面片在视角中可见，则计算面片和视角组合的代价。
* */
void calculate_data_costs(mve::TriangleMesh::ConstPtr mesh, TextureViews * texture_views, Settings const & settings, DataCosts * data_costs);

/**
 * 进行视角选择并将label保存在graph中
 */
void view_selection(DataCosts const & data_costs, UniGraph * graph, Settings const & settings);

/**
  * 使用图形生成纹理块
  */
void generate_texture_patches(UniGraph const & graph,
    mve::TriangleMesh::ConstPtr mesh,
    mve::MeshInfo const & mesh_info,
    TextureViews * texture_views,
    Settings const & settings,
    VertexProjectionInfos * vertex_projection_infos,
    TexturePatches * texture_patches);

/**
  * 全局颜色调整 by Ivanov and Lempitsky
  */
void global_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    mve::MeshInfo const & mesh_info,
    VertexProjectionInfos const & vertex_projection_infos,
    TexturePatches * texture_patches);

void generate_texture_atlases(TexturePatches * texture_patches,
    Settings const & settings, TextureAtlases * texture_atlases);

/**
  * Builds up an model for the mesh by constructing materials and
  * texture atlases form the texture_patches
  */
void build_model(mve::TriangleMesh::ConstPtr mesh,
    TextureAtlases const & texture_atlas, Model * model);

TEX_NAMESPACE_END

#endif /* TEX_TEXTURING_HEADER */