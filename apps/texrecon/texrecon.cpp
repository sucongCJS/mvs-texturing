#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <tbb/task_scheduler_init.h>
#include <omp.h>

#include <util/system.h>
#include <util/file_system.h>
#include <mve/mesh_io_ply.h>

#include "tex/util.h"
#include "tex/texturing.h"
#include "tex/debug.h"
#include "tex/progress_counter.h"

#include "arguments.h"

int main(int argc, char **argv){
    util::system::print_build_timestamp(argv[0]);  // 打印程序名
    // util::system::register_segfault_handler();


    Arguments conf;  // 输入的一系列配置, 参数
    try{
        conf = parse_args(argc, argv);  // 解析参数
    }
    catch(std::invalid_argument & ia){
        std::cerr<<ia.what()<<std::endl;
        std::exit(EXIT_FAILURE);
        
    }

    std::string const out_dir = util::fs::dirname(conf.out_prefix);  // 获取文件路径
    if (!util::fs::dir_exists(out_dir.c_str())) {
        std::cerr << "目标目录不存在!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string const tmp_dir = util::fs::join_path(out_dir, "tmp");
    if (!util::fs::dir_exists(tmp_dir.c_str())) {
        util::fs::mkdir(tmp_dir.c_str());
    }
    else {
        std::cerr
            << "临时目录 \"tmp\" 在目标目录中已存在\n"
            << "无法继续, 需要先将 \"tmp\" 目录删除。\n"
            << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // 设置运行时使用的线程数
    tbb::task_scheduler_init schedule(conf.num_threads > 0 ? conf.num_threads : tbb::task_scheduler_init::automatic);
    if (conf.num_threads > 0) {
        omp_set_dynamic(0);
        omp_set_num_threads(conf.num_threads);
    }

    std::cout<<"加载和准备网格: "<<std::endl;
    mve::TriangleMesh::Ptr mesh;  // 包含顶点等信息
    try{
        mesh = mve::geom::load_ply_mesh(conf.in_mesh);
    }
    catch(std::exception & e){
        std::cerr<<"\t无法加载网格: "<<e.what()<<std::endl;
        remove_tmp(tmp_dir);
        std::exit(EXIT_FAILURE);
    }
    mve::MeshInfo mesh_info(mesh);
    tex::prepare_mesh(&mesh_info, mesh);  // 准备网格以进行纹理处理

    std::cout<<"生成纹理视角: "<<std::endl;
    tex::TextureViews texture_views;
    tex::generate_texture_views(conf.in_scene, &texture_views, tmp_dir);

    std::size_t const num_faces = mesh->get_faces().size() / 3;

    std::cout<<"构建邻接图: "<<std::endl;
    tex::Graph graph(num_faces);  // 邻接图graph
    tex::build_adjacency_graph(mesh, mesh_info, &graph);

    std::cout<<"视角选择: "<<std::endl;
    if (conf.labeling_file.empty()) {
        
        tex::DataCosts data_costs(num_faces, texture_views.size());  // 列为每个面, 行为每个视角, 值为cost
        tex::calculate_data_costs(mesh, &texture_views, conf.settings, &data_costs);
        
        try {
            tex::view_selection(data_costs, &graph, conf.settings);
        }
        catch (std::runtime_error& e) {
            std::cerr << "\t优化失败: " << e.what() << std::endl;
            remove_tmp(tmp_dir);
            std::exit(EXIT_FAILURE);
        }

        /* Write labeling to file. */
        if (conf.write_intermediate_results) {
            std::cout << "\t保存标签信息... " << std::flush;
            std::vector<std::size_t> labeling(graph.num_nodes());
            for (std::size_t i = 0; i < graph.num_nodes(); ++i) {
                labeling[i] = graph.get_label(i);
            }
            vector_to_file(conf.out_prefix + "_labeling.vec", labeling);
            std::cout << "完成." << std::endl;
        }
    }
    else {
        std::cout << "\t从指定文件中加载标签... " << std::flush;

        /* Load labeling from file. */
        std::vector<std::size_t> labeling = vector_from_file<std::size_t>(conf.labeling_file, tmp_dir);
        if (labeling.size() != graph.num_nodes()) {
            std::cerr << "Wrong labeling file for this mesh/scene combination... aborting!" << std::endl;
            remove_tmp(tmp_dir);
            std::exit(EXIT_FAILURE);
        }

        /* Transfer labeling to graph. */
        for (std::size_t i = 0; i < labeling.size(); ++i) {
            const std::size_t label = labeling[i];
            if (label > texture_views.size()){
                std::cerr << "Wrong labeling file for this mesh/scene combination... aborting!" << std::endl;
                remove_tmp(tmp_dir);
                std::exit(EXIT_FAILURE);
            }
            graph.set_label(i, label);
        }
        std::cout << "完成." << std::endl;
    }

    tex::TextureAtlases texture_atlases;
    {
        // 创建纹理面片并进行调整
        tex::TexturePatches texture_patches;
        tex::VertexProjectionInfos vertex_projection_infos;  // 顶点的信息, 包括是哪个面的, 这个面用的纹理块, 这个点投影到纹理块上的二维坐标
        std::cout<<"生成纹理块: "<<std::endl;
        tex::generate_texture_patches(graph, mesh, mesh_info, &texture_views, conf.settings, &vertex_projection_infos, &texture_patches);

        // if(conf.settings.global_seam_leveling){
        //     std::cout<<"全局颜色调整: "<<std::endl;
        //     tex::global_seam_leveling(graph, mesh, mesh_info, vertex_projection_infos, &texture_patches);
        // }

        if(conf.settings.local_seam_leveling){
            std::cout<<"局部颜色调整: "<<std::endl;
            tex::local_seam_leveling(graph, mesh, mesh_info, vertex_projection_infos, &texture_patches);
        }

        std::cout<<"生成纹理集: "<<std::endl;
        tex::generate_texture_atlases(&texture_patches, conf.settings, &texture_atlases);
    }
    
    /* Create and write out obj model. */
    {
        std::cout << "创建 obj模型:" << std::endl;
        tex::Model model;
        tex::build_model(mesh, texture_atlases, &model);

        std::cout << "\t保存模型... " << std::flush;
        tex::Model::save(model, conf.out_prefix);
        std::cout << "完成." << std::endl;
    }

    if(conf.write_view_selection_model){
        texture_atlases.clear();
        std::cout << "生成用于debug的纹理块: " << std::endl;
        {
            tex::TexturePatches texture_patches;
            tex::generate_debug_embeddings(&texture_views);
            std::cout<<texture_views.size()<<std::endl;
            tex::VertexProjectionInfos vertex_projection_infos;
            tex::generate_texture_patches(graph, mesh, mesh_info, &texture_views, conf.settings, &vertex_projection_infos, &texture_patches);
            tex::generate_texture_atlases(&texture_patches, conf.settings, &texture_atlases);
            std::cout<<texture_atlases.size()<<std::endl;
        }

        std::cout << "创建debug的.obj模型: " << std::endl;
        {
            tex::Model model;
            tex::build_model(mesh, texture_atlases, &model);

            std::cout << "\t保存模型... " << std::flush;
            tex::Model::save(model, conf.out_prefix + "_view_selection");
            std::cout << "完成." << std::endl;
        }
    }

    /* 移除临时文件(夹) */
    // remove_tmp(tmp_dir);
    remove_tmp(tmp_dir);

    return EXIT_SUCCESS;
}