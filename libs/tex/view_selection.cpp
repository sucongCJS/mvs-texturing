#include <util/timer.h>

#include "util.h"
#include "texturing.h"
#include "mapmap/full.h"

TEX_NAMESPACE_BEGIN

/**
 * 进行视角选择并将label保存在graph中
 */
void view_selection(DataCosts const & data_costs, UniGraph * graph, Settings const & settings){
    using uint_t = unsigned int;
    using cost_t = float;

    // 设置并行的参数, 自动获取字节长度啥的
    constexpr uint_t simd_w = mapmap::sys_max_simd_width<cost_t>();

    // 构建graph: 相邻的面片添加边
    mapmap::Graph<cost_t> mgraph(graph->num_nodes());  // 单向
    for (std::size_t i = 0; i < graph->num_nodes(); ++i) {
        if (data_costs.col(i).empty()) continue;  // 面i没有被看到

        std::vector<std::size_t> adj_faces = graph->get_adj_nodes(i);
        for (std::size_t j = 0; j < adj_faces.size(); ++j) {
            std::size_t adj_face = adj_faces[j];
            if (data_costs.col(adj_face).empty()) continue;
            
            if (i < adj_face) {  // 防止变成双向的
                mgraph.add_edge(i, adj_face, 1.0f);  // 权重是1
            }
        }
    }
    mgraph.update_components();

    // 创建 标签集
    mapmap::LabelSet<cost_t, simd_w> label_set(graph->num_nodes(), false);  // 与cost大小无关. "面和这个面的labels"的集合, label为data_cost中的顺序
    for (std::size_t i = 0; i < data_costs.cols(); ++i) {  // 遍历所有的面
        DataCosts::Column const & data_costs_for_node = data_costs.col(i);  //  一列数据, 当前面所有view

        std::vector<mapmap::_iv_st<cost_t, simd_w> > labels;  // 保存能看到面的视角(label)
        if (data_costs_for_node.empty()) {
            labels.push_back(0);
        }
        else{
            labels.resize(data_costs_for_node.size());  // 调整为能看到的view的个数
            for(std::size_t j = 0; j < data_costs_for_node.size(); ++j) {  // 遍历能看到的面
                labels[j] = data_costs_for_node[j].first + 1;  // view的id(0表示没有, 所以都要+1)
            }
        }

        label_set.set_label_set_for_node(i, labels);
    }

    std::vector<mapmap::UnaryTable<cost_t, simd_w> > unaries;  // 全都绑到一起: i, label, cost, 
    unaries.reserve(data_costs.cols());
    for (std::size_t i = 0; i < data_costs.cols(); ++i) {  // 遍历每个面
        DataCosts::Column const & data_costs_for_node = data_costs.col(i);  //  一列数据, 当前面所有view

        std::vector<mapmap::_s_t<cost_t, simd_w> > costs;  // 面对各个view的cost, 0~1, 不可见的cost为1, 如果这个面只有一个视角看到, 那么就存一个1, 如果有多个视角看到, 就存各个视角的cost
        if (data_costs_for_node.empty()) {
            costs.push_back(1.0f);
        }
        else{
            costs.resize(data_costs_for_node.size());  // 调整为能看到的view的个数
            for(std::size_t j = 0; j < data_costs_for_node.size(); ++j) {
                float cost = data_costs_for_node[j].second;  // face对每个view的cost
                costs[j] = cost;
            }
        }
        unaries.emplace_back(i, &label_set);  // 和label_set绑定
        unaries.back().set_costs(costs);  // costs的顺序和label_set的顺序一样, 全都绑到一起: i, label, cost
    }
    
    mapmap::PairwisePotts<cost_t, simd_w> pairwise(1.0f);  // Setting the pairwise costs

    // 展示消耗时间的函数
    auto display = [](const mapmap::luint_t time_ms,
            const mapmap::_iv_st<cost_t, simd_w> objective) {
        std::cout << "\t\t" << time_ms / 1000 << "\t" << objective << std::endl;
    };
    
    mapmap::StopWhenReturnsDiminish<cost_t, simd_w> terminate(5, 0.01);  // Changing the termination criteria: stop if after 5 iterations, less than an improvement of 0.01% in objective value has been made

    std::vector<mapmap::_iv_st<cost_t, simd_w> > solution;

    /* Create mapMAP solver object. */
    mapmap::mapMAP<cost_t, simd_w> solver;
    solver.set_graph(&mgraph);
    solver.set_label_set(&label_set);
    for(std::size_t i = 0; i < graph->num_nodes(); ++i)
        solver.set_unary(i, &unaries[i]);
    solver.set_pairwise(&pairwise);
    solver.set_logging_callback(display);
    solver.set_termination_criterion(&terminate);
    
    /* Pass configuration arguments (optional) for solve. */
    mapmap::mapMAP_control ctr;
    ctr.use_multilevel = true;
    ctr.use_spanning_tree = true;
    ctr.use_acyclic = true;
    ctr.spanning_tree_multilevel_after_n_iterations = 5;
    ctr.force_acyclic = true;
    ctr.min_acyclic_iterations = 5;
    ctr.relax_acyclic_maximal = true;
    ctr.tree_algorithm = mapmap::LOCK_FREE_TREE_SAMPLER;

    /* Set false for non-deterministic (but faster) mapMAP execution. */
    ctr.sample_deterministic = true;
    ctr.initial_seed = 548923723;

    std::cout << "\t 优化:\n\t\t时间 \t 能量" << std::endl;
    solver.optimize(solution, ctr);

    /* Label 0 is undefined. */
    std::size_t num_labels = data_costs.rows() + 1;
    std::size_t undefined = 0;
    
    /* Extract resulting labeling from solver. */
    for (std::size_t i = 0; i < graph->num_nodes(); ++i) {
        int label = label_set.label_from_offset(i, solution[i]);
        if (label < 0 || num_labels <= static_cast<std::size_t>(label)) {
            throw std::runtime_error("标签配置错误");
        }
        if (label == 0) undefined += 1;
        graph->set_label(i, static_cast<std::size_t>(label));
    }

    std::cout << '\t' << undefined << "个面片没有被任何视角看到" << std::endl;
}

TEX_NAMESPACE_END
