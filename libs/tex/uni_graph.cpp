#include <limits>
#include <list>

#include "uni_graph.h"

UniGraph::UniGraph(std::size_t nodes) {
    adj_lists.resize(nodes);
    labels.resize(nodes);
    edges = 0;
}

/**
 * @brief 用给定标签的所有subgraph填充subgraphs, subgraph是包含具有<b>相同标签的连通的</b>面片的索引的向量
 * @param label 要寻找的标签
 * @param subgraphs subgraph是一块连通的用了同一label的面片, subgraph与subgraph之间是不连通的
 */
void UniGraph::get_subgraphs(std::size_t label, std::vector<std::vector<std::size_t> > * subgraphs) const {
    std::vector<bool> used(adj_lists.size(), false);  // 面片是否已经被添加到subgraphs, 防止重复添加

    for(std::size_t i = 0; i < adj_lists.size(); ++i) {  // 遍历网格每个面片的相邻列表
        if (labels[i] == label && !used[i]) {  // 如果用的是同一个标签
            subgraphs->push_back(std::vector<std::size_t>());  // 添加一个subgraph
            
            std::list<std::size_t> queue;

            queue.push_back(i);
            used[i] = true;

            while (!queue.empty()) {  // 搜索这个面片相邻的, 以及相邻的相邻面片有没有用的是同一label的
                std::size_t node = queue.front();
                queue.pop_front();

                subgraphs->back().push_back(node);  // 往subgraph里加元素

                // 将具有相同标签的所有used为false的相邻面片添加到队列中
                std::vector<std::size_t> const & adj_list = adj_lists[node];
                for(std::size_t j = 0; j < adj_list.size(); ++j) {
                    std::size_t adj_node = adj_list[j];
                    assert(adj_node < labels.size() && adj_node < used.size());
                    if (labels[adj_node] == label && !used[adj_node]){
                        queue.push_back(adj_node);
                        used[adj_node] = true;
                    }
                }
            }
        }
    }
}
