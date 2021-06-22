#ifndef TEX_UNIGRAPH_HEADER
#define TEX_UNIGRAPH_HEADER

#include <vector>
#include <cassert>
#include <algorithm>

/**
  * 用邻接表实现固定节点数的单向图
  */
class UniGraph {
private:
    std::vector<std::vector<std::size_t> > adj_lists;  // 一个列表里包含 每个面相邻的面的列表, adj_lists.size()=网格的面片数
    std::vector<std::size_t> labels;  // 面的标签
    std::size_t edges;  // 总共的边数, 如果n1, n2之间有边, 只算一条
    
public:
    /*
    * 输入节点数, 单向图, 没有边
    * */
    UniGraph(std::size_t nodes);

    /*
    * 判断n1, n2之间是否有边
    * */
    bool has_edge(std::size_t n1, std::size_t n2) const;

    /*
    * 添加边
    * */
    void add_edge(std::size_t n1, std::size_t n2);

    /*
    * 返回节点数
    * */
    std::size_t num_nodes() const;

    /**
     * Returns the label of node with index n.
     * @warning asserts that the index is valid.
     */
    std::size_t get_label(std::size_t n) const;

    /*
    * 返回边数
    * */
    std::size_t num_edges() const;

    std::vector<std::size_t> const & get_adj_nodes(std::size_t node) const;

    /**
     * 面n的标签设置为label
     */
    void set_label(std::size_t n, std::size_t label);

    /**
     * 用给定标签的所有子图填充subgraphs
     * 子图是包含具有相同标签的连通的图的所有索引的向量
     */
    void get_subgraphs(std::size_t label, std::vector<std::vector<std::size_t> > * subgraphs) const;
};

inline std::size_t UniGraph::get_label(std::size_t n) const {
    assert(n < num_nodes());
    return labels[n];
}

inline bool UniGraph::has_edge(std::size_t n1, std::size_t n2) const{
    assert(n1 < num_nodes() && n2 < num_nodes());
    std::vector<std::size_t> const & adj_list = adj_lists[n1];
    return std::find(adj_list.begin(), adj_list.end(), n2) != adj_list.end();
}

inline void UniGraph::add_edge(std::size_t n1, std::size_t n2){
    assert(n1 < num_nodes() && n2 < num_nodes());
    if(!has_edge(n1, n2)){
        adj_lists[n1].push_back(n2);
        adj_lists[n2].push_back(n1);
        ++edges;
    }
}

inline std::size_t UniGraph::num_nodes() const{
    return adj_lists.size();
}

inline std::size_t UniGraph::num_edges() const{
    return edges;
}

inline std::vector<std::size_t> const & UniGraph::get_adj_nodes(std::size_t node) const {
    assert(node < num_nodes());
    return adj_lists[node];
}

inline void UniGraph::set_label(std::size_t n, std::size_t label) {
    assert(n < num_nodes());
    labels[n] = label;
}

#endif /* TEX_UNIGRAPH_HEADER */
