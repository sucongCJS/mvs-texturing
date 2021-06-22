#ifndef TEX_SPARSETABLE_HEADER
#define TEX_SPARSETABLE_HEADER

#include <vector>

#include <fstream>
#include <cstring>
#include <cerrno>
#include <iostream>

#include "util/file_system.h"
#include "util/exception.h"

/*
*  优化了行和列访问的稀疏表
* */
template<typename C, typename R, typename T>
class SparseTable{
public:
    typedef std::vector<std::pair<R, T>> Column;
    typedef std::vector<std::pair<C, T>> Row;

    std::size_t nnz;

private:
    std::vector<Column> column_wise_data;  // 索引代表第几列, 行用Column的第一个表示, 值用Column的第二个表示
    std::vector<Row> row_wise_data;  // 索引代表第几行, 列用Row的第一个表示, 值用Row的第二个表示

public:
    SparseTable();
    SparseTable(C cols, R rows);

    C cols() const;  // 返回size
    R rows() const;

    Column const & col(C id) const;
    Row const & row(R id) const;

    void set_value(C col, R row, T value);
};

template <typename C, typename R, typename T> C
SparseTable<C, R, T>::cols() const {
    return column_wise_data.size();
}

template <typename C, typename R, typename T> R
SparseTable<C, R, T>::rows() const {
    return row_wise_data.size();
}

template <typename C, typename R, typename T> typename SparseTable<C, R, T>::Column const &
SparseTable<C, R, T>::col(C id) const {
    return column_wise_data[id];
}

template <typename C, typename R, typename T> typename SparseTable<C, R, T>::Row const &
SparseTable<C, R, T>::row(R id) const {
    return row_wise_data[id];
}

template <typename C, typename R, typename T>
SparseTable<C, R, T>::SparseTable(){
    nnz = 0;  // 稀疏表中元素的个数
}

template <typename C, typename R, typename T>
SparseTable<C, R, T>::SparseTable(C cols, R rows){
    column_wise_data.resize(cols);
    row_wise_data.resize(rows);
    nnz = 0;
}

template <typename C, typename R, typename T> void
SparseTable<C, R, T>::set_value(C col, R row, T value){
    column_wise_data[col].push_back(std::pair<R, T>(row, value));
    row_wise_data[row].push_back(std::pair<C, T>(col, value));
    nnz++;
}

#endif /* TEX_SPARSETABLE_HEADER */
