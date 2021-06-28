#ifndef ARGUMENTS_HEADER
#define ARGUMENTS_HEADER

#include "util/arguments.h"
#include "tex/settings.h"

/*
* 命令行中的参数结构
* */
struct Arguments{
    std::string in_scene;  // 输入的3D scene文件
    std::string in_mesh;  // 输入的网格文件
    std::string out_prefix;  // 输出文件的前缀
    std::string labeling_file;  // 需要读取的标签文件

    int num_threads;

    tex::Settings settings;

    bool write_view_selection_model = true;  // 使用同一视角的面片赋予一样的颜色
    bool write_intermediate_results = true;
};

/*
* 
* */
Arguments parse_args(int argc, char **argv);

#endif /* ARGUMENTS_HEADER */
