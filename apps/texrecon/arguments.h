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

    int num_threads;

    tex::Settings settings;
};

/*
* 
* */
Arguments parse_args(int argc, char **argv);

#endif /* ARGUMENTS_HEADER */
