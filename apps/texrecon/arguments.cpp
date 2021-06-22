#include "arguments.h"
#include "util/file_system.h"

Arguments parse_args(int argc, char **argv){
    util::Arguments args;
    args.set_exit_on_error(true);  // 参数解析出错退出
    args.set_nonopt_maxnum(3);  // 必须有三个非选项参数: 
    args.set_nonopt_minnum(3);
    args.set_helptext_indent(34);
    args.set_description("作用: 使用3D场景形式的图片给网格贴纹理");
    args.set_usage("使用方法: " + std::string(argv[0]) + " [选项] scene目录::图片名 mesh文件 输出前缀"
        "\n - scene目录:使用MVE得到的一个包含相机参数, 图片等信息的目录"
        "\n - 图片名: scene目录中图片的名, 不带扩展名"
        "\n - mesh文件: 要对其进行纹理处理的.ply网格文件, 可以使用MVE生成"
        "\n - 输出前缀: 输出文件的路径和名称, 不需要附加对象扩展名, 因为程序会输出多个包含这个前缀的文件(网格、材质文件、纹理文件)"
    );

    args.parse(argc, argv);

    Arguments conf;
    conf.in_scene = args.get_nth_nonopt(0);
    conf.in_mesh = args.get_nth_nonopt(1);
    conf.out_prefix = util::fs::sanitize_path(args.get_nth_nonopt(2));

    conf.num_threads = -1;

    return conf;
}