#include <util/timer.h>
#include <util/tokenizer.h>
#include <util/file_system.h>
#include <mve/image_io.h>
#include <mve/image_tools.h>
#include <mve/bundle_io.h>
#include <mve/scene.h>

#include "texturing.h"
#include "progress_counter.h"

TEX_NAMESPACE_BEGIN

void from_mve_scene(std::string const & scene_dir, std::string const & image_name, TextureViews * texture_views){
    mve::Scene::Ptr scene;

    try{
        scene = mve::Scene::create(scene_dir);
    }
    catch(std::exception & e){
        std::cerr<<"无法打开scene目录: "<<e.what()<<std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::size_t num_views = scene->get_views().size();
    texture_views->reserve(num_views);

    ProgressCounter view_counter("\t加载视角图片: ", num_views);
    for(std::size_t i=0; i<num_views; ++i){  // 遍历每一个view
        view_counter.progress<SIMPLE>();

        mve::View::Ptr view = scene->get_view_by_id(i);
        if(view == NULL){
            view_counter.inc();
            continue;
        }

        if(!view->has_image(image_name, mve::IMAGE_TYPE_UINT8)){
            std::cout<<"警告: 视角"<<view->get_name()<<"没有字节图片"<<image_name<<std::endl;
            continue;
        }

        mve::View::ImageProxy const * image_proxy = view->get_image_proxy(image_name);

        if (image_proxy->channels < 3) {
            std::cerr<<"图片"<<image_name<<"在视角"<<view->get_name()<<"中不是三通道的图片!"<<std::endl;
            exit(EXIT_FAILURE);
        }

        texture_views->push_back(
            TextureView(view->get_id(), view->get_camera(), util::fs::abspath(util::fs::join_path(view->get_directory(), image_proxy->filename))));
        view_counter.inc();
    }
}

void generate_texture_views(std::string const & in_scene, TextureViews * texture_views, std::string const & tmp_dir){
    /* MVE_SCENE::EMBEDDING */
    size_t pos = in_scene.rfind("::");
    if (pos != std::string::npos) {
        std::string scene_dir = in_scene.substr(0, pos);
        std::string image_name = in_scene.substr(pos + 2, in_scene.size());
        from_mve_scene(scene_dir, image_name, texture_views);
    }

    std::sort(texture_views->begin(), texture_views->end(),
        [] (TextureView const & l, TextureView const & r) -> bool {
            return l.get_id() < r.get_id();
        }
    );  // 按id排序

    std::size_t num_views = texture_views->size();
    if(num_views == 0){
        std::cerr<< "没有正确地导入scene目录里的内容"<< std::endl;
        exit(EXIT_FAILURE);
    }

}

TEX_NAMESPACE_END