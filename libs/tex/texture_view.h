#ifndef TEX_TEXTUREVIEW_HEADER
#define TEX_TEXTUREVIEW_HEADER

#include <string>
#include <vector>
#include <cassert>

#include <math/vector.h>
#include <mve/camera.h>
#include <mve/image.h>

#include "settings.h"

TEX_NAMESPACE_BEGIN

/*
* 这个面在一个view中的质量和平均颜色信息
* */
struct FaceProjectionInfo{
    std::uint16_t view_id;  // 当前face对应的view
    float quality;  // 当前face在这个view的quality
    math::Vec3f mean_color;

    bool operator<(FaceProjectionInfo const & other) const {
        return view_id < other.view_id;
    }
};

class TextureView{
private:
    std::size_t id;  // (不变)
    std::string image_file;  // 图片名
    int width;  // 图片宽
    int height;  // 图片高
    math::Matrix3f projection;  // 图片校准矩阵 (将相机坐标中的一点投影到具有“宽度”和“高度”尺寸的图像平面, 得到的是像素中心的位置, 要得到投影后的像素坐标，必须从各个坐标中减去0.5)
    math::Vec3f pos;  // 相机位置(世界坐标)
    math::Vec3f viewdir;  // 视角的方向(世界坐标)
    math::Matrix4f world_to_cam;  // 世界坐标转到视角的坐标

    mve::ByteImage::Ptr image;  // RGBRGB...
    mve::ByteImage::Ptr gradient_magnitude;  // 
    std::vector<bool> validity_mask;

public:
    std::size_t get_id(void) const;  // 返回纹理视角的id
    math::Vec3f get_viewing_direction(void) const;
    int get_width(void) const;
    int get_height(void) const;
    mve::ByteImage::Ptr get_image(void) const;

    // 参数: 视角id, 视角的相机信息, 视角的图片
    TextureView(std::size_t id, mve::CameraInfo const & camera, std::string const & image_file);

    // 返回位置
    math::Vec3f get_pos(void) const;

    // 生成 validity mask
    void generate_validity_mask(void);
    // 加载相应的图片
    void load_image(void);
    // 生成梯度图
    void generate_gradient_magnitude(void);
    // 侵蚀mask一个像素
    void erode_validity_mask(void);

    // 释放 validity mask
    void release_validity_mask(void);
    // 释放梯度图
    void release_gradient_magnitude(void);
    // 释放相应的图片
    void release_image(void);
    
    // 返回给定点投影到视角中的二维像素坐标
    math::Vec2f get_pixel_coords(math::Vec3f const & vertex) const;

    // 返回像素位置在此视角中是否有效
    // 如果像素位置在可见区域内，则该位置有效
    // 如果生成了validity mask，则所有周围（整数坐标）像素在validity mask中都有效。
    bool valid_pixel(math::Vec2f pixel) const;

    // v1,v2,v3为顶点的三角形是否有效, 点从世界坐标转到view坐标是否还在view的范围内, 以及看validity mask
    bool inside(math::Vec3f const & v1, math::Vec3f const & v2, math::Vec3f const & v3) const;
    
    // 
    void get_face_info(math::Vec3f const & v1, math::Vec3f const & v2, math::Vec3f const & v3, FaceProjectionInfo * face_info, Settings const & settings) const;
};

inline std::size_t TextureView::get_id(void) const{ return id; }

inline math::Vec3f TextureView::get_pos(void) const { return pos; }

inline math::Vec3f TextureView::get_viewing_direction(void) const { return viewdir; }

inline int TextureView::get_width(void) const { return width; }

inline int TextureView::get_height(void) const { return height; }

inline mve::ByteImage::Ptr TextureView::get_image(void) const {
    assert(image != NULL);
    return image;
}

inline void TextureView::release_image(void) {
    assert(image != NULL);
    image.reset();
}

inline void TextureView::release_gradient_magnitude(void) {
    assert(gradient_magnitude != NULL);
    gradient_magnitude.reset();
}

inline void TextureView::release_validity_mask(void) {
    assert(validity_mask.size() == static_cast<std::size_t>(width * height));
    validity_mask = std::vector<bool>();
}

inline math::Vec2f TextureView::get_pixel_coords(math::Vec3f const & vertex) const {
    math::Vec3f pixel = projection * world_to_cam.mult(vertex, 1.0f);  // 
    pixel /= pixel[2];  // 齐次归一
    return math::Vec2f(pixel[0] - 0.5f, pixel[1] - 0.5f);  // 从像素中心到右上角
}

inline bool TextureView::inside(math::Vec3f const & v1, math::Vec3f const & v2, math::Vec3f const & v3) const {
    // 获得投影到视角坐标的二维位置
    math::Vec2f p1 = get_pixel_coords(v1);
    math::Vec2f p2 = get_pixel_coords(v2);
    math::Vec2f p3 = get_pixel_coords(v3);
    return valid_pixel(p1) && valid_pixel(p2) && valid_pixel(p3);
}

TEX_NAMESPACE_END

#endif /* TEX_TEXTUREVIEW_HEADER */