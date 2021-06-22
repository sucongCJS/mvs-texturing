#ifndef TEX_TEXTUREPATCH_HEADER
#define TEX_TEXTUREPATCH_HEADER

#include <vector>
#include <cassert>

#include <math/vector.h>
#include <mve/mesh.h>

int const texture_patch_border = 1;

/**
  * 纹理块
  */
class TexturePatch {
public:
    typedef std::shared_ptr<TexturePatch> Ptr;
    typedef std::shared_ptr<const TexturePatch> ConstPtr;
    typedef std::vector<std::size_t> Faces;
    typedef std::vector<math::Vec2f> Texcoords;


private:
    int label;  // 纹理块的标签
    Faces faces;  // 使用了这个纹理的面片的id
    Texcoords texcoords;  // 纹理块面片顶点在纹理上的二维坐标
    mve::FloatImage::Ptr image;  // 纹理图片

public:
    TexturePatch(int _label, std::vector<std::size_t> const & _faces, std::vector<math::Vec2f>  const & _texcoords, mve::FloatImage::Ptr _image);
    TexturePatch(TexturePatch const & texture_patch);

    mve::FloatImage::Ptr get_image(void);
    mve::FloatImage::ConstPtr get_image(void) const;
    int get_label(void) const;
    int get_width(void) const;
    int get_height(void) const;
    int get_size(void) const;

    bool valid_pixel(math::Vec2i pixel) const;
    bool valid_pixel(math::Vec2f pixel) const;

    std::vector<std::size_t> & get_faces(void);
    std::vector<std::size_t> const & get_faces(void) const;
    std::vector<math::Vec2f> & get_texcoords(void);
    std::vector<math::Vec2f> const & get_texcoords(void) const;

    math::Vec3f get_pixel_value(math::Vec2f pixel) const;

    static TexturePatch::Ptr create(TexturePatch::ConstPtr texture_patch);
    static TexturePatch::Ptr create(int label, std::vector<std::size_t> const & faces, std::vector<math::Vec2f> const & texcoords, mve::FloatImage::Ptr image);
};

inline mve::FloatImage::Ptr TexturePatch::get_image(void) {
    return image;
}

inline mve::FloatImage::ConstPtr TexturePatch::get_image(void) const {
    return image;
}

inline int TexturePatch::get_label(void) const {
    return label;
}

inline int TexturePatch::get_width(void) const {
    return image->width();
}

inline int TexturePatch::get_height(void) const {
    return image->height();
}

inline int TexturePatch::get_size(void) const {
    return get_width() * get_height();
}

inline TexturePatch::Ptr TexturePatch::create(TexturePatch::ConstPtr texture_patch) {
    return std::make_shared<TexturePatch>(*texture_patch);
}

inline TexturePatch::Ptr TexturePatch::create(int label, std::vector<std::size_t> const & faces, std::vector<math::Vec2f>  const & texcoords, mve::FloatImage::Ptr image) {
    return std::make_shared<TexturePatch>(label, faces, texcoords, image);
}

inline std::vector<math::Vec2f> const & TexturePatch::get_texcoords(void) const {
    return texcoords;
}

inline std::vector<math::Vec2f> & TexturePatch::get_texcoords(void) {
    return texcoords;
}

inline std::vector<std::size_t> & TexturePatch::get_faces(void) {
    return faces;
}

inline std::vector<std::size_t> const & TexturePatch::get_faces(void) const {
    return faces;
}

#endif /* TEX_TEXTUREPATCH_HEADER */
