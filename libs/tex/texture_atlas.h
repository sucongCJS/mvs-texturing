#ifndef TEX_TEXTUREATLAS_HEADER
#define TEX_TEXTUREATLAS_HEADER

#include <vector>

#include <util/exception.h>
#include <math/vector.h>
#include <mve/mesh.h>
#include <mve/image.h>

#include "texture_patch.h"
#include "rectangular_bin.h"

/*
* 纹理集
* */
class TextureAtlas {
public:
    typedef std::shared_ptr<TextureAtlas> Ptr;
    typedef std::vector<std::size_t> Faces;
    typedef std::vector<math::Vec2f> Texcoords;
    typedef std::vector<std::size_t> TexcoordIds;

private:
    unsigned int const size;  // 纹理集的大小, 纹理集是正方形的
    unsigned int const padding;  // 
    bool finalized;  // 是否已用纹理块填充完

    Faces faces;  // 纹理集中所有纹理块的面片
    Texcoords texcoords;  // 面片顶点在纹理集上的二维坐标
    TexcoordIds texcoord_ids;

    // void apply_edge_padding(void);
    void merge_texcoords(void);

public:
    TextureAtlas(unsigned int size);

    mve::ByteImage::Ptr image;  // 纹理集图片, 包含很多的纹理块
    RectangularBin::Ptr bin;

    static TextureAtlas::Ptr create(unsigned int size);

    Faces const & get_faces(void) const;
    TexcoordIds const & get_texcoord_ids(void) const;
    Texcoords const & get_texcoords(void) const;
    mve::ByteImage::ConstPtr get_image(void) const;

    bool insert(TexturePatch::ConstPtr texture_patch);

    void finalize(void);
};

inline TextureAtlas::Ptr TextureAtlas::create(unsigned int size) {
    return Ptr(new TextureAtlas(size));
}

inline TextureAtlas::Faces const & TextureAtlas::get_faces(void) const {
    return faces;
}

inline TextureAtlas::TexcoordIds const & TextureAtlas::get_texcoord_ids(void) const {
    return texcoord_ids;
}

inline TextureAtlas::Texcoords const & TextureAtlas::get_texcoords(void) const {
    return texcoords;
}

inline mve::ByteImage::ConstPtr TextureAtlas::get_image(void) const {
    if (!finalized) {
        throw util::Exception("Texture atlas not finalized");
    }
    return image;
}

#endif /* TEX_TEXTUREATLAS_HEADER */
