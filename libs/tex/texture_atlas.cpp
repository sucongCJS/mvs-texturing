#include <set>
#include <map>

#include <util/file_system.h>
#include <mve/image_tools.h>
#include <mve/image_io.h>

#include "texture_atlas.h"

TextureAtlas::TextureAtlas(unsigned int size) :
    size(size), padding(size >> 7), finalized(false) {

    bin = RectangularBin::create(size, size);
    image = mve::ByteImage::create(size, size, 3);
    // validity_mask = mve::ByteImage::create(size, size, 1);
}

/**
 * @brief 将src复制到dest
 */
void copy_into(mve::ByteImage::ConstPtr src, int x, int y,
    mve::ByteImage::Ptr dest, int border = 0) {
    
    assert(x >= 0 && x + src->width() + 2 * border <= dest->width());
    assert(y >= 0 && y + src->height() + 2 * border <= dest->height());

    for(int i = 0; i < src->width() + 2 * border; ++i){
        for(int j = 0; j < src->height() + 2 * border; j++){
            int sx = i - border;
            int sy = j - border;

            if(sx < 0 || sx >= src->width() || sy < 0 || sy >= src->height())
                continue;

            for(int c = 0; c < src->channels(); ++c){
                dest->at(x + i, y + j, c) = src->at(sx, sy, c);
            }
        }
    }
}

/**
 * @brief 向纹理集中插入纹理块texture_patch
 */
bool TextureAtlas::insert(TexturePatch::ConstPtr texture_patch) {
    if (finalized) {
        throw util::Exception("此纹理集已经处理完了, 不能再插入纹理块了");
    }

    assert(bin != NULL);

    int const width = texture_patch->get_width() + 2 * padding;
    int const height = texture_patch->get_height() + 2 * padding;
    Rect<int> rect(0, 0, width, height);  // 准备插入到纹理集后的方框
    if (!bin->insert(&rect)) return false;

    mve::ByteImage::Ptr patch_image = mve::image::float_to_byte_image(texture_patch->get_image(), 0.0f, 1.0f);  // 转成[0,255]格式

    copy_into(patch_image, rect.min_x, rect.min_y, image, padding);
    // mve::ByteImage::ConstPtr patch_validity_mask = texture_patch->get_validity_mask();
    // copy_into(patch_validity_mask, rect.min_x, rect.min_y, validity_mask, padding);

    TexturePatch::Faces const & patch_faces = texture_patch->get_faces();
    TexturePatch::Texcoords const & patch_texcoords = texture_patch->get_texcoords();

    // 将纹理块的面片都加到纹理集中
    faces.insert(faces.end(), patch_faces.begin(), patch_faces.end());
    
    // 将纹理块的顶点二维坐标都加到纹理集中
    math::Vec2f offset = math::Vec2f(rect.min_x + padding, rect.min_y + padding);  // 计算纹理块相对纹理坐标的偏移量
    for (std::size_t i = 0; i < patch_faces.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            math::Vec2f rel_texcoord(patch_texcoords[i * 3 + j]);
            math::Vec2f texcoord = rel_texcoord + offset;

            texcoord[0] = texcoord[0] / this->size;  //...
            texcoord[1] = texcoord[1] / this->size;
            
            texcoords.push_back(texcoord);
        }
    }
    return true;
}

// typedef std::set<std::pair<int, int> > PixelSet;
// typedef std::vector<std::pair<int, int> > PixelVector;

// void TextureAtlas::apply_edge_padding(void) {
//     assert(image != NULL);
//     // assert(validity_mask != NULL);

//     const int width = image->width();
//     const int height = image->height();

//     math::Matrix<float, 3, 3> gauss;
//     gauss[0] = 1.0f; gauss[1] = 2.0f; gauss[2] = 1.0f;
//     gauss[3] = 2.0f; gauss[4] = 4.0f; gauss[5] = 2.0f;
//     gauss[6] = 1.0f; gauss[7] = 2.0f; gauss[8] = 1.0f;
//     gauss /= 16.0f;

//     // 计算纹理块边界处的无效像素集
//     PixelSet invalid_border_pixels;
//     // for (int y = 0; y < height; ++y) {
//     //     for (int x = 0; x < width; ++x) {
//     //         // if (validity_mask->at(x, y, 0) == 255) continue;

//     //         // 检查所有无效像素的直接邻域
//     //         for (int j = -1; j <= 1; ++j) {
//     //             for (int i = -1; i <= 1; ++i) {
//     //                 ...
//     //     }
//     // }

//     for (unsigned int n = 0; n <= padding; ++n) {
//         PixelVector new_valid_pixels;
        
//         PixelSet::iterator it = invalid_border_pixels.begin();

//     }
// }

struct VectorCompare {
    bool operator()(math::Vec2f const & lhs, math::Vec2f const & rhs) const {
        return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
    }
};

typedef std::map<math::Vec2f, std::size_t, VectorCompare> TexcoordMap;

void TextureAtlas::merge_texcoords(){
    Texcoords tmp; 
    tmp.swap(this->texcoords);

    TexcoordMap texcoord_map;
    for (math::Vec2f const & texcoord : tmp) {
        TexcoordMap::iterator iter = texcoord_map.find(texcoord);
        if (iter == texcoord_map.end()) {
            std::size_t texcoord_id = this->texcoords.size();
            texcoord_map[texcoord] = texcoord_id;
            this->texcoords.push_back(texcoord);
            this->texcoord_ids.push_back(texcoord_id);
        }
        else{
            this->texcoord_ids.push_back(iter->second);
        }
    }
}

void TextureAtlas::finalize() {
    if (finalized) {
        throw util::Exception("当前纹理集已经处理完毕了");
    }

    this->bin.reset();
    // this->apply_edge_padding();
    // this->validity_mask.reset();
    this->merge_texcoords();

    this->finalized = true;
}
