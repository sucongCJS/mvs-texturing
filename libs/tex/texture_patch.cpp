#include <set>

#include <math/functions.h>
#include <mve/image_color.h>
#include <mve/image_tools.h>
#include <mve/mesh_io_ply.h>

#include "texture_patch.h"

TexturePatch::TexturePatch(int label, std::vector<std::size_t> const & faces,
    std::vector<math::Vec2f>  const & texcoords, mve::FloatImage::Ptr image)
    : label(label), faces(faces), texcoords(texcoords), image(image) {

    // validity_mask = mve::ByteImage::create(get_width(), get_height(), 1);
    // validity_mask->fill(255);
    // blending_mask = mve::ByteImage::create(get_width(), get_height(), 1);
}

TexturePatch::TexturePatch(TexturePatch const & texture_patch) {
    label = texture_patch.label;
    faces = std::vector<std::size_t>(texture_patch.faces);
    texcoords = std::vector<math::Vec2f>(texture_patch.texcoords);
    image = texture_patch.image->duplicate();

    // validity_mask = texture_patch.validity_mask->duplicate();
    // if (texture_patch.blending_mask != NULL) {
    //     blending_mask = texture_patch.blending_mask->duplicate();
    // }
}

bool TexturePatch::valid_pixel(math::Vec2f pixel) const {
    float x = pixel[0];
    float y = pixel[1];

    float const height = static_cast<float>(get_height());
    float const width = static_cast<float>(get_width());
}

math::Vec3f TexturePatch::get_pixel_value(math::Vec2f pixel) const {
    assert(valid_pixel(pixel));

    math::Vec3f color;
    image->linear_at(pixel[0], pixel[1], *color);
    return color;
}
