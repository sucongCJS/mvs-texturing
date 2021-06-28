#include "debug.h"

TEX_NAMESPACE_BEGIN

void generate_debug_colors(std::vector<math::Vec4f> & colors) {
    for (float s = 1.0f; s > 0.0f; s -= 0.4) {
        for (float v = 1.0f; v > 0.0f; v -= 0.3) {
            for (float h = 0.0f; h < 360.0f; h += 30.0f) {
                float c = v * s;
                float x = c * (1.0f - fabs(fmod(h / 60.0f, 2.0f) - 1.0f));
                float m = v - c;

                math::Vec4f color;
                if (0 <= h && h < 60)
                    color = math::Vec4f(c, x, 0.0f, 1.0f);
                if (60 <= h && h < 120)
                    color = math::Vec4f(x, c, 0.0f, 1.0f);
                if (120 <= h && h < 180)
                    color = math::Vec4f(0.0f, c, x, 1.0f);
                if (180 <= h && h < 240)
                    color = math::Vec4f(0.0f, x, c, 1.0f);
                if (240 <= h && h < 300)
                    color = math::Vec4f(x, 0.0f, c, 1.0f);
                if (300 <= h && h < 360)
                    color = math::Vec4f(c, 0.0f, x, 1.0f);

                color = color + math::Vec4f(m, m, m, 0.0f);
                colors.push_back(color);
            }
        }
    }
}

void generate_debug_embeddings(std::vector<TextureView> * texture_views){
    std::vector<math::Vec4f> colors;
    generate_debug_colors(colors);

    for(std::size_t i=0; i<texture_views->size(); ++i){
        math::Vec4f float_color = colors[i % colors.size()];

        TextureView * texture_view = &(texture_views->at(i));  // 必须要用指针或者引用!!!

        // 根据背景颜色确定字体颜色
        float luminance = math::interpolate(float_color[0], float_color[1], float_color[2], 0.30f, 0.59f, 0.11f);  // 转成灰色
        // math::Vec3uc font_color = luminance > 0.5f ? math::Vec3uc(0,0,0) : math::Vec3uc(255,255,255);  // 转成黑白??
        math::Vec3uc font_color = math::Vec3uc(luminance, luminance, luminance);
        
        math::Vec3uc bg_color;
        bg_color[0] = float_color[0] * 255.0f;
        bg_color[1] = float_color[1] * 255.0f;
        bg_color[2] = float_color[2] * 255.0f;

        mve::ByteImage::Ptr image = mve::ByteImage::create(texture_view->get_width(), texture_view->get_height(), 3);
        image->fill_color(*bg_color);

        texture_view->bind_image(image);
    }
}

TEX_NAMESPACE_END
