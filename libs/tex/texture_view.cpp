#include <list>

#include <math/matrix.h>
#include <mve/image_io.h>
#include <mve/image_tools.h>

#include "texture_view.h"
#include "tri.h"

TEX_NAMESPACE_BEGIN

TextureView::TextureView(std::size_t id, mve::CameraInfo const & camera, std::string const & image_file):id(id),image_file(image_file){
    mve::image::ImageHeaders header;  // 图像元数据

    try{
        header = mve::image::load_file_headers(image_file);
    }
    catch(util::Exception e) {
        std::cerr << "无法加载的图像" << image_file << std::endl;
        std::cerr << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    width = header.width;
    height = header.height;

    camera.fill_calibration(*projection, width, height); // 如果header中有图片的宽高, 则做校准(将相机坐标中的一点投影到具有“宽度”和“高度”尺寸的图像平面, 为了获得投影后的像素坐标，必须从坐标中减去0.5), 否则就是常规的投影
    camera.fill_camera_pos(*pos);  // 相机相对世界坐标的位置
    camera.fill_viewing_direction(*viewdir);
    camera.fill_world_to_cam(*world_to_cam);
}

void TextureView::load_image(void){
    if(image != NULL) return;
    image = mve::image::load_file(image_file);
}

void TextureView::generate_gradient_magnitude(void) {
    assert(image != NULL);
    mve::ByteImage::Ptr bw = mve::image::desaturate<std::uint8_t>(image, mve::image::DESATURATE_LUMINANCE);  // 得到灰度图
    gradient_magnitude = mve::image::sobel_edge<std::uint8_t>(bw);  // x, y方向都处理
}

/*
* 从四个角开始, 如果相邻的有黑色的像素, 就不断往外扩张, validity设为false
* */
void TextureView::generate_validity_mask(void){
    assert(image != NULL);
    validity_mask.resize(width * height, true);
    mve::ByteImage::Ptr checked = mve::ByteImage::create(width, height, 1);  // 如果已checked为255, 否则是0; 1表示单通道

    std::list<math::Vec2i> queue;

    // 从四个角开始
    queue.push_back(math::Vec2i(0,0));
        checked->at(0,0,0) = 255;
    queue.push_back(math::Vec2i(0, height-1));
        checked->at(0, height-1, 0) = 255;
    queue.push_back(math::Vec2i(width-1,0));
        checked->at(width-1,0,0) = 255;
    queue.push_back(math::Vec2i(width-1,height-1));
        checked->at(width-1,height-1,0) = 255;

    while(!queue.empty()){
        math::Vec2i pixel = queue.front();
        queue.pop_front();

        int const x = pixel[0];
        int const y = pixel[1];

        int sum = 0;  // 像素三个通道的值加起来
        for(int c=0; c<image->channels(); ++c){
            sum += image->at(x, y, c);
        }

        // std::cout<<sum<<std::endl;

        if(sum == 0){  // 如果像素是黑色的
            validity_mask[x + y * width] = false;

            std::vector<math::Vec2i> neighbours;  // 上下左右邻居
            neighbours.push_back(math::Vec2i(x+1, y));
            neighbours.push_back(math::Vec2i(x, y+1));
            neighbours.push_back(math::Vec2i(x-1, y));
            neighbours.push_back(math::Vec2i(x, y-1));

            for(std::size_t i=0; i<neighbours.size(); ++i){
                math::Vec2i npixel = neighbours[i];
                int const nx = npixel[0];
                int const ny = npixel[1];
                if (0 <= nx && nx < width && 0 <= ny && ny < height) {  // 在范围内
                    if (checked->at(nx, ny, 0) == 0) {
                        queue.push_front(npixel);
                        checked->at(nx, ny, 0) = 255;
                    }
                }
            }
        }
    }
}

/*
* 如果是0, 那么那周围的9个都是false
* */
void TextureView::erode_validity_mask(void) {
    std::vector<bool> eroded_validity_mask(validity_mask);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (x == 0 || x == width - 1 || y == 0 || y == height - 1) {  // 如果是最外一圈的像素
                validity_mask[x + y * width] = false;
                continue;
            }

            if (validity_mask[x + y * width])  // 如果是true
                continue;

            for (int j = -1; j <= 1; ++j) {  // 9个格子
                for (int i = -1; i <= 1; ++i) {
                    int const nx = x + i;
                    int const ny = y + j;
                    eroded_validity_mask[nx + ny * width] = false;
                }
            }
        }
    }

    validity_mask.swap(eroded_validity_mask);
}

bool TextureView::valid_pixel(math::Vec2f pixel) const {
    float const x = pixel[0];
    float const y = pixel[1];

    // 像素坐标不超过边界
    bool valid = (x >= 0.0f && x < static_cast<float>(width - 1) && 
                  y >= 0.0f && y < static_cast<float>(height - 1));

    if (valid && validity_mask.size() == static_cast<std::size_t>(width * height)) {
        // 只有可以正确插值的像素才有效
        float cx = std::max(0.0f, std::min(static_cast<float>(width - 1), x));
        float cy = std::max(0.0f, std::min(static_cast<float>(height - 1), y));
        int const floor_x = static_cast<int>(cx);
        int const floor_y = static_cast<int>(cy);
        int const floor_xp1 = std::min(floor_x + 1, width - 1);
        int const floor_yp1 = std::min(floor_y + 1, height - 1);

        valid = validity_mask[floor_x + floor_y * width] &&
                validity_mask[floor_x + floor_yp1 * width] &&
                validity_mask[floor_xp1 + floor_y * width] &&
                validity_mask[floor_xp1 + floor_yp1 * width];
    }

    return valid;
}

void TextureView::get_face_info(math::Vec3f const & v1, math::Vec3f const & v2, math::Vec3f const & v3, FaceProjectionInfo * face_info, Settings const & settings) const {
    assert(image != NULL);
    assert(settings.data_term != DATA_TERM_GMI || gradient_magnitude != NULL);
    
    // 将网格上的三维坐标转成view上的二维坐标
    math::Vec2f p1 = get_pixel_coords(v1);
    math::Vec2f p2 = get_pixel_coords(v2);
    math::Vec2f p3 = get_pixel_coords(v3);

    assert(valid_pixel(p1) && valid_pixel(p2) && valid_pixel(p3));

    Tri tri(p1, p2, p3);
    float area = tri.get_area();

    if (area < std::numeric_limits<float>::epsilon()) {  // 小于最小非零浮点数
        face_info->quality = 0.0f;
        return;
    }

    std::size_t num_samples = 0;
    math::Vec3d colors(0.0);  // 所有像素三个通道的颜色之和
    double gmi = 0.0;  // 所有像素梯度之和

    bool sampling_necessary = settings.data_term != DATA_TERM_AREA || settings.outlier_removal != OUTLIER_REMOVAL_NONE;  // !!

    if (sampling_necessary && area > 0.5f) {
        // 像素按y升序排
        while (true)
            if(p1[1] <= p2[1])
                if(p2[1] <= p3[1]) break;
                else std::swap(p2, p3);
            else std::swap(p1, p2);

        // bool fast_sampling_possible;


        Rect<float> aabb = tri.get_aabb();
        for (int y = std::floor(aabb.min_y); y < std::ceil(aabb.max_y); ++y) {
            float min_x = aabb.min_x - 0.5f;  // min 如果点在像素中间偏左, 则算上这个像素, 否则不算
            float max_x = aabb.max_x + 0.5f;  // 

            // if(fast_sampling_possible){}

            for (int x = std::floor(min_x + 0.5f); x < std::ceil(max_x - 0.5f); ++x) {
                math::Vec3d color;

                const float cx = static_cast<float>(x) + 0.5f;
                const float cy = static_cast<float>(y) + 0.5f;

                // if (!fast_sampling_possible && !tri.inside(cx, cy)) continue;
                if (!tri.inside(cx, cy)) continue;

                if (settings.outlier_removal != OUTLIER_REMOVAL_NONE) {
                    for (std::size_t c = 0; c < 3; c++){
                         color[c] = static_cast<double>(image->at(x, y, c)) / 255.0;  // 归一化
                    }
                    colors += color;
                }

                // 梯度累加, 如果投影大, gmi会大; 如果梯度大, gmi也会大
                if (settings.data_term == DATA_TERM_GMI) {
                    gmi += static_cast<double>(gradient_magnitude->at(x, y, 0)) / 255.0;
                }
                ++num_samples;
            }
        }
    }

    if (settings.data_term == DATA_TERM_GMI) {
        if (num_samples > 0) {
            gmi = (gmi / num_samples) * area;
        }
        else {
            double gmv1 = static_cast<double>(gradient_magnitude->linear_at(p1[0], p1[1], 0)) / 255.0;
            double gmv2 = static_cast<double>(gradient_magnitude->linear_at(p2[0], p2[1], 0)) / 255.0;
            double gmv3 = static_cast<double>(gradient_magnitude->linear_at(p3[0], p3[1], 0)) / 255.0;
            gmi = ((gmv1 + gmv2 + gmv3) / 3.0) * area;
        }
    }

    if (settings.outlier_removal != OUTLIER_REMOVAL_NONE) {
        if (num_samples > 0) {
            face_info->mean_color = colors / num_samples;
        } else {
            math::Vec3d c1, c2, c3;
            for (std::size_t i = 0; i < 3; ++i) {
                 c1[i] = static_cast<double>(image->linear_at(p1[0], p1[1], i)) / 255.0;
                 c2[i] = static_cast<double>(image->linear_at(p2[0], p2[1], i)) / 255.0;
                 c3[i] = static_cast<double>(image->linear_at(p3[0], p3[1], i)) / 255.0;
            }
            face_info->mean_color = ((c1 + c2 + c3) / 3.0);
        }
    }

    switch (settings.data_term) {
        case DATA_TERM_AREA: face_info->quality = area; break;
        case DATA_TERM_GMI:  face_info->quality = gmi; break;
    }
}

TEX_NAMESPACE_END