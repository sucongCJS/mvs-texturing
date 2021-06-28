#include <set>
#include <list>
#include <iostream>
#include <fstream>

#include <util/timer.h>
#include <mve/image_tools.h>

#include "defines.h"
#include "settings.h"
#include "texture_atlas.h"
#include "texture_patch.h"

#define MAX_TEXTURE_SIZE (8 * 1024)  // 最大的宽高
#define PREF_TEXTURE_SIZE (4 * 1024)
#define MIN_TEXTURE_SIZE (256)  // 纹理集尺寸的下限

TEX_NAMESPACE_BEGIN

/**
 * @brief 计算纹理集大致的大小
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 */
unsigned int calculate_texture_atlas_size(std::list<TexturePatch::ConstPtr> const & texture_patches) {
    unsigned int size = MAX_TEXTURE_SIZE;  // 纹理集大小, 先从最大的开始, 然后不断缩小到合适的

    while (true) {
        unsigned int total_area = 0;  // 总的纹理块面积
        unsigned int padding = size >> 7;  // 随着size的变化, padding也会变化
        unsigned int max_width = 0;  // 所有纹理块中最大的宽
        unsigned int max_height = 0;  // 所有纹理块中最大的高

        // 找最大的纹理块的宽高, 以及计算总的纹理块占用的面积
        for (TexturePatch::ConstPtr texture_patch : texture_patches) {
            unsigned int width = texture_patch->get_width() + 2 * padding;
            unsigned int height = texture_patch->get_height() + 2 * padding;
            
            max_width = std::max(max_width, width);
            max_height = std::max(max_height, height);

            unsigned int area = width * height;
            unsigned int waste = area - texture_patch->get_size();  // padding部分

            // padding部分比原来的纹理块还大, 则跳过
            if (static_cast<double>(waste) / texture_patch->get_size() > 1.0) {  
                break;  // 因为是降序排的, 后面的纹理块也就不用看了
            }

            total_area += area;
        }

        assert(max_width < MAX_TEXTURE_SIZE);
        assert(max_height < MAX_TEXTURE_SIZE);
        // 缩小size
        if (size > PREF_TEXTURE_SIZE && 
            max_width < PREF_TEXTURE_SIZE && max_height < PREF_TEXTURE_SIZE &&
            total_area / (PREF_TEXTURE_SIZE * PREF_TEXTURE_SIZE) < 8) {
            
            size = PREF_TEXTURE_SIZE;
            continue;
        }

        // 纹理集尺寸的下限
        if (size <= MIN_TEXTURE_SIZE) {
            return MIN_TEXTURE_SIZE;
        }

        // 如果size大于总的纹理面积太多, 缩小size的大小
        if (max_height < size / 2 && max_width < size / 2 &&
            static_cast<double>(total_area) / (size * size) < 0.2) {
            size = size / 2;
            continue;
        }

        return size;
    }
}

bool comp(TexturePatch::ConstPtr first, TexturePatch::ConstPtr second) {
    return first->get_size() > second->get_size();
}

/**
 * @brief 
 * @param texture_atlases 所有纹理集的集合, 包含所有的纹理块, 一个texture_atlas只包含大小差不多的纹理块
 */
void generate_texture_atlases(std::vector<TexturePatch::Ptr> * orig_texture_patches, Settings const & settings, std::vector<TextureAtlas::Ptr> * texture_atlases) {

    std::list<TexturePatch::ConstPtr> texture_patches;

    // gamma较正一下
    while (!orig_texture_patches->empty()) {
        TexturePatch::Ptr texture_patch = orig_texture_patches->back();
        orig_texture_patches->pop_back();

        if (settings.tone_mapping != TONE_MAPPING_NONE) {
            mve::image::gamma_correct(texture_patch->get_image(), 1.0f / 2.2f);
        }

        texture_patches.push_back(texture_patch);
    }

    // 通过对纹理块按大小降序排序, 提高把纹理块(patch)装到纹理集(atalas)中的效率
    std::cout << "\t给纹理块排序... " << std::flush;
    texture_patches.sort(comp);
    std::cout << "完成." << std::endl;

    std::size_t const total_num_patches = texture_patches.size();
    std::size_t remaining_patches = texture_patches.size();
    // std::ofstream tty("/dev/tty", std::ios_base::out);

    #pragma omp parallel
    {
    #pragma omp single
    {
        
    while (!texture_patches.empty()) {
        unsigned int texture_atlas_size = calculate_texture_atlas_size(texture_patches);  // 获得最佳的纹理集大小

        texture_atlases->push_back(TextureAtlas::create(texture_atlas_size));
        TextureAtlas::Ptr texture_atlas = texture_atlases->back();
        
        // 尝试将每个纹理块插入到纹理集中
        for (std::list<TexturePatch::ConstPtr>::iterator it = texture_patches.begin(); 
             it != texture_patches.end();) {
            // std::size_t done_patches = total_num_patches - remaining_patches;
            // int precent = static_cast<float>(done_patches) / total_num_patches * 100.0f;
            // if (total_num_patches > 100 && 
            //     done_patches % (total_num_patches / 100) == 0) {
            //     tty << "\r\t Working on atlas " << texture_atlases->size() << " "<< precent << "%..." << std::flush;
            // }

            if (texture_atlas->insert(*it)) {
                it = texture_patches.erase(it);
                remaining_patches -= 1;
            }
            else {
                ++it;
            }
        }
        texture_atlas->finalize();
    }

    std::cout << "\r\t处理" << texture_atlases->size() << "个纹理集 ... 完成." << std::endl;

    #pragma omp taskwait
    std::cout << "\t所有的纹理集生成完毕" << std::endl;

    /* End of single region */
    }
    /* End of parallel region. */
    }
}

TEX_NAMESPACE_END
