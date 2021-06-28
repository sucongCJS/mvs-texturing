#ifndef TEX_DEBUG_HEADER
#define TEX_DEBUG_HEADER

#include <vector>
#include "texturing.h"

TEX_NAMESPACE_BEGIN

/** 用包含不同颜色的图像替换 texture_views 的封装图像 */
void generate_debug_embeddings(std::vector<TextureView> * texture_views);

TEX_NAMESPACE_END

#endif