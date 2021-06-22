#ifndef TEX_TRI_HEADER
#define TEX_TRI_HEADER

#include "math/vector.h"
#include "math/matrix.h"
#include "rect.h"

class Tri {
private:
    math::Vec2f v1;
    math::Vec2f v2;
    math::Vec2f v3;
    
    float detT;  // 
    Rect<float> aabb;

public:
    /** 计算aabb以及预计算重心坐标 */
    Tri(math::Vec2f v1, math::Vec2f v2, math::Vec2f v3);

    /** 判断点x,y是否在重心坐标以内 */
    bool inside(float x, float y) const;

    /** 计算点x,y的重心坐标 */
    math::Vec3f get_barycentric_coords(float x, float y) const;

    /** 返回三角形面积 */
    float get_area(void) const;

    /** 返回三角形的aabb */
    Rect<float> get_aabb(void) const;

};

inline Rect<float> Tri::get_aabb(void) const {
    return aabb;
}

inline math::Vec3f Tri::get_barycentric_coords(float x, float y) const {
    float const alpha = ((v2[1] - v3[1]) * (x - v3[0]) + (v3[0] - v2[0]) * (y - v3[1])) / detT;
    float const beta = ((v3[1] - v1[1]) * (x - v3[0]) + (v1[0] - v3[0]) * (y - v3[1])) / detT;
    float const gamma = 1.0f - alpha - beta;
    return math::Vec3f(alpha, beta, gamma);
}

// alpha, beta, gamma 都要在[0,1]范围内
inline bool Tri::inside(float x, float y) const {
    float const dx = (x - v3[0]);
    float const dy = (y - v3[1]);

    float const alpha = ((v2[1] - v3[1]) * dx + (v3[0] - v2[0]) * dy) / detT;
    if (alpha < 0.0f || alpha > 1.0f)
        return false;

    float const beta = ((v3[1] - v1[1]) * dx + (v1[0] - v3[0]) * dy) / detT;
    if (beta < 0.0f || beta > 1.0f)
        return false;

    if (alpha + beta > 1.0f)
        return false;

    /* else */
    return true;
}

inline float Tri::get_area(void) const {
    math::Vec2f u = v2 - v1;
    math::Vec2f v = v3 - v1;

    return 0.5f * std::abs(u[0] * v[1] - u[1] * v[0]);
}

#endif /* TEX_TRI_HEADER */
