#pragma once
#include <random>

#include "Graphics.h"
#include "Utility.h"

#ifndef GRID
#define GRID 2  // 原来 32，先降到 12，内存 < 50 KB
#endif

/* ====== Noise 工具 ====== */
namespace Noise {
/* 简单 2-D value-noise + 线性插值 */
inline float value(int x, int y, uint32_t seed) {
    uint32_t n = x + y * 57 + seed * 131;
    n = (n << 13) ^ n;
    return 1.f - ((n * (n * n * 15731u + 789221u) + 1376312589u) & 0x7fffffff) / 1073741824.f;
}

inline float smooth(float x, float y, uint32_t seed) {
    int ix = int(std::floor(x));
    int iy = int(std::floor(y));
    float fx = x - ix;
    float fy = y - iy;

    float v00 = value(ix, iy, seed);
    float v10 = value(ix + 1, iy, seed);
    float v01 = value(ix, iy + 1, seed);
    float v11 = value(ix + 1, iy + 1, seed);

    float i1 = Lerp(v00, v10, fx);
    float i2 = Lerp(v01, v11, fx);
    return Lerp(i1, i2, fy);
}

/* 2-层 FBM 就够用了 */
inline float fbm(float x, float y, uint32_t seed) {
    float v = 0.f, amp = 0.5f;
    for (int i = 0; i < 2; ++i) {
        v += smooth(x, y, seed + i * 37u) * amp;
        x *= 2.f;
        y *= 2.f;
        amp *= 0.5f;
    }
    return v;
}
}  // namespace Noise

/* ====== 共享常量 ====== */

constexpr float LIGHT_X = 0.35f;  // 光源 ∈[0,1]
constexpr float LIGHT_Y = 0.75f;

Color OCEAN_LIGHT = CLR(102, 176, 199);
Color OCEAN_MID = CLR(71, 97, 124);
Color OCEAN_DARK = CLR(53, 57, 85);
Color LAND_HIGH = CLR(200, 215, 80);
Color LAND_MID = CLR(120, 180, 63);
Color LAND_LOW = CLR(60, 115, 90);
Color LAND_DARK = CLR(40, 60, 80);
Color CLOUD_A = CLR(225, 240, 255);
Color CLOUD_B = CLR(180, 205, 230);

/* ====== 一个网格面片基类 ====== */
class GridPlane final : public Mesh {
public:
    GridPlane(float scale, Color (*colorFn)(float, float, float))
        : Mesh(scale) {
        /* 1. 生成顶点 */
        vertices = new List<Vec3>();
        vertices->reserve((GRID + 1) * (GRID + 1));
        for (int j = 0; j <= GRID; ++j)
            for (int i = 0; i <= GRID; ++i) {
                float u = (float)i / GRID - 0.5f;
                float v = (float)j / GRID - 0.5f;
                vertices->emplace_back(u, v, 0);
            }

        /* 2. 生成索引 + 三角形着色 */
        indices = new List<int>();
        triangles = new List<Triangle>();
        uint32_t seed = randomSeed();

        for (int j = 0; j < GRID; ++j) {
            for (int i = 0; i < GRID; ++i) {
                int idx = j * (GRID + 1) + i;
                int idxR = idx + 1;
                int idxD = idx + GRID + 1;
                int idxDR = idxD + 1;

                addTri(idx, idxD, idxR, seed, colorFn);
                addTri(idxR, idxD, idxDR, seed, colorFn);
            }
        }
        finalize();
    }

private:
    static uint32_t randomSeed() {
        static std::mt19937 rng{12345u};
        return rng();
    }

    void addTri(int a, int b, int c, uint32_t seed,
                Color (*fn)(float, float, float)) {
        /* 三角质心作采样点 */
        Vec3 p = ((*vertices)[a] + (*vertices)[b] + (*vertices)[c]) / 3.f;
        float u = p.x + 0.5f;  // => [0,1]
        float v = p.y + 0.5f;

        /* 距光源距离，用来做晨昏带 */
        float dLight = sqrtf((u - LIGHT_X) * (u - LIGHT_X) + (v - LIGHT_Y) * (v - LIGHT_Y));

        Color col = fn(u, v, dLight);
        Triangle tri((*vertices)[a], (*vertices)[b], (*vertices)[c]);
        tri.color = col;
        indices->push_back(a);
        indices->push_back(b);
        indices->push_back(c);
        // TODO: optimize
        triangles->push_back(tri);
    }
};

/* ====== 层工厂 ====== */
namespace PlanetLayerFactory {

/* --- 1. 基础海洋层 --- */
inline Color baseColorFn(float u, float v, float dL) {
    /* 外圈裁掉：像素着色时会被 Screen-space 裁圆 */
    float d = sqrtf((u - 0.5f) * (u - 0.5f) + (v - 0.5f) * (v - 0.5f));
    if (d > .5f)
        return Color::black;

    /* 晨昏 */
    if (dL > .6f)
        return OCEAN_DARK;
    if (dL > .4f)
        return OCEAN_MID;
    return OCEAN_LIGHT;
}

inline Mesh* createBaseOceanLayer(Transform& parent) {
    auto* m = new GridPlane(1.f, baseColorFn);
    m->SetParent(&parent);
    return m;
}

/* --- 2. 陆地层 --- */
inline Color landColorFn(float u, float v, float dL) {
    float n = Noise::fbm(u * 8.f, v * 8.f, 777u);
    if (n < 0.5f)
        return Color::black;  // 海 → 透明（留给下层）
    /* 依噪声高度给不同绿度，再叠光照 */
    Color c = (n > 0.8f) ? LAND_HIGH : (n > 0.65f) ? LAND_MID
                                   : (n > 0.55f)   ? LAND_LOW
                                                   : LAND_DARK;

    if (dL > .6f)
        c *= 0.35f;
    else if (dL > .4f)
        c *= 0.6f;
    return c;
}

inline Mesh* createLandMassLayer(Transform& parent, float cutoff) {
    (void)cutoff;                                  // 目前固定 0.5，留接口
    auto* m = new GridPlane(1.001f, landColorFn);  // 稍微放大，盖住海面
    m->SetParent(&parent);
    return m;
}

/* --- 3. 云层 --- */
inline Color cloudColorFn(float u, float v, float dL) {
    float n = Noise::fbm(u * 16.f + 3.1f, v * 16.f - 1.7f, 123u);
    if (n < 0.58f)
        return Color::black;  // 透明
    Color c = (n > 0.8f) ? CLOUD_A : CLOUD_B;
    if (dL > .55f)
        c *= 0.7f;
    return c;
}

inline Mesh* createCloudLayer(Transform& parent) {
    auto* m = new GridPlane(1.005f, cloudColorFn);
    m->SetParent(&parent);
    m->ignoreLighting = true;
    return m;
}

}  // namespace PlanetLayerFactory
