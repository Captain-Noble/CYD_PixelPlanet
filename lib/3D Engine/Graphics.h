#pragma once
#ifndef GRAPHICS_H
#define GRAPHICS_H
#include <Matrix.h>
#include <Utility.h>
#include <math.h>
#include <string.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <unordered_map>
#include <vector>

//------------------------------------------------------------
// 纯软件帧缓冲：修改这里即可改变分辨率
//------------------------------------------------------------
#define DISPLAY_WIDTH_PIX 120
#define DISPLAY_HEIGHT_PIX 120
#define FIELD_OF_VIEW_DEG 20

#define PI 3.1415926

typedef uint16_t Pixel565;  // 5-6-5

inline Pixel565 rgb565(uint8_t r, uint8_t g, uint8_t b) {
    return (r >> 3) << 11 | (g >> 2) << 5 | (b >> 3);
}

struct __attribute__((packed)) PixelRGB {
    uint8_t r, g, b;

    // 让 PixelRGB → Pixel565 可以自动发生
    constexpr operator Pixel565() const noexcept {
        return (Pixel565)((r & 0xF8) << 8) |
               ((g & 0xFC) << 3) |
               (b >> 3);
    }
};

// 行优先的 2-D 显存；所有绘图 API 都写到这里
Pixel565 display[DISPLAY_HEIGHT_PIX][DISPLAY_WIDTH_PIX];

//------------------------------------------------------------
// 小工具：安全写像素 & 线段/三角形栅格化
//------------------------------------------------------------

#pragma region buffer_functions

void drawLineSW(int x0, int y0, int x1, int y1, const PixelRGB& c);  // 见下方实现
void fillTriSW(Vec3 p1, Vec3 p2, Vec3 p3, const PixelRGB& c);        // 见下方实现
inline void putPixel(int x, int y, const PixelRGB& c) {
    if (x < 0 || x >= DISPLAY_WIDTH_PIX || y < 0 || y >= DISPLAY_HEIGHT_PIX)
        return;
    display[y][x] = c;
}

// Bresenham 直线
inline void drawLineSW(int x0, int y0, int x1, int y1, const PixelRGB& c) {
    bool steep = std::abs(y1 - y0) > std::abs(x1 - x0);
    if (steep) {
        std::swap(x0, y0);
        std::swap(x1, y1);
    }
    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0, dy = std::abs(y1 - y0), err = dx / 2, ystep = (y0 < y1 ? 1 : -1);
    for (int x = x0, y = y0; x <= x1; ++x) {
        steep ? putPixel(y, x, c) : putPixel(x, y, c);
        err -= dy;
        if (err < 0) {
            y += ystep;
            err += dx;
        }
    }
}

// 扫描线三角形填充 (面积符号法)
inline void fillTriSW(Vec3 p1, Vec3 p2, Vec3 p3, const PixelRGB& c) {
    auto edge = [](const Vec3& a, const Vec3& b, const Vec3& p) {
        return (p.x - a.x) * (b.y - a.y) - (p.y - a.y) * (b.x - a.x);
    };
    int minx = std::max(0, int(std::floor(std::min({p1.x, p2.x, p3.x}))));
    int maxx = std::min(DISPLAY_WIDTH_PIX - 1, int(std::ceil(std::max({p1.x, p2.x, p3.x}))));
    int miny = std::max(0, int(std::floor(std::min({p1.y, p2.y, p3.y}))));
    int maxy = std::min(DISPLAY_HEIGHT_PIX - 1, int(std::ceil(std::max({p1.y, p2.y, p3.y}))));
    for (int y = miny; y <= maxy; ++y)
        for (int x = minx; x <= maxx; ++x) {
            Vec3 p(x + 0.5f, y + 0.5f, 0);
            float w0 = edge(p2, p3, p), w1 = edge(p3, p1, p), w2 = edge(p1, p2, p);
            if ((w0 >= 0 && w1 >= 0 && w2 >= 0) || (w0 <= 0 && w1 <= 0 && w2 <= 0))
                putPixel(x, y, c);
        }
}
#pragma endregion

using namespace std;

class Graphics;
struct Plane;
struct Point;
struct Line;
struct Triangle;
class Transform;
class Mesh;
class Camera;

struct Direction {
    static Vec3 forward;
    static Vec3 back;
    static Vec3 right;
    static Vec3 left;
    static Vec3 up;
    static Vec3 down;
};
Vec3 Direction::forward = Vec3(0, 0, -1);
Vec3 Direction::back = Vec3(0, 0, 1);
Vec3 Direction::right = Vec3(1, 0, 0);
Vec3 Direction::left = Vec3(-1, 0, 0);
Vec3 Direction::up = Vec3(0, 1, 0);
Vec3 Direction::down = Vec3(0, -1, 0);

struct Color {
    float r;
    float g;
    float b;
    float a;

    static Color black;
    static Color white;
    static Color gray;
    static Color red;
    static Color green;
    static Color blue;
    static Color pink;
    static Color purple;
    static Color yellow;
    static Color turquoise;
    static Color orange;

    Color() {
        this->r = 0.0;
        this->g = 0.0;
        this->b = 0.0;
        this->a = 1;
    }

    Color(float r, float g, float b, float a = 1.0) {
        this->r = r;
        this->g = g;
        this->b = b;
        this->a = a;
    }

    Color(Vec4 vec) {
        this->r = vec.x;
        this->g = vec.y;
        this->b = vec.z;
        this->a = vec.w;
    }

    static Color Random() {
        Color c = Color(Clamp(rand(), 0, 255), Clamp(rand(), 0, 255), Clamp(rand(), 0, 255));
        return c;
    }

    Color operator+(const Color& other) {
        Color color;
        color.r = Clamp(this->r + other.r, 0, 255);
        color.g = Clamp(this->g + other.g, 0, 255);
        color.b = Clamp(this->b + other.b, 0, 255);
        return color;
    }

    Color operator-(const Color& other) {
        Color color;
        color.r = Clamp(this->r - other.r, 0, 255);
        color.g = Clamp(this->g - other.g, 0, 255);
        color.b = Clamp(this->b - other.b, 0, 255);
        return color;
    }

    Color operator+(const Vec3& v3) {
        Color color;
        color.r = Clamp(this->r + v3.x, 0, 255);
        color.g = Clamp(this->g + v3.y, 0, 255);
        color.b = Clamp(this->b + v3.z, 0, 255);
        return color;
    }

    Color operator-(const Vec3& v3) {
        Color color;
        color.r = Clamp(this->r - v3.x, 0, 255);
        color.g = Clamp(this->g - v3.y, 0, 255);
        color.b = Clamp(this->b - v3.z, 0, 255);
        return color;
    }

    Color operator*(const float& scalar) {
        Color color;
        color.r = Clamp(this->r * scalar, 0, 255);
        color.g = Clamp(this->g * scalar, 0, 255);
        color.b = Clamp(this->b * scalar, 0, 255);
        return color;
    }

    Color& operator+=(const Color& other) {
        this->r = Clamp(this->r + other.r, 0, 255);
        this->g = Clamp(this->g + other.g, 0, 255);
        this->b = Clamp(this->b + other.b, 0, 255);
        return *this;
    }

    Color& operator+=(const Vec3& v3) {
        this->r = Clamp(this->r + v3.x, 0, 255);
        this->g = Clamp(this->g + v3.y, 0, 255);
        this->b = Clamp(this->b + v3.z, 0, 255);
        return *this;
    }

    Color& operator-=(const Color& other) {
        this->r = Clamp(this->r - other.r, 0, 255);
        this->g = Clamp(this->g - other.g, 0, 255);
        this->b = Clamp(this->b - other.b, 0, 255);
        return *this;
    }

    Color& operator-=(const Vec3& v3) {
        this->r = Clamp(this->r - v3.x, 0, 255);
        this->g = Clamp(this->g - v3.y, 0, 255);
        this->b = Clamp(this->b - v3.z, 0, 255);
        return *this;
    }

    Color& operator*=(const float scalar) {
        this->r = Clamp(this->r * scalar, 0, 255);
        this->g = Clamp(this->g * scalar, 0, 255);
        this->b = Clamp(this->b * scalar, 0, 255);
        return *this;
    }

    Color& operator/=(const float divisor) {
        if (divisor != 0.0) {
            this->r = Clamp(this->r / divisor, 0, 255);
            this->g = Clamp(this->g / divisor, 0, 255);
            this->b = Clamp(this->b / divisor, 0, 255);
            return *this;
        }
    }
    operator Vec3();
};

#define CLR(r, g, b) Color(float(r), float(g), float(b))

Color Color::black = CLR(0, 0, 0);
Color Color::white = CLR(255, 255, 255);
Color Color::gray = CLR(128, 128, 128);
Color Color::red = CLR(255, 0, 0);
Color Color::green = CLR(0, 255, 0);
Color Color::blue = CLR(0, 0, 255);
Color Color::pink = CLR(255, 0, 255);
Color Color::purple = CLR(128, 0, 128);
Color Color::yellow = CLR(255, 255, 0);
Color Color::turquoise = CLR(0, 255, 255);
Color Color::orange = CLR(255, 158, 0);

struct Material {
    std::string name;
    Color color{Color::white};
    Material(std::string id = "", Color c = Color::white) : name(std::move(id)), color(c) {}
};

class Transform {
public:
    Vec3 localScale{1, 1, 1}, localPosition{0, 0, 0};
    Matrix3x3 localRotation{Matrix3x3::identity};
    Transform *root{this}, *parent{nullptr};

    Transform(const float& scale = 1, const Vec3& position = Vec3(0, 0, 0), const Vec3& rotationEuler = Vec3(0, 0, 0)) {
        this->localScale.x = scale;
        this->localScale.y = scale;
        this->localScale.z = scale;
        this->localPosition = position;
        this->localRotation = YPR(rotationEuler.x, rotationEuler.y, rotationEuler.z);
        this->root = this;
    }

    Transform(const Vec3& scale, const Vec3& position = Vec3(0, 0, 0), const Matrix3x3& rotation = Matrix3x3::identity) {
        this->localScale = scale;
        this->localPosition = position;
        this->localRotation = rotation;
        this->root = this;
    }

    Vec3 Forward() { return Rotation() * Direction::forward; }
    Vec3 Back() { return Rotation() * Direction::back; }
    Vec3 Right() { return Rotation() * Direction::right; }
    Vec3 Left() { return Rotation() * Direction::left; }
    Vec3 Up() { return Rotation() * Direction::up; }
    Vec3 Down() { return Rotation() * Direction::down; }

    Matrix3x3 Rotation();

    Vec3 Position();

    Vec3 Scale();

    void SetParent(Transform* newParent, bool changeOfBasisTransition = true);

    Matrix4x4 LocalScale4x4();

    Matrix4x4 LocalScale4x4Inverse();

    Matrix4x4 LocalRotation4x4();

    Matrix4x4 LocalTranslation4x4();

    Matrix4x4 LocalTranslation4x4Inverse();

    // 1:Scale, 2:Rotate, 3:Translate
    Matrix4x4 TRS();

    // S^-1 * R^-1 * T^-1
    Matrix4x4 TRSInverse();

    // 1:Rotate, 2:Translate
    Matrix4x4 TR();

    // R^-1 * T^-1
    Matrix4x4 TRInverse();
};

/*───────────────────────────  Mesh.h / Graphics.h  ───────────────────────────*/
class Mesh : public Transform, public ManagedObjectPool<Mesh> {
public:
    /*----------- 成员 -----------*/
    static int worldTriangleDrawCount;
    List<Vec3>* vertices{nullptr};
    List<int>* indices{nullptr};
    List<Triangle>* triangles{nullptr};

    Color color = Color::white;
    bool ignoreLighting = false;
    bool forceWireFrame = false;

    /*----------- 构造 / 析构 -----------*/
    Mesh(float scale = 1,
         const Vec3& pos = Vec3(0, 0, 0),
         const Vec3& rotEuler = Vec3(0, 0, 0))
        : Transform(scale, pos, rotEuler),
          ManagedObjectPool<Mesh>(this) {
        /* 注意：此时派生类尚未分配顶点索引，因此
           不要在这里调用 MapVertsToTriangles() ！*/
    }

    Mesh(const Vec3& scale,
         const Vec3& pos = Vec3(0, 0, 0),
         const Matrix3x3& rot = Matrix3x3::identity)
        : Transform(scale, pos, rot),
          ManagedObjectPool<Mesh>(this) {}

    virtual ~Mesh() {
        delete vertices;
        delete indices;
        delete triangles;
    }

    /*----------- 延后初始化 -----------*/
    void finalize() /* ← 新增 */
    {
        if (!vertices || !indices)
            return;

        if (!triangles)
            triangles = new List<Triangle>(indices->size() / 3);

        MapVertsToTriangles();  // 现在才安全
        SetColor(color);
    }

    /*----------- 接口 -----------*/
    bool SetVisibility(bool visible);
    void SetColor(Color&& c) { SetColor(c); }
    void SetColor(Color& c);

    virtual List<Triangle>* MapVertsToTriangles();  // 下面给出实现
    List<Vec3> WorldVertices();
    void TransformTriangles();
};
/*───────────────────────────  END Mesh  ───────────────────────────*/

List<Point>* pointBuffer = new List<Point>();
List<Line>* lineBuffer = new List<Line>();
List<Triangle>* triBuffer = new List<Triangle>();

static float worldScale = 1;
constexpr int screenWidth = DISPLAY_WIDTH_PIX;    // ← 统一到 100
constexpr int screenHeight = DISPLAY_HEIGHT_PIX;  // ← 统一到 100

static inline void ToScreenCoordinates(Vec3& v) {
    static float m[4][4]{
        {0.5f * DISPLAY_WIDTH_PIX, 0, 0, 0.5f * DISPLAY_WIDTH_PIX},
        {0, -0.5f * DISPLAY_HEIGHT_PIX, 0, 0.5f * DISPLAY_HEIGHT_PIX},
        {0, 0, 1, 0},
        {0, 0, 0, 0}};
    v = m * v;
}

//-----------GRAPHICS---------------

class Graphics {
private:
    static PixelRGB _pen;

public:
    // 调试/状态位与原版保持一致
    static bool frustumCulling, backFaceCulling, invertNormals, debugNormals, debugVertices,
        debugAxes, debugBoxCollisions, debugSphereCollisions, debugPlaneCollisions,
        perspective, fillTriangles, displayWireFrames, lighting, vfx, matrixMode;

    static void SetDrawColor(int r, int g, int b) { _pen = {uint8_t(r), uint8_t(g), uint8_t(b)}; }
    static void SetDrawColor(Color c) { SetDrawColor(int(c.r), int(c.g), int(c.b)); }
    static void SetLineWidth(int) {}
    static void SetPointSize(int) {}

    static void DrawPoint(Vec3 p) {
        ToScreenCoordinates(p);
        putPixel(int(p.x + 0.5f), int(p.y + 0.5f), _pen);
    }
    static void DrawLine(Vec3 a, Vec3 b) {
        ToScreenCoordinates(a);
        ToScreenCoordinates(b);
        drawLineSW(int(a.x), int(a.y), int(b.x), int(b.y), _pen);
    }
    static void DrawTriangle(Vec3 p1, Vec3 p2, Vec3 p3) {
        DrawLine(p1, p2);
        DrawLine(p2, p3);
        DrawLine(p3, p1);
    }
    static void DrawTriangleFilled(Vec3 p1, Vec3 p2, Vec3 p3) {
        ToScreenCoordinates(p1);
        ToScreenCoordinates(p2);
        ToScreenCoordinates(p3);
        fillTriSW(p1, p2, p3, _pen);
    }
};

PixelRGB Graphics::_pen{255, 255, 255};
bool Graphics::frustumCulling = true, Graphics::backFaceCulling = true,
     Graphics::invertNormals = false, Graphics::debugNormals = false,
     Graphics::debugVertices = false, Graphics::debugAxes = false,
     Graphics::debugBoxCollisions = false, Graphics::debugSphereCollisions = false,
     Graphics::debugPlaneCollisions = false, Graphics::perspective = true,
     Graphics::fillTriangles = true, Graphics::displayWireFrames = false,
     Graphics::lighting = true, Graphics::vfx = false, Graphics::matrixMode = false;

Vec3 lightSource = .25 * Direction::up + Direction::back * .5;
float nearClippingPlane = -0.1;
float farClippingPlane = -100000.0;
float fieldOfViewDeg = FIELD_OF_VIEW_DEG;
float fov = ToRad(fieldOfViewDeg);
float aspect = (float)screenHeight / (float)screenWidth;

// Perspective Projection Matrix
float persp[4][4] = {
    {aspect * 1 / tan(fov / 2), 0, 0, 0},
    {0, 1 / tan(fov / 2), 0, 0},
    {0, 0, 1, 0},
    {0, 0, -1, 0}};

// Perspective Projection Matrix
float weakPersp[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, -1, 0}};

// Orthographic Projection Matrix
float ortho[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 1}};

float undoAspect[4][4] = {
    {1 / aspect, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 0}};

Matrix4x4 perspectiveProjectionMatrix = persp;
Matrix4x4 weakPerspectiveProjectionMatrix = weakPersp;
Matrix4x4 orthographicProjectionMatrix = ortho;
Matrix4x4 undoAspectRatio = undoAspect;

void FOV(int deg) {
    fieldOfViewDeg = deg;
    fov = ToRad(deg);
    float newPerspectiveProjectionMatrix[4][4] = {
        {aspect * 1 / tan(fov / 2), 0, 0, 0},
        {0, 1 / tan(fov / 2), 0, 0},
        {0, 0, 1, 0},
        {0, 0, -1, 0}};

    perspectiveProjectionMatrix = newPerspectiveProjectionMatrix;
}

Matrix4x4 ProjectionMatrix() {
    if (Graphics::perspective) {
        return perspectiveProjectionMatrix;
    } else {
        return orthographicProjectionMatrix;
    }
}

Matrix4x4 worldToViewMatrix;
Matrix4x4 projectionMatrix;

struct Point {
    Vec3 position;
    Color color;
    int size;

    Point(Vec3 position, Color color = Color::white, int size = 2) {
        this->position = position;
        this->color = color;
        this->size = size;
    }

    void Draw() {
        Graphics::SetPointSize(size);
        Graphics::SetDrawColor(color);
        Graphics::DrawPoint(position);
    }

    static void AddPoint(Point point) {
        pointBuffer->emplace_back(point);
    }

    static void AddWorldPoint(Point point) {
        auto matrix = ProjectionMatrix() * worldToViewMatrix;
        point.position = matrix * point.position;
        pointBuffer->emplace_back(point);
    }
};

struct Line {
    Vec3 from;
    Vec3 to;
    Color color;
    int width;

    Line(Vec3 from, Vec3 to, Color color = Color::white, int width = 2) {
        this->from = from;
        this->to = to;
        this->color = color;
        this->width = width;
    }

    void Draw() {
        Graphics::SetLineWidth(width);
        Graphics::SetDrawColor(color.r, color.g, color.b);
        Graphics::DrawLine(from, to);
    }

    static void AddLine(Line line) {
        lineBuffer->emplace_back(line);
    }

    static void AddWorldLine(Line line) {
        auto matrix = ProjectionMatrix() * worldToViewMatrix;
        line.from = matrix * line.from;
        line.to = matrix * line.to;
        lineBuffer->emplace_back(line);
    }
};

struct Triangle : Plane {
    Vec4 centroid = Vec4();
    Color color = Color::white;
    bool forceWireFrame = false;
    Mesh* mesh = nullptr;

    Triangle() : Plane() {
        centroid = Vec4();
        color = Color::white;
        mesh = nullptr;
    }
    Triangle(Vec3 p1, Vec3 p2, Vec3 p3, Mesh* owner = nullptr) : Plane(p1, p2, p3) {
        color = Color::white;
        centroid = Centroid();
        mesh = owner;
    }

    Vec4 Centroid() {
        centroid = Vec4(
            (verts[0].x + verts[1].x + verts[2].x) / 3.0,
            (verts[0].y + verts[1].y + verts[2].y) / 3.0,
            (verts[0].z + verts[1].z + verts[2].z) / 3.0,
            centroid.w);

        return centroid;
    }

    void Draw() {
        Vec4 p1 = verts[0];
        Vec4 p2 = verts[1];
        Vec4 p3 = verts[2];

        if (Graphics::fillTriangles == false) {
            Graphics::displayWireFrames = true;
        }

        // glColor3ub(255, 255, 255);
        if (Graphics::matrixMode) {
            Graphics::SetDrawColor(0, 255, 0);
        }

        if (Graphics::fillTriangles) {
            Graphics::SetDrawColor(color.r, color.g, color.b);
            Graphics::DrawTriangleFilled(p1, p2, p3);
        }

        bool drawWireFrame = Graphics::displayWireFrames || forceWireFrame || (mesh != nullptr ? mesh->forceWireFrame : false);
        if (drawWireFrame) {
            if (Graphics::fillTriangles) {
                float c = Clamp(1.0 / (0.000001 + (color.r + color.g + color.b) / 3), 0, 255);
                Graphics::SetDrawColor(c, c, c);
            }
            Graphics::DrawLine(p1, p2);
            Graphics::DrawLine(p2, p3);
            Graphics::DrawLine(p3, p1);
        }

        Graphics::SetDrawColor(255, 255, 255);
    }
};

//-------------------------------TRANSFORM---------------------------------------------

Vec3 ExtractPosition(const Matrix4x4& trs) {
    return {trs.m[0][3], trs.m[1][3], trs.m[2][3]};
}

Vec3 ExtractScale(const Matrix4x4& trs) {
    Vec3 c1 = {trs.m[0][0], trs.m[1][0], trs.m[2][0]};
    Vec3 c2 = {trs.m[0][1], trs.m[1][1], trs.m[2][1]};
    Vec3 c3 = {trs.m[0][2], trs.m[1][2], trs.m[2][2]};
    return {c1.Magnitude(), c2.Magnitude(), c3.Magnitude()};
}

// Pass scale parameter if already known since finding the rotation matrix for a TRS matrix
// requires the scale vector, which finding can be expensive.
Matrix3x3 ExtractRotation(const Matrix4x4 trs, Vec3* scale = NULL) {
    static Vec3 v = Vec3::zero;
    Vec3* s = &v;
    if (scale == NULL) {
        *s = ExtractScale(trs);
    } else {
        s = scale;
    }

    float rot[3][3] = {
        {trs.m[0][0] / s->x, trs.m[1][0] / s->y, trs.m[2][0] / s->z},
        {trs.m[0][1] / s->x, trs.m[1][1] / s->y, trs.m[2][1] / s->z},
        {trs.m[0][2] / s->x, trs.m[1][2] / s->y, trs.m[2][2] / s->z},
    };

    return rot;
}

struct TRSInfo {
    Vec3 scale;
    Vec3 position;
    Matrix3x3 rotation;

    TRSInfo(const Vec3& s, const Vec3& pos, const Matrix3x3& rot) {
        scale = s;
        position = pos;
        rotation = rot;
    }
};

TRSInfo ExtractTRS(Matrix4x4& trs) {
    Vec3 scale = ExtractScale(trs);
    return TRSInfo(scale, ExtractPosition(trs), ExtractRotation(trs, &scale));
}

Matrix3x3 Transform::Rotation() {
    if (parent) {  // TRS already checks for parent but checks again here because cheaper to return local position than creating the matrix
        return ExtractRotation(TRS());
    }
    return localRotation;
}

Vec3 Transform::Position() {
    if (parent) {  // TRS already checks for parent but checks again here because cheaper to return local position than creating the matrix
        return ExtractPosition(TRS());
    }

    return localPosition;
}

Vec3 Transform::Scale() {
    if (parent) {
        return ExtractScale(TRS());
    }
    return localScale;
}

void Transform::SetParent(Transform* newParent, bool changeOfBasisTransition) {
    static Matrix4x4 T;
    if (newParent == nullptr) {
        if (changeOfBasisTransition) {
            // Results in seemless unparenting. Assigns (possibly parented) global values to the local values.
            // This way of unparenting will maintain the previous parented position and rotation but in the new reference frame.
            // Also nothing changes if never parented.
            // if already had a parent

            T = this->TR();

            // this->scale = Scale();// Vec3(1, 1, 1);
            float rot[3][3] = {
                {T.m[0][0], T.m[0][1], T.m[0][2]},
                {T.m[1][0], T.m[1][1], T.m[1][2]},
                {T.m[2][0], T.m[2][1], T.m[2][2]}};
            this->localRotation = rot;
            this->localPosition = Vec3(T.m[0][3], T.m[1][3], T.m[2][3]);
        }

        this->parent = NULL;
        this->root = this;
    } else if (newParent != this) {
        if (changeOfBasisTransition) {
            // If you immediatley parent a transform, everything is calculated relative to its immediate parent's reference frame.
            // This would result in a transformation if the original coordinates didn't change.
            // Instead what we want is a change of basis.
            T = newParent->TRInverse() * this->TR();

            float rot[3][3] = {
                {T.m[0][0], T.m[0][1], T.m[0][2]},
                {T.m[1][0], T.m[1][1], T.m[1][2]},
                {T.m[2][0], T.m[2][1], T.m[2][2]}};
            this->localScale = this->Scale();
            this->localRotation = rot;
            this->localPosition = Vec3(T.m[0][3], T.m[1][3], T.m[2][3]);
        }

        this->parent = newParent;
        this->root = newParent->root;
    }
}

Matrix4x4 Transform::LocalScale4x4() {
    float matrix[4][4] =
        {
            {this->localScale.x, 0, 0, 0},
            {0, this->localScale.y, 0, 0},
            {0, 0, this->localScale.z, 0},
            {0, 0, 0, 1}};

    return Matrix4x4(matrix);
}

Matrix4x4 Transform::LocalScale4x4Inverse() {
    float inverse[4][4] =
        {
            {1.0 / this->localScale.x, 0, 0, 0},
            {0, 1.0 / this->localScale.y, 0, 0},
            {0, 0, 1.0 / this->localScale.z, 0},
            {0, 0, 0, 1}};

    return Matrix4x4(inverse);
}

Matrix4x4 Transform::LocalRotation4x4() {
    float matrix[4][4] =
        {
            {this->localRotation.m[0][0], this->localRotation.m[0][1], this->localRotation.m[0][2], 0},
            {this->localRotation.m[1][0], this->localRotation.m[1][1], this->localRotation.m[1][2], 0},
            {this->localRotation.m[2][0], this->localRotation.m[2][1], this->localRotation.m[2][2], 0},
            {0, 0, 0, 1}};

    return Matrix4x4(matrix);
}

Matrix4x4 Transform::LocalTranslation4x4() {
    float matrix[4][4] =
        {
            {1, 0, 0, this->localPosition.x},
            {0, 1, 0, this->localPosition.y},
            {0, 0, 1, this->localPosition.z},
            {0, 0, 0, 1}};

    return Matrix4x4(matrix);
}

Matrix4x4 Transform::LocalTranslation4x4Inverse() {
    float matrix[4][4] =
        {
            {1, 0, 0, -this->localPosition.x},
            {0, 1, 0, -this->localPosition.y},
            {0, 0, 1, -this->localPosition.z},
            {0, 0, 0, 1}};

    return Matrix4x4(matrix);
}

// 1:Scale, 2:Rotate, 3:Translate
Matrix4x4 Transform::TRS() {
    float trs[4][4] =
        {
            {this->localRotation.m[0][0] * localScale.x, this->localRotation.m[0][1] * localScale.y, this->localRotation.m[0][2] * localScale.z, localPosition.x},
            {this->localRotation.m[1][0] * localScale.x, this->localRotation.m[1][1] * localScale.y, this->localRotation.m[1][2] * localScale.z, localPosition.y},
            {this->localRotation.m[2][0] * localScale.x, this->localRotation.m[2][1] * localScale.y, this->localRotation.m[2][2] * localScale.z, localPosition.z},
            {0, 0, 0, 1}};

    if (parent) {
        return parent->TRS() * trs;
    }

    return trs;
}

// S^-1 * R^-1 * T^-1
Matrix4x4 Transform::TRSInverse() {
    if (parent) {
        return LocalScale4x4Inverse() * Matrix4x4::Transpose(LocalRotation4x4()) * LocalTranslation4x4Inverse() * parent->TRSInverse();
    }

    return LocalScale4x4Inverse() * Matrix4x4::Transpose(LocalRotation4x4()) * LocalTranslation4x4Inverse();
}

// 1:Rotate, 2:Translate
Matrix4x4 Transform::TR() {
    float tr[4][4] =
        {
            {this->localRotation.m[0][0], this->localRotation.m[0][1], this->localRotation.m[0][2], localPosition.x},
            {this->localRotation.m[1][0], this->localRotation.m[1][1], this->localRotation.m[1][2], localPosition.y},
            {this->localRotation.m[2][0], this->localRotation.m[2][1], this->localRotation.m[2][2], localPosition.z},
            {0, 0, 0, 1}};
    if (parent) {
        return parent->TR() * tr;
    }

    return tr;
}

// R^-1 * T^-1
Matrix4x4 Transform::TRInverse() {
    if (parent) {
        return Matrix4x4::Transpose(LocalRotation4x4()) * LocalTranslation4x4Inverse() * parent->TRInverse();
    }

    return Matrix4x4::Transpose(LocalRotation4x4()) * LocalTranslation4x4Inverse();
}

//-----------------------------CAMERA-------------------------------------------------
struct CameraSettings {
    static bool outsiderViewPerspective;
    static bool displayReticle;
};
bool CameraSettings::outsiderViewPerspective = false;
bool CameraSettings::displayReticle = true;

class Camera : public Transform {
public:
    static Camera* main;
    static Camera* projector;
    static List<Camera*> cameras;
    static int cameraCount;

    std::string name;

    Camera(const Vec3& position = Vec3(0, 0, 0), const Vec3& rotationEuler = Vec3(0, 0, 0))
        : Transform(1, position, rotationEuler) {
        cameras.emplace(cameras.begin() + cameraCount++, this);
        name = "Camera " + cameraCount;
    }
};
List<Camera*> Camera::cameras = List<Camera*>();
int Camera::cameraCount = 0;
Camera* Camera::projector = new Camera();
Camera* camera1 = new Camera();
Camera* camera2 = new Camera(Vec3(0, 50, 0), Vec3(-90 * PI / 180, 0, 0));
Camera* Camera::main = camera1;

//---------------------------------MESH---------------------------------------------

bool Mesh::SetVisibility(bool visible) {
    if (visible) {
        ManagedObjectPool<Mesh>::addToPool(this);
        return true;
    } else {
        ManagedObjectPool<Mesh>::removeFromPool(this);
        return false;
    }
}
// void Mesh::SetColor(Color&& c) {
//     SetColor(c);
// }

void Mesh::SetColor(Color& c) {
    if (triangles) {
        for (int i = 0; i < triangles->size(); i++) {
            (*triangles)[i].color = c;
        }
    }
}

inline List<Triangle>* Mesh::MapVertsToTriangles() {
    if (!vertices || !indices)  // ← 保护
        return triangles;

    size_t triCnt = indices->size() / 3;
    if (!triangles || triangles->size() < triCnt)
        triangles = new List<Triangle>(triCnt);

    size_t t = 0;
    for (size_t i = 0; i < indices->size(); i += 3, ++t) {
        int p1 = (*indices)[i];
        int p2 = (*indices)[i + 1];
        int p3 = (*indices)[i + 2];

        (*triangles)[t].verts[0] = (*vertices)[p1];
        (*triangles)[t].verts[1] = (*vertices)[p2];
        (*triangles)[t].verts[2] = (*vertices)[p3];
    }
    return triangles;
}

// Convert to world coordinates
List<Vec3> Mesh::WorldVertices() {
    List<Vec3> verts = *vertices;
    Matrix4x4 matrix = TRS();
    for (size_t i = 0; i < verts.size(); i++) {
        verts[i] = (Vec3)(matrix * ((Vec4)verts[i]));
    }

    return verts;
}

void Mesh::TransformTriangles() {
    // Scale/Distance ratio culling
    /* bool tooSmallToSee = scale.SqrMagnitude() / (position - Camera::main->position).SqrMagnitude() < 0.000000125;
    if (tooSmallToSee) {
        return;
    }*/

    Matrix4x4 modelToWorldMatrix = this->TRS();

    // Transform Triangles
    List<Triangle>* tris = MapVertsToTriangles();

    for (int i = 0; i < tris->size(); i++) {
        Triangle tri = (*tris)[i];
        tri.mesh = this;
        Triangle worldSpaceTri = tri;
        Triangle camSpaceTri = tri;
        Triangle projectedTri = tri;
        for (int j = 0; j < 3; j++) {
            // Homogeneous coords (x, y, z, w=1)
            Vec4 vert = tri.verts[j];

            // =================== WORLD SPACE ===================
            // Transform local coords to world-space coords.
            Vec4 worldPoint = modelToWorldMatrix * vert;
            worldSpaceTri.verts[j] = worldPoint;

            // ================ VIEW/CAM/EYE SPACE ================
            // Transform world coordinates to view coordinates.
            Vec4 cameraSpacePoint = worldToViewMatrix * worldPoint;
            camSpaceTri.verts[j] = cameraSpacePoint;

            // ================ SCREEN SPACE ==================
            // Project to screen space (image space)
            Vec4 projectedPoint = projectionMatrix * cameraSpacePoint;
            projectedTri.verts[j] = projectedPoint;
        };

        //------------------- Normal/Frustum Culling (view space)------------------------
        Vec3 p1_c = camSpaceTri.verts[0];
        Vec3 p2_c = camSpaceTri.verts[1];
        Vec3 p3_c = camSpaceTri.verts[2];
        camSpaceTri.Centroid();

        if (Graphics::frustumCulling) {
            bool tooCloseToCamera = (p1_c.z >= nearClippingPlane || p2_c.z >= nearClippingPlane || p3_c.z >= nearClippingPlane || camSpaceTri.centroid.z >= nearClippingPlane);
            if (tooCloseToCamera) {
                continue;
            }

            bool tooFarFromCamera = (p1_c.z <= farClippingPlane || p2_c.z <= farClippingPlane || p3_c.z <= farClippingPlane || camSpaceTri.centroid.z <= farClippingPlane);
            if (tooFarFromCamera) {
                continue;
            }
            /*
            bool behindCamera = DotProduct((Vec3)camSpaceTri.centroid, Direction::forward) <= 0.0;
            if (behindCamera) {
                continue; // Skip triangle if it's out of cam view.
            }*/

            Range xRange = ProjectVertsOntoAxis(projectedTri.verts, 3, Direction::right);
            Range yRange = ProjectVertsOntoAxis(projectedTri.verts, 3, Direction::up);
            if ((xRange.min > 1 || yRange.min > 1) || (xRange.max < -1 || yRange.max < -1)) {
                continue;
            }
        }

        // Calculate triangle suface Normal
        camSpaceTri.Normal();  // camSpaceTri.normal = worldToViewMatrix * modelToWorldMatrix * Vec4(camSpaceTri.normal, 0);

        if (Graphics::invertNormals) {
            camSpaceTri.normal = ((Vec3)camSpaceTri.normal) * -1.0;
        }

        // Back-face Culling - Checks if the triangles backside is facing the camera.
        // Condition makes this setting optional when drawing wireframes alone, but will force culling if triangles are filled.
        if (Graphics::backFaceCulling || Graphics::fillTriangles) {
            Vec3 posRelativeToCam = camSpaceTri.centroid;  // Since camera is (0,0,0) in view space, the displacement vector from camera to centroid IS the centroid itself.
            bool faceInvisibleToCamera = DotProduct(posRelativeToCam, (Vec3)camSpaceTri.normal) >= 0;
            if (faceInvisibleToCamera) {
                continue;  // Skip triangle if it's out of cam view or it's part of the other side of the mesh.
            }
        }

        //------------------------ Lighting (world space)------------------------

        if (Graphics::lighting && Graphics::fillTriangles) {
            if (!ignoreLighting) {
                float amountFacingLight = DotProduct((Vec3)worldSpaceTri.Normal(), lightSource);
                Color colorLit = projectedTri.color * Clamp(amountFacingLight, 0.15, 1);
                projectedTri.color = colorLit;
            }
        }

        if (Graphics::vfx) {
            Vec3 screenLeftSide = Vec3(-1, 0, 0);
            Vec3 screenRightSide = Vec3(1, 0, 0);
            Range range = ProjectVertsOntoAxis(projectedTri.verts, 3, screenRightSide);
            bool leftHalfScreenX = range.min > 0 && range.max > 0;

            if (leftHalfScreenX) {
                projectedTri.color = Color(0, 0, 255);  // std::cout << "Inside" << std::endl;
            } else {
                projectedTri.color = Color::red;
            }
        }
        // ---------- Debugging -----------
        if (Graphics::debugNormals) {
            //---------Draw point at centroid and a line from centroid to normal (view space & projected space)-----------
            Vec2 centroidToNormal_p = projectionMatrix * ((Vec3)camSpaceTri.centroid + camSpaceTri.normal);
            Vec2 centroid_p = projectionMatrix * camSpaceTri.centroid;
            Point::AddPoint(Point(centroid_p));
            Line::AddLine(Line(centroid_p, centroidToNormal_p));
        }

        projectedTri.centroid = projectionMatrix * camSpaceTri.centroid;

        // Nested Projection or Double Projection
        if (CameraSettings::outsiderViewPerspective) {
            Matrix4x4 nestedProjectionMatrix = weakPerspectiveProjectionMatrix * Camera::projector->TRInverse();

            for (size_t k = 0; k < 3; k++) {
                projectedTri.verts[k] = nestedProjectionMatrix * projectedTri.verts[k];
            }
        }

        // Add projected tri
        triBuffer->emplace_back(projectedTri);
    }
}
// List<Mesh*> Mesh::objects = List<Mesh*>(1000);
// int Mesh::meshCount = 0;
int Mesh::worldTriangleDrawCount = 0;

//------------------------------------CUBE MESH------------------------------------------
/*───────────────────────────  CubeMesh  ───────────────────────────*/
class CubeMesh : public Mesh {
public:
    CubeMesh(float scale = 1,
             const Vec3& pos = Vec3(0, 0, 0),
             const Vec3& rot = Vec3(0, 0, 0))
        : Mesh(scale, pos, rot) {
        /* 1. 顶点 */
        vertices = new List<Vec3>{
            // south
            {-0.5, -0.5, 0.5},
            {-0.5, 0.5, 0.5},
            {0.5, 0.5, 0.5},
            {0.5, -0.5, 0.5},
            // north
            {-0.5, -0.5, -0.5},
            {-0.5, 0.5, -0.5},
            {0.5, 0.5, -0.5},
            {0.5, -0.5, -0.5}};

        /* 2. 三角形索引 */
        indices = new List<int>{
            0, 1, 2, 0, 2, 3,  // South
            7, 6, 5, 7, 5, 4,  // North
            3, 2, 6, 3, 6, 7,  // Right
            4, 5, 1, 4, 1, 0,  // Left
            1, 5, 6, 1, 6, 2,  // Top
            3, 7, 4, 3, 4, 0   // Bottom
        };

        /* 3. 创建三角形容器并映射 */
        triangles = new List<Triangle>(indices->size() / 3);
        finalize();  // ← 映射 + 着色
    }
};
/*───────────────────────────  END CubeMesh  ───────────────────────*/

//---------------------------------PLANE---------------------------------------------
/*───────────────────────────  PlaneMesh  ──────────────────────────*/
class PlaneMesh : public Mesh {
public:
    PlaneMesh(float scale = 1,
              const Vec3& pos = Vec3(0, 0, 0),
              const Vec3& rot = Vec3(0, 0, 0))
        : Mesh(scale, pos, rot) {
        vertices = new List<Vec3>{
            {-0.5, 0, 0.5},
            {-0.5, 0, -0.5},
            {0.5, 0, -0.5},
            {0.5, 0, 0.5}};

        indices = new List<int>{0, 1, 2, 0, 2, 3};
        triangles = new List<Triangle>(2);
        finalize();  // ← 同样映射三角形
    }
};
/*───────────────────────────  END PlaneMesh  ──────────────────────*/
/*───────────────────────────  SphereMesh  ─────────────────────────*/
class SphereMesh : public Mesh {
public:
    /// @param subdivisions 细分级数，≥1；1 表示保持立方体拓扑
    SphereMesh(float scale = 1.f,
               const Vec3& pos = Vec3(0, 0, 0),
               const Vec3& rot = Vec3(0, 0, 0),
               int subdivisions = 3)
        : Mesh(scale, pos, rot) {
        /* ---------- 0. 立方体初始顶点 & 索引 ---------- */
        vertices = new List<Vec3>{
            // south (+Z)
            {-0.5f, -0.5f, 0.5f},
            {-0.5f, 0.5f, 0.5f},
            {0.5f, 0.5f, 0.5f},
            {0.5f, -0.5f, 0.5f},
            // north (-Z)
            {-0.5f, -0.5f, -0.5f},
            {-0.5f, 0.5f, -0.5f},
            {0.5f, 0.5f, -0.5f},
            {0.5f, -0.5f, -0.5f}};
        indices = new List<int>{
            0, 1, 2, 0, 2, 3,   // +Z
            7, 6, 5, 7, 5, 4,   // –Z
            3, 2, 6, 3, 6, 7,   // +X
            4, 5, 1, 4, 1, 0,   // –X
            1, 5, 6, 1, 6, 2,   // +Y
            3, 7, 4, 3, 4, 0};  // –Y

        /* ---------- 1. 细分 ---------- */
        for (int s = 1; s < std::max(1, subdivisions); ++s)
            subdivide();  // 每级把所有三角形四分

        /* ---------- 2. 投射到球面 ---------- */
        for (auto& v : *vertices)
            v = v.Normalized() * 0.5f;  // 半径 0.5 的球

        /* ---------- 3. 生成三角形容器并完成映射 ---------- */
        triangles = new List<Triangle>(indices->size() / 3);
        finalize();  // ← 映射 + 着色
    }

private:
    /*── 将当前 indices 所有三角形四分，写回 indices ──*/
    void subdivide() {
        auto* newIndices = new List<int>();
        std::unordered_map<uint64_t, int> midpointCache;

        auto edgeKey = [](int a, int b) -> uint64_t {
            return (static_cast<uint64_t>(std::min(a, b)) << 32) |
                   static_cast<uint32_t>(std::max(a, b));
        };

        auto getMidpoint = [&](int a, int b) -> int {
            uint64_t key = edgeKey(a, b);
            auto it = midpointCache.find(key);
            if (it != midpointCache.end())
                return it->second;

            Vec3 mid = ((*vertices)[a] + (*vertices)[b]) * 0.5f;
            int idx = vertices->size();
            vertices->push_back(mid);
            midpointCache[key] = idx;
            return idx;
        };

        /* 遍历旧三角形，按 (a,b,c) → (a,ab,ca)(ab,b,bc)(ca,bc,c)(ab,bc,ca) 四分 */
        for (size_t i = 0; i < indices->size(); i += 3) {
            int a = (*indices)[i];
            int b = (*indices)[i + 1];
            int c = (*indices)[i + 2];
            int ab = getMidpoint(a, b);
            int bc = getMidpoint(b, c);
            int ca = getMidpoint(c, a);

            newIndices->insert(newIndices->end(),
                               {a, ab, ca,
                                ab, b, bc,
                                ca, bc, c,
                                ab, bc, ca});
        }

        delete indices;
        indices = newIndices;
    }
};
/*───────────────────────────  END SphereMesh  ─────────────────────*/

//------------------------------HELPER FUNCTIONS------------------------------------------------

inline Mesh* LoadMeshFromOBJFile(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "OBJ open fail\n";
        return nullptr;
    }
    std::vector<Vec3> verts;
    std::vector<int> idx;
    std::string tok;
    while (in >> tok) {
        if (tok == "v") {
            float x, y, z;
            in >> x >> y >> z;
            verts.emplace_back(x, y, z);
        } else if (tok == "f") {
            int a, b, c;
            in >> a >> b >> c;
            idx.emplace_back(a - 1);
            idx.emplace_back(b - 1);
            idx.emplace_back(c - 1);
        } else
            std::getline(in, tok);  // 跳过其它行
    }
    auto* m = new Mesh();
    m->vertices = new List<Vec3>(verts.begin(), verts.end());
    m->indices = new List<int>(idx.begin(), idx.end());
    m->triangles = new List<Triangle>(idx.size() / 3);
    m->MapVertsToTriangles();
    return m;
}

// Mesh* LoadMeshFromOBJFile(string objFileName) {
//     if (!SPIFFS.begin(true)) {
//         Serial.println("Failed while mounting SPIFFS.");
//         return (Mesh*)nullptr;
//     }
//     //---------------------------------
//     static string filePath = "/";
//     string mtlFileName = "";
//     // ----------- Read object file -------------
//     List<string> strings;
//     string line;
//     File objFile = SPIFFS.open((filePath + objFileName).c_str(), "r");
//     if (objFile)  // if (objFile.is_open())
//     {
//         while (objFile.available()) {
//             // 1st. Gets the next line.
//             // 2nd. Seperates each word from that line then stores each word into the strings array.
//             String l = objFile.readStringUntil('\n');
//             line = l.c_str();
//             Serial.println(l);
//             string word;
//             stringstream ss(line);
//             while (getline(ss, word, ' ')) {
//                 if (word == "mtllib") {
//                     getline(ss, word, ' ');
//                     mtlFileName = word;
//                     Serial.println(String(mtlFileName.c_str()));
//                 } else {
//                     strings.emplace_back(word);
//                 }
//             }
//         }
//     }
//     objFile.close();
//     // -----------------Construct new mesh-------------------
//     Mesh* mesh = new Mesh();
//     List<Vec3>* verts = new List<Vec3>();
//     List<int>* indices = new List<int>();
//     List<Triangle>* triangles = new List<Triangle>(indices->size() / 3);
//     Material material;
//     for (size_t i = 0; i < strings.size(); i++) {
//         // v = vertex
//         // f = face
//         string objFileSubString = strings[i];
//         // if .obj encounters "usemtl" then the next string will be material id.
//         // if (objFileSubString == "usemtl")
//         // {
//         //     // Try opening Material file to see if it exists
//         //     File mtlFile = SPIFFS.open((filePath + mtlFileName).c_str());
//         //     // check if using material before looking for material key words
//         //     while (mtlFile)
//         //     {
//         //         string mtlID = strings[++i];
//         //         bool mtlIDFound = false;
//         //         //materials->emplace_back(Material(materialID));
//         //         // search .mtl for material properties under the the current materialID
//         //         if (!mtlIDFound)// && mtlFile)
//         //         {
//         //             // 1st. Gets the next line.
//         //             // 2nd. Seperates each word from that line then stores each word into the strings array.
//         //             String l = mtlFile.readStringUntil('\n');
//         //             line = l.c_str();
//         //             Serial.println(l);
//         //             Serial.println("!mtlIDFound");
//         //             string word;
//         //             stringstream ss(line);
//         //             while (!mtlIDFound && getline(ss, word, ' '))
//         //             {
//         //                 Serial.println("!mtlIDFound && getline(ss, word, ' ')");
//         //                 if (word == "newmtl")
//         //                 {
//         //                     getline(ss, word, ' ');
//         //                     if (mtlID == word) {
//         //                         material.name = word;
//         //                     }
//         //                 }
//         //                 else if (mtlID == material.name && word == "Kd") {
//         //                     getline(ss, word, ' ');
//         //                     float r = stof(word);
//         //                     getline(ss, word, ' ');
//         //                     float g = stof(word);
//         //                     getline(ss, word, ' ');
//         //                     float b = stof(word);
//         //                     material.color = Color(255*r, 255*g, 255*b);
//         //                     mtlIDFound = true;
//         //                 }
//         //             }
//         //         }
//         //         mtlFile.close();
//         //     }
//         // }
//         if (false) {
//         } else if (objFileSubString == "v") {
//             float x = stof(strings[++i]);
//             float y = stof(strings[++i]);
//             float z = stof(strings[++i]);
//             verts->emplace_back(Vec3(x, y, z));
//         }
//         // f means the next 3 strings will be the indices for mapping vertices
//         else if (objFileSubString == "f") {
//             int p3Index = stof(strings[++i]) - 1;
//             int p2Index = stof(strings[++i]) - 1;
//             int p1Index = stof(strings[++i]) - 1;
//             indices->emplace_back(p1Index);
//             indices->emplace_back(p2Index);
//             indices->emplace_back(p3Index);
//             Triangle tri = Triangle((*verts)[p1Index], (*verts)[p2Index], (*verts)[p3Index]);
//             tri.color = material.color;
//             triangles->emplace_back(tri);
//         }
//     }
//     mesh->vertices = verts;
//     mesh->indices = indices;
//     mesh->triangles = triangles;
//     return mesh;
// }
/*
Mesh* LoadMeshFromOBJFile(string objFileName)
{
    static string filePath = "./Objects/";
    string mtlFileName = "";
    std::ifstream mtlFile;
    // ----------- Read object file -------------
    List<string> strings;
    string line;
    std::ifstream objFile;
    objFile.open(filePath+objFileName);
    if (objFile.is_open())
    {
        while (objFile) {
            // 1st. Gets the next line.
            // 2nd. Seperates each word from that line then stores each word into the strings array.
            getline(objFile, line);
            string word;
            stringstream ss(line);
            while (getline(ss, word, ' '))
            {
                if (word == "mtllib") {
                    getline(ss, word, ' ');
                    mtlFileName = word;
                }
                else {
                    strings.emplace_back(word);
                }
            }
        }
    }
    objFile.close();

    // -----------------Construct new mesh-------------------
    Mesh* mesh = new Mesh();
    List<Vec3>* verts = new List<Vec3>();
    List<int>* indices = new List<int>();
    List<Triangle>* triangles = new List<Triangle>(indices->size() / 3);
    Material material;
    for (size_t i = 0; i < strings.size(); i++)
    {
        // v = vertex
        // f = face
        string objFileSubString = strings[i];

        // if .obj encounters "usemtl" then the next string will be material id.
        if (objFileSubString == "usemtl")
        {
            // Try opening Material file to see if it exists
            mtlFile.open(filePath + mtlFileName);
            // check if using material before looking for material key words
            if (mtlFile.is_open())
            {
                string mtlID = strings[++i];
                bool mtlIDFound = false;
                //materials->emplace_back(Material(materialID));
                // search .mtl for material properties under the the current materialID
                while (!mtlIDFound && mtlFile)
                {
                    // 1st. Gets the next line.
                    // 2nd. Seperates each word from that line then stores each word into the strings array.
                    getline(mtlFile, line);
                    string word;
                    stringstream ss(line);
                    while (!mtlIDFound && getline(ss, word, ' '))
                    {
                        if (word == "newmtl")
                        {
                            getline(ss, word, ' ');
                            if (mtlID == word) {
                                material.name = word;
                            }
                        }
                        else if (mtlID == material.name && word == "Kd") {
                            getline(ss, word, ' ');
                            float r = stof(word);
                            getline(ss, word, ' ');
                            float g = stof(word);
                            getline(ss, word, ' ');
                            float b = stof(word);
                            material.color = Color(255*r, 255*g, 255*b);

                            mtlIDFound = true;
                        }
                    }
                }
                mtlFile.close();
            }
        }

        else if (objFileSubString == "v") {
            float x = stof(strings[++i]);
            float y = stof(strings[++i]);
            float z = stof(strings[++i]);
            verts->emplace_back(Vec3(x, y, z));
        }
        //f means the next 3 strings will be the indices for mapping vertices
        else if (objFileSubString == "f") {
            int p3Index = stof(strings[++i]) - 1;
            int p2Index = stof(strings[++i]) - 1;
            int p1Index = stof(strings[++i]) - 1;

            indices->emplace_back(p1Index);
            indices->emplace_back(p2Index);
            indices->emplace_back(p3Index);

            Triangle tri = Triangle((*verts)[p1Index], (*verts)[p2Index], (*verts)[p3Index]);
            tri.color = material.color;
            triangles->emplace_back(tri);
        }
    }
    //mtlFile.close();

    mesh->vertices = verts;
    mesh->indices = indices;
    mesh->triangles = triangles;

    return mesh;
}
*/

inline void Draw() {
    worldToViewMatrix = Camera::main->TRInverse();
    projectionMatrix = ProjectionMatrix();
    Matrix4x4 vp = projectionMatrix * worldToViewMatrix;

    for (int i = 0; i < Mesh::count; ++i) {
        if (Graphics::frustumCulling) {
            float sd = (Mesh::objects[i]->root->localPosition - Camera::main->Position()).SqrMagnitude();
            if (sd != 0.0f && Mesh::objects[i]->root->localScale.SqrMagnitude() / sd < 1e-13)
                continue;
        }
        Mesh::objects[i]->TransformTriangles();
    }
    sort(triBuffer->begin(), triBuffer->end(),
         [](const Triangle& a, const Triangle& b) { return a.centroid.w > b.centroid.w; });

    for (auto& t : *triBuffer) t.Draw();
    for (auto& l : *lineBuffer) l.Draw();
    for (auto& p : *pointBuffer) p.Draw();

    pointBuffer->clear();
    lineBuffer->clear();
    triBuffer->clear();
}

#endif /* GRAPHICS_H */
