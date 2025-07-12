// ============================================================================
// main.cpp  --  Software 3D renderer → SimpleFBDisplay (ESP32, 单核版)
// ============================================================================

#include <Arduino.h>
#include <cstring>

#include "Graphics.h"          // 纯软件 3D 渲染头
#include "SimpleFBDisplay.h"

// 显存全局定义（Graphics.h 里是 extern）
extern Pixel565 display[DISPLAY_HEIGHT_PIX][DISPLAY_WIDTH_PIX];

// 场景：一个正方体
SphereMesh *sphere = nullptr;

void setup()
{
    Serial.begin(115200);

    // 初始化 LCD（行缓冲驱动，横屏）
    SimpleFBDisplay::begin(DISPLAY_WIDTH_PIX, DISPLAY_HEIGHT_PIX, /*rot=*/1);

    // 创建一个位于 (0,0,0) 的立方体
    sphere = new SphereMesh(1.0f, Vec3(0, 0, 0), Vec3(0, 0, 0), 3);

    // 把主相机往后移
    Camera::main->localPosition = Vec3(0, 0,1.0f/FIELD_OF_VIEW_DEG * 60);

    // 渲染选项
    Graphics::fillTriangles    = true;
    Graphics::lighting         = true;
    Graphics::displayWireFrames = false;
}

void loop()
{
    /* ── FPS 统计 ──────────────────────────────*/
    static uint32_t lastMicros   = micros();
    static uint32_t lastPrintMs  = millis();
    static uint32_t frameCounter = 0;

    uint32_t nowMicros = micros();
    float dt           = (nowMicros - lastMicros) / 1e6f;
    lastMicros         = nowMicros;

    /* ── 更新场景 ──────────────────────────────*/
    if (sphere)
        sphere->localRotation = Matrix3x3::RotY(1.0f * dt) * sphere->localRotation;

    /* ── 清空软件帧缓冲 ────────────────────────*/
    std::memset(display, 0, sizeof(display));  // 全黑

    /* ── 渲染 ─────────────────────────────────*/
    Draw();  // 往 display[][] 写像素 (RGB565)

    /* ── 整屏 DMA 推送 ────────────────────────*/
    SimpleFBDisplay::pushFrame(reinterpret_cast<uint16_t *>(display));

    /* ── LVGL 心跳 ────────────────────────────*/
    SimpleFBDisplay::tick();

    /* ── FPS 打印 ─────────────────────────────*/
    ++frameCounter;
    uint32_t nowMs = millis();
    if (nowMs - lastPrintMs >= 1000)
    {
        Serial.printf("FPS: %u\n", frameCounter);
        frameCounter = 0;
        lastPrintMs  = nowMs;
    }
}
