#pragma once
#include <Arduino.h>
#include <TFT_eSPI.h>  // 用于真正把像素推到 ILI9341（DMA）
#include <lvgl.h>

#define LV_BUFFER_ROWS 1

/**
 * 超轻量 LVGL 显示驱动（行缓冲刷新）
 * ------------------------------------------------
 * - 只占 w × ROWS × 2 × 2 字节 RAM
 * - 不再给用户暴露 FrameBuffer，而是让用户直接用 LVGL API 画
 */
class SimpleFBDisplay {
public:
    /** 在 setup() 调一次；rot=0~3 同 ILI9341 setRotation */
    static void begin(uint16_t w, uint16_t h, uint8_t rot = 1);

    /** 每帧调用，维持 LVGL 心跳 */
    static void tick();
    // SimpleFBDisplay.h 里新增
    /** 把一幅 16-bit RGB565 帧缓冲整屏推送到 LCD
     *  src 指向 w×h 个 Pixel565（与 begin() 的 w/h 一致）
     */
    static void pushFrame(const uint16_t *src);


private:
    // ------ 屏幕 SPI 驱动 ------
    static TFT_eSPI tft;

    // ------ LVGL draw buffer ------
    // static constexpr uint16_t LV_BUFFER_ROWS = 20;
    static lv_disp_draw_buf_t drawBuf;
    static lv_color_t *bufA;
    static lv_color_t *bufB;

    // ------ LVGL drivers ------
    static lv_disp_drv_t dispDrv;

    // ------ 内部回调 ------
    static void disp_flush(lv_disp_drv_t *disp,
                           const lv_area_t *area,
                           lv_color_t *color_p);

    // ------ 计时 ------
    static uint32_t lastTick;
};
