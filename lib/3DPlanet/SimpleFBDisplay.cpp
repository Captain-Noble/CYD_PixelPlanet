// SimpleFBDisplay.cpp  顶部
#include "SimpleFBDisplay.h"
#include <lvgl.h>                 // 确保先包含

#define LV_TICK_CUSTOM 1

#if LV_TICK_CUSTOM               // 若用户开启了自定义 tick
    extern "C" void lv_tick_inc(uint32_t ms);   // 手动声明一下即可 
#endif

/* ========== 静态成员定义 ========== */
TFT_eSPI SimpleFBDisplay::tft;
lv_disp_draw_buf_t SimpleFBDisplay::drawBuf;
lv_color_t *SimpleFBDisplay::bufA = nullptr;
lv_color_t *SimpleFBDisplay::bufB = nullptr;
lv_disp_drv_t SimpleFBDisplay::dispDrv;
uint32_t SimpleFBDisplay::lastTick = 0;

/* ========== 公有接口 ========== */
void SimpleFBDisplay::begin(uint16_t w, uint16_t h, uint8_t rot) {
    /* 1. 初始化底层 TFT_eSPI */
    tft.begin();
    tft.setRotation(rot);
    tft.initDMA();  // 启用 DMA
    tft.fillScreen(TFT_BLACK);

    /* 2. 准备 LVGL */
    lv_init();

    /* 3. 分配双缓冲：w × ROWS */
    size_t lineBytes = w * LV_BUFFER_ROWS * sizeof(lv_color_t);
    bufA = static_cast<lv_color_t *>(heap_caps_malloc(lineBytes, MALLOC_CAP_DMA));
    bufB = static_cast<lv_color_t *>(heap_caps_malloc(lineBytes, MALLOC_CAP_DMA));

    lv_disp_draw_buf_init(&drawBuf, bufA, bufB, w * LV_BUFFER_ROWS);

    /* 4. 注册显示驱动 */
    lv_disp_drv_init(&dispDrv);
    dispDrv.hor_res = w;
    dispDrv.ver_res = h;
    dispDrv.flush_cb = disp_flush;
    dispDrv.draw_buf = &drawBuf;
    lv_disp_drv_register(&dispDrv);

    lastTick = millis();
}

void SimpleFBDisplay::tick()
{
    uint32_t now = millis();
    #if LV_TICK_CUSTOM == 0          // 官方 tick 口子可直接用
        lv_tick_inc(now - lastTick);
    #endif
    lastTick = now;
    lv_timer_handler();
}


/* ========== LVGL → LCD 刷新回调 ========== */
void SimpleFBDisplay::disp_flush(lv_disp_drv_t *disp,
                                 const lv_area_t *area,
                                 lv_color_t *color_p) {
    uint32_t w = area->x2 - area->x1 + 1;
    uint32_t h = area->y2 - area->y1 + 1;
    uint32_t len = w * h;

    tft.startWrite();
    tft.setAddrWindow(area->x1, area->y1, w, h);
    tft.pushPixelsDMA(reinterpret_cast<uint16_t *>(&color_p->full), len);
    tft.endWrite();

    lv_disp_flush_ready(disp);
}

// SimpleFBDisplay.cpp 里实现（放在 disp_flush() 之后即可）
void SimpleFBDisplay::pushFrame(const uint16_t *src)
{
    const uint16_t w = dispDrv.hor_res;
    const uint16_t h = dispDrv.ver_res;

    tft.startWrite();
    for (uint16_t y = 0; y < h; ++y)
    {
        tft.setAddrWindow(0, y, w, 1);
        // ---- 关键修改行 ----
        tft.pushPixelsDMA(const_cast<uint16_t*>(src + y * w), w);
    }
    tft.endWrite();
}
