#pragma once
#include "Graphics.h"
#include "Planet/EarthPlanet.h"

class PlanetManager {
public:
    enum PlanetType { Earth /*, Mars, Jupiter …*/ };

    explicit PlanetManager(PlanetType t = Earth) {
        loadPlanet(t);
    }

    ~PlanetManager() { unload(); }

    /* 改变行星类型 */
    void loadPlanet(PlanetType t) {
        unload();
        switch (t) {
        case Earth:
            _planet = new EarthPlanet();
            break;
        default:
            _planet = new EarthPlanet();
            break;
        }
    }

    /* 每帧调用：自转 & 动画 */
    void update(float dt) {
        if (!_planet)
            return;
        _planet->update(dt);
    }

private:
    void unload() {
        if (_planet) {
            delete _planet;
            _planet = nullptr;
        }
    }

    Planet* _planet = nullptr;
};
