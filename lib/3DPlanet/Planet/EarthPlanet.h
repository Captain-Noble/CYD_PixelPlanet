#pragma once
#include "../Layer/PlanetLayers.h"
#include "Graphics.h"

/*-------------------------------------------------
 * Planet 基类：只是集合体
 *------------------------------------------------*/
class Planet : public Transform {
public:
    virtual ~Planet() {
        for (auto* m : layers) delete m;
    }
    virtual void update(float dt) {
        /* 简单自转：绕 Y 轴 15°/s */
        // localRotation = Matrix3x3::RotY(15.f * PI / 180.f * dt) * localRotation;
    }

protected:
    std::vector<Mesh*> layers;
};

/*-------------------------------------------------
 * Earth 实现
 *------------------------------------------------*/
class EarthPlanet : public Planet {
public:
    EarthPlanet() {
        using namespace PlanetLayerFactory;
        /* 注意层次顺序：先画深色，后画浅色 */
        layers.push_back(createBaseOceanLayer(*this));
        layers.push_back(createLandMassLayer(*this, 0.55f));
        layers.push_back(createCloudLayer(*this));
    }
};
